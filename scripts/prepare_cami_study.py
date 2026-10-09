#!/usr/bin/env python3
"""Prepare a CAMI2-derived amplicon study from local original genome references."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genome_reference import extract_markers, RANKS

SITES = ('Airways', 'Gastrointestinal', 'Oral', 'Skin', 'Urogenital')
FORWARD = 'GTGYCAGCMGCCGCGGTAA'
REVERSE = 'GGACTACNVGGGTWTCTAAT'


def scan(path):
    path = Path(path)
    hits = []
    for rec in SeqIO.parse(path, 'fasta'):
        for sequence, strand, start, end in extract_markers(str(rec.seq), FORWARD, REVERSE, 230, 310):
            hits.append((sequence, rec.id, strand, start, end))
    return str(path), hashlib.sha256(path.read_bytes()).hexdigest(), hits


def lineage_consensus(lineages):
    result = []
    for i in range(len(RANKS)):
        values = {lineage[i] for lineage in lineages}
        if len(values) != 1 or '' in values:
            break
        result.append(values.pop())
    return result + ['']*(len(RANKS)-len(result))


def prepare(source, output, threads=4):
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists():
        raise ValueError('Choose a new preparation directory; existing inputs are never overwritten')
    compositions = {}; sources = {}; taxonomy = {}; source_files = set(); samples = []
    # Preflight all original genome references before creating output.
    for site in SITES:
        short = source/'MetaGs'/f'CAMI_II_{site}'/'short_read'
        sag = source/'SAGs'/f'CAMI_II_{site}'
        mapping = dict(csv.reader((sag/'genome_to_id.tsv').open(), delimiter='\t'))
        ranks = {}
        with (sag/'genome_taxa_info.tsv').open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                lineage = row['TAXPATHSN'].split('|')[:7]
                ranks[row['_CAMI_genomeID']] = (lineage+['']*7)[:7]
        tables = sorted(short.glob('abundance*.tsv'))
        if not tables: raise ValueError(f'No abundance tables: {short}')
        for table in tables:
            sample = f'{site}_{table.stem.removeprefix("abundance")}'
            samples.append(dict(sample_id=sample, group=site, body_site=site, Type_Group=site))
            values = {}
            with table.open() as handle:
                for gid, abundance in csv.reader(handle, delimiter='\t'):
                    value = float(abundance)
                    if not 0 <= value < float('inf'): raise ValueError(f'Invalid abundance: {table}: {gid}')
                    if not value: continue
                    key = f'{site}:{gid}'
                    path = sag/'fasta'/Path(mapping[gid]).name
                    if not path.is_file(): raise ValueError(f'Missing original genome: {path}')
                    sources[key] = path
                    taxonomy[key] = ranks.get(gid, ['']*7)
                    values[key] = value
            if not values: raise ValueError(f'Empty source composition: {table}')
            compositions[sample] = values
            source_files.add(table)
        source_files.update([sag/'genome_to_id.tsv', sag/'genome_taxa_info.tsv'])
    # Identical accession filenames in site collections select the first local copy;
    # every selected source file is hashed and its exact path is recorded.
    chosen_paths = {}
    for path in sources.values(): chosen_paths.setdefault(path.name, path)
    print(f'Preparing {len(samples)} samples; scanning {len(chosen_paths)} genome FASTAs with {threads} workers.', flush=True)
    extracted = {}
    with ProcessPoolExecutor(max_workers=threads) as pool:
        for i,(path,sha,hits) in enumerate(pool.map(scan, map(str,chosen_paths.values())),1):
            extracted[Path(path).name] = (path,sha,hits)
            if i % 50 == 0: print(f'Genomes scanned: {i}/{len(chosen_paths)}', flush=True)
    for path in set(sources.values()):
        actual, sha, _ = extracted[path.name]
        if str(path) != actual and hashlib.sha256(path.read_bytes()).hexdigest() != sha:
            raise ValueError(f'Conflicting copies of reference accession: {path.name}')
    output.mkdir(parents=True)
    ref = output/'reference'; ref.mkdir()
    sequences = {}; lineages = defaultdict(list); associations = []; audits = []; copies = {}
    for key,path in sorted(sources.items()):
        actual,sha,hits = extracted[path.name]
        counts = Counter()
        for sequence,contig,strand,start,end in hits:
            aid = 'ASV_'+hashlib.sha256(sequence.encode()).hexdigest()[:16]
            sequences[aid] = sequence; counts[aid] += 1
            associations.append(dict(source_genome=key, source_fasta=actual, source_sha256=sha, ASV_ID=aid,
                                     contig=contig, strand=strand, start=start, end=end))
        for aid in counts: lineages[aid].append(taxonomy[key])
        copies[key] = counts
        audits.append(dict(source_genome=key, source_fasta=actual, source_sha256=sha,
                           marker_loci=sum(counts.values()), unique_asvs=len(counts),
                           status='retained' if counts else 'no_exact_primer_bounded_V4'))
    weights = pd.DataFrame(0.0, index=sorted(sequences), columns=list(compositions))
    retained = []
    for sample, values in compositions.items():
        total = sum(values.values()); eligible = sum(v for g,v in values.items() if copies[g])
        for gid,value in values.items():
            for aid,n in copies[gid].items(): weights.loc[aid,sample] += value*n
        if not weights[sample].sum(): raise ValueError(f'{sample}: no eligible V4 loci')
        retained.append(dict(sample_id=sample, source_genomes=len(values), genomes_with_markers=sum(bool(copies[g]) for g in values),
                             retained_source_abundance_percent=100*eligible/total))
    weights /= weights.sum()
    weights.index.name = 'ASV_ID'
    weights.to_csv(ref/'amplicon_abundances.tsv',sep='\t')
    pd.DataFrame(compositions).fillna(0).rename_axis('source_genome').to_csv(ref/'source_genome_abundances.tsv',sep='\t')
    SeqIO.write([SeqRecord(Seq(seq),id=aid,description='') for aid,seq in sorted(sequences.items())], ref/'v4_asv_pool.fasta','fasta')
    taxrows = []
    for aid in sorted(sequences):
        taxrows.append(dict(ASV_ID=aid, **dict(zip(RANKS,lineage_consensus(lineages[aid]))),
                            representative_reference_id='CAMI2_shared_ASV', n_identical_references=len(lineages[aid])))
    pd.DataFrame(taxrows).to_csv(ref/'v4_asv_taxonomy.tsv',sep='\t',index=False)
    pd.DataFrame(associations).to_csv(ref/'genome_marker_loci.tsv',sep='\t',index=False)
    pd.DataFrame(audits).to_csv(ref/'genome_selection.tsv',sep='\t',index=False)
    pd.DataFrame(retained).to_csv(ref/'sample_reference_retention.tsv',sep='\t',index=False)
    metadata = dict(source_type='cami2_amplicon', n_unique_v4=len(sequences), n_samples=len(samples),
                    primers=dict(forward=FORWARD,reverse=REVERSE),
                    abundance_basis='CAMI supplied genome abundance multiplied by exact primer-bounded allele copy count; renormalized over eligible loci',
                    source_files={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source_files)},
                    limitations='Exact primer matches only; missing loci excluded with explicit audit. Shared ASVs collapsed; taxonomy is common lineage. Chemicals are simulated, not CAMI measurements.')
    (ref/'reference_manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (output/'study.yaml').write_text(yaml.safe_dump(dict(study_name='CAMI2_body_site_amplicon_chemistry',samples=samples),sort_keys=False))
    template = Path(__file__).resolve().parents[1]/'config/cami_body_sites.yaml'
    cfg = yaml.safe_load(template.read_text())
    cfg['simulation'].update(n_samples=len(samples), n_asvs=len(sequences))
    cfg['study_design_file'] = str(output/'study.yaml')
    cfg['reference']['output_dir'] = str(ref)
    (output/'config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
    print(f'Ready: {output}/config.yaml; {len(samples)} samples, {len(sequences)} ASVs. No FASTQs generated.',flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('-i','--input',required=True,type=Path);p.add_argument('-o','--output',required=True,type=Path)
    p.add_argument('--threads',type=int,default=4);a=p.parse_args()
    if a.threads<1:p.error('--threads must be positive')
    prepare(a.input,a.output,a.threads)
