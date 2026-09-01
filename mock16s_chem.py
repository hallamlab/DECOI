#!/usr/bin/env python3
"""DECOI: generate a truth-aware V4 16S and chemistry benchmark study."""
from __future__ import annotations

import argparse, csv, gzip, hashlib, html, json, logging, platform, shutil, subprocess, sys, urllib.request, uuid
from pathlib import Path
from typing import TextIO

import numpy as np
import pandas as pd
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

LOG = logging.getLogger("mock16s_chem")
RANKS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
IUPAC_BASES = {"A":"A","C":"C","G":"G","T":"T","R":"AG","Y":"CT","S":"GC","W":"AT","K":"GT","M":"AC","B":"CGT","D":"AGT","H":"ACT","V":"ACG","N":"ACGT"}


def open_text(path: Path, mode: str = "rt") -> TextIO:
    return gzip.open(path, mode) if path.suffix == ".gz" else path.open(mode)


def run(cmd: list[str], capture: bool = False) -> str:
    LOG.info("Running: %s", " ".join(map(str, cmd)))
    p = subprocess.run(cmd, check=True, text=True, capture_output=capture)
    return p.stdout.strip() if capture else ""


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def taxonomy_from_header(description: str) -> list[str]:
    parts = [x.strip() for x in description.lstrip(">").strip().split(";") if x.strip()]
    if len(parts) < 2:
        raise ValueError("Header does not resemble a DADA2 taxonomy header")
    return (parts + [""] * len(RANKS))[:len(RANKS)]


def stable_asv_id(sequence: str) -> str:
    return "ASV_" + hashlib.sha256(sequence.encode()).hexdigest()[:16]


def stable_uuid(namespace: str, value: str) -> str:
    return str(uuid.uuid5(uuid.uuid5(uuid.NAMESPACE_URL, namespace), value))


def resolve_iupac(sequence: str, key: str) -> str:
    digest = hashlib.sha256(key.encode()).digest(); out = []
    for i, char in enumerate(sequence.upper()):
        choices = IUPAC_BASES.get(char)
        if choices is None:
            raise ValueError(f"Unsupported IUPAC base {char!r}")
        out.append(choices[digest[i % len(digest)] % len(choices)])
    return "".join(out)


def iss_amplicon_sequence(insert: str, include_primers: bool, forward: str, reverse: str, key: str) -> str:
    if not include_primers:
        return insert
    fwd = resolve_iupac(forward, key + ":F")
    rev = resolve_iupac(reverse, key + ":R")
    return fwd + insert + str(Seq(rev).reverse_complement())


def prepare_reference(cfg: dict) -> None:
    rcfg = cfg["reference"]; source = Path(rcfg["dada2_silva_fasta"]); outdir = Path(rcfg["output_dir"])
    outdir.mkdir(parents=True, exist_ok=True)
    if not source.exists(): raise FileNotFoundError(source)
    if not shutil.which("cutadapt"): raise RuntimeError("cutadapt is required")
    numbered = outdir / "silva_numbered_full_length.fasta"; raw_tax = outdir / "silva_numbered_taxonomy.tsv"
    n = 0
    with open_text(source) as ih, numbered.open("w") as fo, raw_tax.open("w", newline="") as to:
        w = csv.writer(to, delimiter="\t"); w.writerow(["reference_id", *RANKS, "source_header"])
        for rec in SeqIO.parse(ih, "fasta"):
            try: tax = taxonomy_from_header(rec.description)
            except ValueError: continue
            n += 1; rid = f"SILVA_{n:09d}"
            SeqIO.write(SeqRecord(rec.seq, id=rid, description=""), fo, "fasta"); w.writerow([rid, *tax, rec.description])
    rev_rc = str(Seq(rcfg["reverse_primer"]).reverse_complement())
    trimmed = outdir / "v4_trimmed.fasta"
    run(["cutadapt", "-g", f'{rcfg["forward_primer"]};required...{rev_rc};required', "--discard-untrimmed",
         "--minimum-length", str(rcfg.get("min_length", 230)), "--maximum-length", str(rcfg.get("max_length", 310)),
         "--info-file", str(outdir / "cutadapt.info.tsv"), "-o", str(trimmed), str(numbered)])
    tax = pd.read_csv(raw_tax, sep="\t").set_index("reference_id"); seq_to_refs: dict[str,list[str]] = {}
    for rec in SeqIO.parse(trimmed, "fasta"):
        seq = str(rec.seq).upper().replace("U", "T")
        if set(seq) <= set("ACGT"): seq_to_refs.setdefault(seq, []).append(rec.id)
    with (outdir/"v4_asv_pool.fasta").open("w") as fa, (outdir/"v4_asv_taxonomy.tsv").open("w",newline="") as to, (outdir/"v4_asv_source_map.tsv").open("w",newline="") as mo:
        tw=csv.writer(to,delimiter="\t"); mw=csv.writer(mo,delimiter="\t")
        tw.writerow(["ASV_ID",*RANKS,"representative_reference_id","n_identical_references"]); mw.writerow(["ASV_ID","reference_id"])
        for seq, refs in sorted(seq_to_refs.items()):
            aid=stable_asv_id(seq); rep=refs[0]; row=tax.loc[rep]
            SeqIO.write(SeqRecord(Seq(seq),id=aid,description=""),fa,"fasta")
            tw.writerow([aid,*[row.get(r,"") for r in RANKS],rep,len(refs)])
            for rid in refs: mw.writerow([aid,rid])
    metadata={"source":str(source),"source_sha256":file_sha256(source),"n_source_records":n,"n_unique_v4":len(seq_to_refs),"primers":dict(forward=rcfg["forward_primer"],reverse=rcfg["reverse_primer"])}
    (outdir/"reference_manifest.json").write_text(json.dumps(metadata,indent=2)+"\n")


def load_reference(cfg: dict):
    d=Path(cfg["reference"]["output_dir"])
    tax=pd.read_csv(d/"v4_asv_taxonomy.tsv",sep="\t").set_index("ASV_ID")
    seqs={r.id:str(r.seq).upper() for r in SeqIO.parse(d/"v4_asv_pool.fasta","fasta")}
    return tax,seqs


def load_study_samples(cfg: dict) -> pd.DataFrame:
    path=cfg.get("study_design_file")
    if not path:
        n=int(cfg["simulation"]["n_samples"]); return pd.DataFrame({"sample_id":[f"sample_{i:03d}" for i in range(1,n+1)]})
    study=yaml.safe_load(Path(path).read_text()); rows=[]
    for cohort in study.get("cohorts",[]):
        fixed=dict(cohort.get("metadata",{})); cohort_name=str(cohort["name"])
        metadata_cycle=list(cohort.get("participant_metadata_cycle",[]))
        participant_prefix=str(cohort.get("participant_prefix",cohort_name.upper()))
        for participant_index in range(1,int(cohort["n_participants"])+1):
            participant_id=f"{participant_prefix}_{participant_index:03d}"
            participant_metadata=dict(metadata_cycle[(participant_index-1)%len(metadata_cycle)]) if metadata_cycle else {}
            for type_spec in cohort.get("sample_types",[]):
                if isinstance(type_spec,str): type_spec={"name":type_spec}
                type_name=str(type_spec["name"]); type_code=str(type_spec.get("code",type_name)).replace(" ","_")
                for replicate in range(1,int(type_spec.get("replicates",1))+1):
                    sample_id=f"{participant_id}_{type_code}"
                    if int(type_spec.get("replicates",1))>1: sample_id+=f"_{replicate}"
                    rows.append({"sample_id":sample_id,"group":cohort_name,"Participant_ID":participant_id,
                                 "Case":fixed.get("Case",cohort_name),"Type_Group":type_name,
                                 **fixed,**participant_metadata,**dict(type_spec.get("metadata",{}))})
    for g in study.get("groups",[]):
        fixed=dict(g.get("metadata",{})); name=str(g["name"]); prefix=str(g.get("sample_prefix",name))
        for i in range(1,int(g["n_samples"])+1): rows.append({"sample_id":f"{prefix}_{i:03d}","group":name,**fixed})
    if not rows: raise ValueError("No samples in study design")
    df=pd.DataFrame(rows)
    if df.sample_id.duplicated().any(): raise ValueError("Duplicate sample IDs")
    return df


def choose_asvs(cfg,tax,seqs,rng,n=None,exclude=None):
    filt=cfg["simulation"].get("selection",{}); c=tax.copy(); exclude=set(exclude or [])
    if filt.get("domains"): c=c[c.Domain.isin(filt["domains"])]
    if int(filt.get("min_taxonomy_ranks",0)): c=c[c[RANKS].notna().sum(axis=1)>=int(filt["min_taxonomy_ranks"])]
    ids=[x for x in c.index.intersection(seqs) if x not in exclude]; n=int(n or cfg["simulation"]["n_asvs"])
    if len(ids)<n: raise ValueError(f"Only {len(ids)} reference ASVs available; requested {n}")
    return rng.choice(np.array(ids),size=n,replace=False).tolist()


def run_sparsedossa(cfg,outdir):
    s=cfg["simulation"]; script=Path(__file__).parent/"scripts/simulate_sparsedossa2.R"
    c=outdir/"_sparsedossa_counts.tsv"; r=outdir/"_sparsedossa_relative.tsv"
    run(["Rscript",str(script),
         "--n-samples",str(s["n_samples"]),
         "--n-features",str(s["n_asvs"]),
         "--median-depth",str(s["median_read_depth"]),
         "--template",str(s.get("sparsedossa_template","Stool")),
         "--seed",str(cfg["seed"]),
         "--counts-out",str(c),
         "--rel-out",str(r)])
    return pd.read_csv(c,sep="\t",index_col=0),pd.read_csv(r,sep="\t",index_col=0)


def renormalize_counts(weights: pd.DataFrame, depths: pd.Series, rng) -> pd.DataFrame:
    out=pd.DataFrame(0,index=weights.index,columns=weights.columns,dtype=int)
    for sample in weights.columns:
        w=np.clip(weights[sample].to_numpy(float),0,None); p=w/w.sum() if w.sum() else np.repeat(1/len(w),len(w))
        out[sample]=rng.multinomial(int(depths[sample]),p)
    return out


def apply_group_effects(counts,meta,cfg,rng):
    effects=cfg.get("artifacts",{}).get("group_differential_abundance",{}); truth=[]; truth_columns=["grouping_column","group","ASV_ID","log_fold_change","baseline_prevalence"]
    columns=effects.get("group_columns",["group"])
    if not effects.get("enabled",False): return counts.copy(),pd.DataFrame(columns=truth_columns)
    w=counts.astype(float)+1e-9;n=int(effects.get("asvs_per_group",5));sd=float(effects.get("log_fold_change_sd",1.0));minimum=float(effects.get("min_abs_log_fold_change",0.0));disjoint=bool(effects.get("disjoint_drivers",False));min_prevalence=float(effects.get("candidate_min_prevalence",0.0))
    prevalence=(counts>0).mean(axis=1); mean_abundance=counts.mean(axis=1)
    ranked=pd.DataFrame({"prevalence":prevalence,"mean_abundance":mean_abundance}).sort_values(["prevalence","mean_abundance"],ascending=False).index.tolist()
    preferred=[aid for aid in ranked if prevalence[aid]>=min_prevalence]; candidates=preferred+[aid for aid in ranked if aid not in preferred]; used=set()
    for column in columns:
        if column not in meta: raise ValueError(f"Differential-abundance column not found in metadata: {column}")
        for group in sorted(meta[column].unique()):
            available=[aid for aid in candidates if not disjoint or aid not in used]
            if len(available)<n: raise ValueError(f"Only {len(available)} ASVs available for {n} disjoint drivers in {column}={group}")
            pool=available[:max(n,min(len(available),n*3))];drivers=rng.choice(np.asarray(pool),size=n,replace=False);used.update(drivers)
            betas=minimum+np.abs(rng.normal(0,sd,len(drivers)))
            samples=meta.loc[meta[column]==group,"sample_id"]
            for aid,b in zip(drivers,betas): w.loc[aid,samples]*=np.exp(b);truth.append((column,group,aid,float(b),float(prevalence[aid])))
    return renormalize_counts(w,counts.sum(axis=0),rng),pd.DataFrame(truth,columns=truth_columns)


def apply_microbiome_batch_effects(counts,meta,cfg,rng):
    effects=cfg.get("artifacts",{}).get("microbiome_batch_effect",{}); columns=["batch_column","batch","ASV_ID","log_fold_change","direction"]
    if not effects.get("enabled",False): return counts.copy(),pd.DataFrame(columns=columns)
    column=str(effects.get("column","batch"))
    if column not in meta: raise ValueError(f"Microbiome batch-effect column not found in metadata: {column}")
    batches=sorted(meta[column].dropna().astype(str).unique())
    if len(batches)<2: raise ValueError(f"Microbiome batch effects require at least two values in metadata column '{column}'")
    n=min(int(effects.get("asvs_per_batch",10)),len(counts)); magnitude=float(effects.get("log_fold_change",1.25)); direction=str(effects.get("direction","mixed")).lower()
    if magnitude<=0: raise ValueError("microbiome_batch_effect.log_fold_change must be positive")
    if direction not in {"positive","negative","mixed"}: raise ValueError("microbiome_batch_effect.direction must be positive, negative, or mixed")
    w=counts.astype(float)+1e-9; truth=[]
    for batch in batches:
        drivers=rng.choice(counts.index,size=n,replace=False)
        signs=np.ones(n) if direction=="positive" else (-np.ones(n) if direction=="negative" else rng.choice([-1.0,1.0],size=n))
        samples=meta.loc[meta[column].astype(str)==batch,"sample_id"]
        for aid,sign in zip(drivers,signs):
            effect=float(sign*magnitude); w.loc[aid,samples]*=np.exp(effect); truth.append((column,batch,aid,effect,"higher" if effect>0 else "lower"))
    return renormalize_counts(w,counts.sum(axis=0),rng),pd.DataFrame(truth,columns=columns)


def apply_microbial_network(counts,cfg,rng):
    """Implant disjoint latent abundance modules for network-recovery testing."""
    settings=cfg.get("artifacts",{}).get("microbial_network",{})
    columns=["module","ASV_ID","loading","baseline_prevalence"]
    if not settings.get("enabled",False): return counts.copy(),pd.DataFrame(columns=columns)
    n_modules=int(settings.get("n_modules",4)); per_module=int(settings.get("asvs_per_module",6)); loading=float(settings.get("loading",1.0)); min_prev=float(settings.get("candidate_min_prevalence",0.5))
    if n_modules<1 or per_module<2 or loading<=0: raise ValueError("microbial_network requires n_modules >= 1, asvs_per_module >= 2, and loading > 0")
    prevalence=(counts>0).mean(axis=1); abundance=counts.mean(axis=1)
    candidates=pd.DataFrame({"prevalence":prevalence,"abundance":abundance}).query("prevalence >= @min_prev").sort_values(["prevalence","abundance"],ascending=False).index.tolist()
    needed=n_modules*per_module
    if len(candidates)<needed: raise ValueError(f"microbial_network requires {needed} ASVs with prevalence >= {min_prev}; found {len(candidates)}")
    selected=rng.choice(np.asarray(candidates),size=needed,replace=False)
    w=counts.astype(float)+1e-9; truth=[]
    for module in range(n_modules):
        members=selected[module*per_module:(module+1)*per_module]
        latent=rng.normal(0,1,len(counts.columns))
        member_loadings=rng.uniform(0.85*loading,1.15*loading,len(members))
        for aid,member_loading in zip(members,member_loadings):
            w.loc[aid]*=np.exp(member_loading*latent)
            truth.append((f"M{module+1}",aid,float(member_loading),float(prevalence[aid])))
    return renormalize_counts(w,counts.sum(axis=0),rng),pd.DataFrame(truth,columns=columns)


def apply_pcr_bias(counts,cfg,rng):
    a=cfg.get("artifacts",{}).get("pcr_bias",{})
    if not a.get("enabled",False): return counts.copy(),pd.DataFrame({"ASV_ID":counts.index,"pcr_efficiency":1.0})
    eff=rng.lognormal(float(a.get("log_mean",0)),float(a.get("log_sd",0.35)),len(counts)); w=counts.mul(eff,axis=0)
    return renormalize_counts(w,counts.sum(axis=0),rng),pd.DataFrame({"ASV_ID":counts.index,"pcr_efficiency":eff})


def add_contaminants(counts,chosen,tax,seqs,cfg,rng):
    a=cfg.get("artifacts",{}).get("contaminants",{}); truth=[]; out=counts.copy()
    if not a.get("enabled",False): return out,truth
    ids=choose_asvs(cfg,tax,seqs,rng,n=int(a.get("n_asvs",3)),exclude=chosen)
    prevalence=float(a.get("prevalence",0.5)); mean=float(a.get("mean_reads",50))
    for aid in ids:
        vals=[]
        for _ in out.columns: vals.append(int(rng.poisson(mean)) if rng.random()<prevalence else 0)
        out.loc[aid]=vals; truth.append(aid)
    return out,truth


def prepare_mitochondrial_features(cfg, outdir):
    a=cfg.get("reference_fixtures",{}).get("mitochondria",{})
    empty_tax=pd.DataFrame(columns=[*RANKS,"representative_reference_id"])
    empty_truth=pd.DataFrame(columns=["ASV_ID","source","source_record","description","sequence_length","sequence_sha256"])
    if not a.get("enabled",False): return [],{},empty_tax,empty_truth
    source_fasta=a.get("source_fasta"); source_url=a.get("source_url")
    refs=outdir/"references"; refs.mkdir(parents=True,exist_ok=True)
    if source_fasta:
        source=Path(source_fasta).expanduser()
        if not source.is_absolute(): source=Path(__file__).resolve().parent/source
        source=source.resolve()
        if not source.is_file(): raise FileNotFoundError(f"Mitochondrial source FASTA not found: {source}")
        source_label=str(source)
    elif source_url:
        source=refs/"mitochondria_source.fasta"
        LOG.info("Downloading mitochondrial source: %s",source_url)
        urllib.request.urlretrieve(source_url,source)
        source_label=source_url
    else:
        raise ValueError("Enabled mitochondrial simulation requires source_fasta or source_url")
    with open_text(source) as handle: records=list(SeqIO.parse(handle,"fasta"))
    if not records: raise ValueError(f"No records found in mitochondrial source FASTA: {source}")
    limit=max(1,int(a.get("n_asvs",1))); expected=a.get("expected_sequence_sha256"); ids=[]; seqs={}; tax_rows=[]; truth=[]
    for record in records[:limit]:
        sequence=str(record.seq).upper().replace("U","T")
        if not sequence or set(sequence)-set("ACGT"): raise ValueError(f"Mitochondrial record {record.id} must contain only A/C/G/T")
        digest=hashlib.sha256(sequence.encode()).hexdigest()
        if expected and digest.lower()!=str(expected).lower(): raise ValueError(f"Mitochondrial sequence checksum mismatch for {record.id}: expected {expected}, observed {digest}")
        aid="MITO_"+digest[:16]; ids.append(aid); seqs[aid]=sequence
        tax_rows.append({"ASV_ID":aid,"Domain":"Eukaryota","Phylum":"Mitochondria","Class":"Mammalia","Order":"Primates","Family":"Hominidae","Genus":"Homo","Species":"Homo sapiens","representative_reference_id":record.id})
        truth.append({"ASV_ID":aid,"source":source_label,"source_record":record.id,"description":record.description,"sequence_length":len(sequence),"sequence_sha256":digest})
    return ids,seqs,pd.DataFrame(tax_rows).set_index("ASV_ID"),pd.DataFrame(truth)


def add_mitochondrial_reads(counts, mitochondrial_ids, cfg, rng):
    a=cfg.get("reference_fixtures",{}).get("mitochondria",{}); out=counts.copy(); rows=[]
    prevalence=float(a.get("prevalence",0.5)); mean=float(a.get("mean_reads",50))
    for aid in mitochondrial_ids:
        values=[int(rng.poisson(mean)) if rng.random()<prevalence else 0 for _ in out.columns]
        out.loc[aid]=values; rows.append({"ASV_ID":aid,"configured_prevalence":prevalence,"configured_mean_reads":mean,"samples_present":sum(v>0 for v in values),"total_reads":sum(values)})
    return out,pd.DataFrame(rows)


def make_chimeras(counts,seqs,cfg,rng):
    a=cfg.get("artifacts",{}).get("chimeras",{}); out=counts.copy(); newseq={}; rows=[]
    if not a.get("enabled",False): return out,newseq,pd.DataFrame(columns=["chimera_ASV_ID","parent_1","parent_2","breakpoint"])
    frac=float(a.get("fraction_of_reads",0.02)); n=int(a.get("n_chimeras",5)); parents=list(counts.index)
    if len(parents)<2:return out,newseq,pd.DataFrame()
    for _ in range(n):
        p1,p2=rng.choice(parents,size=2,replace=False); m=min(len(seqs[p1]),len(seqs[p2])); bp=int(rng.integers(max(30,m//3),min(m-30,2*m//3)))
        s=seqs[p1][:bp]+seqs[p2][bp:]; cid="CHIM_"+hashlib.sha256(s.encode()).hexdigest()[:16]
        newseq[cid]=s; transfer=np.floor(out.loc[p1]*frac/2+out.loc[p2]*frac/2).astype(int)
        out.loc[p1]=np.maximum(0,out.loc[p1]-transfer//2); out.loc[p2]=np.maximum(0,out.loc[p2]-(transfer-transfer//2)); out.loc[cid]=transfer
        rows.append((cid,p1,p2,bp))
    return out,newseq,pd.DataFrame(rows,columns=["chimera_ASV_ID","parent_1","parent_2","breakpoint"])


def create_chemistry(rel,compounds,drivers_per,coef_sd,noise_sd,zero_inflation,log_transform,rng,sample_meta=None,batch_sd=0.0,min_abs_coefficient=0.0,driver_min_prevalence=0.0):
    x=np.log(rel.T.to_numpy(float)+1e-8); x-=x.mean(axis=1,keepdims=True); chem=np.zeros((x.shape[0],len(compounds))); truth=[]; batch_rows=[]
    prevalence=(rel>0).mean(axis=1).to_numpy(); eligible=np.flatnonzero(prevalence>=driver_min_prevalence)
    if len(eligible)<drivers_per: raise ValueError(f"Chemistry requires {drivers_per} ASV drivers with prevalence >= {driver_min_prevalence}; found {len(eligible)}")
    for j,compound in enumerate(compounds):
        k=min(drivers_per,len(eligible)); drivers=rng.choice(eligible,size=k,replace=False); raw=rng.normal(0,coef_sd,size=k); beta=np.sign(raw)*np.maximum(np.abs(raw),min_abs_coefficient)
        signal=x[:,drivers]@beta; signal=(signal-signal.mean())/(signal.std()+1e-12); y=signal+rng.normal(0,noise_sd,x.shape[0])
        if sample_meta is not None and batch_sd>0 and "batch" in sample_meta:
            for batch in sorted(sample_meta.batch.unique()):
                effect=float(rng.normal(0,batch_sd)); mask=sample_meta.batch.to_numpy()==batch; y[mask]+=effect; batch_rows.append((compound,batch,effect))
        if log_transform:y=np.exp(y)
        y[rng.random(x.shape[0])<zero_inflation]=0; chem[:,j]=y
        for idx,b in zip(drivers,beta):truth.append((compound,rel.index[idx],float(b),"positive" if b>0 else "negative"))
    return pd.DataFrame(chem,index=rel.columns,columns=compounds),pd.DataFrame(truth,columns=["compound","ASV_ID","coefficient","direction"]),pd.DataFrame(batch_rows,columns=["compound","batch","additive_effect"])


def count_fastq_records(path):
    with open_text(path) as h:n=sum(1 for _ in h)
    if n%4:raise RuntimeError(f"Malformed FASTQ: {path}")
    return n//4


def pad_for_iss_readthrough(molecule: str, cfg: dict) -> tuple[str, int]:
    """Pad short primer-bearing amplicons so ISS MiSeq does not skip them.

    The built-in InSilicoSeq MiSeq model can request reads longer than a V4
    molecule. Real sequencing would continue into adapter/read-through sequence;
    ISS instead skips short templates unless we make the template long enough.
    Returns the padded molecule and the number of synthetic read-through bases
    added.
    """
    minimum_template_length = int(cfg["fastq"].get("minimum_template_length", 350))
    if len(molecule) >= minimum_template_length:
        return molecule, 0

    padding = cfg["fastq"].get(
        "readthrough_padding_sequence",
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    )
    if not padding or set(padding.upper()) - set("ACGT"):
        raise ValueError("fastq.readthrough_padding_sequence must contain only A/C/G/T bases")

    required = minimum_template_length - len(molecule)
    molecule += (padding.upper() * ((required // len(padding)) + 1))[:required]
    return molecule, required


def write_iss_inputs(sample,sample_counts,seqs,cfg,workdir):
    positive=sample_counts[sample_counts.astype(int)>0].astype(int); fa=workdir/f"{sample}.amplicons.fasta"; rc=workdir/f"{sample}.readcounts.tsv"
    raw_lengths=[]; padded_lengths=[]; total_padding=0
    with fa.open("w") as fh,rc.open("w") as rh:
        for aid,n in positive.items():
            molecule=iss_amplicon_sequence(
                seqs[aid],
                bool(cfg["fastq"].get("include_primers",True)),
                cfg["reference"]["forward_primer"],
                cfg["reference"]["reverse_primer"],
                aid,
            )
            raw_lengths.append(len(molecule))
            molecule,added=pad_for_iss_readthrough(molecule,cfg)
            padded_lengths.append(len(molecule)); total_padding+=added
            SeqIO.write(SeqRecord(Seq(molecule),id=aid,description=""),fh,"fasta")

            # InSilicoSeq counts individual reads, not read pairs. For paired-end
            # output, request 2 * n reads so the resulting R1/R2 files contain
            # n records each. This preserves ASV counts as read-pair counts.
            rh.write(f"{aid}\t{int(n)*2}\n")

    expected_pairs=int(positive.sum())
    if raw_lengths:
        LOG.info(
            "%s amplicon lengths before padding: min=%d, median=%d, max=%d; "
            "after padding: min=%d, median=%d, max=%d; added=%d total bases",
            sample,
            min(raw_lengths), int(np.median(raw_lengths)), max(raw_lengths),
            min(padded_lengths), int(np.median(padded_lengths)), max(padded_lengths),
            total_padding,
        )
    else:
        LOG.warning("%s has no positive ASV counts", sample)
    return fa,rc,expected_pairs


def write_fastqs(counts,seqs,cfg,outdir):
    if not shutil.which("iss"):raise RuntimeError("InSilicoSeq executable `iss` is required")
    fq=outdir/"fastq"; inp=outdir/"iss_inputs";fq.mkdir(exist_ok=True);inp.mkdir(exist_ok=True)
    f=cfg["fastq"]; summary=[]
    for i,sample in enumerate(counts.columns):
        fa,rc,expected=write_iss_inputs(sample,counts[sample],seqs,cfg,inp);prefix=fq/sample
        cmd=["iss","generate","--genomes",str(fa),"--readcount_file",str(rc),"--sequence_type","amplicon","--model",str(f.get("model","miseq")),"--cpus",str(f.get("cpus",2)),"--seed",str(int(cfg["seed"])+i),"--output",str(prefix)]
        if f.get("gzip",True):cmd.append("--compress")
        run(cmd);suffix=".fastq.gz" if f.get("gzip",True) else ".fastq";r1=Path(f"{prefix}_R1{suffix}");r2=Path(f"{prefix}_R2{suffix}")
        n1=count_fastq_records(r1);n2=count_fastq_records(r2)
        if n1!=expected or n2!=expected:raise RuntimeError(f"Read-count mismatch {sample}: expected {expected} read pairs, R1={n1}, R2={n2}. Check ISS log for skipped templates.")
        summary.append((sample,expected,n1,n2,r1.name,r2.name))
    pd.DataFrame(summary,columns=["sample_id","expected_pairs","observed_R1","observed_R2","R1","R2"]).to_csv(outdir/"fastq_validation.tsv",sep="\t",index=False)
    with (outdir/"fastq_manifest.tsv").open("w") as manifest:
        manifest.write("sample_id\tfastq_r1\tfastq_r2\n")
        for sample,_,_,_,r1,r2 in summary: manifest.write(f"{sample}\tfastq/{r1}\tfastq/{r2}\n")
    if not f.get("keep_iss_inputs",True):shutil.rmtree(inp)


def write_reference_fixtures(outdir, seqs, contaminants, mitochondrial_ids):
    refs=outdir/"references";refs.mkdir(exist_ok=True);contaminant_ids=list(contaminants)
    if not contaminant_ids: raise ValueError("A complete benchmark requires at least one simulated contaminant")
    with (refs/"contaminants.fasta").open("w") as handle:
        SeqIO.write([SeqRecord(Seq(seqs[x]),id=x,description="simulated_contaminant") for x in contaminant_ids],handle,"fasta")
    with (refs/"mitochondria.fasta").open("w") as handle:
        SeqIO.write([SeqRecord(Seq(seqs[x]),id=x,description="genuine_mitochondrial_source") for x in mitochondrial_ids],handle,"fasta")
    pd.DataFrame([(x,"contaminant") for x in contaminant_ids]+[(x,"mitochondrial") for x in mitochondrial_ids],
                 columns=["ASV_ID","expected_filter"]).to_csv(outdir/"ground_truth_reference_filters.tsv",sep="\t",index=False)


def command_version(exe,args=("--version",)):
    if not shutil.which(exe):return None
    try:return run([exe,*args],capture=True).splitlines()[0]
    except Exception:return "installed (version unavailable)"


def write_report(outdir,manifest):
    counts=pd.read_csv(outdir/"asv_counts_final.tsv",sep="\t",index_col=0); chem=pd.read_csv(outdir/"chemistry.tsv",sep="\t",index_col=0); meta=pd.read_csv(outdir/"sample_metadata.tsv",sep="\t")
    rows=[("Samples",counts.shape[1]),("Final sequence features",counts.shape[0]),("Total read pairs",int(counts.to_numpy().sum())),("Chemical compounds",chem.shape[1]),("Median sample depth",int(counts.sum().median()))]
    table="".join(f"<tr><th>{html.escape(str(k))}</th><td>{html.escape(str(v))}</td></tr>" for k,v in rows)
    body=f"""<!doctype html><html><head><meta charset='utf-8'><title>DECOI report</title><style>body{{font-family:system-ui;max-width:950px;margin:2rem auto;padding:0 1rem}}table{{border-collapse:collapse}}th,td{{border:1px solid #bbb;padding:.45rem;text-align:left}}code{{background:#eee;padding:.1rem .3rem}}</style></head><body><h1>{html.escape(str(manifest['study_name']))}</h1><p>DECOI truth-aware V4 16S and chemistry benchmark.</p><table>{table}</table><h2>Artifact settings</h2><pre>{html.escape(yaml.safe_dump(manifest['artifacts'],sort_keys=False))}</pre><h2>Primary truth files</h2><p><code>asv_counts_biological.tsv</code>, <code>asv_counts_post_pcr.tsv</code>, <code>asv_counts_final.tsv</code>, <code>ground_truth_feature_registry.tsv</code>, <code>ground_truth_microbiome_batch_effects.tsv</code>, <code>ground_truth_mitochondria.tsv</code>, and <code>ground_truth_asv_chem.tsv</code>.</p></body></html>"""
    (outdir/"report.html").write_text(body)


def simulate(cfg,output):
    output.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(int(cfg["seed"]));meta=load_study_samples(cfg);cfg["simulation"]["n_samples"]=len(meta)
    tax,allseq=load_reference(cfg);chosen=choose_asvs(cfg,tax,allseq,rng);counts,rel=run_sparsedossa(cfg,output)
    mitochondrial_ids,mitochondrial_seqs,mitochondrial_tax,mitochondrial_source_truth=prepare_mitochondrial_features(cfg,output)
    allseq.update(mitochondrial_seqs)
    if not mitochondrial_tax.empty: tax=pd.concat([tax,mitochondrial_tax],axis=0)
    counts.index=chosen;rel.index=chosen;counts.columns=meta.sample_id;rel.columns=meta.sample_id;counts=counts.astype(int)
    biological_groups,group_truth=apply_group_effects(counts,meta,cfg,rng)
    biological_network,network_truth=apply_microbial_network(biological_groups,cfg,rng)
    biological,batch_truth=apply_microbiome_batch_effects(biological_network,meta,cfg,rng)
    post_pcr,pcr_truth=apply_pcr_bias(biological,cfg,rng)
    with_mito,mitochondrial_abundance_truth=add_mitochondrial_reads(post_pcr,mitochondrial_ids,cfg,rng)
    with_contam,contaminants=add_contaminants(with_mito,chosen+mitochondrial_ids,tax,allseq,cfg,rng);seqs={x:allseq[x] for x in with_contam.index}
    final,newseq,chim_truth=make_chimeras(with_contam,seqs,cfg,rng);seqs.update(newseq)
    final_rel=final.div(final.sum(axis=0),axis=1).fillna(0)
    biological.to_csv(output/"asv_counts_biological.tsv",sep="\t");post_pcr.to_csv(output/"asv_counts_post_pcr.tsv",sep="\t");final.to_csv(output/"asv_counts_final.tsv",sep="\t");final.to_csv(output/"asv_counts.tsv",sep="\t");final_rel.to_csv(output/"asv_relative_abundance.tsv",sep="\t")
    registry=[]
    for aid in final.index:
        kind="chimera" if aid in newseq else ("contaminant" if aid in contaminants else ("mitochondrial" if aid in mitochondrial_ids else "biological"))
        t=tax.loc[aid].to_dict() if aid in tax.index else {r:"" for r in RANKS}; rep=t.get("representative_reference_id","")
        registry.append({"feature_uuid":stable_uuid("mock16s-chem-feature",aid),"ASV_ID":aid,"feature_type":kind,"sequence_sha256":hashlib.sha256(seqs[aid].encode()).hexdigest(),"v4_sequence":seqs[aid],"silva_representative":rep,**{r:t.get(r,"") for r in RANKS},"sparsedossa_feature":f"feature_{chosen.index(aid)+1:04d}" if aid in chosen else ""})
    reg=pd.DataFrame(registry);reg.to_csv(output/"ground_truth_feature_registry.tsv",sep="\t",index=False)
    with (output/"asv_sequences.fasta").open("w") as h:SeqIO.write([SeqRecord(Seq(seqs[x]),id=x,description="") for x in final.index],h,"fasta")
    reg.set_index("ASV_ID")[[*RANKS,"silva_representative","feature_type"]].to_csv(output/"asv_taxonomy.tsv",sep="\t")
    cc=cfg["chemistry"];batch_sd=float(cfg.get("artifacts",{}).get("chemistry_batch_effect",{}).get("sd",0)) if cfg.get("artifacts",{}).get("chemistry_batch_effect",{}).get("enabled",False) else 0
    chemistry_rel=biological.div(biological.sum(axis=0),axis=1).fillna(0)
    chem,ct,bt=create_chemistry(chemistry_rel,list(cc["compounds"]),int(cc["drivers_per_compound"]),float(cc["coefficient_sd"]),float(cc["noise_sd"]),float(cc.get("zero_inflation",0)),bool(cc.get("log_transform",True)),rng,meta,batch_sd,float(cc.get("min_abs_coefficient",0)),float(cc.get("driver_min_prevalence",0)))
    chem.to_csv(output/"chemistry.tsv",sep="\t",index_label="sample_id");ct.to_csv(output/"ground_truth_asv_chem.tsv",sep="\t",index=False);bt.to_csv(output/"ground_truth_chemistry_batch.tsv",sep="\t",index=False)
    pcr_truth.to_csv(output/"ground_truth_pcr_bias.tsv",sep="\t",index=False);group_truth.to_csv(output/"ground_truth_group_effects.tsv",sep="\t",index=False);network_truth.to_csv(output/"ground_truth_network_modules.tsv",sep="\t",index=False);batch_truth.to_csv(output/"ground_truth_microbiome_batch_effects.tsv",sep="\t",index=False);chim_truth.to_csv(output/"ground_truth_chimeras.tsv",sep="\t",index=False)
    mitochondrial_source_truth.merge(mitochondrial_abundance_truth,on="ASV_ID",how="left").to_csv(output/"ground_truth_mitochondria.tsv",sep="\t",index=False)
    meta=meta.copy();meta["biological_read_depth"]=biological.sum(axis=0).to_numpy();meta["final_read_depth"]=final.sum(axis=0).to_numpy();meta.to_csv(output/"sample_metadata.tsv",sep="\t",index=False)
    write_reference_fixtures(output,seqs,contaminants,mitochondrial_ids)
    write_fastqs(final,seqs,cfg,output)
    study=yaml.safe_load(Path(cfg["study_design_file"]).read_text()) if cfg.get("study_design_file") else {"study_name":"mock_study"}
    manifest={"schema_version":"1.0","software":"DECOI","study_name":study.get("study_name","mock_study"),"seed":cfg["seed"],"config":cfg,"artifacts":cfg.get("artifacts",{}),"dimensions":{"samples":final.shape[1],"biological_asvs":len(chosen),"final_features":final.shape[0],"total_read_pairs":int(final.to_numpy().sum())},"versions":{"python":sys.version.split()[0],"platform":platform.platform(),"insilicoseq":command_version("iss"),"R":command_version("Rscript"),"cutadapt":command_version("cutadapt")},"reference_manifest":json.loads((Path(cfg["reference"]["output_dir"])/"reference_manifest.json").read_text())}
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n");write_report(output,manifest)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--config",required=True,type=Path);p.add_argument("--verbose",action="store_true");p.add_argument("--silva-fasta",type=Path);p.add_argument("--reference-dir",type=Path);p.add_argument("--study-design",type=Path)
    sub=p.add_subparsers(dest="command",required=True);sub.add_parser("prepare-reference");s=sub.add_parser("simulate");s.add_argument("--output",type=Path,default=Path("mock_dataset"));args=p.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,format="%(levelname)s: %(message)s");cfg=yaml.safe_load(args.config.read_text())
    if args.silva_fasta:cfg["reference"]["dada2_silva_fasta"]=str(args.silva_fasta)
    if args.reference_dir:cfg["reference"]["output_dir"]=str(args.reference_dir)
    if args.study_design:cfg["study_design_file"]=str(args.study_design)
    prepare_reference(cfg) if args.command=="prepare-reference" else simulate(cfg,args.output)

if __name__=="__main__":main()
