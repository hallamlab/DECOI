#!/usr/bin/env python3
"""Prepare the CAMI-derived airway/oral study with skin and extraction controls.

Repository study preparation only; use normal `decoi run` for simulation.
"""
import argparse
import json
import shutil
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.prepare_cami_study import prepare
from mock16s_chem import choose_asvs, load_reference

ROOT = Path(__file__).resolve().parents[1]


def merge_settings(target, overrides):
    for key, value in overrides.items():
        if isinstance(value, dict):
            if not isinstance(target.get(key), dict):
                target[key] = {}
            merge_settings(target[key], value)
        else:
            target[key] = value


def clinical_compositions(source, recipe, rng, ref):
    """Construct explicitly synthetic patients with shared latent mixture effects."""
    design = recipe['clinical_design']
    patient_sd, sample_sd = (float(design[k]) for k in ('patient_log_sd','sample_log_sd'))
    if any(not np.isfinite(v) or v < 0 for v in (patient_sd,sample_sd)):
        raise ValueError('Clinical mixture standard deviations must be finite and nonnegative')
    rows, columns, patient_effects = [], {}, {}
    for case, key, prefix in [('Control','control_patients','CTRL'),('Cancer','cancer_patients','CASE')]:
        n = int(design[key])
        if n < 1: raise ValueError('Clinical cohorts require at least one patient')
        for i in range(n):
            patient = f'{prefix}_{i+1:03d}'
            effects = rng.normal(0,patient_sd,len(source))
            patient_effects[patient] = effects
            for site,settings in recipe['sites'].items():
                originals = sorted(c for c in source if c.startswith(site+'_'))
                if not originals: raise ValueError(f'Missing site {site}')
                original = originals[i % len(originals)]
                statuses = (['TumorSide','Contralateral'] if case=='Cancer' else ['Healthy']) if site=='Airways' else [None]
                for status in statuses:
                    sid = f'{patient}_{site}' + (f'_{status}' if status else '')
                    weights = source[original].to_numpy() * np.exp(effects + rng.normal(0,sample_sd,len(source)))
                    columns[sid] = weights/weights.sum()
                    rows.append(dict(sample_id=sid,Sample=sid,Participant_ID=patient,Case=case,
                                     disease_status='case' if case=='Cancer' else 'control',lung_status=status,
                                     group=case,body_site=site,Type_Group=site,study_role=settings['role'],
                                     synthetic_clinical=True,is_biological_control=site=='Skin',
                                     include_primary_comparison=site in recipe['aspire']['comparison_groups'],
                                     source_sample=original,batch=f'plate_{i%2+1}',
                                     DNA_conc=float(rng.lognormal(3,0.4))))
    pd.DataFrame(patient_effects,index=source.index).rename_axis('ASV_ID').to_csv(ref/'patient_log_effects.tsv',sep='\t')
    return rows, columns


def configure(prepared, recipe):
    """Derive a controlled cohort from an unmodified prepared CAMI bundle."""
    prepared = Path(prepared).resolve()
    cfg = yaml.safe_load((prepared/'config.yaml').read_text())
    ref = prepared/'reference'
    source = pd.read_csv(ref/'amplicon_abundances.tsv', sep='\t', index_col=0)
    source.to_csv(ref/'source_amplicon_abundances.tsv', sep='\t')
    rng = np.random.default_rng(int(recipe['seed']))
    sd = float(recipe['composition_log_sd'])
    if not np.isfinite(sd) or sd < 0:
        raise ValueError('composition_log_sd must be nonnegative and finite')
    if recipe.get('clinical_design'):
        rows, columns = clinical_compositions(source, recipe, rng, ref)
    else:
        rows, columns = [], {}
        for site, settings in recipe['sites'].items():
            originals = sorted(c for c in source if c.startswith(site+'_'))
            n = int(settings['samples'])
            if not originals or n < 1:
                raise ValueError(f'Invalid sample count or missing source compositions: {site}')
            for i in range(n):
                original = originals[i % len(originals)]
                sid = f'{site}_{i+1:03d}'
                weights = source[original].to_numpy() * rng.lognormal(0, sd, len(source))
                columns[sid] = weights / weights.sum()
                rows.append(dict(sample_id=sid, Sample=sid, group=site, body_site=site,
                                 Type_Group=site, study_role=settings['role'],
                                 is_biological_control=settings['role']=='biological_control',
                                 include_primary_comparison=site in recipe['aspire']['comparison_groups'],
                                 source_sample=original, source_variant=i//len(originals)+1,
                                 batch=f'plate_{i%2+1}', DNA_conc=float(rng.lognormal(3, 0.4))))
    weights = pd.DataFrame(columns, index=source.index)
    weights = weights.loc[weights.sum(axis=1)>0]
    merge_settings(cfg, recipe.get('simulation_overrides', {}))
    # Generated configs resolve custom paths relative to their input bundle,
    # not the installed DECOI resource directory. Bundle the exact fixture.
    fixture = cfg.get('reference_fixtures', {}).get('mitochondria', {})
    if fixture.get('enabled') and fixture.get('source_fasta'):
        source_fixture = Path(fixture['source_fasta']).expanduser()
        if not source_fixture.is_absolute():
            source_fixture = ROOT / source_fixture
        bundled_fixture = prepared / 'fixtures' / source_fixture.name
        bundled_fixture.parent.mkdir(exist_ok=True)
        if source_fixture.resolve() != bundled_fixture.resolve():
            shutil.copyfile(source_fixture, bundled_fixture)
        fixture['source_fasta'] = str(bundled_fixture)
    cfg['seed'] = int(recipe['seed'])
    cfg['simulation'].update(n_samples=len(rows), n_asvs=len(weights),
                             median_read_depth=int(recipe['amplicon_read_pairs']))
    cfg['wgs']['enabled'] = False
    cfg['artifacts']['contaminants'].update(enabled=True, n_asvs=int(recipe['contaminant_asvs']))
    cfg['artifacts']['extraction_controls'].update(enabled=True, labels=recipe['extraction_blanks'])
    # Select contaminants from reference-derived alleles absent from this cohort.
    # These are synthetic spike-ins, not claims about real reagent contaminants.
    tax, seqs = load_reference(cfg)
    candidates = choose_asvs(cfg, tax, seqs, np.random.default_rng(int(recipe['seed'])),
                            n=int(recipe['contaminant_asvs']), exclude=list(weights.index), role=None)
    tax['reference_role'] = 'unused'
    tax.loc[weights.index, 'reference_role'] = 'biological'
    tax.loc[candidates, 'reference_role'] = 'contaminant'
    tax.to_csv(ref/'v4_asv_taxonomy.tsv', sep='\t')
    tax.loc[candidates].to_csv(ref/'synthetic_contaminants.tsv', sep='\t')
    weights.to_csv(ref/'amplicon_abundances.tsv', sep='\t')
    (prepared/'study.yaml').write_text(yaml.safe_dump(dict(study_name=recipe['study_name'], samples=rows),sort_keys=False))
    (prepared/'config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
    (prepared/'preparation.yaml').write_text(yaml.safe_dump(recipe,sort_keys=False))
    meta = dict(source='source_amplicon_abundances.tsv', seed=recipe['seed'],
                composition_log_sd=None if recipe.get('clinical_design') else sd,
                clinical_design=recipe.get('clinical_design'), biological_libraries=len(rows),
                extraction_blanks=len(recipe['extraction_blanks']),
                synthetic_contaminant_asvs=candidates,
                interpretation='Variants share CAMI source compositions; not independent human participants. Skin supplies the BIO control class for ASPIRE prevalence-based decontamination.')
    (ref/'study_derivation.json').write_text(json.dumps(meta,indent=2)+'\n')
    # Merge this fragment into the ASPIRE config; fill paths after simulation.
    aspire = {
        'core': {
            'table_filter': {'min_sample_reads': 5000, 'min_relative_abundance_pct': 0},
            'control_decontam': {
                'enabled': True, 'metadata': 'REPLACE_WITH_DATASET/sample_metadata.tsv',
                'metadata_sample_col': 'sample_id', 'class_col': 'Type_Group',
                'biological_labels': recipe['aspire']['comparison_groups'],
                'technical_labels': ['Control'], 'bio_control_labels': ['Skin'],
                'positive_labels': ['Positive'], 'min_biological_reads': 5000,
                'technical_enabled': True, 'bio_control_enabled': True,
                'technical_score_threshold': 0.1, 'bio_control_score_threshold': 0.1,
            },
        },
        'standard': {
            'mito': {
                'enabled': True,
                'mito_fasta': 'REPLACE_WITH_DATASET/references/mitochondria.fasta',
                'contaminant_fasta': 'REPLACE_WITH_DATASET/references/contaminants.fasta',
                'mito_db': None, 'biof_db': None,
            },
            'filter_counts': {
                'enabled': True, 'metadata': 'REPLACE_WITH_DATASET/sample_metadata.tsv',
                'sample_id_col': 'sample_id', 'group_col': 'Type_Group',
                'min_relative_abundance_pct': 0.1, 'min_prevalence_fraction': 0.05,
            },
            'metadata_plots': {
                'enabled': True, 'metadata': 'REPLACE_WITH_DATASET/sample_metadata.tsv',
                'sample_col': 'sample_id', 'type_col': 'Type_Group',
                'keep_types': recipe['aspire']['comparison_groups'],
                'input_table': 'filtered', 'group_order': recipe['aspire']['comparison_groups'],
                'run_micro': True, 'run_mito': False,
            },
        },
    }
    (prepared/'aspire_metadata.yaml').write_text(yaml.safe_dump(aspire,sort_keys=False))
    print(f'Ready: {prepared}/config.yaml; {len(rows)} biological libraries plus '
          f"{len(recipe['extraction_blanks'])} extraction blanks. No reads generated.", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-i','--input',type=Path,required=True,help='Local mock_2022 directory')
    p.add_argument('-o','--output',type=Path,required=True,help='New study-input directory')
    p.add_argument('--recipe',type=Path,default=ROOT/'study/cami_patient_showcase.yaml')
    p.add_argument('--threads',type=int,default=4)
    a=p.parse_args()
    if a.threads < 1: p.error('--threads must be positive')
    recipe=yaml.safe_load(a.recipe.read_text())
    prepare(a.input,a.output,a.threads)
    configure(a.output,recipe)


if __name__ == '__main__':
    main()
