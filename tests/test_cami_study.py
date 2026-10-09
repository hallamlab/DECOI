import json
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
import pytest
import mock16s_chem as engine
from scripts.prepare_cami_study import lineage_consensus


def test_shared_asv_taxonomy_stops_at_first_disagreement():
    assert lineage_consensus([['Bacteria','P','C','O','F','G','S'],['Bacteria','P','C','O','F','H','T']]) == ['Bacteria','P','C','O','F','','']


def test_supplied_body_site_composition_does_not_use_sparsedossa(tmp_path,monkeypatch):
    ref=tmp_path/'ref';ref.mkdir()
    (ref/'v4_asv_pool.fasta').write_text('>ASV_a\n'+'A'*250+'\n>ASV_b\n'+'C'*250+'\n')
    (ref/'v4_asv_taxonomy.tsv').write_text('ASV_ID\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\nASV_a\tBacteria\tP\tC\tO\tF\tG\tS\nASV_b\tBacteria\tP\tC\tO\tF\tH\tT\n')
    (ref/'reference_manifest.json').write_text(json.dumps({'source_type':'cami2_amplicon'}))
    (ref/'amplicon_abundances.tsv').write_text('ASV_ID\tAirways_4\tSkin_1\nASV_a\t1\t0\nASV_b\t0\t1\n')
    study=tmp_path/'study.yaml';study.write_text(yaml.safe_dump({'study_name':'test','samples':[{'sample_id':'Airways_4','body_site':'Airways'},{'sample_id':'Skin_1','body_site':'Skin'}]}))
    cfg=yaml.safe_load((Path(__file__).parents[1]/'config/cami_body_sites.yaml').read_text())
    cfg['reference']['output_dir']=str(ref);cfg['study_design_file']=str(study)
    cfg['simulation']['median_read_depth']=100
    cfg['chemistry']['drivers_per_compound']=1
    cfg['artifacts']['pcr_bias']['enabled']=False
    monkeypatch.setattr(engine,'run_sparsedossa',lambda *a:pytest.fail('Must preserve supplied compositions'))
    monkeypatch.setattr(engine,'write_fastqs',lambda *a:None)
    monkeypatch.setattr(engine,'command_version',lambda *a:'test')
    out=tmp_path/'output';engine.simulate(cfg,out)
    counts=pd.read_csv(out/'asv_counts_biological.tsv',sep='\t',index_col=0)
    assert counts.to_numpy().tolist()==[[100,0],[0,100]]
    assert not (out/'wgs').exists()
    assert len(pd.read_csv(out/'chemistry.tsv',sep='\t'))==2
    truth=pd.read_csv(out/'ground_truth_feature_registry.tsv',sep='\t')
    assert truth.sparsedossa_feature.isna().all()


def test_prepare_preserves_shared_alleles_and_copy_weighting(tmp_path,monkeypatch):
    from scripts import prepare_cami_study as cami
    from Bio.Seq import Seq
    class SerialPool:
        def __init__(self, **kw):pass
        def __enter__(self):return self
        def __exit__(self,*a):pass
        def map(self,fn,args):return map(fn,args)
    monkeypatch.setattr(cami,'ProcessPoolExecutor',SerialPool)
    f=engine.resolve_iupac(cami.FORWARD,'f')
    r=str(Seq(engine.resolve_iupac(cami.REVERSE,'r')).reverse_complement())
    locus=f+'A'*250+r
    for site in cami.SITES:
        short=tmp_path/'input/MetaGs'/f'CAMI_II_{site}'/'short_read';short.mkdir(parents=True)
        sag=tmp_path/'input/SAGs'/f'CAMI_II_{site}';(sag/'fasta').mkdir(parents=True)
        (short/'abundance1.tsv').write_text('g1\t1\ng2\t1\n')
        (sag/'genome_to_id.tsv').write_text('g1\t/old/g1.fa\ng2\t/old/g2.fa\n')
        (sag/'genome_taxa_info.tsv').write_text('_CAMI_genomeID\tTAXPATHSN\ng1\tBacteria|P|C|O|F|G|S\ng2\tBacteria|P|C|O|F|H|T\n')
        (sag/'fasta/g1.fa').write_text('>r\n'+locus+'N'*40+locus+'\n')
        (sag/'fasta/g2.fa').write_text('>r\n'+f+'C'*250+r+'\n')
    out=tmp_path/'prepared';cami.prepare(tmp_path/'input',out,1)
    weights=pd.read_csv(out/'reference/amplicon_abundances.tsv',sep='\t',index_col=0)
    assert weights.shape==(2,5)
    assert sorted(weights.iloc[:,0].tolist())==pytest.approx([1/3,2/3])
    assert not yaml.safe_load((out/'config.yaml').read_text())['wgs']['enabled']


def test_controlled_demo_exercises_amplicon_modules(tmp_path, monkeypatch):
    from scripts.prepare_cami_demo import configure
    ref=tmp_path/'reference';ref.mkdir()
    repo=Path(__file__).parents[1]
    cfg=yaml.safe_load((repo/'config/cami_body_sites.yaml').read_text())
    cfg['reference']['output_dir']=str(ref)
    cfg['study_design_file']=str(tmp_path/'study.yaml')
    (tmp_path/'config.yaml').write_text(yaml.safe_dump(cfg))
    ids=[f'ASV_{i}' for i in range(40)]
    rng=np.random.default_rng(42)
    (ref/'v4_asv_pool.fasta').write_text(''.join(f'>{aid}\n'+''.join(rng.choice(list('ACGT'),250))+'\n' for aid in ids))
    tax=pd.DataFrame([dict(ASV_ID=aid, Domain='Bacteria',Phylum='P',Class='C',Order='O',Family='F',Genus='G',Species='S') for aid in ids])
    tax.to_csv(ref/'v4_asv_taxonomy.tsv',sep='\t',index=False)
    (ref/'reference_manifest.json').write_text(json.dumps({'source_type':'cami2_amplicon'}))
    columns={f'{site}_{i}': np.array([1/36]*36+[0]*4) for site in ('Airways','Oral','Skin') for i in range(10)}
    pd.DataFrame(columns,index=ids).rename_axis('ASV_ID').to_csv(ref/'amplicon_abundances.tsv',sep='\t')
    recipe=yaml.safe_load((repo/'study/cami_airway_oral_skin.yaml').read_text())
    configure(tmp_path,recipe)
    cfg=yaml.safe_load((tmp_path/'config.yaml').read_text())
    from decoi.cli import read_config
    resolved=read_config(tmp_path/'config.yaml',repo)
    fixture=Path(resolved['reference_fixtures']['mitochondria']['source_fasta'])
    assert fixture.parent==tmp_path/'fixtures'
    assert fixture.read_bytes()==(repo/'fixtures'/fixture.name).read_bytes()
    study=engine.load_study_samples(cfg)
    assert study.groupby('body_site').size().to_dict()=={'Airways':20,'Oral':20,'Skin':20}
    assert study.loc[study.body_site=='Skin','is_biological_control'].all()
    assert not study.loc[study.body_site=='Skin','include_primary_comparison'].any()
    assert study.groupby('body_site').source_sample.nunique().eq(10).all()
    assert cfg['validation']['run_dada2'] and not cfg['wgs']['enabled']
    monkeypatch.setattr(engine,'write_fastqs',lambda *a:None)
    monkeypatch.setattr(engine,'command_version',lambda *a:'test')
    monkeypatch.setattr(engine,'run_sparsedossa',lambda *a:pytest.fail('Source communities must be retained'))
    out=tmp_path/'output';engine.simulate(cfg,out)
    meta=pd.read_csv(out/'sample_metadata.tsv',sep='\t')
    assert len(meta)==64
    assert meta.is_negative_control.sum()==4
    assert meta.is_biological_control.sum()==20
    assert meta.include_primary_comparison.sum()==40
    assert not meta.loc[meta.body_site=='Skin','is_negative_control'].any()
    assert len(pd.read_csv(out/'chemistry.tsv',sep='\t'))==60
    for suffix in ('group_effects','network_modules','microbiome_batch_effects','pcr_bias','chimeras','mitochondria','extraction_controls','asv_chem','chemistry_batch'):
        assert not pd.read_csv(out/f'ground_truth_{suffix}.tsv',sep='\t').empty, suffix
    registry=pd.read_csv(out/'ground_truth_feature_registry.tsv',sep='\t')
    assert set(registry.feature_type)=={'biological','mitochondrial','contaminant','chimera'}
    assert not (out/'wgs').exists()
    # The same study seed reproduces both blank depths and biological artifacts.
    second=tmp_path/'repeat';engine.simulate(cfg,second)
    assert (out/'asv_counts_final.tsv').read_bytes()==(second/'asv_counts_final.tsv').read_bytes()
    blanks=pd.read_csv(out/'ground_truth_extraction_controls.tsv',sep='\t').set_index('sample_id')
    assert (blanks.loc[['PBS','PBS_twz'],'total_reads']<5000).all()
    assert (blanks.loc[['Negative_96','Negative_man'],'total_reads']>5000).all()
    # Databases contain the exact injected molecules, not just matching names.
    from Bio import SeqIO
    for filename,kind in [('contaminants.fasta','contaminant'),('mitochondria.fasta','mitochondrial')]:
        actual={r.id:str(r.seq) for r in SeqIO.parse(out/'references'/filename,'fasta')}
        expected=registry.loc[registry.feature_type==kind].set_index('ASV_ID').v4_sequence.to_dict()
        assert actual==expected and actual
    aspire=yaml.safe_load((tmp_path/'aspire_metadata.yaml').read_text())
    assert aspire['core']['table_filter']['min_sample_reads']==5000
    assert aspire['core']['control_decontam']['bio_control_labels']==['Skin']
    assert aspire['core']['control_decontam']['technical_enabled']
    assert aspire['core']['control_decontam']['bio_control_enabled']
    assert aspire['standard']['filter_counts']['min_prevalence_fraction']==0.05
    assert aspire['standard']['metadata_plots']['input_table']=='filtered'
    assert 'subtraction_groups' not in aspire['standard']['metadata_plots']
    assert aspire['standard']['mito']['contaminant_fasta'].endswith('references/contaminants.fasta')
    assert aspire['standard']['mito']['mito_fasta'].endswith('references/mitochondria.fasta')


def test_blank_mean_validation_and_scalar_compatibility():
    counts=pd.DataFrame({'sample':[100,10]},index=['biological','contaminant'])
    meta=pd.DataFrame({'sample_id':['sample']})
    settings={'enabled':True,'labels':['blank'],'contaminant_mean_reads':600}
    cfg={'artifacts':{'extraction_controls':settings}}
    first=engine.add_extraction_controls(counts,meta,['contaminant'],cfg,np.random.default_rng(42))[0]
    settings['contaminant_mean_reads']={'blank':600}
    second=engine.add_extraction_controls(counts,meta,['contaminant'],cfg,np.random.default_rng(42))[0]
    pd.testing.assert_frame_equal(first,second)
    for invalid in ({'other':600},{'blank':-1},{'blank':float('nan')}):
        settings['contaminant_mean_reads']=invalid
        with pytest.raises(ValueError):
            engine.add_extraction_controls(counts,meta,['contaminant'],cfg,np.random.default_rng(42))


def test_clinical_cohort_has_paired_airways_and_reproducible_patient_effects(tmp_path):
    from scripts.prepare_cami_demo import clinical_compositions
    repo=Path(__file__).parents[1]
    recipe=yaml.safe_load((repo/'study/cami_patient_showcase.yaml').read_text())
    source=pd.DataFrame({f'{site}_{i}':[0.2,0.8,0] for site in ('Airways','Oral','Skin') for i in range(10)},index=['a','b','c'])
    rows,weights=clinical_compositions(source,recipe,np.random.default_rng(42),tmp_path)
    repeated,other=clinical_compositions(source,recipe,np.random.default_rng(42),tmp_path)
    assert rows==repeated
    pd.testing.assert_frame_equal(pd.DataFrame(weights),pd.DataFrame(other))
    meta=pd.DataFrame(rows)
    assert len(meta)==175 and meta.Participant_ID.nunique()==50
    assert meta.groupby('Case').Participant_ID.nunique().to_dict()=={'Cancer':25,'Control':25}
    assert meta.groupby('body_site').size().to_dict()=={'Airways':75,'Oral':50,'Skin':50}
    for _,patient in meta.groupby('Participant_ID'):
        assert patient.batch.nunique()==1
        assert {'Airways','Oral','Skin'}==set(patient.body_site)
        airway=patient[patient.body_site=='Airways']
        assert set(airway.lung_status)==({'TumorSide','Contralateral'} if patient.Case.iloc[0]=='Cancer' else {'Healthy'})
        assert airway.source_sample.nunique()==1
    assert meta.loc[meta.body_site!='Airways','lung_status'].isna().all()
    assert pd.DataFrame(weights,index=source.index).loc['c'].eq(0).all()
    assert pd.read_csv(tmp_path/'patient_log_effects.tsv',sep='\t').shape==(3,51)


def test_lung_effect_does_not_create_non_lung_group():
    counts=pd.DataFrame({'airway':[100,50],'oral':[100,50]},index=['a','b'])
    meta=pd.DataFrame({'sample_id':['airway','oral'],'lung_status':['Healthy',None]})
    cfg={'artifacts':{'group_differential_abundance':{'enabled':True,'group_columns':['lung_status'],'asvs_per_group':1}}}
    _,truth=engine.apply_group_effects(counts,meta,cfg,np.random.default_rng(42))
    assert truth.group.tolist()==['Healthy']
