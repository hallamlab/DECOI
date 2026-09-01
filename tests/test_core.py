from pathlib import Path
import importlib.util
import hashlib
import numpy as np
import pandas as pd
spec=importlib.util.spec_from_file_location('m',Path(__file__).parents[1]/'mock16s_chem.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_ids():
 assert m.stable_asv_id('ACGT')==m.stable_asv_id('ACGT')
 assert m.stable_uuid('x','a')==m.stable_uuid('x','a')

def test_primers():
 assert m.iss_amplicon_sequence('ACGT',True,'AAAA','CCCC','x')=='AAAAACGTGGGG'
 assert set(m.resolve_iupac('GTGYCAGCMGCCGCGGTAA','x'))<=set('ACGT')

def test_study(tmp_path):
 p=tmp_path/'s.yaml';p.write_text('groups:\n  - name: case\n    sample_prefix: C\n    n_samples: 2\n    metadata:\n      batch: b1\n')
 d=m.load_study_samples({'study_design_file':str(p),'simulation':{'n_samples':1}});assert d.sample_id.tolist()==['C_001','C_002']

def test_cohort_study(tmp_path):
 p=tmp_path/'s.yaml';p.write_text('cohorts:\n  - name: case\n    participant_prefix: P\n    n_participants: 2\n    metadata:\n      Case: Cancer\n    sample_types:\n      - name: Brush\n        code: B\n      - name: BAL\n')
 d=m.load_study_samples({'study_design_file':str(p),'simulation':{'n_samples':1}})
 assert len(d)==4 and d.Participant_ID.nunique()==2
 assert set(d.Type_Group)=={'Brush','BAL'}

def test_participant_metadata_cycle(tmp_path):
 p=tmp_path/'s.yaml';p.write_text('cohorts:\n  - name: case\n    n_participants: 4\n    metadata:\n      Case: Cancer\n    participant_metadata_cycle:\n      - batch: plate_1\n      - batch: plate_2\n    sample_types:\n      - name: Brush\n      - name: BAL\n')
 d=m.load_study_samples({'study_design_file':str(p),'simulation':{'n_samples':1}})
 assert d.groupby('Participant_ID').batch.nunique().eq(1).all()
 assert d.groupby('batch').Participant_ID.nunique().to_dict()=={'plate_1':2,'plate_2':2}

def test_example_study_has_balanced_primary_case_control_comparisons():
 cfg={'study_design_file':str(Path(__file__).parents[1]/'study'/'example_study.yaml'),'simulation':{'n_samples':1}}
 d=m.load_study_samples(cfg)
 primary=d[(d.Case=='Control') | (d.lung_status=='TumorSide')]
 patient_counts=primary.groupby(['Type_Group','Case']).Participant_ID.nunique()
 assert patient_counts.to_dict()=={('BAL','Cancer'):10,('BAL','Control'):10,('Bronchial Brush','Cancer'):10,('Bronchial Brush','Control'):10}
 assert d.groupby('Participant_ID').batch.nunique().eq(1).all()
 assert set(d.groupby('batch').Case.unique().apply(tuple))=={('Control','Cancer')}

def test_renormalize_preserves_depth():
 x=pd.DataFrame({'s1':[4,6],'s2':[2,8]},index=['a','b']); y=m.renormalize_counts(x,x.sum(),np.random.default_rng(1));assert y.sum().tolist()==[10,10]

def test_group_effects_can_use_strong_disjoint_prevalent_drivers():
 counts=pd.DataFrame(np.arange(1,41).reshape(10,4),index=[f'a{i}' for i in range(10)],columns=['s1','s2','s3','s4'])
 meta=pd.DataFrame({'sample_id':counts.columns,'Type':['A','A','B','B'],'Case':['X','Y','X','Y']})
 cfg={'artifacts':{'group_differential_abundance':{'enabled':True,'group_columns':['Type','Case'],'asvs_per_group':2,'log_fold_change_sd':0,'min_abs_log_fold_change':1.5,'candidate_min_prevalence':1,'disjoint_drivers':True}}}
 out,truth=m.apply_group_effects(counts,meta,cfg,np.random.default_rng(3))
 assert len(truth)==8 and truth.ASV_ID.nunique()==8 and (truth.log_fold_change==1.5).all()
 assert (truth.baseline_prevalence==1).all() and out.sum().equals(counts.sum())

def test_pcr_bias_disabled():
 x=pd.DataFrame({'s':[2,3]},index=['a','b']);y,t=m.apply_pcr_bias(x,{'artifacts':{'pcr_bias':{'enabled':False}}},np.random.default_rng(1));pd.testing.assert_frame_equal(x,y);assert (t.pcr_efficiency==1).all()

def test_microbiome_batch_effects_are_reproducible_and_preserve_depth():
 counts=pd.DataFrame({'s1':[100,200,300,400],'s2':[100,200,300,400],'s3':[100,200,300,400],'s4':[100,200,300,400]},index=['a','b','c','d'])
 meta=pd.DataFrame({'sample_id':counts.columns,'batch':['p1','p1','p2','p2']})
 cfg={'artifacts':{'microbiome_batch_effect':{'enabled':True,'column':'batch','asvs_per_batch':2,'log_fold_change':1.25,'direction':'mixed'}}}
 out1,truth1=m.apply_microbiome_batch_effects(counts,meta,cfg,np.random.default_rng(7));out2,truth2=m.apply_microbiome_batch_effects(counts,meta,cfg,np.random.default_rng(7))
 pd.testing.assert_frame_equal(out1,out2);pd.testing.assert_frame_equal(truth1,truth2)
 assert out1.sum().equals(counts.sum()) and len(truth1)==4
 assert set(truth1.direction)<= {'higher','lower'} and set(truth1.batch)=={'p1','p2'}

def test_microbial_network_is_reproducible_and_preserves_depth():
 counts=pd.DataFrame(np.arange(1,97).reshape(24,4),index=[f'a{i}' for i in range(24)],columns=[f's{i}' for i in range(4)])
 cfg={'artifacts':{'microbial_network':{'enabled':True,'n_modules':3,'asvs_per_module':4,'loading':1.0,'candidate_min_prevalence':1.0}}}
 out1,truth1=m.apply_microbial_network(counts,cfg,np.random.default_rng(9));out2,truth2=m.apply_microbial_network(counts,cfg,np.random.default_rng(9))
 pd.testing.assert_frame_equal(out1,out2);pd.testing.assert_frame_equal(truth1,truth2)
 assert out1.sum().equals(counts.sum()) and len(truth1)==12
 assert truth1.ASV_ID.nunique()==12 and truth1.groupby('module').size().eq(4).all()

def test_contaminants():
 tax=pd.DataFrame({'Domain':['Bacteria']*3,'Phylum':['p']*3,'Class':['c']*3,'Order':['o']*3,'Family':['f']*3,'Genus':['g']*3,'Species':['s']*3},index=['a','b','c']);seq={'a':'A','b':'C','c':'G'};counts=pd.DataFrame({'s':[10]},index=['a'])
 cfg={'simulation':{'n_asvs':1,'selection':{'domains':['Bacteria'],'min_taxonomy_ranks':5}},'artifacts':{'contaminants':{'enabled':True,'n_asvs':1,'prevalence':1,'mean_reads':5}}}
 out,ids=m.add_contaminants(counts,['a'],tax,seq,cfg,np.random.default_rng(2));assert len(ids)==1 and len(out)==2

def test_mitochondrial_fixture(tmp_path):
 sequence='ACGTACGT'; source=tmp_path/'mitochondria.fasta';source.write_text('>NC_TEST Homo sapiens mitochondrion\n'+sequence+'\n')
 cfg={'reference_fixtures':{'mitochondria':{'enabled':True,'source_fasta':str(source),'expected_sequence_sha256':hashlib.sha256(sequence.encode()).hexdigest(),'prevalence':1,'mean_reads':5}}}
 ids,seqs,tax,truth=m.prepare_mitochondrial_features(cfg,tmp_path)
 counts=pd.DataFrame({'s1':[10],'s2':[20]},index=['ASV1']);out,abundance=m.add_mitochondrial_reads(counts,ids,cfg,np.random.default_rng(2))
 assert len(ids)==1 and seqs[ids[0]]==sequence and tax.loc[ids[0],'Phylum']=='Mitochondria'
 assert (out.loc[ids[0]]>0).all() and abundance.total_reads.iloc[0]>0 and truth.sequence_sha256.iloc[0]==hashlib.sha256(sequence.encode()).hexdigest()

def test_chimera_creation():
 seq={'a':'A'*150+'C'*100,'b':'G'*150+'T'*100};counts=pd.DataFrame({'s':[100,100]},index=['a','b']);cfg={'artifacts':{'chimeras':{'enabled':True,'n_chimeras':1,'fraction_of_reads':.1}}}
 out,new,truth=m.make_chimeras(counts,seq,cfg,np.random.default_rng(2));assert len(new)==1 and len(truth)==1 and out.sum().iloc[0]==200

def test_chemistry_batch():
 rel=pd.DataFrame([[.7,.2],[.3,.8]],index=['a','b'],columns=['s1','s2']);meta=pd.DataFrame({'sample_id':['s1','s2'],'batch':['x','y']})
 chem,truth,batch=m.create_chemistry(rel,['c'],1,.5,.1,0,True,np.random.default_rng(1),meta,.2);assert chem.shape==(2,1) and len(truth)==1 and len(batch)==2

def test_chemistry_driver_constraints():
 rel=pd.DataFrame([[.6,.5,.4,.5],[.4,.5,.6,.5],[0,0,.01,0]],index=['a','b','rare'],columns=['s1','s2','s3','s4'])
 chem,truth,_=m.create_chemistry(rel,['c'],2,.2,.1,0,True,np.random.default_rng(4),min_abs_coefficient=.9,driver_min_prevalence=.5)
 assert chem.shape==(4,1) and set(truth.ASV_ID)=={'a','b'} and (truth.coefficient.abs()>=.9).all()
