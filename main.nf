nextflow.enable.dsl=2
params.config = params.config ?: "${projectDir}/config/example.yaml"
params.study = params.study ?: "${projectDir}/study/example_study.yaml"
params.outdir = params.outdir ?: "results"
params.silva_url = params.silva_url ?: "https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toSpecies_trainset.fa.gz?download=1"
params.silva_md5 = params.silva_md5 ?: "e0455a1a7de820039d684aad2053e937"
params.silva_fasta = params.silva_fasta ?: null
params.run_dada2 = params.run_dada2 ?: false
params.dada2_threads = params.dada2_threads ?: 8

process DOWNLOAD_SILVA {
  tag "SILVA 138.2"; publishDir "${params.outdir}/reference/source", mode:'copy'; conda "/home/ryan/mambaforge-pypy3/envs/mock16s-chem-v1"
  output: path 'silva_nr99_v138.2_toSpecies_trainset.fa.gz', emit:fasta
  script:
  """
  curl --fail --location --retry 5 '${params.silva_url}' -o silva_nr99_v138.2_toSpecies_trainset.fa.gz
  echo '${params.silva_md5}  silva_nr99_v138.2_toSpecies_trainset.fa.gz' | md5sum --check -
  """
}
process PREPARE_REFERENCE {
  tag "extract SILVA V4"; publishDir "${params.outdir}/reference/v4", mode:'copy'; conda "/home/ryan/mambaforge-pypy3/envs/mock16s-chem-v1"
  input: path silva; path config
  output: path 'reference_v4', emit:reference
  script:
  """
  python '${projectDir}/mock16s_chem.py' --config '${config}' --silva-fasta '${silva}' --reference-dir reference_v4 prepare-reference
  """
}
process SIMULATE_STUDY {
  tag "simulate mock study"; publishDir "${params.outdir}/dataset", mode:'copy'; conda "/home/ryan/mambaforge-pypy3/envs/mock16s-chem-v1"
  input: path reference; path config; path study
  output: path 'mock_dataset', emit:dataset
  script:
  """
  python '${projectDir}/mock16s_chem.py' --config '${config}' --reference-dir '${reference}' --study-design '${study}' simulate --output mock_dataset
  """
}
process VALIDATE_DADA2 {
  tag "DADA2 validation"; publishDir "${params.outdir}/validation", mode:'copy'; conda "/home/ryan/mambaforge-pypy3/envs/mock16s-chem-v1"; cpus params.dada2_threads
  input: path dataset
  output: path 'dada2_validation', emit:validation
  script:
  """
  Rscript '${projectDir}/scripts/validate_dada2.R' \
    --fastq-dir '${dataset}/fastq' --outdir dada2_validation \
    --forward-primer GTGYCAGCMGCCGCGGTAA --reverse-primer GGACTACNVGGGTWTCTAAT \
    --threads ${task.cpus}
  """
}
workflow {
  config_ch=Channel.value(file(params.config,checkIfExists:true)); study_ch=Channel.value(file(params.study,checkIfExists:true))
  if(params.silva_fasta){silva_ch=Channel.value(file(params.silva_fasta,checkIfExists:true))} else {DOWNLOAD_SILVA();silva_ch=DOWNLOAD_SILVA.out.fasta}
  PREPARE_REFERENCE(silva_ch,config_ch); SIMULATE_STUDY(PREPARE_REFERENCE.out.reference,config_ch,study_ch)
  if(params.run_dada2){VALIDATE_DADA2(SIMULATE_STUDY.out.dataset)}
}
