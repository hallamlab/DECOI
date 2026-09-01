#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dada2); library(optparse)})
o <- parse_args(OptionParser(option_list=list(
 make_option('--fastq-dir'), make_option('--outdir'), make_option('--forward-primer'), make_option('--reverse-primer'),
 make_option('--trunc-f',type='integer',default=240), make_option('--trunc-r',type='integer',default=200),
 make_option('--max-ee-f',type='double',default=2), make_option('--max-ee-r',type='double',default=2), make_option('--threads',type='integer',default=4))))
dir.create(o$outdir,recursive=TRUE,showWarnings=FALSE); trimdir=file.path(o$outdir,'primer_trimmed'); filtdir=file.path(o$outdir,'filtered');dir.create(trimdir);dir.create(filtdir)
fnFs=sort(list.files(o$fastq_dir,pattern='_R1.fastq.gz$',full.names=TRUE));fnRs=sort(list.files(o$fastq_dir,pattern='_R2.fastq.gz$',full.names=TRUE));stopifnot(length(fnFs)==length(fnRs),length(fnFs)>0)
samples=sub('_R1.fastq.gz$','',basename(fnFs)); revcomp=function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
tFs=file.path(trimdir,paste0(samples,'_R1.fastq.gz'));tRs=file.path(trimdir,paste0(samples,'_R2.fastq.gz'))
for(i in seq_along(fnFs)) system2('cutadapt',c('-g',o$forward_primer,'-G',o$reverse_primer,'--discard-untrimmed','-o',tFs[i],'-p',tRs[i],fnFs[i],fnRs[i]))
fFs=file.path(filtdir,paste0(samples,'_F_filt.fastq.gz'));fRs=file.path(filtdir,paste0(samples,'_R_filt.fastq.gz'))
filt=filterAndTrim(tFs,fFs,tRs,fRs,truncLen=c(o$trunc_f,o$trunc_r),maxN=0,maxEE=c(o$max_ee_f,o$max_ee_r),truncQ=2,rm.phix=TRUE,compress=TRUE,multithread=o$threads)
errF=learnErrors(fFs,multithread=o$threads);errR=learnErrors(fRs,multithread=o$threads);dF=dada(fFs,err=errF,multithread=o$threads);dR=dada(fRs,err=errR,multithread=o$threads)
mergers=mergePairs(dF,fFs,dR,fRs,verbose=TRUE);seqtab=makeSequenceTable(mergers);seqtab.nochim=removeBimeraDenovo(seqtab,method='consensus',multithread=o$threads,verbose=TRUE)
write.table(seqtab.nochim,file.path(o$outdir,'dada2_sequence_table.tsv'),sep='\t',quote=FALSE,col.names=NA)
track=cbind(filt,sapply(dF,getN),sapply(dR,getN),sapply(mergers,getN),rowSums(seqtab.nochim));colnames(track)=c('input','filtered','denoisedF','denoisedR','merged','nonchim');rownames(track)=samples;write.table(track,file.path(o$outdir,'dada2_tracking.tsv'),sep='\t',quote=FALSE,col.names=NA)
