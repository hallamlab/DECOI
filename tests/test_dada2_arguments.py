"""Exercise the actual R parser without running read processing or DADA2."""
from pathlib import Path
import shutil
import subprocess

import pytest


def test_dada2_hyphenated_options_match_fields_used_by_script():
    rscript = shutil.which('Rscript')
    if not rscript:
        pytest.skip('Requires the supporting R environment')
    script = Path(__file__).parents[1]/'scripts/validate_dada2.R'
    code = r'''
    suppressPackageStartupMessages(library(optparse))
    exprs <- parse(commandArgs(trailingOnly=TRUE)[1])
    parser_expr <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                          identical(x[[2]], as.name("o")), as.list(exprs))
    stopifnot(length(parser_expr)==1L)
    commandArgs <- function(trailingOnly=FALSE) c(
      "--fastq-dir", "reads with spaces", "--outdir", "validation",
      "--forward-primer", "ACGT", "--reverse-primer", "TGCA",
      "--trunc-f", "123", "--trunc-r", "145", "--max-ee-f", "1.5",
      "--max-ee-r", "2.5", "--threads", "3")
    # parse_args evaluates its default args in the optparse namespace.
    parser_expr[[1]][[3]]$args <- commandArgs()
    eval(parser_expr[[1]])
    stopifnot(o$fastq_dir=="reads with spaces", o$outdir=="validation",
              o$forward_primer=="ACGT", o$reverse_primer=="TGCA",
              o$trunc_f==123L, o$trunc_r==145L, o$max_ee_f==1.5,
              o$max_ee_r==2.5, o$threads==3L)
    '''
    subprocess.run([rscript,'-e',code,str(script)],check=True,capture_output=True,text=True)


def test_dada2_tracking_table_counts_denoised_merged_and_empty_results(tmp_path):
    rscript = shutil.which('Rscript')
    if not rscript:
        pytest.skip('Requires the supporting R/DADA2 environment')
    script = Path(__file__).parents[1]/'scripts/validate_dada2.R'
    code = r'''
    suppressPackageStartupMessages(library(dada2))
    args <- commandArgs(trailingOnly=TRUE)
    exprs <- as.list(parse(args[1]))
    is_assignment <- function(x, name) is.call(x) && identical(x[[1]],as.name("<-")) && identical(x[[2]],as.name(name))
    # R's parser also preserves '=' assignments used in the script.
    is_assignment <- function(x, name) is.call(x) && as.character(x[[1]]) %in% c("<-","=") && identical(x[[2]],as.name(name))
    helper <- Filter(function(x) is_assignment(x,"getN"),exprs)
    stopifnot(length(helper)==1L)
    eval(helper[[1]])
    o <- list(outdir=args[2]); samples <- c("sample_A","sample_B")
    filt <- matrix(c(12L,10L,8L,0L),nrow=2,byrow=TRUE)
    dF <- list(structure(list(denoised=c(ACGT=6L,TGCA=3L)),class="dada"),
               structure(list(denoised=setNames(integer(0),character(0))),class="dada"))
    dR <- list(structure(list(denoised=c(ACGT=7L,TGCA=2L)),class="dada"), dF[[2]])
    mergers <- list(data.frame(sequence=c("ACGT","TGCA"),abundance=c(5L,2L)),
                    data.frame(sequence=character(0),abundance=integer(0)))
    seqtab.nochim <- matrix(c(5L,0L),ncol=1,dimnames=list(samples,"ACGT"))
    start <- which(vapply(exprs,function(x) is_assignment(x,"track"),logical(1)))
    stopifnot(length(start)==1L)
    for (i in seq.int(start,length(exprs))) eval(exprs[[i]])
    saved <- read.delim(file.path(o$outdir,"dada2_tracking.tsv"),row.names=1,check.names=FALSE)
    stopifnot(identical(rownames(saved),samples),
              identical(colnames(saved),c("input","filtered","denoisedF","denoisedR","merged","nonchim")),
              all(as.numeric(saved[1,])==c(12,10,9,9,7,5)),
              all(as.numeric(saved[2,])==c(8,0,0,0,0,0)))
    '''
    subprocess.run([rscript,'-e',code,str(script),str(tmp_path)],check=True,capture_output=True,text=True)
