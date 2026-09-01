#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(SparseDOSSA2)
})

opts <- list(
  make_option("--n-samples", type="integer", dest="n_samples"),
  make_option("--n-features", type="integer", dest="n_features"),
  make_option("--median-depth", type="integer", dest="median_depth"),
  make_option("--template", type="character", dest="template", default="Stool"),
  make_option("--seed", type="integer", dest="seed"),
  make_option("--counts-out", type="character", dest="counts_out"),
  make_option("--rel-out", type="character", dest="rel_out")
)
args <- parse_args(OptionParser(option_list=opts))
stopifnot(!is.null(args$counts_out), !is.null(args$rel_out))
set.seed(args$seed)
sim <- SparseDOSSA2(
  template=args$template,
  n_sample=args$n_samples,
  new_features=TRUE,
  n_feature=args$n_features,
  spike_metadata='none',
  median_read_depth=args$median_depth,
  verbose=TRUE
)
counts <- round(sim$simulated_data)
rel <- sim$simulated_matrices$rel
colnames(counts) <- sprintf('sample_%03d', seq_len(ncol(counts)))
colnames(rel) <- colnames(counts)
write.table(counts, args$counts_out, sep='\t', quote=FALSE, col.names=NA)
write.table(rel, args$rel_out, sep='\t', quote=FALSE, col.names=NA)
