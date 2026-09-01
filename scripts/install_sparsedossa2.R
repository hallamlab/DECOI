#!/usr/bin/env Rscript

options(
    repos = c(CRAN = "https://cloud.r-project.org"),
    timeout = 600
)

required <- c(
    "ks",
    "mvtnorm",
    "huge",
    "future.apply",
    "magrittr",
    "truncnorm",
    "igraph",
    "Rmpfr"
)

missing <- required[
    !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0L) {
    stop(
        "Missing preinstalled dependencies: ",
        paste(missing, collapse = ", "),
        "\nThese must come from the Conda environment."
    )
}

if (!requireNamespace("SparseDOSSA2", quietly = TRUE)) {
    remotes::install_github(
        "biobakery/SparseDOSSA2@26a998a",
        dependencies = FALSE,
        upgrade = "never",
        build_vignettes = FALSE,
        quiet = FALSE
    )
}

stopifnot(requireNamespace("SparseDOSSA2", quietly = TRUE))

message(
    "SparseDOSSA2 installed: ",
    as.character(utils::packageVersion("SparseDOSSA2"))
)
