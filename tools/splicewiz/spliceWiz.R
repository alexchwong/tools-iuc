#!/usr/bin/env Rscript

# A command-line interface to SpliceWiz for use with Galaxy

# setup R error handling to go to stderr
options(show.error.messages = FALSE, error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, FALSE)
})

library("getopt")
library("tools")
library()
options(stringAsFactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "help", "h", 0, "logical",
  "cores", "s", 0, "integer",
  "mode", "m", 1, "character",
  "fasta", "G", 1, "character",
  "gtf", "T", 1, "character",
  "genome_type", "g", 1, "character",
  "ref", "r", 1, "character",
  "pbOutput", "P", 1, "character",
  "pbCOV", "C", 1, "character",
  "nxtSE", "N", 1, "character"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

suppressPackageStartupMessages({
    library("SpliceWiz")
})

if(opt$mode == "buildRef") {
    refPath <- file.path(tempdir(), "reference")
    if(dir.exists(refPath)) unlink(refPath, recursive = TRUE)
    species = ""
    if(opt$genome_type %in% c("hg38", "hg19")) {
        species = "Homo sapiens"
    } else if(opt$genome_type %in% c("mm9", "mm10")) {
        species = "Mus musculus"
    }
    buildRef(
        reference_path = refPath,
        fasta = opt$fasta,
        gtf = opt$gtf,
        genome_type = opt$genome_type,
        ontologySpecies = species,
        n_threads = opt$cores
    )
    
    # zip output
    zip(zipfile = opt$ref, files = list.files(
        path = refPath,
        full.names = TRUE, recursive = TRUE
    ))
    
    unlink(refPath, recursive = TRUE)
}