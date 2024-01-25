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
  "cores", "c", 0, "integer",
  "mode", "m", 1, "character",
  "fasta", "G", 1, "character",
  "gtf", "T", 1, "character",
  "genome_type", "g", 1, "character",
  "ref", "r", 1, "character",
  "inbam", "b", 1, "character",
  "pbOutput", "P", 1, "character",
  "pbCOV", "C", 1, "character",
  "nxtNames", "n", 1, "character",
  "nxtNovelSplicing", "s", 1, "character",
  "nxtPackageCOV", "p", 1, "character",
  "nxtSE", "N", 1, "character"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

if ((zip <- Sys.getenv("R_ZIPCMD", "zip")) == "") zip <- "zip"

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
        ontologySpecies = species
    )
    
    # tar output
    tmptar <- paste0(tempfile(), ".tar")
    
    setwd(refPath) # to archive all files under current directory
    filestodo <- list.files(".", recursive = TRUE, full.names = TRUE)
    utils::tar(tarfile = tmptar, files = filestodo, compression = "none")
    setwd(tempdir())
    
    file.copy(from = tmptar, to = opt$ref, overwrite = TRUE)
    
    file.remove(tmptar)
    unlink(refPath, recursive = TRUE)
} else if(opt$mode == "processBAM") {
    
    # extract input ref
    refPath <- file.path(tempdir(), "reference")
    if(dir.exists(refPath)) unlink(refPath, recursive = TRUE)
    dir.create(refPath)
    untar(opt$ref, exdir = refPath)

    outPath <- file.path(tempdir(), "pb_outs")
    if(dir.exists(outPath)) unlink(outPath, recursive = TRUE)
    dir.create(outPath)
    
    bams <- unlist(strsplit(opt$inbam, split = ",", fixed = TRUE))
    pbouts <- unlist(strsplit(opt$pbOutput, split = ",", fixed = TRUE))
    pbcovs <- unlist(strsplit(opt$pbCOV, split = ",", fixed = TRUE))
    sampleNames <- paste0("sample_", as.character(seq_len(length(bams))))
    
    processBAM(
        bamfiles = bams,
        sample_names = sampleNames,
        reference_path = refPath,
        output_path = outPath,
        n_threads = opt$cores
    )
    # expect output as 
        # "[tempdir]/bams/sample_1.txt.gz"
        # "[tempdir]/bams/sample_1.cov"
    
    for(i in seq_len(length(bams))) {
        file.copy(
            from = file.path(outPath, paste0("sample_", as.character(i), ".txt.gz")), 
            to = pbouts[i], overwrite = TRUE
        )
        file.copy(
            from = file.path(outPath, paste0("sample_", as.character(i), ".cov")), 
            to = pbcovs[i], overwrite = TRUE
        )    
    }
    unlink(refPath, recursive = TRUE)
    unlink(outPath, recursive = TRUE)
} else if(opt$mode == "collateData") {
    # extract input ref
    refPath <- file.path(tempdir(), "reference")
    if(dir.exists(refPath)) unlink(refPath, recursive = TRUE)
    dir.create(refPath)
    untar(opt$ref, exdir = refPath)

    outPath <- file.path(tempdir(), "cd_outs")
    if(dir.exists(outPath)) unlink(outPath, recursive = TRUE)
    dir.create(outPath)

    pbouts <- unlist(strsplit(opt$pbOutput, split = ",", fixed = TRUE))
    pbcovs <- unlist(strsplit(opt$pbCOV, split = ",", fixed = TRUE))
    nxtNames <- unlist(strsplit(opt$nxtNames, split = ",", fixed = TRUE))

    # rename pbouts and pbcovs to assume nxtNames
    
    inPath <- file.path(tempdir(), "cd_ins")
    if(dir.exists(inPath)) unlink(inPath, recursive = TRUE)
    dir.create(inPath)
    for(i in seq_len(length(nxtNames))) {
        sampleName <- nxtNames[i]
        file.copy(
            from = pbouts[i], 
            to = file.path(inPath, paste0(sampleName, ".txt.gz")), 
            overwrite = TRUE
        )
        file.copy(
            from = pbcovs[i], 
            to = file.path(inPath, paste0(sampleName, ".cov")), 
            overwrite = TRUE
        )    
    }
    
    expr <- findSpliceWizOutput(inPath)
    doNovelSplicing <- (opt$nxtNovelSplicing == "yes")
    doPackageCov <- (opt$nxtPackageCOV == "yes")
    collateData(
        expr,
        reference_path = refPath,
        output_path = outPath,
        packageCOVfiles = FALSE,
        novelSplicing = doNovelSplicing,
        n_threads = opt$cores
    )
    
    # Not sure why collateData's packageCOVfiles doesn't work
    # - do it manually
    if(doPackageCov) {
        rds <- readRDS(file.path(outPath, "colData.Rds"))
        rds$df.files$cov_file <- expr$cov_file[
            match(rds$df.files$sample, expr$sample)
        ]
        
        dirCOV <- file.path(outPath, "COV")
        if(!dir.exists(dirCOV)) dir.create(dirCOV)
        breakCopyCOV <- FALSE
        for(i in seq_len(nrow(rds$df.files))) {
            cov <- rds$df.files$cov_file[i]
            if(isCOV(cov)) {      
                covOut <- file.path(dirCOV, basename(cov))
                file.copy(cov, covOut, overwrite = TRUE)
                if(file.exists(covOut)) {
                    rds$df.files$cov_file[i] <- covOut
                    message(cov, " is copied over as ", covOut)
                } else {
                    message(cov, " failed to be transferred - aborting")
                    breakCopyCOV <- TRUE
                    break
                }
            } else {
                message(cov, " is NOT a COV file")
                breakCopyCOV <- TRUE
                break
            }
        }
        if(!breakCopyCOV) {
            saveRDS(rds, file.path(outPath, "colData.Rds"))
        } else {
            unlink(dirCOV, recursive = TRUE)
        }
    }
    
    # tar output
    tmptar <- paste0(tempfile(), ".tar")
    
    setwd(outPath) # to archive all files under current directory
    filestodo <- list.files(".", recursive = TRUE, full.names = TRUE)
    print("Files to be packaged:")
    print(filestodo)
    utils::tar(tarfile = tmptar, files = filestodo, compression = "none")
    setwd(tempdir())
    
    file.copy(from = tmptar, to = opt$nxtSE, overwrite = TRUE)
    file.remove(tmptar)
    
    unlink(refPath, recursive = TRUE)
    unlink(inPath, recursive = TRUE)
    unlink(outPath, recursive = TRUE)
}