#!/usr/bin/env Rscript
# Helper script to create QC report
# Stephen Zhang, 2018
# Usage:
# create_qc.r [bamFilePath] [sampleName] [genomeBioStrings] [genomeTxDb]
#             [Rmd_script] [outputDir]
args <- commandArgs(trailingOnly = T)
bamFilePath <- args[1]
sampleName <- args[2]
genomeBioStrings <- args[3]
genomeTxDb <- args[4]
rmdPath <- args[5]
outputDir <- args[6]

library(rmarkdown)
render(rmdPath, output_format = "html_document",
       output_dir = outputDir,
       intermediates_dir = outputDir,
       knit_root_dir = outputDir,
       clean = T,
       params = list(
         bamFilePath = bamFilePath,
         sampleName = sampleName,
         genomeBioStrings = genomeBioStrings,
         genomeTxDb = genomeTxDb
       ))
