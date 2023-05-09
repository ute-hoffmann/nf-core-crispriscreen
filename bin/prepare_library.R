#!/usr/bin/env Rscript

# Purpose:
# R script to import a *.fasta file and convert it to an *.saf genome
# annotation file. This script is part of the nextflow pipeline
# 'nf-core-crispriscreen' for processing of CRISPRi library data.
# The pipeline can be found at https://github.com/MPUSP/nf-core-crispriscreen.
#
# Date: 2022-04-11
# Author: Michael Jahn, PhD
# Affilation: Science For Life Laboratory (KTH), Stockholm, Sweden

# input parameters
args <- commandArgs(trailingOnly = TRUE)
fasta_file <- args[1] # path to fasta file, mandatory
gene_controls <- args[2] # pattern for control barcodes, default: "" aka empty string

# read and process fasta file
## Install Biostrings if not available
# add library dir for additional packages
homedir <- system("echo ${HOME}", intern = TRUE)
libdir <- .libPaths()
if (any(grepl(homedir, libdir))) {
  pkdir <- grep(homedir, libdir, value = TRUE)
} else {
  pkdir <- paste0(homedir, "/.R")
  system(paste0("mkdir ", pkdir))
  .libPaths(new = pkdir)
}
list_bioc_packages <- c("Biostrings")
list_to_install <- setdiff(list_bioc_packages, rownames(installed.packages()))
if (length(list_to_install)) {
  message(paste0("Missing package(s) ", paste(list_to_install, collapse = ", "), " are installed to '", pkdir, "'."))
  install.packages(pkgs = "BiocManager", lib = pkdir, repos = "https://cloud.r-project.org")
  BiocManager::install(
    pkgs = list_to_install, lib = pkdir,
    update = FALSE, ask = FALSE
  )
}

library(Biostrings)
## parse fasta file with Biostrings function to ensure error handling / ensuring that really fasta format etc.
fasta_df = readDNAStringSet(fasta_file)
saf_df <- data.frame(
  GeneID = names(fasta_df),
  Sequence = paste(fasta_df)
)

# check for duplications
stopifnot(!any(duplicated(saf_df$GeneID)))
stopifnot(!any(duplicated(saf_df$Sequence)))

# add remaining columns
saf_df$Chr <- saf_df$GeneID
saf_df$Start <- 1
saf_df$End <- sapply(saf_df$Sequence, nchar)
saf_df$Strand <- "*"
saf_df <- saf_df[c("GeneID", "Chr", "Start", "End", "Strand", "Sequence")]

# replace file ending by .saf (everything after last . in file name is interpreted as file ending) - if there is no file ending, append ".saf"
if (length(strsplit(fasta_file, "\\.")[[1]][-1])) {
  saf_name <- gsub(strsplit(fasta_file, "\\.")[[1]][-1], "saf", basename(fasta_file))
} else {
  saf_name <- paste(fasta_file, ".saf", sep="")
}
write.table(
    x = saf_df, file = saf_name, sep = "\t", row.names = FALSE,
    col.names = FALSE, quote = FALSE
)

# optionally identify controls and save them for reference
if (gene_controls != "") {
    ctrl_hits <- grepl(gene_controls, saf_df$GeneID)
    if (!sum(ctrl_hits)) {
        stop(paste0("No barcode name matches supplied pattern '", gene_controls, "'."))
    }
    saf_df_ctrl <- subset(saf_df, ctrl_hits)
    ctrl_name <- gsub(".fasta$", "_controls.tsv", basename(fasta_file))
    write.table(
        x = saf_df_ctrl, file = ctrl_name, sep = "\t", row.names = FALSE,
        col.names = FALSE, quote = FALSE
    )
}
