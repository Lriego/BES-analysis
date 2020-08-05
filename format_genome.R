# -------------------------------------------------------------------------------
# R script to format FASTA genome headers to look like ContigName|siezXXXXX
# "########"
#
# Example:
# Rscript format_genome.R genome.fasta out_name.fasta
# Copyright Cintia Gómez-Muñoz. 2020
# cintia.gomez<at>ipicyt.edu.mx
# -------------------------------------------------------------------------------

# Load required libraries
library(Biostrings)

# Extract arguments form command line

args = commandArgs(trailingOnly=TRUE)

Fas <- readDNAStringSet(filepath = args[1])

OutFile <- args[2]

# Extract data

Nm <- names(Fas)

Sz <- width(Fas)

FasEd <- Fas

names(FasEd) <- paste0(Nm, "|size", Sz)

writeXStringSet(FasEd, filepath = OutFile)

