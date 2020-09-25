# -------------------------------------------------------------------------------
# R script to generate a grpahic summary from long paired-end sequences as decribed in: 
# "########"
#
# Example:
# Rscript graphic_summary.R example_data/790v2_Leon_blast.summary 4004 example_data/790v2_Leon_ed_unp_short.bed graphic_summary.tiff
# Copyright Cintia Gómez-Muñoz. 2019
# cintia.gomez<at>ipicyt.edu.mx
# -------------------------------------------------------------------------------


# Load required libraries

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(RColorBrewer)


# Extract arguments form command line

args = commandArgs(trailingOnly=TRUE)
BesDf <- read.table(args[1], header = TRUE, sep = "\t")
BesNo <- as.numeric(args[2])
unpScanDf <-read.table(args[3], header = FALSE, sep = "\t")
outName <- args[4]

# Analyse data

# Total BESs aligned

algBesDf <- data.frame("BESs group" = c("Aligned BESs", "Not aligned BESs"), "Value" = c(BesDf[9,2], BesNo-BesDf[9,2]))

algPer <- round((algBesDf$Value*100)/BesNo, 2)

algPlot <- ggplot(algBesDf, aes(x = "", y = algBesDf$Value, fill = algBesDf$BESs.group)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0, direction = -1) +
  scale_fill_manual(values = c("#4575B4", "#91BFDB")) +
  guides(fill=guide_legend(title=NULL)) +
  geom_text_repel(label = c(" ", paste(algPer[2], "%")), nudge_x = c(0, 0.8))+
  geom_text(label = c(paste(algPer[1], "%"), " "), position = position_stack(vjust = 0.5, reverse = TRUE)) +
  theme_void()


# BESs type orientation

BesNoTot <- BesDf[c(1:5, 8),]

BesNoTot$BES.type <- factor(BesNoTot$BES.type, levels = rev(BesNoTot$BES.type))

BesNoPer <- round((BesNoTot[,2]*100)/BesNo, 2)

BesPlot <- ggplot(BesNoTot, aes(x = "", y = BesNoPer, fill = BesNoTot$BES.type)) +
  geom_bar(width = 0.5, stat = "identity") +
  scale_fill_brewer(palette="RdYlBu") +
  guides(fill=guide_legend(title=NULL)) +
  labs(x = "BESs type", y = "Percentage (%)") +
  geom_text(label = ifelse(BesNoPer > 1, paste(BesNoPer, "%"), ""), position = position_stack(vjust = 0.5))+
  theme_classic()


# Unpaired-end BESs

unpDf <- BesDf[c(7,6),]

unpPer <- round(((unpDf[,2]*100) / BesNo), 2)

unpPlot <- ggplot(unpDf, aes(x = unpDf$BES.type, y = unpPer, fill = unpDf$BES.type)) +
  geom_bar(width = 0.5, stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = c("#4575B4", "#FC8D59")) +
  guides(fill=guide_legend(title=NULL)) +
  labs(x = "Unpaired-end BESs", y = "Percentage (%)") +
  geom_text(label = paste(unpPer, "%"), position = position_stack(vjust = 0.5, reverse = TRUE))+
  theme_classic()


# Consistent unpaired-end BESs

colnames(unpScanDf) <- c("Contig", "Start", "End", "Partner")

unpTab <- ggtexttable(head(unpScanDf, 10), rows = NULL)


# Combine plots

combined_p <- ggarrange(algPlot, BesPlot, unpPlot, unpTab, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)


# Save image

tiff(file = outName, width=300, height=200, units="mm", res=600)
combined_p
dev.off()