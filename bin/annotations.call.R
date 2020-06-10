#! /usr/bin/R

argsIn <- commandArgs(trailingOnly = TRUE)
source(argsIn[1])

##function to annotate GRanges from biomaRt
GRanges_anno_func("biomaRt_all.extensent.RData")

##Pfam annotation
##annotation of all transcripts with P(rotein)ID(entifier)
EnsDb.Hsapiens.v86_func("EnsDb.Hsapiens.v86.txs.PFAM.db.RData")
