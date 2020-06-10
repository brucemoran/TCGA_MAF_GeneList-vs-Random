#! /usr/bin/R

##read VCFs into Granges

libs <- c("ensembldb", "EnsDb.Hsapiens.v86", "ensemblVEP", "customProDB", "tidyverse", "biomaRt", "magrittr", "reshape2", "ggridges", "PFAM.db")
libsLoaded <- lapply(libs,function(l){suppressMessages(library(l, character.only = TRUE))})

strSplitFun <- function(input,sepn){
  lapply(input,function(f){strsplit(f,sepn)[[1]]})
}
strSplitVec <- function(inVec,sepn){
  sapply(seq_along(inVec),function(f){strsplit(inVec[f],sepn)[[1]]})
}

##function to annotate GRanges from biomaRt
vcfParseAnnoGR <- function(vcfIn, annoIn){

  ##read inputs
  print("Reading VCF input...")
  grVcf <- granges(readVcf(file=vcfIn))
  gr <- suppressWarnings(InputVcf(vcfIn))
  chrTest <- length(grep("chr",levels(seqnames(gr[[1]]))[1]))>0

  if(!class(annoIn)[1] == "GRanges"){
    if(chrTest==TRUE){
      annoIn$seqnames <- paste0("chr", annoIn$seqnames)
    }
    biomartAll <- GRanges(annoIn)
  }
  if(class(annoIn)[1] == "GRanges"){
    if(chrTest==TRUE){
      annoin <- as_tibble(annoIn)
      annoin$seqnames <- paste0("chr", annoin$seqnames)
      biomartAll <- GRanges(annoin)
    } else {
      biomartAll <- annoIn
    }
  }

  ##f(ind) o(ver)l(aps) with annotation
  fol <- as_tibble(as.data.frame(findOverlaps(grVcf, biomartAll)))

  print("Annotating...")
  grVcfAnno <- do.call(rbind, lapply(unique(fol$queryHits), function(ff){
    shits <- fol[fol$queryHits==ff,]$subjectHits
     do.call(cbind, apply(as.data.frame(values(biomartAll[shits])), 2, function(f){
       data.frame(paste(gsub(" ","",f),collapse=","))
     }))
   }
  ))
  names(grVcfAnno) <- names(values(biomartAll))

  ##take all grVcf in fol
  grVcf <- grVcf[unique(fol$queryHits)]
  values(grVcf) <- cbind(values(grVcf), DataFrame(grVcfAnno))

  ##apply samples GT into mcols of grVcf

  grT <- DataFrame(do.call(cbind, lapply(seq_along(names(gr)), function(f){
    grt <- gr[[names(gr)[f]]]
    grta <- grt[unique(fol$queryHits)]$GT
    as.data.frame(grta)
  })))
  names(grT) <- names(gr)
  values(grVcf) <- cbind(values(grVcf), grT)
  return(list(grVcf, biomartAll))
}

##return vector of mutations per annotated genes
variantsPerGene <- function(vcfIn, geneVecIn, rownameTag){

  extGeneCol <- grep("external_gene_name",colnames(values(vcfIn)))
  rowTagCols <- grep(rownameTag,colnames(values(vcfIn)))
  genoVcf <- vcfIn[,c(extGeneCol, rowTagCols)]
  ugeneVec <- unique(genoVcf$external_gene_name)
  geneVec <- genoVcf$external_gene_name
  allTab <- table(unlist(strsplit(geneVec,",")))
  varPerGene <- allTab[names(allTab) %in% geneVecIn]
  return(varPerGene)

}
