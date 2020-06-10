#! /bin/R

##script reads sorted VCFs, annotates to gene level, subsets to geneLists
##random sampling of same size as geneLists for Wilcoxon tests
##test: divergent level of mutations in genes

#input args
argsIn <- commandArgs(trailingOnly = TRUE)
source(argsIn[1])
TCGAID <- argsIn[2]
VCFIN <- argsIn[3]
GENELISTSET <- paste(dir(pattern=".txt$", argsIn[4], full.names=TRUE), collapse=",")

GENELIST <- unlist(strSplitVec(GENELISTSET, ","))
GLDF <- t(as.data.frame(strSplitVec(GENELIST,"/")))
GENELISTNAMES <- strSplitVec(GLDF[,ncol(GLDF)],"\\.")[1,]
VCFNAME <- rev(strSplitVec(VCFIN, "/"))[1]

##output dirs
dir.create("results", showWarnings=FALSE)

##parse, annotate VCF into Grange
vcfGr <- vcfParseAnnoGR(VCFIN)
extEnt <- tibble(entrezgene_id = unlist(strSplitVec(vcfGr$entrezgene_id, ",")),
                 external_gene_name = unlist(strSplitVec(vcfGr$external_gene_name, ",")))
annoMart <- useMart(biomart="ensembl",
                    dataset="hsapiens_gene_ensembl",
                    host="www.ensembl.org")
annoGenenameEnsEnt <- as_tibble(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand', 'external_gene_name', 'ensembl_gene_id', 'entrezgene_id'), mart = annoMart))

##parse out geneLists
genesFound <- kbVpgGlList <- glVcfGrList <- as.list(GENELISTNAMES)

for (x in seq_along(GENELIST)){
  genelist <- suppressMessages(read_tsv(GENELIST[x]))
  #print(x)

  ##grep out positions supporting gene_name
  ##NB these need to be comma'd, and regexed to allow overlap with other genes
  ##these are then accounted for and summed further on
  ggoList <- as.list(1:4)
  ggoList[[1]] <- paste0("^",genelist[[2]],"$")
  ggoList[[2]] <- paste0("^",genelist[[2]],",")
  ggoList[[3]] <- paste0(",",genelist[[2]],"$")
  ggoList[[4]] <- paste0(",",genelist[[2]],",")
  ggo <- c()
  for(gg in 1:4){
    #print(paste0("Working on :",x))
    ggo <- c(ggo, grep(paste(ggoList[[gg]],collapse="|"), vcfGr$external_gene_name, perl=TRUE))
  }
  glVcfGr <- unique(vcfGr[ggo])
  vpgGl <- variantsPerGene(glVcfGr, geneVecIn=genelist$Gene_Name, rownameTag="TCGA")
  genesFound[[x]] <- vpgGl
  names(genesFound)[x] <- paste0(length(vpgGl), "/", length(genelist$Gene_Name))

  ##adjust to gene sizes
  annovpgGl <- annoGenenameEnsEnt %>%
               dplyr::filter(external_gene_name %in% names(vpgGl)) %>%
               dplyr::filter(chromosome_name %in% c(1:22,"X","Y")) %>%
               dplyr::mutate(length = end_position - start_position) %>%
               group_by(external_gene_name) %>%
               dplyr::slice(which.max(length)) %>%
               arrange(external_gene_name)

  glVcfGrList[[x]] <- glVcfGr
  kbVpgGlList[[x]] <- vpgGl / (as.vector(annovpgGl$length)/1000)

  if(length(kbVpgGlList[[x]])!=0){
    oSNV_tb <- left_join(tibble(Gene_Name=names(vpgGl), Total_SNV=vpgGl),
                         tibble(Gene_Name=names(kbVpgGlList[[x]]), SNV_per_KB=kbVpgGlList[[x]]))
    write_tsv(oSNV_tb, paste0("results/", TCGAID,".", GENELISTNAMES[x], ".txt"))
  } else {
    write_tsv(tibble(Gene_Name=NA, Total_SNV=NA, SNV_per_KB=NA)[-1,], paste0("results/", TCGAID,".", GENELISTNAMES[x], ".txt"))
  }
}

##add full VcfGr to list
glVcfGrList[["fullVcfGr"]] <- vcfGr
names(glVcfGrList) <- c(GENELISTNAMES, "fullVcfGr")

##save output
assignName <- paste0(TCGAID,".glVcfGrList")
assign(assignName,value=glVcfGrList)
saveFile <- paste0(TCGAID,".glVcfGrList.RData")
save(list=assignName, file=paste0(saveFile))

##Pfam annotation
##annotation of all transcripts with P(rotein)ID(entifier)
ANNOFILE <- "EnsDb.Hsapiens.v86.txs.PFAM.db.RData"
EnsDb.Hsapiens.v86_func(ANNOFILE)
