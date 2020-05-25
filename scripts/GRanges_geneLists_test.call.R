#! /bin/R

##script reads sorted VCFs, annotates to gene level, subsets to geneLists
##random sampling of same size as geneLists for Wilcoxon tests
##test: divergent level of mutations in genes

#input args
argsIn <- commandArgs(trailingOnly = TRUE)
RLIBPATH <- argsIn[1]
BASEDIR <- argsIn[2]
WORKDIR <- argsIn[3]
OUTDIR <- argsIn[4]
VCFIN <- dir(pattern="tumour.vcf$", WORKDIR, full.names=TRUE)
GENELISTSET <- paste(dir(pattern=".txt$", argsIn[5], full.names=TRUE), collapse=",")
RUNTEST <- argsIn[6]

if(length(argsIn)<6){
  RUNTEST <- "run"
}
.libPaths(RLIBPATH)
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
setwd(BASEDIR)
source("scripts/vcfParseBiomartAnno2GR.func.R")

GENELIST <- unlist(strSplitVec(GENELISTSET, ","))
GLDF <- t(as.data.frame(strSplitVec(GENELIST,"/")))
GENELISTNAMES <- strSplitVec(GLDF[,ncol(GLDF)],"\\.")[1,]
VCFNAME <- rev(strSplitVec(VCFIN, "/"))[1]

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
  print(genelist)

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
    ggo <- c(ggo,grep(paste(ggoList[[gg]],collapse="|"), vcfGr$external_gene_name, perl=TRUE))
  }
  glVcfGr <- unique(vcfGr[ggo])
  vpgGl <- variantsPerGene(glVcfGr, geneVecIn=genelist$Gene_Name, rownameTag="TCGA")
  genesFound[[x]] <- vpgGl
  names(genesFound)[x] <- paste0(length(vpgGl), "/", length(genelist$Gene_Name))

  ##adjust to gene sizes
  annovpgGl <- annoGenenameEnsEnt %>% filter(external_gene_name %in% names(vpgGl)) %>% filter(chromosome_name %in% c(1:22,"X","Y")) %>% mutate(length = end_position - start_position) %>% group_by(external_gene_name) %>% slice(which.max(length)) %>% arrange(external_gene_name)
  glVcfGrList[[x]] <- glVcfGr
  kbVpgGlList[[x]] <- vpgGl / (as.vector(annovpgGl$length)/1000)
}

##add full VcfGr to list
glVcfGrList[["fullVcfGr"]] <- vcfGr
names(glVcfGrList) <- c(GENELISTNAMES, "fullVcfGr")

##save output
assignName <- paste0(VCFNAME,".glVcfGrList")
assign(assignName,value=glVcfGrList)
saveFile <- paste0(OUTDIR,"/",VCFNAME,".glVcfGrList.RData")
save(list=assignName, file=paste0(saveFile))

uniqGenes <- sort(unique(grep(",",vcfGr$external_gene_name, invert=TRUE, value=TRUE)))

#if(RUNTEST == "run"){
  genePVecList <- lapply(seq_along(GENELIST), function(x){
    print(paste0("Working on: ", GENELIST[x]))

    pVec <- c()

    pb <- txtProgressBar(min=1, max=100, initial=1, style=3)
    for (p in 1:100){
      setTxtProgressBar(pb,p)
      geneVec100 <- sample(uniqGenes,size=length(genesFound[[x]]), replace=FALSE)

      ##subset by three possible matches:
      ##1:gene alone; 2:gene with comma; 3:comma with gene
      ggoList100 <- as.list(1:4)
      ggoList100[[1]] <- paste0("^",geneVec100,"$")
      ggoList100[[2]] <- paste0("^",geneVec100,",")
      ggoList100[[3]] <- paste0(",",geneVec100,"$")
      ggoList100[[4]] <- paste0(",",geneVec100,",")
      ggo100 <- c()
      for(g in 1:4){
        #print(paste0("Working on :",x))
        ggo100 <- c(ggo100,grep(paste(ggoList100[[g]],collapse="|"), vcfGr$external_gene_name, perl=TRUE))
      }
      foundGenes100 <- unique(vcfGr$external_gene_name[ggo100])
      glVcfGr100 <- unique(vcfGr[ggo100])
      vpgGl100 <- variantsPerGene(glVcfGr100, geneVecIn=geneVec100, rownameTag="TCGA")
      vpgGl100 <- vpgGl100[names(vpgGl100) %in% geneVec100]

      #per Kb
      annovpgGl100 <- annoGenenameEnsEnt %>% filter(external_gene_name %in% names(vpgGl100)) %>% filter(chromosome_name %in% c(1:22,"X","Y")) %>% mutate(length = end_position - start_position) %>% group_by(external_gene_name) %>% slice(which.max(length)) %>% arrange(external_gene_name)

      vpgGl100 <- vpgGl100[names(vpgGl100) %in% as.vector(annovpgGl100$external_gene_name)]
      kbvpgGl100 <- vpgGl100 / (as.vector(annovpgGl100$length)/1000)

      wto <- wilcox.test(as.vector(kbVpgGlList[[x]]), as.vector(kbvpgGl100), exact=FALSE)
      pVec <- c(pVec, wto$p.value)
    }
    close(pb)
    #pVec <- p.adjust(pVec,method="BH")
    return(list(pVec, glVcfGr100))
  })
  names(genePVecList) <- paste(GENELISTNAMES, names(genesFound))

  ##plot all pVec density to visualise which sets have more significant results
  ##first make melt DF of all pVecs and associated genelist
  genePVecMlt <- do.call(rbind, lapply(seq_along(genePVecList), function(f){
                  if(!is.null(genePVecList[[f]][[1]])){
                    data.frame(p.value=genePVecList[[f]][[1]],
                               genelist=names(genePVecList)[f])
                             }
                           }))
  genePVecMlt$p.adj <- p.adjust(genePVecMlt$p.value, method="BH")

  ##ggplot
  setwd(OUTDIR)

  pp <- ggplot(genePVecMlt,aes(x=log2(p.value), y=genelist, fill=genelist)) +
  geom_density_ridges() +
  geom_vline(xintercept = log2(0.01), colour="red", linetype = "dashed") +
  ggtitle(VCFNAME) +
  guides(fill=FALSE)
  ggsave(pp, file=paste0(VCFIN,".genelists.geom_density-pvalue.nm.pdf"), dpi=1200)

  qq <- ggplot(genePVecMlt,aes(x=log2(p.adj), y=genelist, fill=genelist)) +
  geom_density_ridges() +
  geom_vline(xintercept = log2(0.01), colour="red", linetype = "dashed") +
  ggtitle(paste0(VCFNAME, " - BH p.value adjusted")) +
  guides(fill=FALSE)
  ggsave(qq, file=paste0(VCFIN,".genelists.geom_density-padjBH.nm.pdf"), dpi=1200)

  ##save data
  assignName <- paste0(VCFNAME,".pVecList")
  assign(assignName,value=genePVecList)
  saveFile <- paste0(OUTDIR,"/",VCFNAME,".nmList.RData")
  save(list=assignName, file=paste0(saveFile))
#}
