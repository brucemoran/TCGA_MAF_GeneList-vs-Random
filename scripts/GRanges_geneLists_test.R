#! /bin/R

##script reads sorted VCFs, annotates to gene level, subsets to geneLists
##random sampling of same size as geneLists for Wilcoxon tests
##test: divergent level of mutations in genes

#input args
args <- commandArgs(trailingOnly = TRUE)
RLIBPATH <- args[1]
WORKDIR <- OUTDIR <- args[2]
VCFIN <- dir(pattern="tumour.vcf$", WORKDIR, full.names=TRUE)
GENELISTSET <- args[3] #comma-sep

.libPaths(RLIBPATH)
source("scripts/vcfParseBiomartAnno2GR.func.R")

GENELIST <- unlist(strSplitVec(GENELISTSET, ","))
GENELISTNAMES <- strSplitVec(strSplitVec(GENELIST,"/")[3,],"\\.")[1,]
VCFNAME <- rev(strSplitVec(VCFIN, "/"))[1]

##parse, annotate VCF into Grange
vcfGr <- vcfParseAnnoGR(VCFIN)
extEnt <- tibble(entrezgene = unlist(strSplitVec(vcfGr$entrezgene, ",")),
                 external_gene_name = unlist(strSplitVec(vcfGr$external_gene_name, ",")))
		 stspGene <- strSplitVec(vcfGr$external_gene_name, ",")

##parse out geneLists
genePVecList <- lapply(GENELIST, function(genelist){
  ggenelist <- suppressMessages(read_tsv(genelist))
    print(paste0("Working on: ", genelist, " (n = ", dim(ggenelist)[1], ")"))
      glVcfGr <- unique(vcfGr[sort(unlist(lapply(ggenelist[[grep("Gene_Name",colnames(ggenelist))]],function(gene){
          grep(gene, stspGene)})))])
	    vpgGl <- variantsPerGene(glVcfGr, rownameTag="TCGA")

  ##test vs. 100 genelists of same size
    pVec <- c()
      sanitVcfGr <- vcfGr[vcfGr$external_gene_name %in% unique(unlist(stspGene)),]
        pb <- txtProgressBar(min=1, max=100, initial=1, style=3)
	  for (x in 1:100){
	      setTxtProgressBar(pb,x)
	          geneVec100 <- unique(sanitVcfGr[floor(runif(length(unique(glVcfGr$external_gene_name)),
		                                 min=1,
						                                max=length(sanitVcfGr)))]$external_gene_name)
										    ##subset by three possible matches:
										        ##1:gene alone; 2:gene with comma; 3:comma with gene
											    matchVec <- unique(unlist(lapply(geneVec100,function(gene){
											          o1 <- paste0("^", gene, "$")
												        o2 <- paste0(",",gene,",")
													      o22 <- paste0(",",gene,"$")
													            o3 <- paste0("^", gene, ",")
														          return(c(o1,o2,o22,o3))
															      })))
															          glVcfGr100 <- unique(vcfGr[unlist(lapply(matchVec,function(mgene){grep(mgene,vcfGr$external_gene_name)}))])
																      vpgGl100 <- variantsPerGene(glVcfGr100, rownameTag="TCGA")

    wto <- wilcox.test(c(unlist(vpgGl)), c(unlist(vpgGl100)), exact=FALSE)
        pVec <- c(pVec, wto$p.value)
	  }
	    close(pb)
	      return(pVec)
	      })
	      names(genePVecList) <- GENELISTNAMES

##plot all pVec density to visualise which sets have more significant results
##first make melt DF of all pVecs and associated genelist
genePVecMlt <- melt(genePVecList)
colnames(genePVecMlt) <- c("p.value", "genelist")
##ggplot
pp <- ggplot(genePVecMlt,aes(x=log2(p.value),fill=genelist)) +
geom_density(position="stack", alpha=0.6) +
geom_vline(xintercept = log2(0.01), colour="red", linetype = "dashed") +
ggtitle(VCFNAME)

ggsave(pp, file=paste0(VCFIN,"genelists.geom_density-pvalue.pdf"), dpi=1200)
