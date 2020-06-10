#! /usr/bin/env R

##annotate protein domains
##determine if certain domains are more highly mutated: hotspots
##takes 3 arguments
##1 <- functions file
##2 <- the *.glVcfGrList.RData object output from GRanges_geneLists_test.R script
##3 <- the path to genelists

argsIn <- commandArgs(trailingOnly = TRUE)
source(argsIn[1])
INPUTFILE <- argsIn[2]
OUTDIR <- paste0("/PID_output")

dir.create(OUTDIR, showWarnings=F, recursive = T)
OUTPUTFILE <- gsub("glVcfGrList", "PIDTableList", INPUTFILE)
RDATANAME <- rev(strsplit(INPUTFILE, "/")[[1]])[1]
GLVCFGRLIST <- gsub(".RData", "", RDATANAME)
suppressMessages(load(INPUTFILE))
GLVCFGRLIST1 <- suppressMessages(get(GLVCFGRLIST))

GENELISTSET <- paste(dir(pattern=".txt$", argsIn[3], full.names=TRUE), collapse=",")
GENELIST <- unlist(str_split(GENELISTSET, ","))
GENELISTNAMES <- unlist(lapply(GENELIST, function(f){
  gsub(".txt", "", rev(unlist(str_split(f, "/")))[1])
}))
BOOTSTRAPLIST <- lapply(GENELISTNAMES, function(f){
  dir(pattern=f, "bootstraps", full.names=TRUE)
})
names(BOOTSTRAPLIST) <- GENELISTNAMES

##Pfam annotation
##annotation of all transcripts with P(rotein)ID(entifier)
load("EnsDb.Hsapiens.v86.txs.PFAM.db.RData")

#load("results/TCGA.BRCA.mutect/results/TCGA.BRCA.mutect/TCGA.BRCA.mutect.somatic.tumour.vcf.glVcfGrList.RData")
grPIDgeneList <- lapply(seq_along(GLVCFGRLIST1), function(g){

  if(names(GLVCFGRLIST1)[g] != "fullVcfGr"){

    print(paste0("Working on: ", names(GLVCFGRLIST1)[g]))
    gr <- GLVCFGRLIST1[[g]][,c()]

    if(length(gr)!=0){

      ##annotate PID
      gr_txs_pid <- sqlvl_pid_gr_func(gr, txs_pid)

      ################
      ## per gene ##
      ################
      ##domains within a gene may have increased mutatbility
      ##here this is tested
      ##multiple annotations sources per mutation exist
      ##need to check per source as may influence result

      ##lapply(unique(gr_txs_pid$protein_domain_source) for all sources
      s <- "pfam"
      writeLines(paste0("\n", s))
      pb <- txtProgressBar(min=0, max=length(unique(gr_txs_pid$gene_id)),  style = 3)
      gg <- do.call(rbind, lapply(seq_along(unique(gr_txs_pid$gene_id)), function(f){
              setTxtProgressBar(pb, value=f, title = NULL, label = NULL)
              ##get unique SNVs
              gr_txs_pid_snvs <- gr_txs_pid %>%
                                 plyranges::filter(gene_id %in% unique(gr_txs_pid$gene_id)[f])
              ##create a mean domain width for multiple annotations
              ##here we remove protein annotation to reduce the set to domain-gene-level
              ##count number of variants per each PID, dividing by width of domain
              per_gene_func(gr_txs_pid_snvs, s)
            }))

      per_gene_domain_func(gg, GENELISTNAMES, g)
    }
  }
})
names(grPIDgeneList) <- names(GLVCFGRLIST1)

grPIDbootList <- lapply(seq_along(BOOTSTRAPLIST), function(b){

    print(paste0("Working on: ", names(BOOTSTRAPLIST)[b]))
    bootList <- lapply(seq_along(BOOTSTRAPLIST[[b]]), function(ff){
      genesb <- suppressMessages(read_tsv(BOOTSTRAPLIST[[b]][ff])) %>%
                left_join(., annoGenenameEnsEnt, by=c("Gene_Name" = "external_gene_name")) %>%
                dplyr::select(ensembl_gene_id) %>%
                unlist()
      genesb_txs_pid <- txs_pid %>% plyranges::filter(gene_id %in% genesb)
      grb_txs_pid <- sqlvl_pid_gr_func(GLVCFGRLIST1$fullVcfGr, genesb_txs_pid) %>%
                     plyranges::select(-contains("TCGA"), -c("paramRangeID","REF","ALT","QUAL","FILTER","external_gene_name","ensembl_gene_id","entrezgene_id"))

      if(length(unique(grb_txs_pid$gene_id)) != 0){

        ################
        ## per gene ##
        ################
        ##domains within a gene may have increased mutatbility
        ##here this is tested
        ##multiple annotations sources per mutation exist
        ##need to check per source as may influence result
        pbb <- txtProgressBar(min=0, max=length(unique(grb_txs_pid$gene_id)),  style = 3)
        ggb <- do.call(rbind, lapply(seq_along(unique(grb_txs_pid$gene_id)), function(f){
                setTxtProgressBar(pbb, value=f, title = NULL, label = NULL)
                ##get unique SNVs
                grb_txs_pid_snvs <- grb_txs_pid %>%
                                   plyranges::filter(gene_id %in% unique(grb_txs_pid$gene_id)[f])
                ##create a mean domain width for multiple annotations
                ##here we remove protein annotation to reduce the set to domain-gene-level
                ##count number of variants per each PID, dividing by width of domain
                per_gene_func(grb_txs_pid_snvs, s="pfam")
              }))
        print(paste0(names(BOOTSTRAPLIST)[b], "_", ff))
        per_gene_domain_func(ggb, names(BOOTSTRAPLIST)[b], b)
      }
  })
  return(bootList)
})

##make ggridges plot
names(grPIDbootList) <- names(BOOTSTRAPLIST)
genePIDridgeList <- lapply(names(grPIDbootList), function(f){

  if(!is.null(grPIDgeneList[[f]])){
    print(f)
    gl <- data.frame(pg_var_per_width=grPIDgeneList[[f]][[1]][,"pg_var_per_width"],
                 Geneset=f,
                 Dataset="Genelist", stringsAsFactors=F)
    bl <- do.call(rbind, lapply(seq_along(grPIDbootList[[f]]), function(b){
      if(!is.null(grPIDbootList[[f]][[b]][[1]])){
      data.frame(pg_var_per_width=grPIDbootList[[f]][[b]][[1]]$pg_var_per_width,
                 Geneset=f,
                 Dataset=paste0("Bootlist_", b), stringsAsFactors=F)
      }
    }))

    ##wilcox.test each bootstrap vs. genelist in turn
    pvalList <- unlist(lapply(unique(bl$Dataset), function(bf){
      wbl <- bl %>% dplyr::filter(Dataset %in% bf)
      wt <- wilcox.test(wbl$pg_var_per_width, gl$pg_var_per_width, correct=TRUE, exact=FALSE)
      wt$p.value
    }))
    bl$Dataset <- "Bootlist"
    tt <- table(pvalList>0.01)["FALSE"]
    if(is.na(tt)){
      tt <- 0
    }
    return(list(rbind(gl, bl), unname((100-tt[1])/100), pvalList))
  }
})
names(genePIDridgeList) <- names(grPIDbootList)

##summary of pvalues
pval_summary <- tibble(Geneset = names(genePIDridgeList), Bootstrap_p = do.call(rbind, lapply(names(genePIDridgeList), function(f){
  if(is.null(genePIDridgeList[[f]][[2]])){
    return(1)
  } else {
    unlist(genePIDridgeList[[f]][[2]])
  }
}))) %>% dplyr::rename("Bootstrap_p" = 2)

##construct full set
genePIDridgeMlt <- as_tibble(do.call(rbind, lapply(genePIDridgeList, function(f){
  f[[1]]
}))) %>% dplyr::mutate(Dataset = fct_rev(as.factor(Dataset))) %>%
         left_join(., pval_summary)
colnames(genePIDridgeMlt) <- c("pg_var_per_width", "Geneset", "Dataset",  "pval_summary")
x_text <- max(log2(genePIDridgeMlt$pg_var_per_width))+0.5
xlims <- c(min(log2(genePIDridgeMlt$pg_var_per_width))-0.1,
           max(log2(genePIDridgeMlt$pg_var_per_width))+1.5)
pp <- ggplot(genePIDridgeMlt, aes(x=log2(pg_var_per_width), y=Geneset)) +
      geom_density_ridges(
        aes(fill = paste(Geneset, Dataset)),
        alpha = .6, color = "white"
      ) +
      scale_y_discrete(expand = c(0.01, 0)) +
      scale_fill_cyclical(
        breaks = c("CIT Genelist", "CIT Bootlist"),
        labels = c(`CIT Genelist` = "Genelist", `CIT Bootlist` = "Bootstrap"),
        values = c("grey40", "forestgreen"),
        name = "Dataset", guide = "legend"
      ) +
      geom_text(
        data = pval_summary,
        aes(y = Geneset, x = x_text, label = paste0("p = ", Bootstrap_p)),
        nudge_y = 0.08
      ) +
      xlim(xlims)
ggsave(pp, file=paste0(OUTPUTFILE, ".genelist-vs-bootstraps.geom_density_ridges.pdf"), dpi=1200)

save(grPIDgeneList, grPIDbootList, genePIDridgeList,  file=OUTPUTFILE)
