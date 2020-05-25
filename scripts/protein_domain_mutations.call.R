#! /usr/bin/env R

##annotate protein domains
##determine if certain domains are more highly mutated: hotspots
##takes 3 arguments
##1 <- BASEDIR, base directory of project
##2 <- the *.glVcfGrList.RData object output from GRanges_geneLists_test.R script
##3 <- the path to genelists

argsIn <- commandArgs(trailingOnly = TRUE)
BASEDIR <- argsIn[1]
source(paste0(BASEDIR,"/scripts/protein_domain_mutations.func.R"))

INPUTFILE <- argsIn[2]
INPUTDIR <- dirname(INPUTFILE)
OUTDIR <- paste0(dirname(INPUTFILE), "/PID_output")
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
  dir(pattern=paste0(f,"_"), paste0(INPUTDIR,"/bootstraps"), full.names=TRUE)
})
names(BOOTSTRAPLIST) <- GENELISTNAMES

##Pfam annotation
##annotation of all transcripts with P(rotein)ID(entifier)
ANNOFILE <- paste0(BASEDIR,"/data/EnsDb.Hsapiens.v86.txs.PFAM.db.RData")
EnsDb.Hsapiens.v86_func(ANNOFILE)
load(ANNOFILE)

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
              per_gene_func(gr_txs_pid_snvs)
            }))

      per_gene_domain_func(gg, OUTDIR, GENELISTNAMES, g)
    }
  }
})

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
              per_gene_func(grb_txs_pid_snvs)
            }))
      print(paste0(names(BOOTSTRAPLIST)[b], "_", ff))
      per_gene_domain_func(ggb)
  })
  return(bootList)
})

##make ggridges plot
names(grPIDbootList) <- names(grPIDgeneList) <- names(BOOTSTRAPLIST)
genePIDridgeList <- lapply(seq_along(grPIDbootList), function(f){

                    gl <- data.frame(pg_var_per_width=grPIDgeneList[[f]][[1]][,"pg_var_per_width"],
                                 Geneset=names(grPIDbootList)[f],
                                 Dataset="Genelist", stringsAsFactors=F)
                    bl <- do.call(rbind, lapply(seq_along(grPIDbootList[[f]]), function(b){ data.frame(pg_var_per_width=grPIDbootList[[f]][[b]][[1]][,"pg_var_per_width"],
                                 Geneset=names(grPIDbootList)[f],
                                 Dataset=paste0("Bootlist_",b), stringsAsFactors=F)
                    }))
                    ##wilcox.test each bootstrap vs. genelist in turn
                    pvalList <- unlist(lapply(unique(bl$Dataset), function(bf){
                      wbl <- bl %>% dplyr::filter(Dataset %in% bf)
                      wt <- wilcox.test(wbl$pg_var_per_width, gl$pg_var_per_width, correct=TRUE, exact=FALSE)
                      wt$p.value
                    }))
                    bl$Dataset <- "Bootlist"
                    tt <- table(pvalList>0.01)["FALSE"]
                    return(list(rbind(gl, bl), unname((100-tt[1])/100), pvalList))
})
names(genePIDridgeList) <- names(grPIDbootList)

##summary of pvalues
pval_summary <- tibble(Geneset = names(genePIDridgeList), Bootstrap_p = do.call(rbind, lapply(names(genePIDridgeList), function(f){
  unlist(genePIDridgeList[[f]][[2]])
}))) %>% dplyr::rename("Bootstrap_p" = 2)

##construct full set
genePIDridgeMlt <- as_tibble(do.call(rbind, lapply(genePIDridgeList, function(f){
  f[[1]]
}))) %>% dplyr::mutate(Dataset = fct_rev(as.factor(Dataset))) %>%
         left_join()

pp <- ggplot(genePIDridgeMlt,aes(x=log2(pg_var_per_width), y=Geneset)) +
      geom_density_ridges(
        aes(fill = paste(Geneset, Dataset)),
        alpha = .6, color = "white"
      ) +
      scale_y_discrete(expand = c(0.01, 0)) +
      scale_x_continuous(expand = c(-0.1, 0)) +
      scale_fill_cyclical(
        breaks = c("CIT Genelist", "CIT Bootlist"),
        labels = c(`CIT Genelist` = "Genelist", `CIT Bootlist` = "Bootstrap"),
        values = c("grey40", "forestgreen"),
        name = "Dataset", guide = "legend"
      ) +
      geom_text(
        data = pval_summary,
        aes(y = Geneset, x = 3, label = paste0("p = ", Bootstrap_p)),
        nudge_y = 0.08
      )
ggsave(pp, file=paste0(OUTPUTFILE, ".genelist-vs-bootstraps.geom_density_ridges.pdf"), dpi=1200)

save(grPIDgeneList, grPIDbootList, genePIDridgeList,  file=OUTPUTFILE)
