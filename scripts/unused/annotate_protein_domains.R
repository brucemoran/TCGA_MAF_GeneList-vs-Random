#! /usr/bin/env R

##annotate resulting data structure from 'protein_domain_mutations.R'
##here we use the PFAM.db library with keys mapped to multiple annotations
##all annotations are then mapped back to PFAM descriptions
##this generates a tibble of all PFAM description and associated IDs in the db

libs <- c("PFAM.db", "tidyverse", "plyranges", "biomaRt")
libsLoaded <- lapply(libs,function(l){suppressMessages(library(l, character.only = TRUE))})
source("scripts/vcfParseBiomartAnno2GR.func.R")

##two arguments, the grPIDVcfList and the path to genelists
argsIn <- commandArgs(trailingOnly = TRUE)

##grPIDList is a list
##[[1]] is GRanges object for geneset
##[[2]] is occurrence of domains in the geneset, per annotation source
##[[3]] is occurrence of domains across all genes -> deprecated

##annotate each domain, counting occurrence per gene also
##output a tibble of protein_domain_id, protein_domain_source, gene_count, genes
#f<-grPIDList[[1]]
grPIDAnnoDomList <- lapply(grPIDList, function(f){
    gr <- f[[1]] %>% plyranges::filter(protein_domain_source %in% "pfam") %>%
                     plyranges::select(2, 3, 4, 5, gene_id) %>%
                     unique() %>%
                     mutate(domain_width = prot_dom_end - prot_dom_start)
    ##mutations per pid within gene
    count_gene_pid <- as_tibble(mcols(gr)) %>%
                      dplyr::count(protein_domain_id, gene_id)
    ##mutations per pid total
    count_gene_pid_n <- count_gene_pid %>%
                        dplyr::group_by(protein_domain_id) %>%
                        summarise(nn = sum(n)) %>%
                        left_join(count_gene_pid, .) %>%
                        dplyr::rename(gene_pid_n=n, pid_n=nn) %>%
                        left_join(., unique(as_tibble(mcols(gr))[c("protein_domain_id", "domain_width")])) %>%
                        dplyr::group_by(protein_domain_id) %>%
                        dplyr::mutate(domain_width_mean = mean(domain_width)) %>%
                        ungroup() %>%
                        dplyr::select(-domain_width) %>%
                        distinct() %>%
                        dplyr::mutate(variants_per_bp = gene_pid_n/domain_width_mean) %>%
                        dplyr::arrange(desc(variants_per_bp))

    left_join(count_gene_pid_n, pfam_tib, by=c("protein_domain_id" = "pfam")) %>%
    left_join(., annoGenenameEnsEnt, by=c("gene_id" = "ensembl_gene_id")) %>%
    dplyr::select(1, pfam_de, variants_per_bp, gene_id, external_gene_name, everything(), -pfam_prosite, -pfam_smart, -pfam_prints) %>% distinct()
})
names(grPIDAnnoDomList) <- GENELISTNAMES

##write output to CSV
dir.create(paste0(OUTDIR,"/per_domain_gene"), showWarnings=F, recursive = T)
lapply(seq_along(grPIDAnnoDomList), function(wo){
  nmo1 <- GENELISTNAMES[[wo]]
  if(dim(grPIDAnnoDomList[[wo]])[1] != 0){
      write_csv(grPIDAnnoDomList[[wo]], path=paste0(OUTDIR, "/per_domain_gene/per_domain_gene.", nmo1, ".pfam.csv"))
    }
})

##per_gene annotation
grPIDAnnoGeneList <- lapply(grPIDList, function(f){
  lapply(unlist(names(f[[2]])), function(ff){
    df <- as.data.frame(f[[2]][[ff]])
    if(ff != "pfam"){
      ff <- paste0("pfam_", ff)
    }
    if(ff %in% colnames(pfam_tib)){
      if(length(df)>1){
      colnames(df) <- c(ff, "variants_per_bp")
      left_join(df, pfam_tib, by=ff) %>% dplyr::select(!!ff, pfam, pfam_de, variants_per_bp) %>% dplyr::arrange(desc(mutations_per_bp)) %>%
      distinct()
    }}
  })
})

dir.create(paste0(OUTDIR,"/per_domain"), showWarnings=F, recursive = T)
lapply(seq_along(grPIDAnnoGeneList), function(wo){
  nmo1 <- GENELISTNAMES[[wo]]
  lapply(seq_along(grPIDAnnoGeneList[[wo]]), function(po){
    nmo2 <- names(grPIDList[[wo]][[2]])[po]
    if(length(grPIDAnnoGeneList[[wo]][[po]]) != 0){
        write_csv(unique(grPIDAnnoGeneList[[wo]][[po]]), path=paste0(OUTDIR, "/per_domain/per_domain.", nmo1, ".", nmo2, ".csv"))
    }
  })
})
