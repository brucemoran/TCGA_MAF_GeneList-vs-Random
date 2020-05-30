#! /usr/bin/env R

##annotate protein domains
##determine if certain domains are more highly mutated: hotspots
libs <- c("ensembldb", "EnsDb.Hsapiens.v86", "plyranges", "PFAM.db", "tidyverse", "biomaRt", "ggridges")
libsLoaded <- lapply(libs,function(l){suppressMessages(library(l, character.only = TRUE))})

##Pfam annotation
##annotation of all transcripts with P(rotein)ID(entifier)
EnsDb.Hsapiens.v86_func <- function(ANNOFILE){
  if(!file.exists(ANNOFILE)){
    print(paste0("Did not find: ", ANNOFILE, " ...generating..."))
    edb <- EnsDb.Hsapiens.v86
    txs <- transcripts(edb,
               columns = c("protein_id", "protein_domain_id", "protein_domain_source", "prot_dom_start", "prot_dom_end", "uniprot_id", "tx_biotype", "gene_id"))
    txs_pid <- txs[!is.na(txs$protein_domain_id),]

    pfam_de <- PFAMDE
    pfam_prosite <- PFAMPROSITEPROFILE
    pfam_smart <- PFAMSMART
    pfam_prints <- PFAMPRINTS
    pfam_vec <- c("pfam_de", "pfam_prosite", "pfam_smart", "pfam_prints")
    pfam_list <- as.list(pfam_vec)

    ##map keys, tibble
    pfam_tibs <- lapply(seq_along(pfam_vec), function(f){
      mapped_keys <- mappedkeys(get(pfam_vec[f]))
      mapped_keys <- as.list(get(pfam_vec[f])[mapped_keys])
      do.call(rbind, lapply(seq_along(mapped_keys), function(ff){
        tibble(pfam=names(mapped_keys)[ff], pp=unname(unlist(mapped_keys[[ff]])))
      })) %>% rename(!!pfam_vec[f]:=pp)
    })

    ##join
    pfam_tib <- pfam_tibs[[1]]
    for(ft in 2:length(pfam_tibs)){
      pfam_tib <- left_join(pfam_tib, pfam_tibs[[ft]])
    }

    ##gene_annotation
    biomartCacheClear()
    annoMart <- useMart(biomart="ensembl",
                        dataset="hsapiens_gene_ensembl")
    annoGenenameEnsEnt <- as_tibble(getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), mart = annoMart))

    save(txs_pid, pfam_tib, annoGenenameEnsEnt, file=ANNOFILE)
  }
  if(file.exists(ANNOFILE)){
    print(paste0("Found: ", ANNOFILE, " ...loading..."))
  }
}

##rename seqnames chr -> no_chr
##join, make a table of the most represented domains
sqlvl_pid_gr_func <- function(gr, txs_pid){
    seqlevelsStyle(gr) <- "Ensembl"
    wntsqlvl <- c(1:22, "X","Y")
    sqlvlsgr <- wntsqlvl[wntsqlvl %in% seqlevels(gr)]
    suppressWarnings(seqinfo(gr, pruning.mode="coarse", new2old=c(1:length(sqlvlsgr))) <- seqinfo(txs_pid)[sqlvlsgr])
    join_overlap_intersect(gr, txs_pid)
}

##domains within a gene may have increased mutatbility
##here this is tested
##multiple annotations sources per mutation exist
##need to check per source as may influence result
per_gene_func <- function(gr_txs_pid_snvs, s){
  if(length(gr_txs_pid_snvs)!=0){
  tb_tmp <- as_tibble(gr_txs_pid_snvs) %>%
            dplyr::mutate(variant = names(gr_txs_pid_snvs)) %>%
            dplyr::filter(protein_domain_source %in% s) %>%
            dplyr::select(-1, -2, -3, -4, -5, -protein_domain_source, -protein_id, -uniprot_id, -tx_id, -tx_biotype) %>%
            dplyr::mutate(domain_width = prot_dom_end - prot_dom_start) %>%
            dplyr::select(-prot_dom_start, -prot_dom_end) %>%
            distinct() %>%
            group_by(protein_domain_id) %>%
            dplyr::mutate(domain_width_mean = round(mean(domain_width))) %>%
            dplyr::select(-domain_width) %>%
            ungroup() %>%
            distinct() %>%
            dplyr::add_count(protein_domain_id) %>%
            dplyr::rename(pg_var_count = n) %>%
            dplyr::select(-variant) %>%
            dplyr::arrange(protein_domain_id) %>%
            dplyr::mutate(pg_var_per_width = pg_var_count / domain_width_mean) %>%
            distinct()
          }
  if(length(gr_txs_pid_snvs)==0){
    tb_tmp <- tibble(protein_domain_id="",
                     gene_id="",
                     domain_width_mean="",
                     pg_var_count="",
                     pg_var_per_width="")[-1,]
  }
  return(tb_tmp)
}

##output per-gene and per-domain results
per_gene_domain_func <- function(gg, GENELISTNAMES, g){
  if(dim(gg)[1] > 0){
    per_gene <- gg %>%
                distinct() %>%
                left_join(., pfam_tib, by=c("protein_domain_id" = "pfam")) %>%
                left_join(., annoGenenameEnsEnt, by=c("gene_id" = "ensembl_gene_id")) %>%
                dplyr::select(protein_domain_id, ensembl_gene_id=gene_id, external_gene_name, pg_var_count, pg_var_per_width, pfam_de) %>%
                dplyr::arrange(desc(pg_var_per_width))

    per_domain <- do.call(rbind, lapply(unique(gg$protein_domain_id), function(d){
                  gg %>% dplyr::filter(protein_domain_id %in% d) %>%
                         dplyr::select(-gene_id) %>%
                         dplyr::mutate(domain_width_mean_total = round(mean(domain_width_mean)),
                                pid_var_count_total = sum(pg_var_count)) %>%
                         dplyr::select(-domain_width_mean) %>%
                         dplyr::distinct() %>%
                         dplyr::mutate(pid_var_per_width = pid_var_count_total / domain_width_mean_total) %>%
                         dplyr::select(protein_domain_id, pid_var_count_total, pid_var_per_width)
                  })) %>% distinct %>%
                          dplyr::arrange(desc(pid_var_per_width)) %>%
                          left_join(., pfam_tib, by=c("protein_domain_id" = "pfam")) %>%
                          dplyr:::select(1, 2, 3, 4)


    lapply(c("per_gene", "per_domain"), function(d){
      dir.create(d, showWarnings = F, recursive = T)
    })
    write_csv(per_gene, path=paste0("per_gene/per_gene.", GENELISTNAMES[g], ".pfam.csv"))
    write_csv(per_domain, path=paste0("per_domain/per_domain.", GENELISTNAMES[g], ".pfam.csv"))

    olist <- list(per_gene, per_domain)
    names(olist) <- c("per_gene", "per_domain")
    return(olist)
  }
}
