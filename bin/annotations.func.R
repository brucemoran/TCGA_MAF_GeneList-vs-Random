#! /usr/bin/R

##annotations
libs <- c("ensembldb", "EnsDb.Hsapiens.v86", "ensemblVEP", "customProDB", "tidyverse", "biomaRt", "magrittr", "reshape2", "ggridges", "PFAM.db")
libsLoaded <- lapply(libs,function(l){suppressMessages(library(l, character.only = TRUE))})

strSplitFun <- function(input,sepn){
  lapply(input,function(f){strsplit(f,sepn)[[1]]})
}
strSplitVec <- function(inVec,sepn){
  sapply(seq_along(inVec),function(f){strsplit(inVec[f],sepn)[[1]]})
}

##function to annotate GRanges from biomaRt
GRanges_anno_func <- function(GRANNOFILE){
  if(!file.exists(GRANNOFILE)){
    print("Generating BiomaRt annotation data...")
    biomartCacheClear()
    annoMart <- useMart(biomart="ensembl",
                        dataset="hsapiens_gene_ensembl",
                        host="www.ensembl.org")
    atts <-  c('chromosome_name','start_position','end_position','strand','external_gene_name','ensembl_gene_id','entrezgene_id')
    biomartAll <- getBM(attributes=atts, mart = annoMart)
    biomartAll <- as_tibble(biomartAll)
    colnames(biomartAll) <- c("seqnames",
                                      "start",
                                      "end",
                                      "strand",
                                      "external_gene_name",
                                      "ensembl_gene_id",
                                      "entrezgene_id")

    biomartAll %<>% dplyr::mutate(strand = unlist(lapply(strand, function(f){
                            if(f == 1){return("+")}
                            if(f == -1){return("-")}
                          })))
  save(biomartAll, file=GRANNOFILE)
  }
  if(file.exists(GRANNOFILE)){
    print(paste0("Found: ", GRANNOFILE, " ...loading..."))
  }
}
GRanges_anno_func("biomaRt_all.extensent.RData")

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
