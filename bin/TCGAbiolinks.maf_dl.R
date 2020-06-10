#! /usr/bin/env R
library(tidyverse)
library(TCGAbiolinks)

tcga_ids <- as_tibble(getGDCprojects()) %>%
            dplyr::select(id) %>%
            dplyr::filter(grepl("TCGA", id)) %>%
            dplyr::mutate(id = gsub("TCGA-", "", id)) %>%
            unlist() %>% unname()

##each time run, some disease types were missed, so rerun twice to get all 32
GDCquery_Maf(tcga_ids, pipelines="mutect2")
GDCquery_Maf(tcga_ids, pipelines="mutect2")
GDCquery_Maf(tcga_ids, pipelines="mutect2")
