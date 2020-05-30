#! /usr/bin/env R
library(tidyverse)
library(TCGAbiolinks)

tcga_ids <- as_tibble(getGDCprojects()) %>%
            dplyr::select(id) %>%
            dplyr::filter(grepl("TCGA", id)) %>%
            dplyr::mutate(id = gsub("TCGA-", "", id)) %>%
            unlist() %>% unname()
GDCquery_Maf(tcga_ids, pipelines="mutect2")
