knitr::opts_chunk$set(echo = TRUE)
library(dplyr)
source('functions.R')

col_direction <- c("Chr"="#0c8cf4","Soluble"="#8ecc7a","NC"="#cdcdcd")
col_enrichment <- c("IM"="#4198ef","IO"="#94c5ed","MO"="#88b6ed","NC"="#cdcdcd","Chr"="#0c8cf4") # BLUE 
#col_enrichment <- c("IM"="#dc585e","IO"="#ebb4b2","MO"="#eacfb5") #orangne
col_fraction <- c("S"="#8ecc7a","P"="#0c8cf4")

results_A549_I <- data.table::fread("result/deseq_merged/DESeqResult_A549_I.gz")
results_A549_M <- data.table::fread("result/deseq_merged/DESeqResult_A549_M.gz")
results_HeLa_I <- data.table::fread("result/deseq_merged/DESeqResult_HeLa_I.gz")
results_HeLa_M <- data.table::fread("result/deseq_merged/DESeqResult_HeLa_M.gz")

cheRNA_list <- readRDS("result/cheRNA_list.rds")
Com_Specifc_CheRNA <- readRDS(file =  "result/Com_Specifc_CheRNA.rds")

gene.list.table <- data.table::fread("resource/merged.expressed.hg38.bed",sep = "\t")


A549_cheRNA <- list(IO = setdiff(cheRNA_list$A549_I,cheRNA_list$A549_M) %>% unique(),
                    MO = setdiff(cheRNA_list$A549_M,cheRNA_list$A549_I) %>% unique(),
                    IM = intersect(cheRNA_list$A549_I,cheRNA_list$A549_M) %>% unique()
)

HeLa_cheRNA <- list(IO = setdiff(cheRNA_list$HelaS3_I,cheRNA_list$HelaS3_M) %>% unique(),
                    MO = setdiff(cheRNA_list$HelaS3_M,cheRNA_list$HelaS3_I) %>% unique(),
                    IM = intersect(cheRNA_list$HelaS3_I,cheRNA_list$HelaS3_M) %>% unique()
)


deresults <- list(A549_I=results_A549_I,
                  A549_M=results_A549_M,
                  HeLa_I=results_HeLa_I,
                  HeLa_M=results_HeLa_M)


A549_cheRNA.bed <- data.table::fread("result/A549.cheRNA.bed",sep = "\t")
HeLa_cheRNA.bed <- data.table::fread("result/HeLa.cheRNA.bed",sep = "\t")
