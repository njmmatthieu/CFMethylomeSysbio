library(tibble)
library(dplyr)

# 1. Loading the CF network

CF_PPI_network.lcc.node_type.interactions <- 
  read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/Meta-analysis article/CFnetwork/data/kegg_diff_pathways_network/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_interactions_df.txt",
             sep = "\t",
             header = T,
             check.names = F)

CF_PPI_network.lcc.node_type.nodes <- 
  read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/Meta-analysis article/CFnetwork/data/kegg_diff_pathways_network/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_nodes_df.txt",
             sep = "\t",
             header = T,
             check.names = F)

# 2. Differentially methylated CpG sites between CF and NCF dataframe

dm_CF_NCF_sup_table3_1267CpG <- 
  read.table(file = "/Users/matthieu/ownCloud/Suite/CFMethylomeSysbio/dm_CpG_CF_NCF_Magalhaes_2018_sup_table_3.tsv",
             skip = 1,
             quote = "\"",
             sep = "\t",
             header = T,
             check.names = F)

dm_CF_NCF_sup_table3_1267CpG <- 
  dm_CF_NCF_sup_table3_1267CpG[order(dm_CF_NCF_sup_table3_1267CpG$`p-value`),]
dm_CF_NCF_sup_table3_1267CpG$CpG_rank <- 1:dim(dm_CF_NCF_sup_table3_1267CpG)[1]

## 2.1 Keeping only CpG sites in the body of the genes

dm_CF_NCF_sup_table3_1267CpG_onlyGenes <- 
  dm_CF_NCF_sup_table3_1267CpG[which(dm_CF_NCF_sup_table3_1267CpG$`Genomic location`=="Body"),]

dm_CF_NCF_sup_table3_1267CpG_onlyGenes$Gene_rank <- 1:dim(dm_CF_NCF_sup_table3_1267CpG_onlyGenes)[1]

## 2.2 Searching for dm genes in the CF network

dm_CF_NCF_sup_table3_1267CpG_Genes.list <- unique(dm_CF_NCF_sup_table3_1267CpG_onlyGenes$Gene)

dm_CF_NCF_Genes_in_CF_network.list <- 
  dm_CF_NCF_sup_table3_1267CpG_Genes.list[which(dm_CF_NCF_sup_table3_1267CpG_Genes.list %in% CF_PPI_network.lcc.node_type.nodes$Symbol)]

dm_CF_NCF_Genes_in_CF_network.CF_network_df <- 
  CF_PPI_network.lcc.node_type.nodes[which(CF_PPI_network.lcc.node_type.nodes$Symbol %in% dm_CF_NCF_Genes_in_CF_network.list),]

dm_CF_NCF_Genes_in_CF_network.dm_CF_NCF_df <-
  dm_CF_NCF_sup_table3_1267CpG[which(dm_CF_NCF_sup_table3_1267CpG$Gene %in% dm_CF_NCF_Genes_in_CF_network.list),]

# 3. Differentially methylated CpG sites between severe and mild dataframe

dm_severe_mild_sup_table4_189CpG <- 
  read.table(file = "/Users/matthieu/ownCloud/Suite/CFMethylomeSysbio/dm_CpG_CF_severe_mild_Magalhaes_2018_sup_table_4.tsv",
             skip = 1,
             quote = "\"",
             sep = "\t",
             header = T,
             check.names = F)

dm_severe_mild_sup_table4_189CpG <- 
  dm_severe_mild_sup_table4_189CpG[order(dm_severe_mild_sup_table4_189CpG$`p-value`),]
dm_severe_mild_sup_table4_189CpG$CpG_rank <- 1:dim(dm_severe_mild_sup_table4_189CpG)[1]

## 3.1 Keeping only CpG sites in the body of the genes

dm_severe_mild_sup_table4_189CpG_onlyGenes <- 
  dm_severe_mild_sup_table4_189CpG[which(dm_severe_mild_sup_table4_189CpG$`Genomic location`=="Body"),]

dm_severe_mild_sup_table4_189CpG_onlyGenes$Gene_rank <- 1:dim(dm_severe_mild_sup_table4_189CpG_onlyGenes)[1]

# 3.2 Searching for dm genes in the CF network

dm_severe_mild_sup_table4_189CpG_Genes.list <- unique(dm_severe_mild_sup_table4_189CpG_onlyGenes$Gene)

dm_CF_NCF_Genes_in_CF_network.list <- 
  dm_severe_mild_sup_table4_189CpG_Genes.list[which(dm_severe_mild_sup_table4_189CpG_Genes.list %in% CF_PPI_network.lcc.node_type.nodes$Symbol)]

dm_CF_NCF_Genes_in_CF_network.CF_network_df <- 
  CF_PPI_network.lcc.node_type.nodes[which(CF_PPI_network.lcc.node_type.nodes$Symbol %in% dm_CF_NCF_Genes_in_CF_network.list),]

dm_CF_NCF_Genes_in_CF_network.dm_CF_NCF_df <-
  dm_severe_mild_sup_table4_189CpG[which(dm_severe_mild_sup_table4_189CpG$Gene %in% dm_CF_NCF_Genes_in_CF_network.list),]
