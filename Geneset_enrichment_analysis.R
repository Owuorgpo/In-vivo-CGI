library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(ggplot2)

#Load list of genes
gse <- read_csv("Data/GSE/Final DEBRA for KEGG1.csv")

#Load custom KEGG functional database
dt <- read_csv("mtu KEGG with all genes.csv")
save(dt, file = "mtb_kegg.RData")
# Load the custom KEGG dataset
load("mtb_kegg.RData")

# Function to perform custom KEGG enrichment analysis
custom_enrichKEGG <- function(genes, custom_kegg) {
  enriched <- enricher(
    genes, 
    TERM2GENE = custom_kegg[, c("PathwayID", "GeneID")], 
    TERM2NAME = custom_kegg[, c("PathwayID", "PathwayName")],
    pvalueCutoff = 0.05
  )
  return(enriched)
}

library(clusterProfiler)
gse1 <- gse %>% #filter(sspB == 6 | sspB == 10) %>% 
  dplyr::select(strains, LFC_PZA, padj_PZA, sspB) %>% 
  filter(!grepl("H37RV_", strains)) %>% 
  filter(LFC_PZA < -1 & padj_PZA < 0.01) %>% group_by(strains) %>% 
  dplyr::select(strains) %>% unique()

gene_list <- gse1$strains

# Perform enrichment analysis
enrichment_result <- custom_enrichKEGG(gene_list, custom_kegg = dt)

data <- as.data.frame(enrichment_result@result)
data[c('Genes', 'Hits')] <- str_split_fixed(data$GeneRatio, '/', 2)
data[c('Background', 'nRv')] <- str_split_fixed(data$BgRatio, '/', 2)
data1 <- data %>% dplyr::select(ID, Description, geneID, Count,Hits, Background,nRv) %>% 
  mutate(nRv = 3906) %>% 
  mutate(Background = as.numeric(Background)) %>% 
  mutate(Hits = as.numeric(Hits)) %>% 
  mutate(Expected = (Background/3906)*Hits) %>% 
  mutate(Enrichment = Count/Expected)

# Create tables
## This makes 2x2 tables for each gene
fet <-as.data.frame(fisher.tables <-apply(data1, 1, 
                                          function(x) {
                                            tbl <- matrix(as.numeric(x[4:7]), ncol=2, byrow=T)
                                            fisher.test(tbl, alternative="two.sided")$p.value
                                          })) %>% 
  dplyr::rename(pval = `fisher.tables <- apply(data1, 1, function(x) {     tbl <- matrix(as.numeric(x[4:7]), ncol = 2, byrow = T)     fisher.test(tbl, alternative = "two.sided")$p.value })`)


fet1 <- fet %>%   mutate(p.adj = p.adjust(pval, "fdr"))
ID <- rownames(fet1)
rownames(fet1) <- NULL
fet2 <- cbind(ID,fet1)

data2 <- data1 %>% dplyr::select(-Hits, -Background, -nRv, -Expected) %>% left_join(fet2) %>% 
  mutate(p.adj = round(p.adj, 6)) %>% 
  filter(p.adj < 0.05) %>% 
  mutate(Count  = as.numeric(Count)) %>% 
  dplyr::rename(FDR = p.adj)

data_ranked <- data2[order(data2$Enrichment, data2$Description), ]

k <- ggplot(data_ranked, aes(x = reorder(Description, Enrichment), y = Enrichment, fill = FDR)) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = Count), color = "black", size = 3, fill = "white", position = position_stack(vjust = 0.5)) +
  coord_flip() +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(legend.position="right") +
  labs(x = "", y = "Fold enrichment", fill = "FDR") +
  theme(axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12))

dev.off()
