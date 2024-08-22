library(tidyr)
library(tidyverse)
library(dplyr)
library(DEBRA)
library(DESeq2)

#Load counts matrix
counts <- read_tsv("Data/Counts/ssb6 counts.txt", col_names = T) %>% 
  remove_rownames %>% column_to_rownames(var="Gene")
# PCA plot
pca <- counts[,12:30]
threshold <- 50
# Set the minimum number of columns with values below the threshold
min_below_threshold <- 4

# Filter rows based on the specified conditions
pca1 <- pca[(rowSums(pca > threshold) >= min_below_threshold), ]


condition <-factor(c(rep("No_Drug", 4), rep("INH", 4), rep("RIF", 4), rep("PZA", 3), rep("HRZE", 4)))
coldata<-data.frame(row.names = colnames(pca), condition)
ddsfile<-DESeqDataSetFromMatrix(countData = pca, 
                                colData = coldata,
                                design = ~1)

colData(ddsfile)$condition <- relevel(colData(ddsfile)$condition, "No_Drug")
ddsfile <- ddsfile[ rowSums(counts(ddsfile)) >1, ]
dds_rld <- rlog(ddsfile, blind = FALSE, fitType = "local")

pca6 <- plotPCA(dds_rld, intgroup = c("condition")) + geom_point(size = 4) + coord_fixed(ratio = 1.5) + scale_color_manual(values = c("blue", "orange", "red", "darkviolet", "green")) + 
  theme(plot.title = element_text(face = "bold", colour = "black", size = 16)) +
  theme(axis.title.x = element_text(face = "bold", size = 14)) +
  theme(axis.title.y = element_text(face = "bold", size = 14)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  
  theme(legend.position = "bottom")

#DEBRA anaylsis

rif <- counts[, c("Nodrug_1", "Nodrug_2", "Nodrug_3", "Nodrug_4", "RIF_1", "RIF_2", "RIF_3", "RIF_4")] 

threshold <- 50
# Set the minimum number of columns with values below the threshold
min_below_threshold <- 4

# Filter rows based on the specified conditions
rif1 <- rif[(rowSums(rif > threshold) >= min_below_threshold), ]

control_names = c("Nodrug_1", "Nodrug_2", "Nodrug_3", "Nodrug_4") #Define the control group
condition_names = c("RIF_1", "RIF_2", "RIF_3", "RIF_4") #Define the case group
x <- DEBRA(rif1, control_names=control_names, condition_names=condition_names, 
           method=c("DESeq2(Wald)"), 
           modified = F, 
           shrunkLFC = T,
           trended = T)
results=resultsDRB(x)
Rresults <- results %>% tibble::rownames_to_column("strains") 
