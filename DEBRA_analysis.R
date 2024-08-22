library(tidyr)
library(tidyverse)
library(dplyr)
library(DEBRA)

#Load counts matrix
counts <- as.matrix(read.table("Data/Counts/ssb6 counts.txt", header = T))

rif <- ssb6[, c("No_drug", "No_drug_1", "No_drug_2", "No_drug_3", "rif", "rif_1", "rif_2", "rif_3")] 

threshold <- 50
# Set the minimum number of columns with values below the threshold
min_below_threshold <- 4

# Filter rows based on the specified conditions
rif1 <- rif[(rowSums(rif > threshold) >= min_below_threshold), ]

control_names = c("No_drug", "No_drug_1", "No_drug_2", "No_drug_3") #Define the control group
condition_names = c("rif", "rif_1", "rif_2", "rif_3") #Define the case group
x <- DEBRA(rif1, control_names=control_names, condition_names=condition_names, 
           method=c("DESeq2(Wald)"), 
           modified = F, 
           shrunkLFC = T,
           trended = T)
results=resultsDRB(x)
Rresults <- results %>% tibble::rownames_to_column("strains") 
