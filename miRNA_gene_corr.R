##============================##
##  miRNA - Gene Interaction 
##============================##
## install required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("multiMiR")

## Loading required packages
library(multiMiR) # for miRNA extraction

## Load data (demo data in 'multiMiR' package for bladder cancer)
load(url("http://multimir.org/bladder.rda"))

##===============================================##
## Retrieve interactions between miRNAs and genes
##===============================================##
## Search all tables & top 10% predictions
## Hypothesis:: interactions between these miRNAs and genes whose expression changed at opposite directions may play a role in disease metastasis.
miRNA_gene <- get_multimir(org = "hsa",
                         mirna = DE.miRNA.up, # DE.miRNA.up contains 9 up-regulated miRNAs
                         target = DE.entrez.dn, # DE.entrez.dn has 47 down-regulated genes
                         table = "all",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 10,
                         use.tibble = TRUE)

table(miRNA_gene@data$type)

## validated interactions
result <- select(miRNA_gene, keys = "validated", keytype = "type", 
                 columns = columns(miRNA_gene))

## non-redundant list
unique_pairs <- result[!duplicated(result[, 
                            c("mature_mirna_id", "target_entrez")]), ]
result
##O/P:: 85 unique miRNA-gene pairs that have been validated


##========================================##
## miRNAs are associated with the disease
##========================================##
mykeytype <- "disease_drug"

mykeys <- keys(miRNA_gene, keytype = mykeytype)
mykeys <- mykeys[grep("bladder", mykeys, ignore.case = TRUE)]

result <- select(miRNA_gene, keytype = "disease_drug", keys = mykeys,
                 columns = columns(miRNA_gene))
result
##O/P:: 2 miRNAs are associated with bladder cancer in miR2Disease and PhenomiR.

predicted <- select(miRNA_gene, keytype = "type", keys = "predicted", 
                    columns = columns(miRNA_gene))
length(unique(predicted$mature_mirna_id))

length(unique(predicted$target_entrez))

unique.pairs <- 
  unique(data.frame(miRNA.ID = as.character(predicted$mature_mirna_id),
                    target.Entrez = as.character(predicted$target_entrez)))
nrow(unique.pairs)
head(unique.pairs)
## Results from each of the predicted databases are ordered by their scores from best to worst.
miRNA_gene.split <- split(predicted, predicted$database)

## get miRNA-gene and interaction score
mykeys <- keys(miRNA_gene_org)[1:4]  ## miRNA_gene_org from miRNA_DEA.R
head(select(miRNA_gene_org, keys = mykeys, 
            columns = c("database", "target_entrez", "score")))