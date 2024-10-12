##=======================================##
##  miRNA Differential Expression Analysis 
##=======================================##
## install required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("multiMiR")

## Loading required packages
library(edgeR)    # for DEG analysis
library(multiMiR) # for miRNA extraction

## Load data (demo data is within 'multiMiR' package)
counts_file  <- system.file("extdata", "counts_table.Rds", package = "multiMiR")
strains_file <- system.file("extdata", "strains_factor.Rds", package = "multiMiR")
counts_table   <- readRDS(counts_file)
strains_factor <- readRDS(strains_file)

## Standard edegR differential expression analysis
design <- model.matrix(~ strains_factor)

## Using trended dispersions
deg <- DGEList(counts = counts_table)
deg <- calcNormFactors(deg)
deg$samples$strains <- strains_factor
deg <- estimateGLMCommonDisp(deg, design)
deg <- estimateGLMTrendedDisp(deg, design)
deg <- estimateGLMTagwiseDisp(deg, design)

## Fit GLM model for strain/batch effect
fit <- glmFit(deg, design)
lrt <- glmLRT(fit)

## Table of unadjusted p-values (PValue) and FDR values
p_val_DE_edgeR <- topTags(lrt, adjust.method = 'BH', n = Inf)

## Getting top differentially expressed miRNA's
top_miRNAs <- rownames(p_val_DE_edgeR$table)[1:10]

## Plug miRNA's into multiMiR and getting validated targets
multimir_results <- get_multimir(org     = 'mmu',
                                 mirna   = top_miRNAs,
                                 table   = 'validated',
                                 summary = TRUE)
multimir_results


##======================================================##
## Retrieve all validated target genes of a given miRNA 
##======================================================##
# search validated interactions in human
miRNA_gene <- get_multimir(mirna = 'hsa-miR-18a-3p', summary = TRUE)
names(miRNA_hsa)

# Check which types of associations were returned
table(miRNA_gene@data$type)

# Detailed information of the validated miRNA-target interaction
head(miRNA_gene@data)

# Which interactions are supported by Luciferase assay?
miRNA_gene@data[grep("Luciferase", example1@data[, "experiment"]), ]
# filter miRNA for that target validated by Luciferase assay maximum time (here, KRAS)
miRNA_gene@summary[miRNA_hsa@summary[,"target_symbol"] == "KRAS",]

##=======================================================================##
## DRUG ## Retrieve miRNA-target interactions associated with drug/disease
##=======================================================================##
# find which miRNAs and their target genes are associated with Cisplatin, a chemotherapy drug used in several cancers.
miRNA_drug <- get_multimir(disease.drug = 'cisplatin', table = 'disease.drug')
head(miRNA_drug@data)


##==========================================##
## Select miRNAs predicted to target a gene
##==========================================##
miRNA_gene_org <- get_multimir(org     = "hsa",
                               target  = "FABP4",
                               table   = "predicted",
                               summary = TRUE,
                               predicted.cutoff = 35,
                               predicted.cutoff.type = "p",
                               predicted.site = "all")
table(miRNA_gene_org@data$type)

## examine how many predictions each of the databases has.
apply(miRNA_gene_org@summary[, 6:13], 2, function(x) sum(x > 0))


##============================================================##
## Select miRNA(s) predicted to target most genes of interest
##============================================================##
miRNA_allgenes <- get_multimir(org = 'hsa',
                               target = c('AKT2', 'CERS6', 'S1PR3', 'SULF2'),
                               table = 'predicted',
                               summary = TRUE,
                               predicted.cutoff.type = 'n',
                               predicted.cutoff = 500000) # searching the top 500,000 predictions in each external database

## count the number of target genes for each miRNA.
miRNA_allgenes.counts <- addmargins(table(miRNA_allgenes@summary[, 2:3]))
miRNA_allgenes.counts <- miRNA_allgenes.counts[-nrow(miRNA_allgenes.counts), ]
miRNA_allgenes.counts <- miRNA_allgenes.counts[order(miRNA_allgenes.counts[, 4], decreasing = TRUE), ] # 4 indicates miRNA targeting atleast 4 genes, but missing in 1 gene
head(miRNA_allgenes.counts)
