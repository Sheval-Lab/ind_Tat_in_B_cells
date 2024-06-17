## Load libraries --------------------------------------------------------------
library(tidyverse)
library(DESeq2)
source("scripts/functions/extract_results.R")

input_dir <- "results/tables"
output_dir <- "results/tables/01_DGEA"
fig_dir <- "results/figures/01_DGEA"


## Load gene annotation --------------------------------------------------------
gene_annotation_file <- "data/metadata/GRCh38.p10_ALL.annotation.IDs.txt"
gene_annotation <- read_tsv(gene_annotation_file, 
                            col_names = c("geneID", "gene_name", "gene_type"))


## Prepare data for DGEA -------------------------------------------------------
### Load counts data 
counts <- read_tsv(file.path(input_dir, "counts_raw.tsv"))

### Prepare matrix of counts
counts_mtx <- counts %>% 
  select_if(is.numeric) %>% 
  as.matrix()
rownames(counts_mtx) <- counts$geneID

### Prepare design dataframe
design <- data.frame(
  sample = str_remove(colnames(counts_mtx), ".\\d$"),
  row.names = colnames(counts_mtx))

### Create DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = counts_mtx, 
  colData = design, 
  design = ~sample)

### Pre-filter low count genes
### Keep genes with 10+ reads in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

### Calculate size factors and save normalized counts
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
dput(sizeFactors(dds), file.path(input_dir, "size_factors.txt"))

counts_norm <- counts(dds, normalized=TRUE) %>% 
  as.data.frame() %>% 
  mutate(geneID = row.names(.)) %>% 
  select(geneID, everything())

write_tsv(counts_norm, file.path(input_dir, "counts_norm.tsv"))


## Run DGEA --------------------------------------------------------------------
dds <- DESeq(dds)

saveRDS(dds, file.path(output_dir, "dds.rds"))


## Extract DGEA results --------------------------------------------------------
### Create contrasts list
contrasts_list <- tibble(
  numerator = design$sample,
  denominator = design$sample) %>% 
  tidyr::expand(numerator, denominator)

contrasts_list <- contrasts_list[c(13, 10, 15, 5),]


### Extract results
log2FC_threshold <- log2(1.5)

#### Tat.stable/Tat.0h vs control
walk2(contrasts_list$numerator, contrasts_list$denominator, 
      safely(extract_results), dds = dds, gene_annotation = gene_annotation, log2FC_threshold, output_dir)


#### Tat.16h vs Tat.0h
dds$sample <- relevel(dds$sample, ref = "Tat.0h")
dds <- DESeq(dds)

walk2(contrasts_list$numerator, contrasts_list$denominator, 
      safely(extract_results), dds = dds, gene_annotation = gene_annotation, log2FC_threshold, output_dir)


#### Tat.stable vs Tat.16h
dds$sample <- relevel(dds$sample, ref = "Tat.16h")
dds <- DESeq(dds)

walk2(contrasts_list$numerator, contrasts_list$denominator, 
      safely(extract_results), dds = dds, gene_annotation = gene_annotation, log2FC_threshold, output_dir)

