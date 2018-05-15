
#1
setwd('/Users/martin/Desktop/Proteomics_VT18/0__PrepareData/1_Parsed/')
proteins <- read.delim('donors.body.full.wide.csv', sep = ',', header = TRUE, row.names = 1)
raw_proteins <- proteins

proteins2 <- raw_proteins

#proteins <- na.omit(proteins)

cor_mat_prot_logged <- cor(log2(proteins + 0.5))

heatmaply(cor_mat_prot_logged,
          main = "Sample-sample correlation, log2 counts",
          plot_method = 'plotly')

#2
center_raw_mat <-
  cor_mat_prot_logged - apply(cor_mat_prot_logged, 1, median)

raw_max <- max(abs(center_raw_mat), na.rm = TRUE)

raw_limits <- c(-raw_max, raw_max)

heatmaply(
  center_raw_mat,
  fontsize_col = 7.5,
  col = cool_warm(100),
  main = 'Centred log2 read counts',
  limits = raw_limits,
  plot_method = 'plotly'
)

# TRANSPOSED
# heatmaply(t(center_raw_mat),
#           fontsize_col = 7.5,
#           col = cool_warm(100),
#           main = 'Centred log2 read counts',
#           limits = raw_limits,
#           plot_method = 'plotly')

heatmaply_cor(cor(center_raw_mat),
              main = 'Sample-sample correlation based on centred, log2 counts',
              plot_method = 'plotly')

voomed_proteins <- as.matrix(voom(proteins))

center_voom_mat <-
  voomed_proteins - apply(voomed_proteins, 1, median)

voom_max <- max(abs(center_voom_mat))
voom_limits <- c(-voom_max, voom_max)

heatmaply(center_voom_mat,
          fontsize_col = 7.5,
          col = cool_warm(50),
          limits = voom_limits,
          main = 'Normalised, centred log2 CPM',
          plot_method = 'plotly')

## Non centered heatmaps

log_raw_mat <- log2(raw_proteins + 0.5)

heatmaply(log_raw_mat,
          showticklabels = c(TRUE, TRUE),
          fontsize_col = 7.5,
          col = gplots::bluered(50),
          main = "Pre-normalisation log2 counts",
          plot_method = 'plotly')

# TRANSPOSED
# heatmaply(t(log_raw_mat),
#           showticklabels = c(TRUE, TRUE),
#           fontsize_col = 7.5,
#           col = gplots::bluered(50),
#           main = "Pre-normalisation log2 counts",
#           plot_method = 'plotly')

heatmaply_cor(cor(log_raw_mat),
              main = 'Sample-sample correlation based on log2-transformed protein abundances',
              plot_method = 'plotly')

heatmaply(voomed_proteins,
        showticklabels = c(TRUE, TRUE),
        fontsize_col = 7.5,
        col = gplots::bluered(50),
        main = 'Normalised log2 CPM',
        plot_method = 'plotly')

heatmaply_cor(cor(voomed_proteins),
              showticklabels = c(FALSE, FALSE),
              main = 'Sample-sample correlation based on normalised protein abundances',
              plot_method = 'plotly')

### Kotsyfakis et al.

library(dendextend)

row.names(proteins) <- NULL
y <- as.matrix(rpkm50)
x <- as.matrix(proteins)
x <- x[which(rowSums(x) > 0),]
rc <- rainbow(nrow(x), start=0, end=.3)
cc <- rainbow(ncol(x), start=0, end=.3)


data("rpkm50")
head(rpkm50)

hr <- hclust(as.dist(1-cor(t(x), method="spearman")), method="complete")
hc <- hclust(as.dist(1-cor(x, method="spearman")), method="complete")

library(gplots)
heatmap.2(
  x,
  col = bluered(75),
  Colv = as.dendrogram(hc),
  Rowv = as.dendrogram(hr),
  scale = "row",
  key = T,
  keysize = 1.5,
  density.info = "none",
  trace = "none",
  cexCol = 0.9,
  cexRow = 0.9,
  labRow = NA,
  dendrogram = "both"
)

heatmaply(
  x,
  col = bluered(75),
  Colv = as.dendrogram(hc),
  Rowv = as.dendrogram(hr),
  scale = "row",
  key = T,
  keysize = 1.5,
  density.info = "none",
  trace = "none",
  cexCol = 0.9,
  cexRow = 0.9,
  labRow = NA,
  dendrogram = "both"
)


# Decisions for clustering (highlight_branches helps identify
# the topology of the dendrogram, it colors each branch based 
# on its height):

DATA <- raw_proteins
d <- dist(sqrt(DATA))

library(dendextend)
dend_expend(d)[[3]]


dend_row <- d %>% hclust(method = "single") %>% as.dendrogram 
dend_row %>% highlight_branches %>% plot

dend_k <- find_k(dend_row)
plot(dend_k)

# logged
DATA <- raw_proteins+0.5
d <- dist(log2(DATA))
dend_expend(d)[[3]]

dend_row <- d %>% hclust(method = "average") %>% as.dendrogram 
dend_row %>% highlight_branches %>% plot

dend_k <- find_k(dend_row)
plot(dend_k)

Rowv <- dend_row %>% color_branches(k = 2)
heatmap.2(
  as.matrix(log2(DATA)),
  Colv = NULL,
  Rowv = Rowv,
  trace = "none",
  col = viridis(200),
  margins = c(3, 9)
)

Rowv <-
  dend_row %>% color_branches(k = 2) %>% seriate_dendrogram(x = d)
heatmap.2(as.matrix(log2(DATA)),
  Colv = NULL,
  Rowv = Rowv,
  trace = "none",
  col = viridis(200),
  margins = c(3, 9)
)

heatmaply(
  as.matrix(log2(DATA)),
  Colv = NULL,
  hclust_method = "average",
  fontsize_row = 8,
  fontsize_col = 6,
  k_row = NA,
  margins = c(60, 170, 70, 40),
  xlab = "Sample",
  ylab = "Protein",
  main = "Heatmap visualization of protein abundances",
  plot_method = "plotly", row_dend_left = TRUE
)


dend <- d %>% find_dend %>% seriate_dendrogram(., d)

heatmaply(
  as.matrix(log2(DATA)),
  limits = c(0, 8),
  Rowv = dend,
  margins = c(85, 40),
  grid_gap = 0.2,
  k_row = 21
)

### and ends here


# Omit NAs from dataset
# Proteins that is not present in all (current) samples are removed
# Make less strict? >= 1/3 (python)
proteins_noNA <- na.omit(proteins)

heatmaply(proteins_noNA)
proteins_noNAlog <- log10(proteins_noNA)
heatmaply(proteins_noNAlog)




heatmaply(proteins_noNAlog, k_row = 1, k_col = 21)


library('biomaRt')

## Annotate using biomaRt
mart <- useDataset('hsapiens_gene_ensembl', useMart('ensembl'))
genes <- rownames(data)
symbols <- getBM(filters= 'ensembl_gene_id', 
                 attributes= c('ensembl_gene_id','hgnc_symbol'), 
                 values = genes,
                 mart= mart)

# Installation

source("https://bioconductor.org/biocLite.R")
biocLite("SummarizedExperiment")     ## R version 2.15 or later
biocLite("biomaRt") 
biocLite("limma")
biocLite("voom")
biocLite("TCGAbiolinks")
install.packages("tidyverse")
install.packages('glmnet')


devtools::install_github("talgalili/heatmaplyExamples", build_vignettes = TRUE)

library(heatmaply)
library(heatmaplyExamples)



# PART TWO
# https://cdn.rawgit.com/talgalili/heatmaplyExamples/master/inst/doc/biological_data_2.html

###

library(heatmaply)
library(heatmaplyExamples)

pam50_genes <- intersect(pam50_genes, rownames(raw_expression))
raw_pam50_expression <- raw_expression[pam50_genes, ]
voomed_pam50_expression <- voomed_expression[pam50_genes, ]

center_raw_mat <- cor_mat_raw_logged - 
  apply(cor_mat_raw_logged, 1, median)

raw_max <- max(abs(center_raw_mat), na.rm=TRUE)
raw_limits <- c(-raw_max, raw_max)
