library(dplyr)
library(ggplot2)
library(GEOquery)
library(REBayes)
library(AnnotationDbi)
library(hgu133plus2.db)
library(matrixStats)
library(Biobase)
library(affy)
library(hgu133plus2cdf)
library(R.utils)

source('../../EPB main/EPB.R')
source('../../helper functions/Distribution Visualization helper.R')

# CEL Data Download

## 1) Get raw CEL files
#base_dir = "GSE9101_raw"
#dir.create(base_dir, showWarnings = FALSE)
#getGEOSuppFiles("GSE9101", baseDir = base_dir)
#untar(file.path(base_dir, "GSE9101", "GSE9101_RAW.tar"), exdir = file.path(base_dir, "CEL"))

## 2) Gunzip all files ended with *.CEL.gz
#cel_dir = file.path(base_dir, "CEL")
#gz_files = list.files(cel_dir, pattern="\\.CEL\\.gz$", full.names=TRUE)
#if (length(gz_files)) {
  #sapply(gz_files, R.utils::gunzip, overwrite=TRUE)
#}
#cel_files = list.files(cel_dir, pattern="\\.CEL$", full.names=TRUE)

# Data Preprocessing

## 1) RMA (log2)  *GPL570 (HG-U133 Plus 2.0) requires hgu133plus2cdf
#abatch = ReadAffy(filenames = cel_files)
#eset_rma = rma(abatch)                 # log2 RMA
#X_rma = exprs(eset_rma)                # probes x samples (log2 RMA)

## 2) Rename samples to GSM IDs and add phenotype
#gset = getGEO("GSE9101", GSEMatrix=TRUE, getGPL=FALSE)
#pd = pData(gset[[1]])[, c("geo_accession","title"), drop=FALSE]
#gsm_from_file = sub(".*(GSM[0-9]+).*", "\\1", basename(sampleNames(eset_rma)))
#m = match(gsm_from_file, pd$geo_accession)
#pd_ordered = pd[m, c("geo_accession","title"), drop = FALSE]
#rownames(pd_ordered) = pd_ordered$geo_accession         
#sampleNames(eset_rma) = pd_ordered$geo_accession 
#pData(eset_rma) = pd_ordered

## 3) Save log2 RMA probe expression matrix as csv
#write.csv(exprs(eset_rma), "./GSE9101 Matrix/GSE9101_RMA_log2_probelevel.csv") # probe expression matrix (probes x samples)

#X = as.matrix(read.csv("./GSE9101 Matrix/GSE9101_RMA_log2_probelevel.csv", row.names = 1, check.names = FALSE)) #dim(X) = 54675 x 12
#mode(X) = "numeric"   # ensure numeric

## 4) Drop Affy control probes
#keep_ctrl = !grepl("^AFFX", rownames(X))
#X = X[keep_ctrl, , drop=FALSE]

## 5) Remove low-expression probes
#cutoff = quantile(X, 0.05)                         # 5th percentile of GSE9101 intensities
#min_samples = ceiling(0.5 * ncol(X))               # expressed in >= half of samples
#keep_expr = rowSums(X > cutoff) >= min_samples
#expr_filtered = X[keep_expr, , drop=FALSE]         # probe-level filtered matrix

## 6) Merge probes into genes (GPL570)

### 1) Create annotation
#ann = AnnotationDbi::select(hgu133plus2.db, keys = rownames(expr_filtered), columns  = c("SYMBOL"), keytype  = "PROBEID")
#ann = ann[!duplicated(ann$PROBEID), ]
#expr_filtered = expr_filtered[ann$PROBEID, , drop=FALSE]
#ann = ann[match(rownames(expr_filtered), ann$PROBEID), ]

### 2) Keep only probes with a valid SYMBOL
#has_symbol = !is.na(ann$SYMBOL) & ann$SYMBOL != ""
#E = expr_filtered[has_symbol, , drop=FALSE]
#sym = ann$SYMBOL[has_symbol]

### 3) Average probes per gene
#gene_sum    = rowsum(E, group = sym)
#gene_counts = as.integer(table(sym)[rownames(gene_sum)])
#gene_expr_pre   = sweep(gene_sum, 1, gene_counts, "/")   # rows = genes, cols = samples

## 7) Remove zero variance genes
#X1_pre = gene_expr_pre[, 1:3] # all Control
#X2_pre = gene_expr_pre[, 4:12] # all cells stimulated with lipoprotein
#info_pre = information_extractor(X1_pre, X2_pre)
#keep_idx = (info_pre$S1_list != 0) & (info_pre$S2_list > 0)
#expr = gene_expr_pre[keep_idx, ]

## 8) Save final gene expression matrix as csv
#write.csv(expr, "./GSE9101 Matrix/GSE9101_gene_expression_matrix.csv", row.names = TRUE)

# DE analysis
expr = as.matrix(read.csv("./GSE9101 Matrix/GSE9101_gene_expression_matrix.csv", row.names = 1, check.names = FALSE)) # expression matrix (probes x samples)
mode(expr) = "numeric"   # ensure numeric
#dim(expr) = 20989 x 12

K1 = 3
K2 = 9
X1 = expr[, 1:3] # filtered Control
X2 = expr[, 4:12] # filtered cells stimulated with lipoprotein

# Run VREPB, DVEPB, Welch, B_F and VE_test

info = information_extractor(X1, X2)

alpha = 0.01
NPMLE_1D_parameter = c(1000, 0, 1.0)
NPMLE_2D_parameter = c(80, 80, 0.01, 1.0)
algorithm_list = c(1,2,3,4,5)

#result = solver(X1, X2, NPMLE_1D_parameter, NPMLE_2D_parameter, algorithm_list)

# Histogram of p-values for each method

#hist_VREPB = plotter_pvalue_histogram(result$VREPB)
#ggsave("./Result/pvalue_hist_VREPB_GSE9101.jpg", hist_VREPB, width = 9, height = 6.5, dpi = 300)

#hist_DVEPB = plotter_pvalue_histogram(result$DVEPB)
#ggsave("./Result/pvalue_hist_DVEPB_GSE9101.jpg", hist_DVEPB, width = 9, height = 6.5, dpi = 300)

#hist_Welch = plotter_pvalue_histogram(result$Welch)
#ggsave("./Result/pvalue_hist_Welch_GSE9101.jpg", hist_Welch, width = 9, height = 6.5, dpi = 300)

#hist_EV = plotter_pvalue_histogram(result$EV_test)
#ggsave("./Result/pvalue_hist_EV_GSE9101.jpg", hist_EV, width = 9, height = 6.5, dpi = 300)

#hist_BF = plotter_pvalue_histogram(result$B_F)
#ggsave("./Result/pvalue_hist_BF_GSE9101.jpg", hist_BF, width = 9, height = 6.5, dpi = 300)

# Number of discoveries for each method
#length(my_BH(result$VREPB, alpha)) # 0.1: 5240 / 0.01: 2322
#length(my_BH(result$DVEPB, alpha)) # 0.1: 5596 / 0.01: 2940
#length(my_BH(result$Welch, alpha)) # 0.1: 4528 / 0.01: 1614
#length(my_BH(result$B_F, alpha)) # 0.1: 413 / 0.01: 0
#length(my_BH(result$EV_test, alpha)) # 0.1: 4721 / 0.01: 2180

# Generate the H plot for GSE9101

#df_2D_GSE9101 = data.frame(x = result$'2D_grid'[,1], y = result$'2D_grid'[,2], prob = result$'2D_mass')
#write.csv(df_2D_GSE9101, './Result/2D_GSE9101.csv')

df_2D_GSE9101 = read.csv('./Result/2D_GSE9101.csv')

Total_mass = 0.999999
plot_H_GSE9101 = plotter_2D_ALL(Total_mass, df_2D_GSE9101, '')
#ggsave(filename = './Result/H(var1, var2) of GSE9101.jpg', plot = plot_H_GSE9101, width = 12, height = 10, dpi = 500)

# Generate the G plot for GSE9101

#df_1D_GSE9101 = data.frame(x = result$'1D_grid', prob = result$'1D_mass')
#write.csv(df_1D_GSE9101, './Result/1D_GSE9101.csv')

df_1D_GSE9101 = read.csv('./Result/1D_GSE9101.csv')
df_1D_GSE9101_with_ecdf = arrange(df_1D_GSE9101, x) %>% mutate(ecdf = cumsum(prob))
df_ecdf = data.frame(log_lambda = log(df_1D_GSE9101_with_ecdf$x), ecdf = df_1D_GSE9101_with_ecdf$ecdf)

ecdf_p = ggplot(df_ecdf, aes(x = log_lambda, y = ecdf)) +
  geom_line(size = 1.5, color = 'blue') +
  xlim(-2, 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = expression(log(lambda[i])),
    y = expression(hat(G)(lambda[i]))
  ) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none")

#ggsave("./Result/GSE9101 ecdf of G.jpg", ecdf_p, width = 8, height = 6, dpi = 400)

plot_G_GSE9101 = plotter_1D(df_1D_GSE9101, '')
#ggsave(filename = './Result/G(lambda) of GSE9101.jpg', plot = plot_G_GSE9101, width = 8, height = 6, dpi = 500)

