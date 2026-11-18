library(GEOquery)
library(Biobase)
library(affy)
library(hgu133plus2cdf)
library(hgu133plus2.db)

source('../../EPB main/EPB.R')

##----Data Download----

## 1) Get raw CEL files
base_dir = "GSE9101_raw"
dir.create(base_dir, showWarnings = FALSE)
getGEOSuppFiles("GSE9101", baseDir = base_dir)
untar(file.path(base_dir, "GSE9101", "GSE9101_RAW.tar"), exdir = file.path(base_dir, "CEL"))

## 2) Gunzip all files ended with *.CEL.gz
cel_dir = file.path(base_dir, "CEL")
gz_files = list.files(cel_dir, pattern="\\.CEL\\.gz$", full.names=TRUE)
if (length(gz_files)) {
  sapply(gz_files, R.utils::gunzip, overwrite=TRUE)
}
cel_files = sort(list.files(cel_dir, pattern="\\.CEL$", full.names=TRUE))

# Data Preprocessing

## 1) RMA (log2)  *GPL570 (HG-U133 Plus 2.0) requires hgu133plus2cdf
abatch = ReadAffy(filenames = cel_files)
eset_rma = rma(abatch)                 # log2 RMA

## 2) Rename samples to GSM IDs and add phenotype
gset = getGEO("GSE9101", GSEMatrix=TRUE, getGPL=FALSE)
pd = pData(gset[[1]])[, c("geo_accession","title"), drop=FALSE]
gsm_from_file = sub(".*(GSM[0-9]+).*", "\\1", basename(sampleNames(eset_rma)))
m = match(gsm_from_file, pd$geo_accession)
pd_ordered = pd[m, c("geo_accession","title"), drop = FALSE]
rownames(pd_ordered) = pd_ordered$geo_accession
sampleNames(eset_rma) = pd_ordered$geo_accession
pData(eset_rma) = pd_ordered

## 3) Save log2 RMA probe expression matrix as csv
dir.create("GSE9101 Matrix", recursive=TRUE, showWarnings=FALSE)
write.csv(exprs(eset_rma), "./GSE9101 Matrix/GSE9101_RMA_log2_probelevel.csv") # probe expression matrix (probes x samples)

X = as.matrix(read.csv("./GSE9101 Matrix/GSE9101_RMA_log2_probelevel.csv", row.names = 1, check.names = FALSE)) #dim(X) = 54675 x 12
mode(X) = "numeric"   # ensure numeric

## 4) Drop Affy control probes
keep_ctrl = !grepl("^AFFX", rownames(X))
X = X[keep_ctrl, , drop=FALSE]

## 5) Remove low-expression probes
cutoff = quantile(X, 0.05)                         # 5th percentile of GSE9101 intensities
min_samples = ceiling(0.5 * ncol(X))               # expressed in >= half of samples
keep_expr = rowSums(X > cutoff) >= min_samples
expr_filtered = X[keep_expr, , drop=FALSE]         # probe-level filtered matrix

## 6) Merge probes into genes (GPL570)

annotation_source = "database"

if (annotation_source == "database") {
  ann = AnnotationDbi::select(hgu133plus2.db, keys=rownames(expr_filtered), columns="SYMBOL", keytype="PROBEID")
  ann = ann[!duplicated(ann$PROBEID), ]
  has_symbol = !is.na(ann$SYMBOL) & ann$SYMBOL != ""
  ann_used = ann[has_symbol, c("PROBEID", "SYMBOL"), drop=FALSE]
  write.csv(ann_used, "./GSE9101 Matrix/GSE9101_probe_to_symbol_mapping.csv", row.names=FALSE)
} else if (annotation_source == "saved") {
  ann_used = read.csv("./GSE9101 Matrix/GSE9101_probe_to_symbol_mapping.csv", stringsAsFactors=FALSE, check.names=FALSE)
} else {
  stop("annotation_source must be 'database' or 'saved'.")
}

if (anyDuplicated(ann_used$PROBEID)) stop("Duplicated PROBEID values in ann_used.")
if (anyNA(ann_used$PROBEID) || anyNA(ann_used$SYMBOL) || any(ann_used$SYMBOL == "")) stop("Invalid PROBEID or SYMBOL values in ann_used.")

probe_index = match(ann_used$PROBEID, rownames(expr_filtered))
if (anyNA(probe_index)) stop("Some probes in ann_used are absent from expr_filtered.")

E = expr_filtered[probe_index, , drop=FALSE]
sym = ann_used$SYMBOL
stopifnot(identical(rownames(E), ann_used$PROBEID))

### Average probes per gene
gene_sum = rowsum(E, group=sym)
gene_counts = as.integer(table(sym)[rownames(gene_sum)])
gene_expr_pre = sweep(gene_sum, 1, gene_counts, "/")

## 7) Remove zero variance genes
X1_pre = gene_expr_pre[, 1:3] # all Control
X2_pre = gene_expr_pre[, 4:12] # all cells stimulated with lipoprotein
info_pre = information_extractor(X1_pre, X2_pre)
keep_idx = (info_pre$S1_list > 0) & (info_pre$S2_list > 0)
expr = gene_expr_pre[keep_idx, ]

## 8) Save final gene expression matrix as csv
write.csv(expr, "./GSE9101 Matrix/GSE9101_gene_expression_matrix.csv", row.names = TRUE)

#----DE analysis----
expr = as.matrix(read.csv("./GSE9101 Matrix/GSE9101_gene_expression_matrix.csv", row.names = 1, check.names = FALSE)) # expression matrix (probes x samples)
mode(expr) = "numeric"   # ensure numeric
#dim(expr) = 20980 * 12

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

result = solver(X1, X2, NPMLE_1D_parameter, NPMLE_2D_parameter, algorithm_list)


# Number of significant discoveries for each method
significance_summary = data.frame(
  method = c("VREPB", "DVEPB", "Welch", "B_F", "EV"),
  n_significant = c(
    length(my_BH(result$VREPB, alpha)), # 0.01: 2325
    length(my_BH(result$DVEPB, alpha)), # 0.01: 2941
    length(my_BH(result$Welch, alpha)), # 0.01: 1614
    length(my_BH(result$B_F, alpha)),   # 0.01: 0
    length(my_BH(result$EV, alpha))     # 0.01: 2180
  )
)
print(significance_summary, row.names = FALSE)

##----Saving results----
dir.create("Result", recursive=TRUE, showWarnings=FALSE)

df_2D_GSE9101 = data.frame(x = result$'2D_grid'[,1], y = result$'2D_grid'[,2], prob = result$'2D_mass')
write.csv(df_2D_GSE9101, './Result/2D_GSE9101.csv')

df_1D_GSE9101 = data.frame(x = result$'1D_grid', prob = result$'1D_mass')
write.csv(df_1D_GSE9101, './Result/1D_GSE9101.csv')

