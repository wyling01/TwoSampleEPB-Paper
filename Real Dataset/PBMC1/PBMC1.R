library(limma)
library(edgeR)
library(SingleCellExperiment)

source('../../helper functions/voomByGroup.R')
source("../../EPB main/EPB.R")
source('../../EPB main/voom_EPB.R')

# #---- Download data----
data_commit = "38c8abfba1d329dcef9d099bc0351f8f293ac307"
dir.create("PBMC1 Data", recursive=TRUE, showWarnings=FALSE)
download.file(
  paste0(
    "https://raw.githubusercontent.com/YOU-k/voomByGroup/",
    data_commit,
    "/data/pbmc_pb.rds"
  ),
  destfile = "./PBMC1 Data/pbmc_pb.rds",
  mode = "wb"
)
pbmc_pb = readRDS("./PBMC1 Data/pbmc_pb.rds")

cell_type = "CD56dim CD16+ NK"
sample_index = c(which(pbmc_pb$group_id == "Moderate"), which(pbmc_pb$group_id == "HC"))

counts = as.matrix(assay(pbmc_pb, cell_type)[, sample_index, drop=FALSE])
group = as.character(pbmc_pb$group_id[sample_index])
names(group) = colnames(counts)

keep_gene = rowSums(counts) > 50
counts = counts[keep_gene, , drop=FALSE]
stopifnot(identical(dim(counts), c(8462L, 8L)))

saveRDS(counts, "./PBMC1 Data/counts_ModHC.rds")
saveRDS(group, "./PBMC1 Data/group_id_ModHC.rds")

# #----Load PBMC1 data----

counts = readRDS("./PBMC1 Data/counts_ModHC.rds")
group = readRDS("./PBMC1 Data/group_id_ModHC.rds")

## Create DGE list
y = DGEList(counts=counts, group=group)

## Calculate TMM normalization factors
y = calcNormFactors(y)

## Create design matrix
group_assignment = y$samples$group
design = model.matrix(~ 0 + group_assignment)

## NPMLE Parameters

alpha = 0.1
VR_parameter = c(1000, 0, 1.0)
DV_parameter = c(80, 80, 0.01, 1.0)

# #----EPB with weight info from voom----

v = voom(y, design = design, plot = FALSE)

is_A = group_assignment == "Moderate"
is_B = group_assignment == "HC"

Y_A = v$E[, is_A, drop=FALSE]
Y_B = v$E[, is_B, drop=FALSE]
w_A = v$weights[, is_A, drop=FALSE]
w_B = v$weights[, is_B, drop=FALSE]


algorithm_list = c(1, 2, 3, 4, 5)
result = solver_weighted(Y_A, Y_B, w_A, w_B,VR_parameter,
                         DV_parameter,algorithm_list)

# ----Saving results----
dir.create("Result", recursive=TRUE, showWarnings=FALSE)
df_2D_PBMC1 = data.frame(x=result$`2D_grid`[, 1],y=result$`2D_grid`[, 2],prob=result$`2D_mass`)
write.csv(df_2D_PBMC1, "./Result/2D_PBMC_voom.csv", row.names=FALSE)

df_1D_PBMC1 = data.frame(x=result$`1D_grid`,prob=result$`1D_mass`)
write.csv(df_1D_PBMC1, "./Result/1D_PBMC_voom.csv", row.names=FALSE)


# ----EPB with weight info from voombygroup----
vbg = voomByGroup(y, design=design, group=group_assignment,
                  plot="none", print=FALSE)

Y_A_vbg = vbg$E[, is_A, drop=FALSE]
Y_B_vbg = vbg$E[, is_B, drop=FALSE]
w_A_vbg = vbg$weights[, is_A, drop=FALSE]
w_B_vbg = vbg$weights[, is_B, drop=FALSE]

result_vbg = solver_weighted(
  Y_A_vbg, Y_B_vbg, w_A_vbg, w_B_vbg,
  VR_parameter, DV_parameter, algorithm_list
)

# ----Saving results----

df_2D_PBMC1_vbg = data.frame(x=result_vbg$`2D_grid`[, 1],y=result_vbg$`2D_grid`[, 2],
                             prob=result_vbg$`2D_mass`)
write.csv(df_2D_PBMC1_vbg, "./Result/2D_PBMC_vbg.csv", row.names=FALSE)

df_1D_PBMC1_vbg = data.frame(x=result_vbg$`1D_grid`,prob=result_vbg$`1D_mass`)
write.csv(df_1D_PBMC1_vbg, "./Result/1D_PBMC_vbg.csv", row.names=FALSE)

# Number of significant discoveries for each voom weighting method
significance_summary = data.frame(
  voom_method = rep(c("Voom", "VoombyGroup"), each=5),
  method = rep(c("VREPB", "DVEPB", "Welch", "B_F", "EV"), times=2),
  n_significant = c(
    length(my_BH(result$VREPB, alpha)),       # 7
    length(my_BH(result$DVEPB, alpha)),       # 286
    length(my_BH(result$Welch, alpha)),       # 0
    length(my_BH(result$B_F, alpha)),         # 0
    length(my_BH(result$EV, alpha)),          # 4
    length(my_BH(result_vbg$VREPB, alpha)),   # 7
    length(my_BH(result_vbg$DVEPB, alpha)),   # 339
    length(my_BH(result_vbg$Welch, alpha)),   # 0
    length(my_BH(result_vbg$B_F, alpha)),     # 0
    length(my_BH(result_vbg$EV, alpha))       # 7
  )
)
print(significance_summary, row.names=FALSE)
