library(MAnorm2)

source('../../EPB main/EPB.R')

# #----Data download----
data_commit = "5fd6bdb1a19f5f4cd000a9ed4c804dcddcc443b5"
base_url = paste0(
  "https://raw.githubusercontent.com/tushiqi/MAnorm2/",
  data_commit,
  "/utility/code-MAnorm2Paper/data/"
)
dir.create("chipseq_data", recursive=TRUE, showWarnings=FALSE)

download.file(
  paste0(base_url, "GM12891_cmp_GM12892.H3K4me3_hieMA.xls"),
  "./chipseq_data/GM12891_cmp_GM12892.H3K4me3_hieMA.xls",
  mode="wb"
)

download.file(
  paste0(base_url, "GM12891_cmp_GM12892.H3K4me3_profile_bins.xls"),
  "./chipseq_data/GM12891_cmp_GM12892.H3K4me3_profile_bins.xls",
  mode="wb"
)

# #----Data Preparation----

profile_bins <- read.table("./chipseq_data/GM12891_cmp_GM12892.H3K4me3_profile_bins.xls",
                           header = TRUE, stringsAsFactors = FALSE)
n1 <- 3
n <- 6
normalized <- read.table("./chipseq_data/GM12891_cmp_GM12892.H3K4me3_hieMA.xls",
                         header = TRUE, stringsAsFactors = FALSE)

rownames(normalized) <- paste(profile_bins$chrom, profile_bins$start, profile_bins$end,
                              sep = "_")
design <- cbind(cond1 = 1, cond2_minus_cond1 = c(rep(0, n1), rep(1, n - n1)))
rownames(design) <- names(normalized)


# #----DE analysis----
expr <- as.matrix(normalized)
K1 = 3
K2 = 3
X1 = expr[, 1:3, drop=FALSE] # GM12891
X2 = expr[, 4:6, drop=FALSE] # GM12892

info = information_extractor(X1, X2)

alpha = 0.01
NPMLE_1D_parameter = c(5000, 0, 1.0)
NPMLE_2D_parameter = c(80, 80, 0.01, 1.0)
algorithm_list = c(1,2,3,4,5)

result = solver(X1, X2, NPMLE_1D_parameter, NPMLE_2D_parameter, algorithm_list)

## ----Saving results----
dir.create("Result", recursive=TRUE, showWarnings=FALSE)

df_2D_chiseq = data.frame(x = result$'2D_grid'[,1], y = result$'2D_grid'[,2], prob = result$'2D_mass')
write.csv(df_2D_chiseq, './Result/2D_chipseq.csv')

df_1D_chipseq = data.frame(x = result$'1D_grid', prob = result$'1D_mass')
write.csv(df_1D_chipseq, './Result/1D_chipseq.csv')


#EV Ratio Test
ev_test_ratio <- function(x1, x2, lambda) {
  stopifnot(lambda > 0)
  
  n1 <- length(x1); n2 <- length(x2)
  m1 <- mean(x1);   m2 <- mean(x2)
  s1 <- var(x1);    s2 <- var(x2)
  
  df <- n1 + n2 - 2
  Sp2 <- ((n1 - 1) * s1 + (n2 - 1) * (lambda * s2)) / df
  se <- sqrt(Sp2 * (1/n1 + 1/(lambda * n2)))
  tstat <- (m1 - m2) / se
  
  2 * pt(-abs(tstat), df = df)
}

P_value_EV_ratio = function (m, X1, X2,lambda) {
  
  P_value_list_EV_ratio = rep(0, m)
  for (i in 1:m) {
    X1_i = X1[i, ]
    X2_i = X2[i, ]
    P_value_list_EV_ratio[i] =  ev_test_ratio(X1_i, X2_i, lambda)
  }
  
  return (P_value_list_EV_ratio)
}

norm2 <- profile_bins
norm2[4:(3 + n)] <- normalized

# Follow the original workflow for differential analysis.
conds2 <- list(GM12891 = bioCond(norm2[4:(3 + n1)], norm2[4:(3 + n1) + n], name = "GM12891"),
               GM12892 = bioCond(norm2[(4 + n1):(3 + n)], norm2[(4 + n1):(3 + n) + n], name = "GM12892"))
conds2 <- fitMeanVarCurve(conds2, method = "parametric", init.coef = c(0.1, 10))
lambda <- conds2$GM12891$fit.info$ratio.var/conds2$GM12892$fit.info$ratio.var
print(c(EVratio = lambda))

result$EV_ratio  <- P_value_EV_ratio(info$m,X1,X2,lambda)

# Number of significant discoveries for each method
significance_summary = data.frame(
  method = c("VREPB", "DVEPB", "Welch", "B_F", "EV", "EVratio"),
  n_significant = c(
    length(my_BH(result$VREPB, alpha)),   # 2966
    length(my_BH(result$DVEPB, alpha)),   # 9342
    length(my_BH(result$Welch, alpha)),   # 0
    length(my_BH(result$B_F, alpha)),     # 0
    length(my_BH(result$EV, alpha)),      # 2938
    length(my_BH(result$EV_ratio, alpha)) # 2968
  )
)
print(significance_summary, row.names = FALSE)

