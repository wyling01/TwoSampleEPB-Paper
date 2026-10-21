library(ggplot2)
library(limma)
library(dplyr)
library(oligo)
library(edgeR)
library(readxl)
library(MAnorm2)

source('../../EPB main/EPB.R')
source('../../helper functions/Distribution Visualization helper.R')

#### Data Preparation ###

profile_bins <- read.table("./chipseq_data/GM12891_cmp_GM12892.H3K4me3_profile_bins.xls",
                           header = TRUE, stringsAsFactors = FALSE)
n1 <- 3
n <- 6
normalized <- read.table("./chipseq_data/GM12891_cmp_GM12892.H3K4me3_hieMA.xls",
                         header = TRUE, stringsAsFactors = FALSE)

head(normalized)
rownames(normalized) <- paste(profile_bins$chrom, profile_bins$start, profile_bins$end,
                              sep = "_")
design <- cbind(cond1 = 1, cond2_minus_cond1 = c(rep(0, n1), rep(1, n - n1)))
rownames(design) <- names(normalized)


#### DE analysis ####
expr <- as.matrix(normalized)
K1 = 3
K2 = 3
X1 = expr[, 1:3] # filtered Control
X2 = expr[, 4:6]


#### Run VREPB, DVEPB, Welch, B_F and VE_test ###

info = information_extractor(X1, X2)

alpha = 0.01
NPMLE_1D_parameter = c(5000, 0.01, 0.99)
NPMLE_2D_parameter = c(50, 50, 0.001, 1.0)
algorithm_list = c(1,2,3,4,5)

result = solver(X1, X2, NPMLE_1D_parameter, NPMLE_2D_parameter, algorithm_list)

#### Histogram of p-values for each method ####

hist_VREPB = plotter_pvalue_histogram(result$VREPB, nbins = 100,ymax=15000)
ggsave("./Result/pvalue_hist_VREPB_chipseq.jpg", hist_VREPB, width = 9, height = 6.5, dpi = 300)

hist_DVEPB = plotter_pvalue_histogram(result$DVEPB,nbins = 100,ymax=15000)
ggsave("./Result/pvalue_hist_DVEPB_chipseq.jpg", hist_DVEPB, width = 9, height = 6.5, dpi = 300)

hist_Welch = plotter_pvalue_histogram(result$Welch,nbins = 100,ymax=15000)
ggsave("./Result/pvalue_hist_Welch_chipseq.jpg", hist_Welch, width = 9, height = 6.5, dpi = 300)

hist_EV = plotter_pvalue_histogram(result$EV_test,nbins = 100,ymax=15000)
ggsave("./Result/pvalue_hist_EV_chipseq.jpg", hist_EV, width = 9, height = 6.5, dpi = 300)

hist_BF = plotter_pvalue_histogram(result$B_F,nbins = 100,ymax=15000)
ggsave("./Result/pvalue_hist_BF_chipseq.jpg", hist_BF, width = 9, height = 6.5, dpi = 300)

# Number of discoveries for each method
length(my_BH(result$VREPB, alpha)) # 0.01: 2966
length(my_BH(result$DVEPB, alpha)) # 0.01: 9345
length(my_BH(result$Welch, alpha)) # 0.01: 0
length(my_BH(result$B_F, alpha)) # 0.01: 0
length(my_BH(result$EV_test, alpha)) # 0.01: 2938

#### Generate the H plot for Chipseq ####

#df_2D_chiseq = data.frame(x = result$'2D_grid'[,1], y = result$'2D_grid'[,2], prob = result$'2D_mass')
#write.csv(df_2D_chiseq, './Result/2D_chipseq.csv')

df_2D_chiseq = read.csv('./Result/2D_chiseq.csv')

Total_mass = 0.999999
plot_H_chiseq = plotter_2D_ALL(Total_mass, df_2D_chiseq, '')
#ggsave(filename = './Result/H(var1, var2) of chipseq.jpg', plot = plot_H_chiseq, width = 12, height = 10, dpi = 500)

#### Generate the G plot for Chipseq ####

#df_1D_chipseq = data.frame(x = result$'1D_grid', prob = result$'1D_mass')
#write.csv(df_1D_chipseq, './Result/1D_chipseq.csv')

df_1D_chipseq = read.csv('./Result/1D_chipseq.csv')
df_1D_chipseq_with_ecdf = arrange(df_1D_chipseq, x) %>% mutate(ecdf = cumsum(prob))
df_ecdf = data.frame(log_lambda = log(df_1D_chipseq_with_ecdf$x), ecdf = df_1D_chipseq_with_ecdf$ecdf)

ecdf_p = ggplot(df_ecdf, aes(x = log_lambda, y = ecdf)) +
  geom_line(size = 1.5, color = 'blue') +
  xlim(0.20, 0.30) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = expression(log(lambda[i])),
    y = expression(hat(G)(lambda[i]))
  ) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none")

#ggsave("./Result/chipseq ecdf_of_G.jpg", ecdf_p, width = 8, height = 6, dpi = 400)

plot_G_chipseq = plotter_1D(df_1D_chipseq, '', c(0.20,0.30))
#ggsave(filename = './Result/G(lambda) of chipseq.jpg', plot = plot_G_chipseq, width = 8, height = 6, dpi = 500)


#### EV Ratio Test####
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

P_value_EV_ratio_test = function (m, X1, X2,lambda) {
  
  P_value_list_EV_ratio_test = rep(0, m)
  for (i in 1:m) {
    X1_i = X1[i, ]
    X2_i = X2[i, ]
    P_value_list_EV_ratio_test[i] =  ev_test_ratio(X1_i, X2_i, lambda)
  }
  
  return (P_value_list_EV_ratio_test)
}

norm2 <- profile_bins
norm2[4:(3 + n)] <- normalized

# Follow the original workflow for differential analysis.
conds2 <- list(GM12891 = bioCond(norm2[4:(3 + n1)], norm2[4:(3 + n1) + n], name = "GM12891"),
               GM12892 = bioCond(norm2[(4 + n1):(3 + n)], norm2[(4 + n1):(3 + n) + n], name = "GM12892"))
conds2 <- fitMeanVarCurve(conds2, method = "parametric", init.coef = c(0.1, 10))
res2 <- diffTest(conds2[[1]], conds2[[2]])
lambda <- conds2$GM12891$fit.info$ratio.var/conds2$GM12892$fit.info$ratio.var

result$EV_ratio  <- P_value_EV_ratio_test(info$m,X1,X2,lambda)
length(my_BH(result$EV_ratio, alpha)) # 2968

hist_ev_ratio = plotter_pvalue_histogram(result$EV_ratio,nbins = 100,ymax=15000)
#ggsave("./Result/pvalue_hist_EV_ratio_chipseq.jpg", hist_ev_ratio, width = 9, height = 6.5, dpi = 300)


x_rng <- c(0.1, 0.4)

df_step <- data.frame(
  x = c(x_rng[1], log(lambda), x_rng[2]),
  y = c(0, 1, 1)
)
ecdf_p_combined <- ggplot(df_ecdf) +
  geom_line(aes(log_lambda, ecdf, color = "VREPB"), linewidth = 1.5) +
  geom_step(
    data = df_step,
    aes(x, y, color = "EV-ratio-test"),
    direction = "hv",          # horizontal then vertical jump at log_lambda0
    linewidth = 1.2,
    linetype = "solid"
  ) +
  coord_cartesian(xlim = c(0.20,0.30), ylim = c(0, 1)) +
  scale_color_manual(
    name = NULL,
    values = c("VREPB" = "blue", "EV-ratio-test" = "black"),
    breaks = c( "VREPB","EV-ratio-test")
  ) +
  labs(
    x = expression(log(lambda[i])),
    y = expression(hat(G)(lambda[i]))
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = c(0.36, 1),
    legend.justification = c("right", "top"),
    legend.text.align = 0,
    legend.key.height = unit(2, "lines")
  )
ggsave("./Result/chipseq ecdf_of_G_combined.jpg", ecdf_p_combined, width = 8, height = 6, dpi = 400)

