library(limma)
library(edgeR)
library(ggplot2)
library(tidyverse)

source('../../helper functions/voomByGroup.R')
source('../../EPB main/voom_EPB.R')
source('../../helper functions/Distribution Visualization helper.R')

# load PBMC1 data

counts = readRDS("./PBMC1 Data/counts.rds")
group = readRDS("./PBMC1 Data/group_id.rds")

# Data Filtering

high_expression_counts = counts[rowSums(counts)>50,]
filter = (group == 'Moderate' | group == 'HC') # 5 Mod vs 3 HC
filtered_counts = high_expression_counts[, filter]
filtered_group = group[filter]

# Create DGE list

y = DGEList(filtered_counts)
y$samples$group = filtered_group

# Calculate Normalizing Factors

y = calcNormFactors(y)

# Create the design matrix based on group information.

group_assignment = y$samples$group
design = model.matrix(~ 0 + group_assignment)

## Normalize library size R

R = y$samples$lib.size * y$samples$norm.factors
Y = t(log2((t(filtered_counts) + 0.5)/(R+1) * 1e6))
Y = normalizeBetweenArrays(Y, method='none')
Y1 = Y[, c(1,2,3,4,8)]
Y2 = Y[, c(5,6,7)]

# NPMLE Parameters

alpha = 0.1
VR_parameter = c(1000, 0, 1.0)
DV_parameter = c(80, 80, 0.01, 1.0)

# EPB with weight info from voom

v = voom(y, design = design, plot = FALSE)

w_A = v$weights[, c(1,2,3,4,8)]
w_B = v$weights[, c(5,6,7)]

Y_A = Y1
Y_B = Y2

information = information_extractor(Y_A, Y_B, w_A, w_B)

## VREPB

VR_NPMLE_result = VR_NPMLE(information,B = VR_parameter[1], lower_quantile = VR_parameter[2], upper_quantile = VR_parameter[3])
p_value_VREPB = P_value_VREPB(information, VR_NPMLE_result$grid,VR_NPMLE_result$mass)
length(my_BH(p_value_VREPB, alpha))  # 7

hist_VREPB_voom = plotter_pvalue_histogram(p_value_VREPB, ymax = 1000)
#ggsave("./Mod vs HC Result/pvalue_hist_VREPB_voom.jpg", hist_VREPB_voom, width = 9, height = 6.5, dpi = 300)

df_1D_weighted_voom = data.frame(x = VR_NPMLE_result$grid, prob = VR_NPMLE_result$mass)
#write.csv(df_1D_weighted_voom, './Mod vs HC Result/1D_weighted_voom_ModHC.csv')

## DVEPB

DV_NPMLE_result = DV_NPMLE(S1_list = information$S1_list, S2_list = information$S2_list, B1 = DV_parameter[1], B2 = DV_parameter[2], m = information$m, n1 = information$n1, n2 = information$n2, lower_quantile = DV_parameter[3], upper_quantile = DV_parameter[4])
p_value_DVEPB = P_value_DVEPB(information, DV_NPMLE_result$grid, DV_NPMLE_result$mass)
length(my_BH(p_value_DVEPB, alpha)) # 260

hist_DVEPB_voom = plotter_pvalue_histogram(p_value_DVEPB, ymax = 1000)
#ggsave("./Mod vs HC Result/pvalue_hist_DVEPB_voom.jpg", hist_DVEPB_voom, width = 9, height = 6.5, dpi = 300)

df_2D_weighted_voom = data.frame(x = DV_NPMLE_result$grid[, 1], y = DV_NPMLE_result$grid[, 2], prob = DV_NPMLE_result$mass)
#write.csv(df_2D_weighted_voom, './Mod vs HC Result/2D_weighted_voom_ModHC.csv')

## Welch

p_value_Welch = P_value_Welch(information)
length(my_BH(p_value_Welch, alpha)) # 0

hist_Welch_voom = plotter_pvalue_histogram(p_value_Welch, ymax = 1000)
#ggsave("./Mod vs HC Result/pvalue_hist_Welch_voom.jpg", hist_Welch_voom, width = 9, height = 6.5, dpi = 300)

## EV-test

p_value_EV = P_value_EV_test(information)
length(my_BH(p_value_EV, alpha)) # 4

hist_EV_voom = plotter_pvalue_histogram(p_value_EV, ymax = 1000)
#ggsave("./Mod vs HC Result/pvalue_hist_EV_voom.jpg", hist_EV_voom, width = 9, height = 6.5, dpi = 300)

## B-F

p_value_BF = P_value_BF_test(information)
length(my_BH(p_value_BF, alpha)) # 0

hist_BF_voom = plotter_pvalue_histogram(p_value_BF, ymax = 1000)
#ggsave("./Mod vs HC Result/pvalue_hist_BF_voom.jpg", hist_BF_voom, width = 9, height = 6.5, dpi = 300)

# EPB with weight info from voombygroup

vbg = voomByGroup(y, design = design, group = group_assignment, plot = FALSE)

w_vbg = vbg$weights
w_A = w_vbg[, c(1,2,3,4,8)]
w_B = w_vbg[, c(5,6,7)]

Y_A = Y1
Y_B = Y2

information_vbg = information_extractor(Y_A, Y_B, w_A, w_B)

## VREPB

VR_NPMLE_result_vbg = VR_NPMLE(information_vbg,B = VR_parameter[1], lower_quantile = VR_parameter[2], upper_quantile = VR_parameter[3])
p_value_VREPB_vbg = P_value_VREPB(information_vbg, VR_NPMLE_result_vbg$grid,VR_NPMLE_result_vbg$mass)
length(my_BH(p_value_VREPB_vbg, alpha))  # 7

hist_VREPB_vbg = plotter_pvalue_histogram(p_value_VREPB_vbg, ymax = 1000)
#ggsave("./Mod vs HC Result/pvalue_hist_VREPB_vbg.jpg", hist_VREPB_vbg, width = 9, height = 6.5, dpi = 300)

df_1D_weighted_vbg = data.frame(x = VR_NPMLE_result_vbg$grid, prob = VR_NPMLE_result_vbg$mass)
#write.csv(df_1D_weighted_vbg, './Mod vs HC Result/1D_weighted_vbg_ModHC.csv')

## DVEPB

DV_NPMLE_result_vbg = DV_NPMLE(S1_list = information_vbg$S1_list, S2_list = information_vbg$S2_list, B1 = DV_parameter[1], B2 = DV_parameter[2], m = information_vbg$m, n1 = information_vbg$n1, n2 = information_vbg$n2, lower_quantile = DV_parameter[3], upper_quantile = DV_parameter[4])
p_value_DVEPB_vbg = P_value_DVEPB(information_vbg, DV_NPMLE_result_vbg$grid, DV_NPMLE_result_vbg$mass)
length(my_BH(p_value_DVEPB_vbg, alpha)) # 316

hist_DVEPB_vbg = plotter_pvalue_histogram(p_value_DVEPB_vbg, ymax = 1000)
#ggsave("./Mod vs HC Result/pvalue_hist_DVEPB_vbg.jpg", hist_DVEPB_vbg, width = 9, height = 6.5, dpi = 300)

df_2D_weighted_vbg = data.frame(x = DV_NPMLE_result_vbg$grid[, 1], y = DV_NPMLE_result_vbg$grid[, 2], prob = DV_NPMLE_result_vbg$mass)
#write.csv(df_2D_weighted_vbg, './Mod vs HC Result/2D_weighted_vbg_ModHC.csv')

## Welch

p_value_Welch_vbg = P_value_Welch(information_vbg)
length(my_BH(p_value_Welch_vbg, alpha)) # 0

hist_Welch_vbg = plotter_pvalue_histogram(p_value_Welch_vbg, ymax = 1000)
#ggsave("./Mod vs HC Result/pvalue_hist_Welch_vbg.jpg", hist_Welch_vbg, width = 9, height = 6.5, dpi = 300)

##  EV-test

p_value_EV_vbg = P_value_EV_test(information_vbg)
length(my_BH(p_value_EV_vbg, alpha)) # 7

hist_EV_vbg = plotter_pvalue_histogram(p_value_EV_vbg, ymax = 1000)
#ggsave("./Mod vs HC Result/pvalue_hist_EV_vbg.jpg", hist_EV_vbg, width = 9, height = 6.5, dpi = 300)

## B-F

p_value_BF_vbg = P_value_BF_test(information_vbg)
length(my_BH(p_value_BF_vbg, alpha)) # 0

hist_BF_vbg = plotter_pvalue_histogram(p_value_BF_vbg, ymax = 1000)
#ggsave("./Mod vs HC Result/pvalue_hist_BF_vbg.jpg", hist_BF_vbg, width = 9, height = 6.5, dpi = 300)