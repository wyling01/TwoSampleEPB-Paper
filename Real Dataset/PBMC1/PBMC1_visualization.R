library(dplyr)
library(ggplot2)

source('../../helper functions/Distribution Visualization helper.R')

# Voom

## H

df_2D_weighted_voom = read.csv('./Mod vs HC Result/2D_weighted_voom_ModHC.csv')
Total_mass = 0.99999
plot_G_voom = plotter_2D(Total_mass, df_2D_weighted_voom, '')
#ggsave(filename = './Mod vs HC Result/H(var1, var2) with Voom.jpg', plot = plot_G_voom, width = 12, height = 10, dpi = 500)

# VoombyGroup

## H

df_2D_weighted_vbg = read.csv('./Mod vs HC Result/2D_weighted_vbg_ModHC.csv')
plot_G_vbg = plotter_2D(Total_mass, df_2D_weighted_vbg, '')
#ggsave(filename = './Mod vs HC Result/H(var1, var2) with Voombygroup.jpg', plot = plot_G_vbg, width = 12, height = 10, dpi = 500)

# Voom vs VoombyGroup

df_1D_weighted_voom = read.csv('./Mod vs HC Result/1D_weighted_voom_ModHC.csv')
df_1D_weighted_vbg = read.csv('./Mod vs HC Result/1D_weighted_vbg_ModHC.csv')


df_lambda = data.frame(x =  df_1D_weighted_voom$x, prob =  df_1D_weighted_voom$prob)
df_lambda_vbg = data.frame(x = df_1D_weighted_vbg$x, prob = df_1D_weighted_vbg$prob)

df_lambda_with_ecdf = arrange(df_lambda, x) %>% mutate(ecdf = cumsum(prob))
df_lambda_with_ecdf_vbg = arrange(df_lambda_vbg, x) %>% mutate(ecdf = cumsum(prob))

df1 = data.frame(log_lambda = log(df_lambda_with_ecdf$x), ecdf = df_lambda_with_ecdf$ecdf, method = "Voom")
df2 = data.frame(log_lambda = log(df_lambda_with_ecdf_vbg$x), ecdf = df_lambda_with_ecdf_vbg$ecdf, method = "VoombyGroup")

df_plot = rbind(df1, df2)

method_labels = c("Voom" = expression(hat(G)~"with Voom"), "VoombyGroup" = expression(hat(G)~"with VoombyGroup"))

p = ggplot(df_plot, aes(x = log_lambda, y = ecdf, color = method)) +
  geom_line(size = 1.5) +
  scale_color_manual(
    name = NULL,
    values = c("Voom" = "#0072B2", "VoombyGroup" = "black"),
    labels = method_labels
  ) +
  xlim(-2, 2) +
  labs(
    x = expression(log(lambda[i])),
    y = expression(hat(G)(lambda[i]))
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = c(0.5, 1),
    legend.justification = c("right", "top"),
    legend.text.align = 0
  )

#ggsave("./Mod vs HC Result/PBMC1 ecdf of G.jpg", p, width = 8, height = 6, dpi = 400)