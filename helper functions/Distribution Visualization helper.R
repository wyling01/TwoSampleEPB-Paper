find_threshold = function(mass, Total_mass = 0.98) {
  sorted_mass = sort(mass, decreasing = TRUE)
  sum = 0
  for (i in c(1:length(mass))) {
    sum = sum + sorted_mass[i]
    if (sum >= Total_mass) {
      return (sorted_mass[i+1])
    }
  }
}

plotter_2D_ALL = function(Total_mass, df_2D, title) {
  threshold = find_threshold(df_2D$prob, Total_mass = Total_mass)
  size_breaks = seq(min(df_2D$prob[df_2D$prob > threshold]), max(df_2D$prob), length.out = 10)
  plot_2D = ggplot(df_2D, aes(x = log(x), y = log(y))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 1) +
    geom_point(aes(size = ifelse(prob < threshold, NA, prob)), color = "blue", alpha = 0.8) +
    scale_size_continuous(range = c(0.2, 5), breaks = size_breaks) +
    theme_minimal() +
    #xlim(-5, 2) +
    labs(title = title,
         x = expression(log({sigma[iA]}^2)),
         y = expression(log({sigma[iB]}^2))) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title = element_text(size = 35),
          axis.text = element_text(size = 28),
          legend.position = 'none')
  return (plot_2D)
}

plotter_2D = function(Total_mass, df_2D, title) {
  threshold = find_threshold(df_2D$prob, Total_mass = Total_mass)
  size_breaks = seq(min(df_2D$prob[df_2D$prob > threshold]), max(df_2D$prob), length.out = 10)
  plot_2D = ggplot(df_2D, aes(x = log(x), y = log(y))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", size = 1) +
    geom_point(aes(size = ifelse(prob < threshold, NA, prob)), color = "blue", alpha = 0.8) +
    scale_size_continuous(range = c(0.2, 5), breaks = size_breaks) +
    theme_minimal() +
    xlim(-2.5, 3) +
    labs(title = title,
         x = expression(log({sigma[iA]}^2)),
         y = expression(log({sigma[iB]}^2))) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title = element_text(size = 35),
          axis.text = element_text(size = 28),
          legend.position = 'none')
  return (plot_2D)
}

plotter_1D = function(df_1D, title = '', x_range = NULL) {
  plot_1D = ggplot(df_1D, aes(x = log(x), y = 0)) +
    geom_segment(aes(xend = log(x), yend = prob), linewidth = 1, color = "blue") +
    scale_y_continuous(name = "Density", limits = c(0, max(df_1D$prob))) +
    theme_minimal() +
    labs(title = title, x = expression(log(lambda[i]))) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 15))
  
  if (!is.null(x_range)) {
    plot_1D = plot_1D + coord_cartesian(xlim = x_range)
  }
  
  return(plot_1D)
}

plotter_pvalue_histogram = function(p_values, nbins = 50, ymax = 5500) {
  
  brks = seq(0, 1, length.out = nbins)
  
  ggplot(data.frame(p = p_values), aes(x = p)) +
    geom_histogram(breaks = brks, closed = "right", fill = "grey80", color = "black", alpha = 0.7) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, ymax)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "p-value", y = "Frequency") +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text  = element_text(size = 14)
    )
}