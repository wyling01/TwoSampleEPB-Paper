library(REBayes)
library(Rmosek)
library(asht)

my_BH = function(P, alpha) {
  n = length(P)
  sorted_index = order(P, decreasing = FALSE)
  k_list = which(P[sorted_index][1:n] <= (1:n) * alpha / n)
  if (length(k_list) == 0) {
    return (c())
  }
  else {
    k = max(k_list)
    discovery = sorted_index[1:k]
    return (discovery)
  }
}

# P value calculator for VREPB

p_tau_j_given_lambda = function(n1, n2, tau_j, lambda) {
  out = (1 / lambda) * (1 / beta((n1-1)/2, (n2-1)/2)) * ((n1-1)/(n2-1))^((n1-1)/2) * (tau_j/lambda) ^ ((n1-3)/2) * ((n1-1)/(n2-1)*tau_j/lambda + 1)^(-(n1+n2-2)/2)
}

p_given_lambda = function(n1, n2, Z1, Z2, S1, S2, lambda) {
  
  upper = (Z1 - Z2) * sqrt(n1 + n2 - 2)
  lower = sqrt((n1-1) * (n2 * lambda + n1) / (n1 * n2 * lambda) * S1 + (n2-1) * (n2 * lambda + n1) / (n1 * n2) * S2)
  
  test_stat = upper/lower
  
  p = pt(q = abs(test_stat), df = n1 + n2 - 2, lower.tail = FALSE) * 2
  
  out = p
}

mass_given_tau_j_VR = function(n1, n2, grid, mass, tau_j) {
  f_tau_j_given_grid = p_tau_j_given_lambda(n1, n2, tau_j, grid)
  f_tau_j_grid = f_tau_j_given_grid * mass
  post_mass = f_tau_j_grid / sum(f_tau_j_grid)
  return(post_mass)
}

p_value_VREPB_j = function(n1, n2, Z1, Z2, S1, S2, grid, mass, tau_j) {
  P_value_joint = p_given_lambda(n1, n2, Z1, Z2, S1, S2, grid)
  post_mass_j = mass_given_tau_j_VR(n1, n2, grid, mass, tau_j)
  P_value_j = sum(P_value_joint * post_mass_j)
  return (P_value_j)
}

VR_NPMLE = function(info, B, lower_quantile, upper_quantile) {
  m = info$m
  n1 = info$n1
  n2 = info$n2
  S1_list = info$S1_list
  S2_list = info$S2_list
  
  tau = S1_list/S2_list

  lower = quantile(tau, lower_quantile)
  upper = quantile(tau, upper_quantile)
  log_u = seq(log(lower), log(upper), length = B)
  u = exp(log_u)
  d = rep(1,B)
  w = rep(1, m)/m
  A = outer(tau, u, FUN = p_tau_j_given_lambda, n1 = n1, n2 = n2)
  row_max = apply(A, 1, max)
  stopifnot(all(is.finite(row_max)), all(row_max > 0))
  A_scaled = sweep(A, 1, row_max, "/")
  
  result = KWPrimal(A_scaled, d, w)
  mass = result$f/sum(result$f)
  
  output = list('grid' = u, 'mass' = mass)
  
}

P_value_VREPB = function(info, grid, mass) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  
  tau = S1_list/S2_list
  
  P_value_list_VREPB = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_VREPB[i] = p_value_VREPB_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i], grid, mass, tau[i])
  }
  
  return(P_value_list_VREPB)
}

# P value calculator for DVEPB

p_s_j_given_sigma2 = function(n1, n2, s1_j, s2_j, var1, var2) {
  out = ((n1-1)/var1) * 1 / (2^((n1-1)/2) * gamma((n1-1)/2)) * ((n1-1) * s1_j/var1)^((n1-3)/2) * exp((-1/2) * (n1-1) * s1_j/var1) * ((n2-1)/var2) * 1 / (2^((n2-1)/2) * gamma((n2-1)/2)) * ((n2-1) * s2_j/var2)^((n2-3)/2) * exp((-1/2) * (n2-1) * s2_j/var2)
}

p_given_var1_var2 = function(n1, n2, Z1, Z2, var1, var2) {
  
  upper = (Z1 - Z2)
  lower = sqrt(var1/n1 + var2/n2)
  
  test_stat = as.numeric(unlist(upper/lower))
  p = pnorm(abs(test_stat), lower.tail = FALSE) * 2
  
  out = p
}

mass_given_tau_j_DV = function(var1, var2, mass, s1_j, s2_j, n1, n2) {
  f_s_j_given_grid = p_s_j_given_sigma2(n1, n2, s1_j, s2_j, var1, var2)
  f_s_j_grid = f_s_j_given_grid * mass
  post_mass = f_s_j_grid / sum(f_s_j_grid)
  return(post_mass)
}

p_value_DVEPB_j = function(n1, n2, Z1, Z2, grid, mass, s1_j, s2_j) {
  P_value_joint_list = p_given_var1_var2(n1, n2, Z1, Z2, grid$var1, grid$var2)
  post_mass_j = mass_given_tau_j_DV(grid$var1, grid$var2, mass, s1_j, s2_j, n1, n2)
  P_value_j = sum(P_value_joint_list * post_mass_j)
  return (P_value_j)
}

DV_NPMLE = function(info, B1, B2, lower_quantile, upper_quantile) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  
  lower1 = quantile(S1_list, lower_quantile)
  upper1 = quantile(S1_list, upper_quantile)
  log_u1 = seq(log(lower1), log(upper1), length = B1)
  u1 = exp(log_u1)
  
  lower2 = quantile(S2_list, lower_quantile)
  upper2 = quantile(S2_list, upper_quantile)
  log_u2 = seq(log(lower2), log(upper2), length = B2)
  u2 = exp(log_u2)
  
  var_df = data.frame('var1' = rep(u1, each = B2), 'var2' = rep(u2, B1))
  
  d = rep(1, B1 * B2)
  w = rep(1, m)/m
  
  A = matrix(0, nrow = m, ncol = B1 * B2)
  for (i in 1:(B1*B2)) {
    A[, i] = p_s_j_given_sigma2(n1, n2, S1_list, S2_list, var_df[i, 1], var_df[i, 2])
  }
  
  row_max = apply(A, 1, max)
  stopifnot(all(is.finite(row_max)), all(row_max > 0))
  A_scaled = sweep(A, 1, row_max, "/")
  
  result = KWPrimal(A_scaled, d, w)
  mass =  result$f/sum(result$f)
  
  output = list('grid' = var_df, 'mass' = mass)
  
  return (output)
  
}

P_value_DVEPB = function(info, grid, mass) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  
  P_value_list_DVEPB = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_DVEPB[i] = p_value_DVEPB_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], grid, mass, S1_list[i], S2_list[i])
  }
  
  return (P_value_list_DVEPB)
}

# P value calculation for Welch Approximation

Welch_approximation_p = function(n1, n2, Z1, Z2, S1, S2) {
  dfw = (S1/n1 + S2/n2)^2/(1 / (n1-1) * (S1/n1)^2 + 1 / (n2-1) * (S2/n2)^2)
  se2 = S1/n1 + S2/n2
  tw = (Z1 - Z2)/sqrt(se2)
  p = pt(q = abs(tw), df = dfw, lower.tail = FALSE) * 2
  return (p)
}

P_value_Welch = function (info) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  
  P_value_list_Welch = rep(0, m)
  for (i in 1:m) {
    P_value_list_Welch[i] = Welch_approximation_p(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i])
  }
  
  return (P_value_list_Welch)
}

# P value calculation for Beherens-Fisher Test

P_value_BF = function(m, X1, X2) {
  
  P_list_BF = rep(0, m)
  
  for (i in c(1:m)) {
    X1_i = X1[i, ]
    X2_i = X2[i, ]
    P_i = bfTest(X1_i, X2_i)$p.value
    P_list_BF[i] = P_i
  }
  
  return (P_list_BF)
}

# P value calculation for Equal Variance t-test

P_value_EV_test = function (m, X1, X2) {
  
  P_value_list_EV_test = rep(0, m)
  for (i in 1:m) {
    X1_i = X1[i, ]
    X2_i = X2[i, ]
    P_value_list_EV_test[i] = t.test(X1_i, X2_i, var.equal = TRUE)$p.value
  }
  
  return (P_value_list_EV_test)
}

# extract sufficient statistics

information_extractor = function(X1, X2) {
  
  n1 = ncol(X1)
  n2 = ncol(X2)
  m = nrow(X1)
  
  Z1_list = rowMeans(X1)
  Z2_list = rowMeans(X2)
  S1_list = apply(X1,1,var)
  S2_list = apply(X2,1,var)
  
  information = list('n1' = n1, 'n2' = n2, 'm' = m, 'Z1_list' = Z1_list, 'Z2_list' = Z2_list, 'S1_list' = S1_list, 'S2_list' = S2_list)
  
  return (information)
}

# solver

solver = function(X1, X2, VR_parameter, DV_parameter, algorithm_list) {
  
  algorithm_name = c('VREPB', 'DVEPB', 'Welch', 'B_F', 'EV_test')
  
  print('Solver Starts')
  
  information = information_extractor(X1, X2)
  
  P_list_VR = NA
  P_list_DV = NA
  P_list_Welch = NA
  P_list_BF = NA
  P_list_EV = NA
  
  VR_NPMLE_result = NA
  DV_NPMLE_result = NA
  
  for (code in algorithm_list) {
    
    if (code == 1) {
      print(paste('start of', algorithm_name[code]))
      VR_NPMLE_result = VR_NPMLE(information,B = VR_parameter[1], lower_quantile = VR_parameter[2], upper_quantile = VR_parameter[3])
      P_list_VR = P_value_VREPB(information, VR_NPMLE_result$grid,VR_NPMLE_result$mass)
    }
    
    if (code == 2) {
      print(paste('start of', algorithm_name[code]))
      DV_NPMLE_result = DV_NPMLE(information, B1 = DV_parameter[1], B2 = DV_parameter[2], lower_quantile = DV_parameter[3], upper_quantile = DV_parameter[4])
      P_list_DV = P_value_DVEPB(information, DV_NPMLE_result$grid, DV_NPMLE_result$mass)
    }
    
    if (code == 3) {
      print(paste('start of', algorithm_name[code]))
      P_list_Welch = P_value_Welch(information)
    }
    
    if (code == 4) {
      print(paste('start of', algorithm_name[code]))
      P_list_BF = P_value_BF(information$m, X1, X2)
    }
    
    if (code == 5) {
      print(paste('start of', algorithm_name[code]))
      P_list_EV = P_value_EV_test(information$m, X1, X2)
    }
    
  }
  
  print('Solver Ends')
  print('')
  
  output = list('VREPB' = P_list_VR, 'DVEPB' = P_list_DV, 'Welch' = P_list_Welch, 'B_F' = P_list_BF, 'EV_test' = P_list_EV)
  
  if (is.list(VR_NPMLE_result)) {
    output$'1D_grid' = VR_NPMLE_result$grid
    output$'1D_mass' = VR_NPMLE_result$mass
  }
  
  if (is.list(DV_NPMLE_result)) {
    output$'2D_grid' = DV_NPMLE_result$grid
    output$'2D_mass' = DV_NPMLE_result$mass
  }
  
  return(output)
}
