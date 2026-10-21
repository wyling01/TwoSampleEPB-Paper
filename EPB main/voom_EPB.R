library(REBayes)
library(Rmosek)
library(asht)

# P value calculator for weighted VREPB

p_tau_j_given_lambda = function(n1, n2, tau_j, lambda) {
  out = (1 / lambda) * (1 / beta((n1-1)/2, (n2-1)/2)) * ((n1-1)/(n2-1))^((n1-1)/2) * (tau_j/lambda) ^ ((n1-3)/2) * ((n1-1)/(n2-1)*tau_j/lambda + 1)^(-(n1+n2-2)/2)
}

p_given_lambda = function(n1, n2, Z1, Z2, S1, S2, W1, W2, lambda) {
  
  nA = sum(W1)
  nB = sum(W2)
  
  tau = S1/S2
  
  c = tau / (tau + nA/nB)
  gamma = lambda / (lambda + nA/nB)
  
  phi = sqrt((((1-c) / (1-gamma)) * (n2-1) + (c / gamma) * (n1-1)) / (n1 + n2 - 2))
  
  left = (Z1 - Z2) / sqrt(S1/nA + S2/nB)
  
  test_stat = left/phi
  
  p = pt(q = abs(test_stat), df = n1 + n2 - 2, lower.tail = FALSE) * 2
  
  return(p)
}

mass_given_tau_j_VR = function(n1, n2, grid, mass, tau_j) {
  f_tau_j_given_grid = p_tau_j_given_lambda(n1, n2, tau_j, grid)
  f_tau_j_grid = f_tau_j_given_grid * mass
  post_mass = f_tau_j_grid / sum(f_tau_j_grid)
  return(post_mass)
}

p_value_VREPB_j = function(n1, n2, Z1, Z2, S1, S2, W1, W2, grid, mass, tau_j) {
  P_value_joint = p_given_lambda(n1, n2, Z1, Z2, S1, S2, W1, W2, grid)
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
  W1_matrix = info$W1_matrix
  W2_matrix = info$W2_matrix
  
  tau = S1_list/S2_list
  
  P_value_list_VREPB = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_VREPB[i] = p_value_VREPB_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i], W1 = W1_matrix[i, ], W2 = W2_matrix[i, ], grid, mass, tau[i])
  }
  
  return(P_value_list_VREPB)
}

# P value calculator for DVEPB

p_s_j_given_sigma2 = function(n1, n2, s1_j, s2_j, var1, var2) {
  out = ((n1-1)/var1) * 1 / (2^((n1-1)/2) * gamma((n1-1)/2)) * ((n1-1) * s1_j/var1)^((n1-3)/2) * exp((-1/2) * (n1-1) * s1_j/var1) * ((n2-1)/var2) * 1 / (2^((n2-1)/2) * gamma((n2-1)/2)) * ((n2-1) * s2_j/var2)^((n2-3)/2) * exp((-1/2) * (n2-1) * s2_j/var2)
}

p_given_var1_var2 = function(n1, n2, Z1, Z2, var1, var2, W1, W2) {
  
  nA = sum(W1)
  nB = sum(W2)
  
  upper = (Z1 - Z2)
  lower = sqrt(var1/nA + var2/nB)
  
  test_stat = as.numeric(unlist(upper/lower))
  p = pnorm(abs(test_stat), lower.tail = FALSE) * 2
  
  return(p)
}

mass_given_tau_j_DV = function(var1, var2, mass, s1_j, s2_j, n1, n2) {
  f_s_j_given_grid = p_s_j_given_sigma2(n1, n2, s1_j, s2_j, var1, var2)
  f_s_j_grid = f_s_j_given_grid * mass
  post_mass = f_s_j_grid / sum(f_s_j_grid)
  return(post_mass)
}

p_value_DVEPB_j = function(n1, n2, Z1, Z2, grid, mass, s1_j, s2_j, W1, W2) {
  P_value_joint_list = p_given_var1_var2(n1, n2, Z1, Z2, grid$var1, grid$var2, W1, W2)
  post_mass_j = mass_given_tau_j_DV(grid$var1, grid$var2, mass, s1_j, s2_j, n1, n2)
  P_value_j = sum(P_value_joint_list * post_mass_j)
  return (P_value_j)
}

DV_NPMLE = function(S1_list, S2_list, B1, B2, m, n1, n2, lower_quantile, upper_quantile) {
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
  mass = result$f/sum(result$f)
  
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
  W1_matrix = info$W1_matrix
  W2_matrix = info$W2_matrix
  
  P_value_list_DVEPB = rep(0, m)
  for (i in c(1:m)) {
    P_value_list_DVEPB[i] = p_value_DVEPB_j(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], grid, mass, S1_list[i], S2_list[i], W1_matrix[i, ], W2_matrix[i, ])
  }
  
  return (P_value_list_DVEPB)
}

# P value calculation for Welch Approximation

Welch_approximation_p = function(n1, n2, Z1, Z2, S1, S2, W1, W2) {
  nA = sum(W1)
  nB = sum(W2)
  
  dfw = (S1/nA + S2/nB)^2/(1 / (n1-1) * (S1/nA)^2 + 1 / (n2-1) * (S2/nB)^2)
  se2 = S1/nA + S2/nB
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
  W1_matrix = info$W1_matrix
  W2_matrix = info$W2_matrix
  
  P_value_list_Welch = rep(0, m)
  for (i in 1:m) {
    P_value_list_Welch[i] = Welch_approximation_p(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i], W1 = W1_matrix[i, ], W2 = W2_matrix[i, ])
  }
  
  return (P_value_list_Welch)
}

# P value calculation for Beherens-Fisher Test

P_value_BF = function(m, X1, X2) {
  
  P_list = rep(0, m)
  
  for (i in c(1:m)) {
    X1_i = X1[i, ]
    X2_i = X2[i, ]
    P_i = bfTest(X1_i, X2_i)$p.value
    P_list[i] = P_i
  }
  
  return (P_list)
}

# P value calculation for Equal Variance t-test

EV_test_p = function(n1, n2, Z1, Z2, S1, S2, W1, W2) {
  nA = sum(W1)
  nB = sum(W2)
  
  Spool = sqrt(((n1-1) * S1 + (n2-1) * S2) / (n1 + n2 - 2))
  se = Spool * sqrt(1/nA + 1/nB)
  t = (Z1 - Z2)/se
  p = pt(q = abs(t), df = n1 + n2 - 2, lower.tail = FALSE) * 2
  return (p)
}

P_value_EV_test = function (info) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  W1_matrix = info$W1_matrix
  W2_matrix = info$W2_matrix
  
  P_value_list_EV_test = rep(0, m)
  for (i in 1:m) {
    P_value_list_EV_test[i] = EV_test_p(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i], W1 = W1_matrix[i, ], W2 = W2_matrix[i, ])
  }
  
  return (P_value_list_EV_test)
}

# P value calculation for Beherens-Fisher test

BF_test_p = function(n1, n2, Z1, Z2, S1, S2, W1, W2) {
  nA = sum(W1)
  nB = sum(W2)
  
  ifunc = function(x){
    pt( (Z1 - Z2 + (sqrt(S1) / sqrt(nA)) * x) / (sqrt(S2) / sqrt(nB)), n2 - 1) * dt(x, n1 - 1)
  }
  
  p.L = integrate(ifunc, -Inf, Inf)$value
  
  p = min(p.L, 1 - p.L) * 2
  return (p)
}

P_value_BF_test = function (info) {
  
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list 
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  W1_matrix = info$W1_matrix
  W2_matrix = info$W2_matrix
  
  P_value_list_BF_test = rep(0, m)
  for (i in 1:m) {
    P_value_list_BF_test[i] = BF_test_p(n1 = n1, n2 = n2, Z1 = Z1_list[i], Z2 = Z2_list[i], S1 = S1_list[i], S2 = S2_list[i], W1 = W1_matrix[i, ], W2 = W2_matrix[i, ])
  }
  
  return (P_value_list_BF_test)
}

# Extract sufficient statistics

information_extractor = function(X1, X2, W1, W2) {
  
  n1 = ncol(X1)
  n2 = ncol(X2)
  m = nrow(X1)
  
  Z1_list = rowSums(X1 * W1) / rowSums(W1)
  Z2_list = rowSums(X2 * W2) / rowSums(W2)
  S1_list = 1 / (n1-1) * rowSums(W1 * (X1 - Z1_list)^2)
  S2_list = 1 / (n2-1) * rowSums(W2 * (X2 - Z2_list)^2)
  
  information = list('n1' = n1, 'n2' = n2, 'm' = m, 'Z1_list' = Z1_list, 'Z2_list' = Z2_list, 'S1_list' = S1_list, 'S2_list' = S2_list, 'W1_matrix' = W1, 'W2_matrix' = W2)
  
  return (information)
}

# BH

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
