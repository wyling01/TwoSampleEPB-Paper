# Source EPB.R before this file to load the shared NPMLE and posterior-mass functions.

# P value calculator for weighted VREPB

p_given_lambda_weighted = function(n1, n2, Z1, Z2, S1, S2, W1, W2, lambda) {
  nA = sum(W1)
  nB = sum(W2)
  tau = S1/S2
  c = tau/(tau + nA/nB)
  gamma = lambda/(lambda + nA/nB)
  phi = sqrt((((1-c)/(1-gamma)) * (n2-1) + (c/gamma) * (n1-1))/(n1+n2-2))
  test_stat = ((Z1-Z2)/sqrt(S1/nA + S2/nB))/phi
  p = pt(abs(test_stat), df = n1+n2-2, lower.tail = FALSE) * 2
  return(p)
}

p_value_VREPB_j_weighted = function(n1, n2, Z1, Z2, S1, S2, W1, W2, grid, mass, tau_j) {
  P_value_joint = p_given_lambda_weighted(n1, n2, Z1, Z2, S1, S2, W1, W2, grid)
  post_mass_j = mass_given_tau_j_VR(n1, n2, grid, mass, tau_j)
  P_value_j = sum(P_value_joint * post_mass_j)
  return(P_value_j)
}

P_value_VREPB_weighted = function(info, grid, mass) {
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
  for (i in 1:m) P_value_list_VREPB[i] = p_value_VREPB_j_weighted(n1, n2, Z1_list[i], Z2_list[i], S1_list[i], S2_list[i], W1_matrix[i, ], W2_matrix[i, ], grid, mass, tau[i])
  return(P_value_list_VREPB)
}

# P value calculator for weighted DVEPB

p_given_var1_var2_weighted = function(n1, n2, Z1, Z2, var1, var2, W1, W2) {
  nA = sum(W1)
  nB = sum(W2)
  test_stat = as.numeric(unlist((Z1-Z2)/sqrt(var1/nA + var2/nB)))
  p = pnorm(abs(test_stat), lower.tail = FALSE) * 2
  return(p)
}

p_value_DVEPB_j_weighted = function(n1, n2, Z1, Z2, grid, mass, s1_j, s2_j, W1, W2) {
  P_value_joint_list = p_given_var1_var2_weighted(n1, n2, Z1, Z2, grid$var1, grid$var2, W1, W2)
  post_mass_j = mass_given_tau_j_DV(grid$var1, grid$var2, mass, s1_j, s2_j, n1, n2)
  P_value_j = sum(P_value_joint_list * post_mass_j)
  return(P_value_j)
}

P_value_DVEPB_weighted = function(info, grid, mass) {
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
  for (i in 1:m) P_value_list_DVEPB[i] = p_value_DVEPB_j_weighted(n1, n2, Z1_list[i], Z2_list[i], grid, mass, S1_list[i], S2_list[i], W1_matrix[i, ], W2_matrix[i, ])
  return(P_value_list_DVEPB)
}

# P value calculation for weighted Welch approximation

Welch_approximation_p_weighted = function(n1, n2, Z1, Z2, S1, S2, W1, W2) {
  nA = sum(W1)
  nB = sum(W2)
  dfw = (S1/nA + S2/nB)^2/((S1/nA)^2/(n1-1) + (S2/nB)^2/(n2-1))
  tw = (Z1-Z2)/sqrt(S1/nA + S2/nB)
  p = pt(abs(tw), df = dfw, lower.tail = FALSE) * 2
  return(p)
}

P_value_Welch_weighted = function(info) {
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
  for (i in 1:m) P_value_list_Welch[i] = Welch_approximation_p_weighted(n1, n2, Z1_list[i], Z2_list[i], S1_list[i], S2_list[i], W1_matrix[i, ], W2_matrix[i, ])
  return(P_value_list_Welch)
}

# P value calculation for weighted equal-variance t-test

EV_p_weighted = function(n1, n2, Z1, Z2, S1, S2, W1, W2) {
  nA = sum(W1)
  nB = sum(W2)
  Spool = sqrt(((n1-1) * S1 + (n2-1) * S2)/(n1+n2-2))
  t = (Z1-Z2)/(Spool * sqrt(1/nA + 1/nB))
  p = pt(abs(t), df = n1+n2-2, lower.tail = FALSE) * 2
  return(p)
}

P_value_EV_weighted = function(info) {
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  W1_matrix = info$W1_matrix
  W2_matrix = info$W2_matrix
  P_value_list_EV = rep(0, m)
  for (i in 1:m) P_value_list_EV[i] = EV_p_weighted(n1, n2, Z1_list[i], Z2_list[i], S1_list[i], S2_list[i], W1_matrix[i, ], W2_matrix[i, ])
  return(P_value_list_EV)
}

# P value calculation for weighted Behrens-Fisher test

BF_p_weighted = function(n1, n2, Z1, Z2, S1, S2, W1, W2) {
  nA = sum(W1)
  nB = sum(W2)
  ifunc = function(x) pt((Z1-Z2 + sqrt(S1/nA) * x)/sqrt(S2/nB), n2-1) * dt(x, n1-1)
  p.L = integrate(ifunc, -Inf, Inf)$value
  p = min(p.L, 1-p.L) * 2
  return(p)
}

P_value_BF_weighted = function(info) {
  n1 = info$n1
  n2 = info$n2
  m = info$m
  Z1_list = info$Z1_list
  Z2_list = info$Z2_list
  S1_list = info$S1_list
  S2_list = info$S2_list
  W1_matrix = info$W1_matrix
  W2_matrix = info$W2_matrix
  P_value_list_BF = rep(0, m)
  for (i in 1:m) P_value_list_BF[i] = BF_p_weighted(n1, n2, Z1_list[i], Z2_list[i], S1_list[i], S2_list[i], W1_matrix[i, ], W2_matrix[i, ])
  return(P_value_list_BF)
}

# Extract weighted sufficient statistics

information_extractor_weighted = function(X1, X2, W1, W2) {
  n1 = ncol(X1)
  n2 = ncol(X2)
  m = nrow(X1)
  Z1_list = rowSums(X1 * W1)/rowSums(W1)
  Z2_list = rowSums(X2 * W2)/rowSums(W2)
  S1_list = rowSums(W1 * (X1-Z1_list)^2)/(n1-1)
  S2_list = rowSums(W2 * (X2-Z2_list)^2)/(n2-1)
  information = list('n1' = n1, 'n2' = n2, 'm' = m, 'Z1_list' = Z1_list, 'Z2_list' = Z2_list, 'S1_list' = S1_list, 'S2_list' = S2_list, 'W1_matrix' = W1, 'W2_matrix' = W2)
  return(information)
}

# Solver
solver_weighted = function(X1, X2, W1, W2, VR_parameter, DV_parameter, algorithm_list) {
  algorithm_name = c("VREPB", "DVEPB", "Welch", "B_F", "EV")
  print("Weighted Solver Starts")
  
  information = information_extractor_weighted(X1, X2, W1, W2)
  
  P_list_VR = NA
  P_list_DV = NA
  P_list_Welch = NA
  P_list_BF = NA
  P_list_EV = NA
  VR_NPMLE_result = NA
  DV_NPMLE_result = NA
  
  for (code in algorithm_list) {
    if (code == 1) {
      print(paste("start of", algorithm_name[code]))
      VR_NPMLE_result = VR_NPMLE(information, B=VR_parameter[1], lower_quantile=VR_parameter[2], upper_quantile=VR_parameter[3])
      P_list_VR = P_value_VREPB_weighted(information, VR_NPMLE_result$grid, VR_NPMLE_result$mass)
    }
    if (code == 2) {
      print(paste("start of", algorithm_name[code]))
      DV_NPMLE_result = DV_NPMLE(information, B1=DV_parameter[1], B2=DV_parameter[2], lower_quantile=DV_parameter[3], upper_quantile=DV_parameter[4])
      P_list_DV = P_value_DVEPB_weighted(information, DV_NPMLE_result$grid, DV_NPMLE_result$mass)
    }
    if (code == 3) {
      print(paste("start of", algorithm_name[code]))
      P_list_Welch = P_value_Welch_weighted(information)
    }
    if (code == 4) {
      print(paste("start of", algorithm_name[code]))
      P_list_BF = P_value_BF_weighted(information)
    }
    if (code == 5) {
      print(paste("start of", algorithm_name[code]))
      P_list_EV = P_value_EV_weighted(information)
    }
  }
  
  print("Weighted Solver Ends")
  print("")
  
  output = list("VREPB"=P_list_VR, "DVEPB"=P_list_DV, "Welch"=P_list_Welch, "B_F"=P_list_BF, "EV"=P_list_EV)
  if (is.list(VR_NPMLE_result)) {
    output$"1D_grid" = VR_NPMLE_result$grid
    output$"1D_mass" = VR_NPMLE_result$mass
  }
  if (is.list(DV_NPMLE_result)) {
    output$"2D_grid" = DV_NPMLE_result$grid
    output$"2D_mass" = DV_NPMLE_result$mass
  }
  return(output)
}
