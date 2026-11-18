source('../EPB main/EPB.R')
library(future)
library(future.apply)
library(progressr)
# ----Power, FDR----

Power <- function(discovery, flag_list) {
  alternative_count <- sum(flag_list == 1)
  sum(flag_list[discovery] == 1) / alternative_count
}

FDP <- function(discovery, flag_list) {
  if (length(discovery) == 0) {
    return(0)
  }
  
  sum(flag_list[discovery] == 0) / length(discovery)
}

# ----Data Generator----

data_generator = function(n1, n2, data_generation_parameter, var_struct) {
  
  m = data_generation_parameter$m
  mu1 = data_generation_parameter$mu1
  mu2 = data_generation_parameter$mu2
  mean_var2 = data_generation_parameter$mean_var2
  var_var2 = data_generation_parameter$var_var2
  pi0 = data_generation_parameter$pi0
  mu0 = data_generation_parameter$mu0
  
  null_count = as.integer(pi0 * m)
  dis_count = m - null_count
  flag_list = c(rep(0, null_count), rep(1, dis_count)) # 0 for null, 1 for alternative
  
  var2 = abs(rnorm(m, mean = mean_var2, sd = sqrt(var_var2)))
  
  if (var_struct == 0) {
    # unequal variance simulation with scaled F distribution
    k = data_generation_parameter$k
    d1 = data_generation_parameter$d1
    d2 = data_generation_parameter$d2
    lambda = k * rf(m, df1 = d1, df2 = d2)
  }
  else if (var_struct == 1) {
    # equal variance
    lambda = rep(1, m)
  }
  else if (var_struct == 2) {
    # diffused uninformative distribution
    a = data_generation_parameter$a
    b = data_generation_parameter$b
    lambda = exp(runif(m, min = a, max = b))
  }
  
  var1 = lambda * var2
  
  X1 = matrix(0, m, n1)
  X2 = matrix(0, m, n2)
  
  for (j in 1:m) {
    if (j <= null_count) {
      X1[j, ] = rnorm(n1, mean = mu0, sd = sqrt(var1[j]))
      X2[j, ] = rnorm(n2, mean = mu0, sd = sqrt(var2[j]))
    }
    else {
      mu1_j = rnorm(1, mean = 0, sd = sqrt(mu1 * var1[j]))
      X1[j, ] = rnorm(n1, mean = mu1_j, sd = sqrt(var1[j]))
      X2[j, ] = rnorm(n2, mean = mu2, sd = sqrt(var2[j]))
    }
  }
  
  output = list(
    'X1' = X1,
    'X2' = X2,
    'flag_list' = flag_list,
    'lambda' = lambda,
    'var1' = var1,
    'var2' = var2
  )
  
  return(output)
}

# ----Solver----
simulation_solver = function(X1, X2, flag_list, VR_parameter, DV_parameter, alpha, algorithm_list) {
  
  P_value_result = solver(X1,X2,VR_parameter,DV_parameter,algorithm_list,verbose = FALSE)
  
  algorithm_name = c('VREPB', 'DVEPB', 'Welch', 'B_F', 'EV')
  
  Power_result = rep(NA_real_, 5)
  FDP_result = rep(NA_real_, 5)
  
  for (code in algorithm_list) {
    P_list = P_value_result[[algorithm_name[code]]]
    discovery = my_BH(P_list, alpha)
    Power_result[code] = Power(discovery, flag_list)
    FDP_result[code] = FDP(discovery, flag_list)
  }
  
  output = list(
    'Power' = Power_result,
    'FDP' = FDP_result
  )
  
  return(output)
}

# ----Simulator----
dir_name = function(var_struct) {
  if (var_struct == 0) var_structure = 'unequal'
  else if (var_struct == 1) var_structure = 'equal'
  else if (var_struct == 2) var_structure = 'diffuse'
  return(file.path('Simulation_result', var_structure))
}

file_name = function(n1, n2) {
  return(paste0('KA_', n1, '_KB_', n2, '_result.csv'))
}

simulator = function(n1, n2, seed, data_generation_parameter, VR_parameter, DV_parameter, alpha, rounds, algorithm_list, var_struct, num_workers = 1L) {
  
  dir = dir_name(var_struct)
  file = file_name(n1, n2)
  
  if (!dir.exists(dir)) {
    print('Create Result Directory')
    dir.create(dir, recursive = TRUE)
  }
  
  print('Simulation Start')
  print(paste('Result file:', file.path(dir, file)))
  print(paste('Number of workers:', num_workers))
  
  run_one = function(r) {
    set.seed(seed + r)
    output_r = data_generator(n1, n2, data_generation_parameter, var_struct)
    result_r = simulation_solver(output_r$X1, output_r$X2, output_r$flag_list, VR_parameter, DV_parameter, alpha, algorithm_list)
    return(result_r)
  }
  
  if (num_workers == 1L) {
    pb = txtProgressBar(min = 0, max = rounds, style = 3)
    round_results = tryCatch(
      lapply(seq_len(rounds), function(r) {
        result_r = run_one(r)
        setTxtProgressBar(pb, r)
        return(result_r)
      }),
      finally = close(pb)
    )
  } else {
    old_plan = future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = num_workers)
    
    old_progressr_enable = getOption('progressr.enable')
    old_progressr_handlers = progressr::handlers()
    on.exit({
      options(progressr.enable = old_progressr_enable)
      progressr::handlers(old_progressr_handlers)
    }, add = TRUE)
    
    options(progressr.enable = TRUE)
    progressr::handlers(progressr::handler_txtprogressbar(clear = FALSE))
    
    round_results = progressr::with_progress({
      progress = progressr::progressor(steps = rounds)
      progress(amount = 0)
      future.apply::future_lapply(
        seq_len(rounds),
        function(r) {
          result_r = run_one(r)
          progress(message = paste('Round', r))
          return(result_r)
        },
        future.seed = NULL,
        future.chunk.size = 1,
        future.packages = c('REBayes', 'Rmosek', 'asht')
      )
    })
  }
  
  Power_of_algorithms = do.call(cbind, lapply(round_results, function(x) x$Power))
  FDP_of_algorithms = do.call(cbind, lapply(round_results, function(x) x$FDP))
  
  result_df = data.frame(
    'vrepb_fdp' = FDP_of_algorithms[1, ], 'vrepb_power' = Power_of_algorithms[1, ],
    'dvepb_fdp' = FDP_of_algorithms[2, ], 'dvepb_power' = Power_of_algorithms[2, ],
    'welch_fdp' = FDP_of_algorithms[3, ], 'welch_power' = Power_of_algorithms[3, ],
    'b_f_fdp' = FDP_of_algorithms[4, ], 'b_f_power' = Power_of_algorithms[4, ],
    'ev_fdp' = FDP_of_algorithms[5, ], 'ev_power' = Power_of_algorithms[5, ]
  )
  
  write.csv(result_df, file = file.path(dir, file), row.names = FALSE)
  
  print('Simulation Over')
  print(paste('Result saved to:', file.path(dir, file)))
  
  return(result_df)
}


args = commandArgs(trailingOnly = TRUE)

get_argument = function(name, default = NULL) {
  position = match(paste0('--', name), args)
  if (is.na(position)) {
    if (!is.null(default)) return(default)
    stop(paste0('Missing argument: --', name))
  }
  if (position == length(args)) stop(paste0('Missing value for argument: --', name))
  return(args[position + 1])
}

var_struct = as.integer(get_argument('var_struct')) # 0: unequal; 1: equal; 2: diffuse
KA = as.integer(get_argument('KA'))
KB = as.integer(get_argument('KB'))
rounds = as.integer(get_argument('rounds', 50))
num_workers = as.integer(get_argument('num_workers'))

seed = 2
alpha = 0.1
VR_parameter = c(1000, 0, 1.0)
DV_parameter = c(80, 80, 0.01, 1)
algorithm_list = 1:5

data_generation_parameter = data.frame('k' = 2, 'd1' = 8, 'd2' = 12, 'm' = 5000, 'mu1' = 12, 'mu2' = 0, 'mean_var2' = 6, 'var_var2' = 4, 'pi0' = 0.9, 'mu0' = 0)
if (var_struct == 2) data_generation_parameter = data.frame('a' = -5, 'b' = 5, 'm' = 5000, 'mu1' = 12, 'mu2' = 0, 'mean_var2' = 6, 'var_var2' = 4, 'pi0' = 0.9, 'mu0' = 0)

invisible(simulator(KA, KB, seed, data_generation_parameter, VR_parameter, DV_parameter, alpha, rounds, algorithm_list, var_struct, num_workers))
