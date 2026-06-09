############################################################
## Unbalanced Fixed-Block RCBD Simulation
## Naive, Adjusted Blinded, and Unblinded Variance Estimators
############################################################

rm(list = ls())
set.seed(20260606)

############################################################
## 1. Helper functions
############################################################

power_rcbd_unbalanced <- function(N, df1, df2, lambda, alpha = 0.05) {
  if (df2 <= 0) return(NA_real_)
  fcrit <- qf(1 - alpha, df1 = df1, df2 = df2)
  1 - pf(fcrit, df1 = df1, df2 = df2, ncp = lambda)
}

make_unbalanced_fixed_data <- function(g, b, n_i, mu, tau, beta, sigma2) {
  dat <- do.call(
    rbind,
    lapply(1:g, function(i) {
      do.call(
        rbind,
        lapply(1:b, function(j) {
          data.frame(
            trt = i,
            block = j,
            rep = 1:n_i[i]
          )
        })
      )
    })
  )
  
  dat$y <- mu +
    tau[dat$trt] +
    beta[dat$block] +
    rnorm(nrow(dat), mean = 0, sd = sqrt(sigma2))
  
  dat
}

make_design_matrix <- function(dat) {
  model.matrix(~ factor(trt) + factor(block), data = dat)
}

contrast_matrix_treatment <- function(g, b) {
  ## With lm parameterization: intercept + treatment dummies 2:g + block dummies 2:b
  ## Test all treatment dummy coefficients = 0
  p <- 1 + (g - 1) + (b - 1)
  C <- matrix(0, nrow = g - 1, ncol = p)
  C[, 2:g] <- diag(g - 1)
  C
}

lambda_unbalanced_fixed <- function(dat, tau, beta, sigma2, g, b) {
  X <- make_design_matrix(dat)
  C <- contrast_matrix_treatment(g, b)
  
  psi <- c(
    mu = 0,
    tau[2:g] - tau[1],
    beta[2:b] - beta[1]
  )
  
  middle <- C %*% solve(t(X) %*% X) %*% t(C)
  num <- t(C %*% psi) %*% solve(middle) %*% (C %*% psi)
  
  as.numeric(num / sigma2)
}

find_sample_size_unbalanced_fixed <- function(g, b, ratio, tau, beta, sigma2,
                                              alpha = 0.05,
                                              target_power = 0.80,
                                              n1_min = 1,
                                              n1_max = 500) {
  df1 <- g - 1
  
  for (n1 in n1_min:n1_max) {
    n_i <- ratio * n1
    
    dat_plan <- make_unbalanced_fixed_data(
      g = g,
      b = b,
      n_i = n_i,
      mu = 0,
      tau = tau,
      beta = beta,
      sigma2 = sigma2
    )
    
    X <- make_design_matrix(dat_plan)
    N <- nrow(dat_plan)
    df2 <- N - qr(X)$rank
    
    lambda <- lambda_unbalanced_fixed(
      dat = dat_plan,
      tau = tau,
      beta = beta,
      sigma2 = sigma2,
      g = g,
      b = b
    )
    
    pow <- power_rcbd_unbalanced(
      N = N,
      df1 = df1,
      df2 = df2,
      lambda = lambda,
      alpha = alpha
    )
    
    if (!is.na(pow) && pow >= target_power) {
      return(list(
        n1 = n1,
        n_i = n_i,
        N = N,
        df2 = df2,
        lambda = lambda,
        power = pow
      ))
    }
  }
  
  return(NULL)
}

compute_ssto <- function(y) {
  sum((y - mean(y))^2)
}

compute_sse_unblinded <- function(dat) {
  fit <- lm(y ~ factor(trt) + factor(block), data = dat)
  sum(resid(fit)^2)
}

############################################################
## 2. Simulation settings
############################################################

g <- 4
b <- 6

alpha <- 0.05
target_power <- 0.80
Nsim <- 5000

mu <- 10

tau <- c(-1.5, -0.5, 0.5, 1.5)
beta <- c(-2, -1, 0, 0.5, 1, 1.5)

ratio <- c(1, 1, 3, 3)

sigma2_vec <- c(2, 4, 6)
theta_vec <- c(0.25, 0.50, 0.75)

############################################################
## 3. Main simulation
############################################################

results <- data.frame()

for (sigma2_true in sigma2_vec) {
  
  plan <- find_sample_size_unbalanced_fixed(
    g = g,
    b = b,
    ratio = ratio,
    tau = tau,
    beta = beta,
    sigma2 = sigma2_true,
    alpha = alpha,
    target_power = target_power
  )
  
  n_i_plan <- plan$n_i
  N_plan <- plan$N
  
  for (theta in theta_vec) {
    
    n_i_star <- ceiling(theta * n_i_plan)
    N_star <- b * sum(n_i_star)
    
    sim_store <- data.frame(
      method = character(),
      sigma_hat = numeric(),
      N_new = numeric(),
      N_final = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (s in 1:Nsim) {
      
      dat <- make_unbalanced_fixed_data(
        g = g,
        b = b,
        n_i = n_i_star,
        mu = mu,
        tau = tau,
        beta = beta,
        sigma2 = sigma2_true
      )
      
      SSTO_star <- compute_ssto(dat$y)
      X_star <- make_design_matrix(dat)
      rank_X_star <- qr(X_star)$rank
      
      ####################################################
      ## Naive blinded estimator
      ####################################################
      
      sigma2_NB <- SSTO_star / (N_star - 1)
      
      ####################################################
      ## Unblinded estimator
      ####################################################
      
      SSE_star <- compute_sse_unblinded(dat)
      sigma2_UB <- SSE_star / (N_star - rank_X_star)
      
      ####################################################
      ## Adjusted blinded estimator: unbalanced fixed-block
      ####################################################
      
      n_ij_star <- table(dat$trt, dat$block)
      
      eta_plan <- matrix(0, nrow = g, ncol = b)
      for (i in 1:g) {
        for (j in 1:b) {
          eta_plan[i, j] <- tau[i] + beta[j]
        }
      }
      
      eta_bar_plan <- sum(n_ij_star * eta_plan) / N_star
      
      Q_struct_UF <- sum(
        n_ij_star *
          (eta_plan - eta_bar_plan)^2
      )
      
      sigma2_AB_UF <- (SSTO_star - Q_struct_UF) / (N_star - 1)
      sigma2_AB_UF <- max(sigma2_AB_UF, 1e-8)
      
      ####################################################
      ## Sample size re-estimation
      ####################################################
      
      est_list <- list(
        NB = sigma2_NB,
        AB_UF = sigma2_AB_UF,
        UB = sigma2_UB
      )
      
      for (method_name in names(est_list)) {
        
        sigma2_hat <- est_list[[method_name]]
        
        new_plan <- find_sample_size_unbalanced_fixed(
          g = g,
          b = b,
          ratio = ratio,
          tau = tau,
          beta = beta,
          sigma2 = sigma2_hat,
          alpha = alpha,
          target_power = target_power
        )
        
        N_new <- new_plan$N
        N_final <- max(N_star, N_new)
        
        sim_store <- rbind(
          sim_store,
          data.frame(
            method = method_name,
            sigma_hat = sigma2_hat,
            N_new = N_new,
            N_final = N_final
          )
        )
      }
    }
    
    summary_tab <- aggregate(
      cbind(sigma_hat, N_new, N_final) ~ method,
      data = sim_store,
      FUN = mean
    )
    
    summary_tab$bias <- summary_tab$sigma_hat - sigma2_true
    
    mse_tab <- aggregate(
      sigma_hat ~ method,
      data = sim_store,
      FUN = function(x) mean((x - sigma2_true)^2)
    )
    
    names(mse_tab)[2] <- "mse"
    
    summary_tab <- merge(summary_tab, mse_tab, by = "method")
    
    summary_tab$sigma2_true <- sigma2_true
    summary_tab$theta <- theta
    summary_tab$N_plan <- N_plan
    summary_tab$N_star <- N_star
    summary_tab$n_i_plan <- paste(n_i_plan, collapse = ":")
    summary_tab$n_i_star <- paste(n_i_star, collapse = ":")
    
    results <- rbind(results, summary_tab)
  }
}

############################################################
## 4. Display and save results
############################################################

results <- results[
  order(results$sigma2_true, results$theta, results$method),
]

print(results, row.names = FALSE)

write.csv(
  results,
  file = "unbalanced_fixed_rcbd_simulation_results.csv",
  row.names = FALSE
)
