############################################################
## Unbalanced Random-Block RCBD Simulation
## Naive, Adjusted Blinded, and Unblinded Variance Estimators
############################################################

rm(list = ls())
set.seed(20260606)

############################################################
## 1. Helper functions
############################################################

power_rcbd_unbalanced <- function(df1, df2, lambda, alpha = 0.05) {
  if (df2 <= 0) return(NA_real_)
  fcrit <- qf(1 - alpha, df1 = df1, df2 = df2)
  1 - pf(fcrit, df1 = df1, df2 = df2, ncp = lambda)
}

make_unbalanced_random_data <- function(g, b, n_i, mu, tau,
                                        sigma2, sigmaB2) {
  B <- rnorm(b, mean = 0, sd = sqrt(sigmaB2))
  
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
  
  dat$B <- B[dat$block]
  dat$y <- mu +
    tau[dat$trt] +
    dat$B +
    rnorm(nrow(dat), mean = 0, sd = sqrt(sigma2))
  
  dat
}

make_design_matrix <- function(dat) {
  model.matrix(~ factor(trt) + factor(block), data = dat)
}

contrast_matrix_treatment <- function(g, b) {
  p <- 1 + (g - 1) + (b - 1)
  C <- matrix(0, nrow = g - 1, ncol = p)
  C[, 2:g] <- diag(g - 1)
  C
}

lambda_unbalanced <- function(dat, tau, g, b, sigma2) {
  X <- make_design_matrix(dat)
  C <- contrast_matrix_treatment(g, b)
  
  psi <- c(
    intercept = 0,
    tau[2:g] - tau[1],
    rep(0, b - 1)
  )
  
  middle <- C %*% solve(t(X) %*% X) %*% t(C)
  num <- t(C %*% psi) %*% solve(middle) %*% (C %*% psi)
  
  as.numeric(num / sigma2)
}

find_sample_size_unbalanced <- function(g, b, ratio, tau, sigma2,
                                        alpha = 0.05,
                                        target_power = 0.80,
                                        n1_min = 1,
                                        n1_max = 500) {
  df1 <- g - 1
  
  for (n1 in n1_min:n1_max) {
    n_i <- ratio * n1
    
    dat_plan <- make_unbalanced_random_data(
      g = g,
      b = b,
      n_i = n_i,
      mu = 0,
      tau = tau,
      sigma2 = sigma2,
      sigmaB2 = 0
    )
    
    X <- make_design_matrix(dat_plan)
    df2 <- nrow(dat_plan) - qr(X)$rank
    
    lambda <- lambda_unbalanced(
      dat = dat_plan,
      tau = tau,
      g = g,
      b = b,
      sigma2 = sigma2
    )
    
    pow <- power_rcbd_unbalanced(
      df1 = df1,
      df2 = df2,
      lambda = lambda,
      alpha = alpha
    )
    
    if (!is.na(pow) && pow >= target_power) {
      return(list(
        n1 = n1,
        n_i = n_i,
        N = nrow(dat_plan),
        df2 = df2,
        lambda = lambda,
        power = pow
      ))
    }
  }
  
  NULL
}

compute_ssto <- function(y) {
  sum((y - mean(y))^2)
}

compute_sse_unblinded <- function(dat) {
  fit <- lm(y ~ factor(trt) + factor(block), data = dat)
  sum(resid(fit)^2)
}

estimate_sigmaB_unbalanced <- function(dat) {
  block_sizes <- as.numeric(table(dat$block))
  b <- length(block_sizes)
  N_star <- nrow(dat)
  
  block_means <- aggregate(y ~ block, data = dat, mean)
  grand_mean <- mean(dat$y)
  
  SS_B <- sum(
    block_sizes *
      (block_means$y - grand_mean)^2
  )
  
  SS_W <- sum(
    (dat$y - ave(dat$y, dat$block))^2
  )
  
  MS_B <- SS_B / (b - 1)
  MS_W <- SS_W / (N_star - b)
  
  m0 <- (N_star - sum(block_sizes^2) / N_star) / (b - 1)
  
  sigmaB2_hat <- (MS_B - MS_W) / m0
  
  max(sigmaB2_hat, 0)
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

ratio <- c(1, 1, 3, 3)

sigmaB2_true <- 4
sigma2_vec <- c(2, 4, 6)
theta_vec <- c(0.25, 0.50, 0.75)

############################################################
## 3. Main simulation
############################################################

results <- data.frame()

for (sigma2_true in sigma2_vec) {
  
  plan <- find_sample_size_unbalanced(
    g = g,
    b = b,
    ratio = ratio,
    tau = tau,
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
      sigmaB_hat = numeric(),
      N_new = numeric(),
      N_final = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (s in 1:Nsim) {
      
      dat <- make_unbalanced_random_data(
        g = g,
        b = b,
        n_i = n_i_star,
        mu = mu,
        tau = tau,
        sigma2 = sigma2_true,
        sigmaB2 = sigmaB2_true
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
      ## Blinded block variance estimator
      ####################################################
      
      sigmaB2_hat <- estimate_sigmaB_unbalanced(dat)
      
      ####################################################
      ## Adjusted blinded estimator: unbalanced random-block
      ####################################################
      
      n_i_dot_star <- as.numeric(table(dat$trt))
      n_dot_j_star <- as.numeric(table(dat$block))
      
      tau_bar_plan <- sum(n_i_dot_star * tau) / N_star
      
      Q_Trt_UR <- sum(
        n_i_dot_star *
          (tau - tau_bar_plan)^2
      )
      
      Q_Block_UR <- sigmaB2_hat *
        (
          N_star -
            sum(n_dot_j_star^2) / N_star
        )
      
      sigma2_AB_UR <- (SSTO_star - Q_Trt_UR - Q_Block_UR) /
        (N_star - 1)
      
      sigma2_AB_UR <- max(sigma2_AB_UR, 1e-8)
      
      ####################################################
      ## Sample size re-estimation
      ####################################################
      
      est_list <- list(
        NB = sigma2_NB,
        AB_UR = sigma2_AB_UR,
        UB = sigma2_UB
      )
      
      for (method_name in names(est_list)) {
        
        sigma2_hat <- est_list[[method_name]]
        
        new_plan <- find_sample_size_unbalanced(
          g = g,
          b = b,
          ratio = ratio,
          tau = tau,
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
            sigmaB_hat = sigmaB2_hat,
            N_new = N_new,
            N_final = N_final
          )
        )
      }
    }
    
    summary_tab <- aggregate(
      cbind(sigma_hat, sigmaB_hat, N_new, N_final) ~ method,
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
    summary_tab$sigmaB2_true <- sigmaB2_true
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
  file = "unbalanced_random_rcbd_simulation_results.csv",
  row.names = FALSE
)
