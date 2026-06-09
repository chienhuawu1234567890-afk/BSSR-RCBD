############################################################
## Balanced Fixed-Block RCBD Simulation
## Blinded SSR with Naive, Adjusted Blinded, and Unblinded
############################################################

rm(list = ls())

set.seed(20260606)

############################################################
## 1. Basic functions
############################################################

power_rcbd <- function(N, g, b, delta2, sigma2, alpha = 0.05) {
  df1 <- g - 1
  df2 <- N - g - b + 1
  
  if (df2 <= 0) return(NA_real_)
  
  lambda <- N * delta2 / sigma2
  fcrit <- qf(1 - alpha, df1 = df1, df2 = df2)
  
  1 - pf(fcrit, df1 = df1, df2 = df2, ncp = lambda)
}

find_sample_size <- function(g, b, delta2, sigma2,
                             alpha = 0.05,
                             target_power = 0.80,
                             N_min = NULL,
                             N_max = 2000) {
  if (is.null(N_min)) {
    N_min <- g * b
  }
  
  for (N in N_min:N_max) {
    pow <- power_rcbd(N, g, b, delta2, sigma2, alpha)
    if (!is.na(pow) && pow >= target_power) {
      return(N)
    }
  }
  
  return(NA_integer_)
}

generate_balanced_fixed_rcbd <- function(g, b, n0, mu, tau, beta, sigma2) {
  dat <- expand.grid(
    trt = 1:g,
    block = 1:b,
    rep = 1:n0
  )
  
  dat$y <- mu +
    tau[dat$trt] +
    beta[dat$block] +
    rnorm(nrow(dat), mean = 0, sd = sqrt(sigma2))
  
  dat
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

sigma2_vec <- c(2, 4, 6)
theta_vec <- c(0.25, 0.50, 0.75)

delta2 <- mean(tau^2)

############################################################
## 3. Main simulation
############################################################

results <- data.frame()

for (sigma2_true in sigma2_vec) {
  
  N_plan <- find_sample_size(
    g = g,
    b = b,
    delta2 = delta2,
    sigma2 = sigma2_true,
    alpha = alpha,
    target_power = target_power
  )
  
  n0_plan <- ceiling(N_plan / (g * b))
  N_plan_balanced <- g * b * n0_plan
  
  for (theta in theta_vec) {
    
    n0_star <- ceiling(theta * n0_plan)
    N_star <- g * b * n0_star
    
    sim_store <- data.frame(
      method = character(),
      sigma_hat = numeric(),
      N_new = numeric(),
      N_final = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (s in 1:Nsim) {
      
      dat <- generate_balanced_fixed_rcbd(
        g = g,
        b = b,
        n0 = n0_star,
        mu = mu,
        tau = tau,
        beta = beta,
        sigma2 = sigma2_true
      )
      
      SSTO_star <- compute_ssto(dat$y)
      
      ####################################################
      ## Naive blinded estimator
      ####################################################
      
      sigma2_NB <- SSTO_star / (N_star - 1)
      
      ####################################################
      ## Unblinded estimator
      ####################################################
      
      SSE_star <- compute_sse_unblinded(dat)
      sigma2_UB <- SSE_star / (N_star - g - b + 1)
      
      ####################################################
      ## Adjusted blinded estimator: balanced fixed-block
      ####################################################
      
      Q_struct_BF <- b * n0_star * sum(tau^2) +
        g * n0_star * sum(beta^2)
      
      sigma2_AB_BF <- (SSTO_star - Q_struct_BF) / (N_star - 1)
      
      ## Optional truncation to avoid negative variance estimates
      sigma2_AB_BF <- max(sigma2_AB_BF, 1e-8)
      
      ####################################################
      ## Sample size re-estimation for each estimator
      ####################################################
      
      est_list <- list(
        NB = sigma2_NB,
        AB_BF = sigma2_AB_BF,
        UB = sigma2_UB
      )
      
      for (method_name in names(est_list)) {
        
        sigma2_hat <- est_list[[method_name]]
        
        N_new <- find_sample_size(
          g = g,
          b = b,
          delta2 = delta2,
          sigma2 = sigma2_hat,
          alpha = alpha,
          target_power = target_power,
          N_min = g * b
        )
        
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
    summary_tab$N_plan <- N_plan_balanced
    summary_tab$n0_plan <- n0_plan
    summary_tab$n0_star <- n0_star
    summary_tab$N_star <- N_star
    
    results <- rbind(results, summary_tab)
  }
}

############################################################
## 4. Display results
############################################################

results <- results[
  order(results$sigma2_true, results$theta, results$method),
]

print(results, row.names = FALSE)

############################################################
## 5. Optional: save results
############################################################

write.csv(
  results,
  file = "balanced_fixed_rcbd_simulation_results.csv",
  row.names = FALSE
)
