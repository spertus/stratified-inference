# conservative inference on stratified populations with bounded RVs
library(tidyverse)
library(cumstats)
library(confseq)
library(lpSolve)


############# confidence sequences and p-values for unstratified samples ##########
shuffle <- function(x){sample(x, size = length(x), replace = FALSE)}
#from appendix of https://arxiv.org/pdf/2008.08536.pdf 
#an efficient implementation of KMART (techincally, KMART is when prior_alpha = prior_beta = 1)
#tests whether mean(population) > 1/2
#inputs:
#population: a length-N vector with all elements in [0,1], the population
#prior_alpha: the alpha parameter for the prior (mixing) distribution
#prior_beta: the beta parameter for the same
#outputs:
#a length-N vector which is the p-value sequence
kmart_p_value_sequence <- function(population, prior_alpha = 1, prior_beta = 1){
  population <- shuffle(population)
  N <- length(population)
  Y_n <- cumsum(population)
  denominator <- 1 + exp(log(2) * (1:N) + lbeta(Y_n + prior_alpha, 1:N - Y_n + prior_beta) - lbeta(prior_alpha, prior_beta) + log(1 - pbeta(q = 1/2, Y_n + prior_alpha, 1:N - Y_n + prior_beta)) - log(1 - pbeta(q = 1/2, shape1 = prior_alpha, shape2 = prior_beta)))
  a_star <- 1 / denominator
  a_star
}


#the hedged p-value from https://arxiv.org/abs/2010.09686 
#inputs:
#mu_0: the population mean under the null hypothesis
#population: the population we are estimating the mean of assuming SRSWR
#theta: theta parameter determining how much weight to put on the larger/smaller alternative
#log: whether to return the p-value on the log scale or not
#shuffle: is the population already in random order or should it be shuffled (should usually be TRUE)
#last: return the entire sequence of p-values or just the last p-value?
#outputs: 
#a finite-sample, sequentially valid p-value for the hypothesis that the mean(population) == m
hedged_pvalue <- function(population, mu_0 = 1/2, theta = 1, log = FALSE, shuffle = TRUE, last = FALSE){
  if(shuffle){
    population <- shuffle(population)
  }
  N <- length(population)
  alpha <- .05
  mu_hat <- cummean(population)
  lagged_mu_hat <- lag(mu_hat)
  lagged_mu_hat[1] <- 1/2
  v_n <- (population - lagged_mu_hat)^2
  lagged_sigma_hat <- lag(sqrt(cummean((population - mu_hat)^2)))
  lagged_sigma_hat[1:2] <- 1/4
  lambda_sequence_plus <- pmin(sqrt(2 * log(2/alpha) / (log(1:N) * 1:N * lagged_sigma_hat)), .9 / mu_0)
  lambda_sequence_minus <- pmin(sqrt(2 * log(2/alpha) / (log(1:N) * 1:N * lagged_sigma_hat)), .9 / (1-mu_0))
  if(log){
    K_plus <- cumsum(log(1 + lambda_sequence_plus * (population - mu_0)))
    K_minus <- cumsum(log(1 - lambda_sequence_minus * (population - mu_0)))
    K_plusminus <- pmax(log(theta) + K_plus, log(1 - theta) + K_minus)
    p_value <- -K_plusminus
  } else{
    K_plus <- cumprod(1 + lambda_sequence_plus * (population - mu_0))
    K_minus <- cumprod(1 - lambda_sequence_minus * (population - mu_0))
    K_plusminus <- pmax(theta * K_plus, (1 - theta) * K_minus)
    p_value <- 1 / K_plusminus
  }
  
  if(last){
    p_value[length(p_value)]
  } else{
    p_value
  }
}

#hedged confidence interval from p-values
#inputs:
#population: a vector of values representing the (true, fixed) population distribution 
#n: the desired sample size
#alpha: the desired confidence level
#outputs:
#a length-2 vector of lower and upper confidence bounds (respectively)
hedged_CI <- function(population, n, alpha = .05){
  shuffled_pop <- sample(population, n, replace = TRUE)
  interval_center <- optimize(function(x){hedged_pvalue(mu_0 = x, population = shuffled_pop, theta = 0.5, log = TRUE, shuffle = FALSE, last = TRUE)}, interval = c(0,1), maximum = TRUE)$maximum
  min_p <- hedged_pvalue(mu_0 = 0, population = shuffled_pop, theta = 0.5, shuffle = FALSE, last = TRUE)
  max_p <- hedged_pvalue(mu_0 = 1, population = shuffled_pop, theta = 0.5, shuffle = FALSE, last = TRUE)
  if(min_p > .05){
    LB <- 0
  } else{
    LB <- uniroot(function(x){hedged_pvalue(mu_0 = x, population = shuffled_pop, theta = 0.5, shuffle = FALSE, last = TRUE) - alpha}, lower = 0, upper = interval_center)$root
  }
  
  if(max_p > .05){
    UB <- 1
  } else{
    UB <- uniroot(function(x){hedged_pvalue(mu_0 = x, population = shuffled_pop, theta = 0.5, shuffle = FALSE, last = TRUE) - alpha}, lower = interval_center, upper = 1)$root
  }
  c(LB, UB)
}



#prior-posterior / mixture likelihood martingale. The mixing distribution is conjugate, yielding a closed-form posterior
#tests whether mean(population) > mu_0
#inputs:
#population: a length-N vector with all elements in [0,1], the population
#prior_alpha: the alpha parameter for the prior (mixing) distribution
#prior_beta: the beta parameter for the same
#mu_0: the population mean under the null hypothesis
#outputs:
#a length-N vector which is the p-value sequence
beta_binomial_ppm <- function(population, mu_0, prior_alpha, prior_beta, last = FALSE){
  population <- shuffle(population)
  N <- length(population)
  Y_n <- cumsum(population)
  log_posterior <- dbeta(x = mu_0, shape1 = prior_alpha + Y_n, shape2 = prior_beta + 1:N - Y_n, log = TRUE)
  log_prior <- dbeta(x = mu_0, shape1 = prior_alpha, shape2 = prior_beta, log = TRUE)
  log_ppr <- log_prior - log_posterior
  if(last){
    p_value <- min(1, exp(-log_ppr))
  } else{
    p_value <- pmin(1, exp(-log_ppr))
  }
  p_value
}


#binomial test for samples with replacement (not sequentially valid)
#tests whether mean(population) > mu_0
#inputs:
#population: a length-N vector with all elements in {0,1}, the population
#mu_0: the hypothesized mean under the null
#outputs:
#a p-value based on the binomial likelihood
binom_pvalue <- function(population, mu_0 = .5){
  population <- shuffle(population)
  N <- length(population)
  Y_n <- cumsum(population)
  p_value <- pbinom(q = Y_n - 1, size = 1:N, prob = mu_0, lower.tail = FALSE)
  p_value
}

#return the stopping time when the sequence "p_values" falls below "risk_limit"
#if the sequence never falls below risk_limit, return Inf
get_stopping_time <- function(p_values, risk_limit = 0.05){
  if(all(p_values >= risk_limit)){
    Inf
  } else{
    min(which(p_values < risk_limit))
  }
}

########## stratified sampling with empirical Bernstein p-values ##########
# get values necessary for EB within each strata
empirical_Vn <- function(sample){
  means = cummean(sample)
  # lag by 1, with first index being 1/2
  lagged_means <- c(1/2, means[-length(means)])
  # cumulative variance
  cumsum((sample - lagged_means)^2)
}

get_psi_E <- function(lambda){
  (-log(1 - lambda) - lambda)
}
get_psi_H <- function(lambda){
  lambda^2 / 8
}
get_mu_hat_n <- function(Y){
  (1/2 + cumsum(Y)) / (1:length(Y) + 1)
}

get_sigma_hat_squared <- function(Y){
  (1/4 + cumsum(Y - get_mu_hat_n(Y))^2) / (1:length(Y) + 1)
}



########## Exponential martingales for stratified samples ############
#function to get samples from a population (or a stratum of a population). Used in get_stratified_pvalue()
#inputs:
  #pop: the population as a vector of values in [0,1]
  #n: the sample size
  #alpha: the targeted significance level, used to compute optimal lambda
  #method: the martingale method to be used, determines how lambda is computed
#outputs:
  #a list with Y (samples), lambda (tuning parameters for some martingales), psi (normalizing function for some martingales), and V_n (accumulating variance measure for some martingales). May be null depending on method
get_statistics <- function(pop, n, alpha = 0.05, method){
  Y <- sample(pop, n, replace = TRUE)
  if(method == "empirical_bernstein"){
    mu_hat <- get_mu_hat_n(Y)
    mu_hat_lag <- lag(mu_hat, 1)
    mu_hat_lag[1] <- 1/2
    sigma_hat_squared <- (1/4 + cumsum(Y - mu_hat)^2) / (1:n + 1) 
    sigma_hat_squared_lag <- lag(sigma_hat_squared, 1)
    sigma_hat_squared_lag[1] <- 1/4
    
    #optimized for fixed sample size n
    lambda <- pmin(sqrt(2 * log(2/alpha) / (n * sigma_hat_squared_lag)), .9)
    psi <- get_psi_E(lambda)
    V_n <- (Y - mu_hat_lag)^2
  } else if(method == "hoeffding"){
    lambda <- rep(sqrt(8 * log(1/alpha) / n), n)
    psi <- get_psi_H(lambda)
    V_n <- rep(1, n)
  } else if(method == "beta-binomial"){
    lambda <- NULL
    psi <- NULL
    V_n <- NULL
  } else if(method == "hedged"){
    mu_hat <- get_mu_hat_n(Y)
    mu_hat_lag <- lag(mu_hat, 1)
    mu_hat_lag[1] <- 1/2
    sigma_hat_squared <- (1/4 + cumsum(Y - mu_hat)^2) / (1:n + 1) 
    sigma_hat_squared_lag <- lag(sigma_hat_squared, 1)
    sigma_hat_squared_lag[1] <- 1/4
    psi <- NULL
    V_n <- NULL
    #this lambda is not thresholded (b/c threshold depends on mu_0), it has to be thresholded later
    lambda <- abs(sqrt(2 * log(2/alpha) / (n * sigma_hat_squared_lag)))
  }
  list(Y = Y, "lambda" = lambda, "psi" = psi, V_n = V_n)
}



#function to compute a p-value after n samples given a population and strata 
#inputs: 
  #population: a 1-D vector, the population to sample from
  #strata: a 1-D vector of stratum indicators with length equal to length(population)
  #mu_0: the hypothesized overall null mean to be tested
  #n: a vector of length length(unique(strata)), specifying how many samples to draw (at random with replacement) from each stratum
  #method: the martingale to be used, one of hoeffding, empirical_bernstein, hedged, or beta-binomial
  #alpha: the targeted significance level, used in optimizing martingales. Note that resulting p-value is valid for any significance level.
#output:
  #a p-value for the null hypothesis: mean(population) = mu_0
get_stratified_pvalue <- function(population, strata, mu_0, n, method = "hoeffding", alpha = .05){
  if(length(population) != length(strata)){
    stop("Strata do not cover population (length(strata) != length(population))")
  }
  strata_names <- unique(strata)
  strata_sizes <- as.numeric(table(strata))
  K <- length(strata_names)
  a <- prop.table(table(strata))
  statistics_strata <- list()
  for(k in 1:K){
    statistics_strata[[k]] <- get_statistics(population[strata == strata_names[k]], n = n[k], method = method, alpha = alpha)
  }
  
  
  if(method %in% c("hoeffding","empirical_bernstein")){
    lambda_k <- statistics_strata %>%
      map(function(x){x$lambda}) %>%
      reduce(c)
    Y <- statistics_strata %>%
      map(function(x){x$Y}) %>%
      reduce(c)
    constants <- statistics_strata %>%
      map(function(x){sum(x$V_n * x$psi - x$lambda * x$Y)}) %>%
      reduce(c)
    coefficients <- statistics_strata %>%
      map(function(x){sum(x$lambda)}) %>%
      reduce(c) 
    
    constraint_matrix <- rbind(
      #linear equations
      cbind(diag(coefficients, nrow = K, ncol = K), -diag(K)),
      #zero part
      cbind(matrix(0, nrow = K, ncol = K), -diag(K)),
      #constraint on sum of within-stratum means
      cbind(t(a), matrix(0, nrow = 1, ncol = K)),
      #less than 1
      cbind(diag(K), matrix(0, nrow = K, ncol = K)),
      #positivity
      cbind(diag(K), matrix(0, nrow = K, ncol = K))
    )
    #call to lpSolve
    mu_star <- lp(
      direction = "max",
      objective.in = -c(rep(0,K), rep(1, K)),
      const.mat = constraint_matrix,
      const.dir = c(rep("<=", 2*K), "==", rep("<=", K), rep(">=", K)),
      const.rhs = c(-constants, rep(0, K), mu_0, rep(1, K), rep(0, K))
      )$solution[1:K]
    combined_test_stat <- -2 * (sum(pmin(0, constants + coefficients * mu_star)))
    p_value <- pchisq(combined_test_stat, df = 2*K, lower.tail = FALSE)
  
  } else if(method == "hedged"){
    lambda_k <- statistics_strata %>%
      map(function(x){x$lambda}) %>%
      reduce(c)
    Y <- statistics_strata %>%
      map(function(x){x$Y}) %>%
      reduce(c)
    if(K == 1){
      p_value <- min(1, 1 / prod(1 + pmin(lambda_k, .75/mu_0) * (Y - mu_0)))
      } else if(K == 2){
      solution <- optimize(
        f = function(mu_01){
          mu_02 <- (mu_0 - a[1] * mu_01) / a[2]
          2 * max(0, sum(log(1 + pmin(statistics_strata[[1]]$lambda, .75/mu_01) * (statistics_strata[[1]]$Y - mu_01)))) + max(0, sum(log(1 + pmin(statistics_strata[[2]]$lambda, .75/mu_02) * (statistics_strata[[2]]$Y - mu_02))))
        },
        interval = c(0,mu_0/a[1]))
      combined_test_stat <- solution$objective
      p_value <- pchisq(q = combined_test_stat, df = 4, lower.tail = FALSE)
      } else{
      stop("Hedged martingale only works for two or fewer strata right now.")
      }
  } else if(method == "beta-binomial"){
      Y_k <- statistics_strata %>%
        map(function(x){sum(x$Y)}) %>%
        reduce(c)
      C_k <- statistics_strata %>%
        map(function(x){
          #note that alpha = beta = 1 currently, which is a flat prior/mixing distribution
          lgamma(1 + sum(x$Y)) + lgamma(1 + length(x$Y) - sum(x$Y)) - lgamma(2 + length(x$Y))
        }) %>%
        reduce(c)
      
      #Lagrangian is used to enforce equality constraint on weighted sum of null means
      #to enforce constraint that 0 <= mu_0k <= 1
      # constraint_matrix <- rbind(
      #   #positivity
      #   cbind(diag(K), 0),
      #   #less than 1
      #   cbind(-diag(K), 0)
      # )
      # solution <- constrOptim(
      #   theta = rep(mu_0/sum(a), K+1),
      #   f = function(x){
      #     -2 * sum(Y_k * log(x[1:K]) + (n - Y_k) * log(1 - x[1:K]) - C_k) - x[K+1] * (sum(a * x[1:K]) - mu_0)
      #   },
      #   grad = function(x){
      #     c(-2 * (Y_k / x[1:K] + (n - Y_k)/(1 - x[1:K])) - x[K+1] * a, sum(a * x[1:K]) - mu_0)
      #   },
      #   ui = constraint_matrix,
      #   ci = c(rep(0, K), rep(-1, K))
      # )
      
      if(K == 1){
        p_value <- exp(Y_k * log(mu_0) + (n - Y_k) * log(1 - mu_0) - C_k)
      } else if(K == 2){
        solution <- optimize(
          f = function(mu_01){
            mu_02 <- (mu_0 - a[1] * mu_01) / a[2]
            mu_0_vec <- c(mu_01, mu_02)
            #-2 * sum(pmin(0, Y_k * log(mu_0_vec) + (n - Y_k) * log(1 - mu_0_vec) - C_k))
            -2 * sum(Y_k * log(mu_0_vec) + (n - Y_k) * log(1 - mu_0_vec) - C_k)
          },
          interval = c(0,mu_0/a[1])
        )
        # mu_01 <- seq(0,mu_0/a[1], length.out = 1000)
        # mu_02 <- (mu_0 - a[1] * mu_01) / a[2]
        # obj <- rep(NA, length(mu_01))
        # for(k in 1:length(mu_01)){ 
        #   mu_0_vec <- c(mu_01[k], mu_02[k])
        #   obj[k] <- -2 * sum(pmin(0, Y_k * log(mu_0_vec) + (n - Y_k) * log(1 - mu_0_vec) - C_k))
        # }
        combined_test_stat <- solution$objective
        p_value <- pchisq(q = combined_test_stat, df = 4, lower.tail = FALSE)
      } else{
        stop("For beta-binomial, only fewer than 2 strata work right now.")
      }
  }
  
  #list("psi_k" = psi_k, "V_nk" = V_nk, "lambda_k" = lambda_k, "Y" = Y, "mu_star" = mu_star, "strata" = strata, "p_value" = p_value)
  p_value 
}



