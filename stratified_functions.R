# conservative inference on stratified populations with bounded RVs
library(tidyverse)
library(confseq)
library(lpSolve)


############# confidence sequences and p-values for unstratified samples ##########
#function to shuffle a vector
shuffle <- function(x){sample(x, size = length(x), replace = FALSE)}

#function to compute a cumulative variance quickly
#see: https://stackoverflow.com/questions/52459711/how-to-find-cumulative-variance-or-standard-deviation-in-r
cumvar <- function (x, sd = FALSE) {
  x <- x - x[sample.int(length(x), 1)]  ## see Remark 2 below
  n <- seq_along(x)
  v <- (cumsum(x ^ 2) - cumsum(x) ^ 2 / n) / (n - 1)
  if (sd) v <- sqrt(v)
  v
}

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

########## helper functions ##########
get_empirical_Vn <- function(sample){
  means = cummean(sample)
  # lag by 1, with first index being 1/2
  lagged_means <- c(1/2, means[-length(means)])
  (sample - lagged_means)^2
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

#helper function for strata sample size allocation, turns doubles to integers with the sum preserved
#input: 
#n_strata: the number of samples to take from each stratum, which could be non-integer, but sums to overall sample size
#output:
#a rounded vector with the sum equal to the sum of the original vector
round_strata_sizes <- function(n_strata){
  rounded_n_strata <- floor(n_strata)
  indices <- tail(order(n_strata-rounded_n_strata), round(sum(n_strata)) - sum(rounded_n_strata))
  rounded_n_strata[indices] <- rounded_n_strata[indices] + 1
  rounded_n_strata
}




########## Sequentiall valid risk-measures for stratified samples ############
#function to get samples from a population (or a stratum of a population). Used in get_stratified_pvalue()
#inputs:
  #pop: the population as a vector of values in [0,1]
  #n: the sample size
  #alpha: the targeted significance level, used to compute optimal lambda
  #method: the martingale method to be used, determines how lambda is computed
  #alpha_pars: a list of additional parameters with elements eta_0, d, and epsilon; needed for alpha
#outputs:
  #a list with Y (samples), lambda (tuning parameters for some martingales), psi (normalizing function for some martingales), and V_n (accumulating variance measure for some martingales). May be null depending on method
get_statistics <- function(pop, n, alpha = 0.05, method, pars){
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
    epsilon <- NULL
  } else if(method == "hoeffding"){
    lambda <- rep(sqrt(8 * log(1/alpha) / n), n)
    psi <- get_psi_H(lambda)
    V_n <- rep(1, n)
    epsilon <- NULL
  } else if(method == "beta-binomial"){
    lambda <- NULL
    psi <- NULL
    V_n <- NULL
    epsilon <- NULL
  } else if(method == "hedged"){
    mu_hat <- get_mu_hat_n(Y)
    mu_hat_lag <- lag(mu_hat, 1)
    mu_hat_lag[1] <- 1/2
    sigma_hat_squared <- (1/4 + cumsum(Y - mu_hat)^2) / (1:n + 1) 
    sigma_hat_squared_lag <- lag(sigma_hat_squared, 1)
    sigma_hat_squared_lag[1] <- 1/4
    psi <- NULL
    V_n <- NULL
    epsilon <- NULL
    #this lambda is not thresholded (b/c threshold depends on mu_0), it has to be thresholded later
    lambda <- sqrt(2 * log(2/alpha) / (n * sigma_hat_squared_lag))
  }
  else if(method == "alpha"){
    d <- pars$d
    eta_0 <- pars$eta_0
    running_sum <- lag(cumsum(Y), 1)
    running_sum[1] <- 1/2
    eta <- (d * eta_0 + running_sum) / (d + 1:n - 1)
    psi <- NULL
    V_n <- NULL
    #what is called eta in ALPHA is called lambda here, for compatability reasons
    lambda <- eta
    epsilon <- .001 / sqrt(d + 1:n - 1)
  }
  list(Y = Y, "lambda" = lambda, "psi" = psi, V_n = V_n, epsilon = epsilon)
}



#function to compute a p-value after n samples given a population and strata 
#inputs: 
  #population: a 1-D vector, the population to sample from
  #strata: a 1-D vector of stratum indicators with length equal to length(population)
  #mu_0: the hypothesized overall null mean to be tested
  #n: a vector of length length(unique(strata)), specifying how many samples to draw (at random with replacement) from each stratum
  #method: the martingale to be used, one of hoeffding, empirical_bernstein, hedged, or beta-binomial
  #alpha: the targeted significance level, used in optimizing martingales. Note that resulting p-value is valid for any significance level.
  #alpha_pars: additional parameters for hedged martingale or alpha
#output:
  #a p-value for the null hypothesis: mean(population) = mu_0
get_stratified_pvalue <- function(population, strata, mu_0, n, method, pool = "fisher", alpha = .05, bounds = c(0,1), pars = NULL){
  if(length(bounds) == 2){
    population <- (population - bounds[1]) / diff(bounds)
    mu_0 <- (mu_0 - bounds[1]) / diff(bounds)
    if(!is.null(pars$eta_0)){
      #null mean for alpha needs to be rescaled to reflect bounds as well
      pars$eta_0 <- (pars$eta_0 - bounds[1]) / diff(bounds)
    }
  }
  if(length(population) != length(strata)){
    stop("Strata do not cover population (length(strata) != length(population))")
  }
  strata_names <- unique(strata)
  strata_sizes <- as.numeric(table(strata))
  K <- length(strata_names)
  a <- prop.table(table(strata))
  statistics_strata <- list()
  for(k in 1:K){
    statistics_strata[[k]] <- get_statistics(population[strata == strata_names[k]], n = n[k], method = method, alpha = alpha, pars = list(d = pars$d, eta_0 = pars$eta_0[k]))
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
    
    if(pool == "fisher"){
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
    } else if(pool == "martingale"){
      mu_star <- lp(
        direction = "max",
        objective.in = coefficients,
        const.mat = constraint_matrix[-(1:(2*K)),1:K],
        const.dir = c("==", rep("<=", K), rep(">=", K)),
        const.rhs = c(mu_0, rep(1, K), rep(0, K))
      )$solution[1:K]
      log_pvalue <- sum(constants + coefficients * mu_star)
      p_value <- min(1, exp(log_pvalue))
    } else{
      stop("Input a valid pooling method to combine: either fisher or intersection")
    }
    
  
  } else if(method == "hedged"){
    lambda_k <- statistics_strata %>%
      map(function(x){x$lambda}) %>%
      reduce(c)
    Y <- statistics_strata %>%
      map(function(x){x$Y}) %>%
      reduce(c)
    if(K == 1){
      p_value <- min(1, 1 / prod(1 + pmin(lambda_k, .75/mu_0) * (Y - mu_0)))
      } else if(K == 2 & pool == "fisher"){
      solution <- optimize(
        f = function(mu_01){
          mu_02 <- (mu_0 - a[1] * mu_01) / a[2]
          2 * max(0, sum(log(1 + pmin(statistics_strata[[1]]$lambda, .75/mu_01) * (statistics_strata[[1]]$Y - mu_01)))) + max(0, sum(log(1 + pmin(statistics_strata[[2]]$lambda, .75/mu_02) * (statistics_strata[[2]]$Y - mu_02))))
        },
        interval = c(0,mu_0/a[1]))
      combined_test_stat <- solution$objective
      p_value <- pchisq(q = combined_test_stat, df = 4, lower.tail = FALSE)
      } else if(K == 2 & pool == "martingale"){
        solution <- optimize(
          f = function(mu_01){
            mu_02 <- (mu_0 - a[1] * mu_01) / a[2]
            sum(log(1 + pmin(statistics_strata[[1]]$lambda, .75/mu_01) * (statistics_strata[[1]]$Y - mu_01))) + sum(log(1 + pmin(statistics_strata[[2]]$lambda, .75/mu_02) * (statistics_strata[[2]]$Y - mu_02)))
          },
          interval = c(0,mu_0/a[1]))
        log_pvalue <- -solution$objective
        p_value <- min(1, exp(log_pvalue))
      } else{
      stop("Hedged martingale only works for two or fewer strata.")
      }
  } else if(method == "alpha"){
    lambda_k <- statistics_strata %>%
      map(function(x){x$lambda}) %>%
      reduce(c)
    Y <- statistics_strata %>%
      map(function(x){x$Y}) %>%
      reduce(c)
    if(K == 1){
      eta <- pmin(pmax(statistics_strata[[1]]$lambda, mu_0 + statistics_strata[[1]]$epsilon), 1)
      p_value <- min(1, 1 / prod(1 + (eta / mu_0 - 1)/(1 - mu_0) * (statistics_strata[[1]]$Y - mu_0)))
    } else if(K == 2 & pool == "fisher"){
      epsilon_1 <- statistics_strata[[1]]$epsilon
      epsilon_2 <- statistics_strata[[2]]$epsilon
      solution <- optimize(
          f = function(mu_01){
            mu_02 <- (mu_0 - a[1] * mu_01) / a[2]
            eta_1 <- pmin(pmax(statistics_strata[[1]]$lambda, mu_01 + epsilon_1), 1)
            eta_2 <- pmin(pmax(statistics_strata[[2]]$lambda, mu_02 + epsilon_2), 1)
            2 * (max(0, sum(log(1 + (eta_1 / mu_01 - 1)/(1 - mu_01) * (statistics_strata[[1]]$Y - mu_01)))) + max(0, sum(log(1 + (eta_2 / mu_02 - 1)/(1 - mu_02) * (statistics_strata[[2]]$Y - mu_02)))))
          },
          interval = c(0,mu_0/a[1]))
        combined_test_stat <- solution$objective
        p_value <- pchisq(q = combined_test_stat, df = 4, lower.tail = FALSE)
      } else if(K == 2 & pool == "martingale"){
        epsilon_1 <- statistics_strata[[1]]$epsilon
        epsilon_2 <- statistics_strata[[2]]$epsilon
        solution <- optimize(
          f = function(mu_01){
            mu_02 <- (mu_0 - a[1] * mu_01) / a[2]
            eta_1 <- pmin(pmax(statistics_strata[[1]]$lambda, mu_01 + epsilon_1), 1)
            eta_2 <- pmin(pmax(statistics_strata[[2]]$lambda, mu_02 + epsilon_2), 1)
            sum(log(1 + (eta_1 / mu_01 - 1)/(1 - mu_01) * (statistics_strata[[1]]$Y - mu_01))) + sum(log(1 + (eta_2 / mu_02 - 1)/(1 - mu_02) * (statistics_strata[[2]]$Y - mu_02)))
          },
          interval = c(0,mu_0/a[1]))
        log_pvalue <- -solution$objective
        p_value <- min(1, exp(log_pvalue))
      } else{
        stop("ALPHA only works for two or fewer strata.")
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
      
      if(K == 1){
        p_value <- exp(Y_k * log(mu_0) + (n - Y_k) * log(1 - mu_0) - C_k)
      } else if(K == 2 & pool == "fisher"){
        solution <- optimize(
          f = function(mu_01){
            mu_02 <- (mu_0 - a[1] * mu_01) / a[2]
            mu_0_vec <- c(mu_01, mu_02)
            -2 * sum(Y_k * log(mu_0_vec) + (n - Y_k) * log(1 - mu_0_vec) - C_k)
          },
          interval = c(0,min(1,mu_0/a[1]))
        )
        combined_test_stat <- solution$objective
        p_value <- pchisq(q = combined_test_stat, df = 4, lower.tail = FALSE)
      } else{
        stop("For beta-binomial, only fewer than 2 strata and fisher pooling work.")
      }
  }
  #list("psi_k" = psi_k, "V_nk" = V_nk, "lambda_k" = lambda_k, "Y" = Y, "mu_star" = mu_star, "strata" = strata, "p_value" = p_value)
  p_value 
}


#a function to help choose eta for ALPHA
#inputs:
  #audit_type: a string, either "polling" or "comparison"
  #diluted_margin: a scalar or vector in [0,1]^K; the diluted margin in each strata 1:K
  #overstatement_rate: a scalar or vector in [0,1]^K, only needed for comparison audits, the rate of overstatement errors (i.e., a guess at \bar{\omega}) 
#outputs:
  #an initial guess at the population mean (eta_0), which is a parameter for alpha
choose_eta_0 <- function(audit_type = "comparison", diluted_margin, overstatement_rate = .001){
  if(audit_type == "polling"){
    eta_0 <- (1 + diluted_margin)/2
  } else if(audit_type == "comparison"){
    eta_0 <- (1 - overstatement_rate) / (2 - diluted_margin)
  } else{
    stop("input a valid audit_type, either polling or comparison")
  }
  eta_0
}




#a function to simulate ballot-level comparison risk-limiting audits
#inputs:
  #population: a length N vector of population values, which will be sampled
  #strata: a length N vector of K unique integers {1,...,K} recording the stratum each ballot belongs to. 
  #sample_sizes: a matrix with K columns describing how many samples to draw from each stratum, rows can be different, simulations will be run at each sample size by iterating over rows
  #n_sims: a scalar integer, the number of simulations to run
  #alpha: a scalar float in (0,1), the risk limit
  #bounds: a length 2 vector or a Kx2 matrix; the known minimum and maximum values in the population or within each stratum; may be wider but not narrower than the population within strata
run_stratified_simulation <- function(population, strata, sample_sizes, mu_0 = 0.5, n_sims = 300, alpha = .05, bounds = c(0,1), pars){
  if(min(population) < bounds[1] | max(population) > bounds[2]){
    stop("range(population) must be contained within bounds")
  }
  
  # power_ppm_unstrat <- matrix(NA, nrow = nrow(sample_sizes), ncol = 1)
  # power_ppm_fisher <- matrix(NA, nrow = nrow(sample_sizes), ncol = 1)
  power_eb_unstrat <- matrix(NA, nrow = nrow(sample_sizes), ncol = 1)
  power_eb_fisher <- matrix(NA, nrow = nrow(sample_sizes), ncol = 1)
  power_eb_martingale <- matrix(NA, nrow = nrow(sample_sizes), ncol = 1)
  power_hedged_unstrat <- matrix(NA, nrow = nrow(sample_sizes), ncol = 1)
  power_hedged_fisher <- matrix(NA, nrow = nrow(sample_sizes), ncol = 1)
  power_hedged_martingale <- matrix(NA, nrow = nrow(sample_sizes), ncol = 1)
  power_alpha_unstrat <- matrix(NA, nrow = nrow(sample_sizes), ncol = 1)
  power_alpha_fisher <- matrix(NA, nrow = nrow(sample_sizes), ncol = 1)
  power_alpha_martingale <- matrix(NA, nrow = nrow(sample_sizes), ncol = 1)
  
  #need to average the eta_0 for the unstratified version of ALPHA
  unstrat_pars <- pars
  unstrat_pars$eta_0 <- mean(pars$eta_0)
  
  for(i in 1:nrow(sample_sizes)){
    replicates_eb_unstrat <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = rep(1,length(population)), mu_0 = mu_0, n = sum(sample_sizes[i,]), method = "empirical_bernstein", bounds = bounds))
    replicates_eb_fisher <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = strata, mu_0 = mu_0, n = sample_sizes[i,], method = "empirical_bernstein", pool = "fisher", bounds = bounds))
    replicates_eb_martingale <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = strata, mu_0 = mu_0, n = sample_sizes[i,], method = "empirical_bernstein", pool = "martingale", bounds = bounds))
    replicates_hedged_unstrat <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = rep(1,length(population)), mu_0 = mu_0, n = sum(sample_sizes[i,]), method = "hedged", bounds = bounds))
    replicates_hedged_fisher <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = strata, mu_0 = mu_0, n = sample_sizes[i,], method = "hedged", pool = "fisher", bounds = bounds))
    replicates_hedged_martingale <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = strata, mu_0 = mu_0, n = sample_sizes[i,], method = "hedged", pool = "martingale", bounds = bounds))
    replicates_alpha_unstrat <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = rep(1, length(population)), mu_0 = mu_0, n = sum(sample_sizes[i,]), method = "alpha", bounds = bounds, pars = unstrat_pars))
    replicates_alpha_fisher <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = strata, mu_0 = mu_0, n = sample_sizes[i,], method = "alpha", pool = "fisher", bounds = bounds, pars = pars))
    replicates_alpha_martingale <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = strata, mu_0 = mu_0, n = sample_sizes[i,], method = "alpha", pool = "martingale", bounds = bounds, pars = pars))
    
    power_alpha_unstrat[i,1] <- mean(replicates_alpha_unstrat < alpha)
    power_alpha_fisher[i,1] <- mean(replicates_alpha_fisher < alpha)
    power_alpha_martingale[i,1] <- mean(replicates_alpha_martingale < alpha)
    power_eb_unstrat[i,1] <- mean(replicates_eb_unstrat < alpha)
    power_eb_fisher[i,1] <- mean(replicates_eb_fisher < alpha)
    power_eb_martingale[i,1] <- mean(replicates_eb_martingale < alpha)
    power_hedged_unstrat[i,1] <- mean(replicates_hedged_unstrat < alpha)
    power_hedged_fisher[i,1] <- mean(replicates_hedged_fisher < alpha)
    power_hedged_martingale[i,1] <- mean(replicates_hedged_martingale < alpha)
  }
  
  # colnames(power_ppm_unstrat) <- "power_ppm_unstrat"
  # colnames(power_ppm_fisher) <- "power_ppm_fisher"
  
  colnames(power_alpha_fisher) <- "power_alpha_fisher"
  colnames(power_alpha_martingale) <- "power_alpha_martingale"
  colnames(power_alpha_unstrat) <- "power_alpha_unstrat"
  colnames(power_eb_fisher) <- "power_eb_fisher"
  colnames(power_eb_martingale) <- "power_eb_martingale"
  colnames(power_eb_unstrat) <- "power_eb_unstrat"
  colnames(power_hedged_fisher) <- "power_hedged_fisher"
  colnames(power_hedged_martingale) <- "power_hedged_martingale"
  colnames(power_hedged_unstrat) <- "power_hedged_unstrat"
  
  
  power_frame <- data.frame(
    power_alpha_fisher,
    power_alpha_martingale,
    power_alpha_unstrat,
    power_eb_fisher,
    power_eb_martingale,
    power_eb_unstrat, 
    power_hedged_fisher,
    power_hedged_martingale,
    power_hedged_unstrat,
    # power_ppm_unstrat, 
    # power_ppm_fisher, 
    "total_sample_size" = rowSums(sample_sizes)) %>%
    pivot_longer(cols = starts_with("power"), names_to = "design", values_to = "power", names_prefix = "power_") 
  power_frame
} 




# compute a P-value from two strata using ALPHA with a particular stratumwise sample allocation scheme
#inputs:
  #stratum_1: vector of doubles in [0,u[1]], population values for the first stratum
  #stratum_2: vector of doubles in [0,u[2]], population values for the second stratum
  #replace: boolean, TRUE if sampling is with replacement, FALSE if simple random sampling
  #u: length-2 double, the upper bound on the population within each stratum
  #d: length-2 double, the d parameter for ALPHA in each stratum; can be left undefined
  #eta_0: length-2 double, the eta_0 parameter for ALPHA in each stratum; can be left undefined
  #rule: string, the allocation rule to be used
  #resolution: integer, the resolution of the grid of mu_0, defaults to the equivalent of 1 vote from largest stratum
  #mu_0: double, hypothesized global null, defaults to 1/2 which is the usual global null for ALPHA
#output:
  #a dataframe with a sequence of P-values and indicators for whether samples were taken from each stratum
get_two_strata_alpha <- function(stratum_1, stratum_2, replace, u, d = NULL, eta_0 = NULL, rule = "equal", resolution = NULL, mu_0 = 0.5, combine = "product"){
  N <- c(length(stratum_1), length(stratum_2))
  if(is.null(resolution)){resolution <- 2*max(N)}
  w <- N / sum(N)
  shuffled_1 <- sample(stratum_1, size = N[1], replace = replace)
  shuffled_2 <- sample(stratum_2, size = N[2], replace = replace)
  S_1 <- c(0, cumsum(shuffled_1)[-length(shuffled_1)])
  S_2 <- c(0, cumsum(shuffled_2)[-length(shuffled_2)])
  
  #epsilon is defined so that eta is always one assorter value (1 invalid vote) above the null mean 
  epsilon_1 <- 1/(2*N[1])
  epsilon_2 <- 1/(2*N[2])
  
  mu_01 <- matrix(seq(epsilon_1, mu_0/w[1] - epsilon_1, length.out = resolution), ncol = resolution, nrow = N[1], byrow = TRUE)
  mu_02 <- matrix((mu_0 - w[1] * mu_01[1,]) / w[2], ncol = ncol(mu_01), nrow = N[2], byrow = TRUE)
  
  if(replace == FALSE){
    m_1 <- (N[1] * mu_01 - S_1) / (N[1] - 1:N[1] + 1)
    m_2 <- (N[2] * mu_02 - S_2) / (N[2] - 1:N[2] + 1)
  } else{
    m_1 <- mu_01
    m_2 <- mu_02
  }
  
  #eta as in ALPHA, or just the running mean if parameters are undefined
  if(!is.null(d) & !is.null(eta_0)){
    eta_1 <- matrix((d[1]*eta_0[1] + c(0, lag(cumsum(shuffled_1), 1)[2:N[1]])) / (d[1] + 1:N[1] - 1), nrow = N[1], ncol = ncol(m_1))
    eta_2 <- matrix((d[2]*eta_0[2] + c(0, lag(cumsum(shuffled_2), 1)[2:N[2]])) / (d[2] + 1:N[2] - 1), nrow = N[2], ncol = ncol(m_2))
  } else{
    eta_1 <- matrix(c(0, lag(cummean(shuffled_1), 1)[2:N[1]]), nrow = N[1], ncol = ncol(m_1))
    eta_2 <- matrix(c(0, lag(cummean(shuffled_2), 1)[2:N[2]]), nrow = N[2], ncol = ncol(m_2))
  }
  
  
  
  #truncation for eta
  eta_1 <- pmin(matrix(u[1]*(1-.Machine$double.eps), nrow = N[1], ncol = ncol(m_1)), pmax(eta_1, m_1 + epsilon_1))
  eta_2 <- pmin(matrix(u[2]*(1-.Machine$double.eps), nrow = N[2], ncol = ncol(m_2)), pmax(eta_2, m_2 + epsilon_2))
  
  #there can be terms equal to 0, because eta_1 is allowed to be as high as u. This is a problem, we may end up in a situation where we have 0 * Inf. The martingales can and should always be bounded away from 0.
  #logging may help too
  terms_1 <- (shuffled_1 / m_1) * (eta_1 - m_1) / (u[1] - m_1) + (u[1] - eta_1) / (u[1] - m_1)
  terms_2 <- (shuffled_2 / m_2) * (eta_2 - m_2) / (u[2] - m_2) + (u[2] - eta_2) / (u[2] - m_2)
  
  terms_1 <- rbind(terms_1, matrix(1, ncol = ncol(terms_1), nrow = max(0, N[2] - N[1])))
  terms_2 <- rbind(terms_2, matrix(1, ncol = ncol(terms_2), nrow = max(0, N[1] - N[2])))
  
  terms_1[m_1 < 0] <- Inf
  terms_2[m_2 < 0] <- Inf
  terms_1[m_1 > u[1]] <- 0
  terms_2[m_2 > u[2]] <- 0
  
  
  #is the running max a valid rule? The running max is a valid rule for a single martingale, but it might not be when the combination is done.
  if(rule == "equal"){
    allocation <- function(x){x}
    terms_1_stopped <- apply(terms_1, 2, allocation)
    terms_2_stopped <- apply(terms_2, 2, allocation)
  } else if(rule == "hard_threshold"){
    allocation <- function(x){
      mart <- cumprod(x)
      n <- length(x)
      crossed <- min(which(mart < .95 & 1:n > 30))
      if(is.finite(crossed)){
        x[crossed:n] <- 1 
      }
     x 
    }
    terms_1_stopped <- apply(terms_1, 2, allocation)
    terms_2_stopped <- apply(terms_2, 2, allocation)
  } else if(rule == "gaussian_ucb"){
    allocation <- function(x){
      log_x <- log(x)
      n <- length(x)
      running_variance <- c(1/2, cumvar(log_x)[2:n])
      #bound <- sqrt(4 * running_variance * log(n) / 1:n)
      bound <- qnorm(.95) * sqrt(running_variance / 1:n)
      ucb <- cummean(log_x) + bound
      crossed <- min(which(ucb < 0 & (1:n) > 30))
      if(is.finite(crossed)){
        x[crossed:n] <- 1 
      }
      x 
    } 
    terms_1_stopped <- apply(terms_1, 2, allocation)
    terms_2_stopped <- apply(terms_2, 2, allocation)
  } else if(rule == "alpha_ucb"){
      #compute terms for lower-sided hypothesis test
      ucb_eta_1 <- matrix(c(0, lag(cummean(shuffled_1), 1)[2:N[1]]), nrow = N[1], ncol = ncol(m_1))
      ucb_eta_2 <- matrix(c(0, lag(cummean(shuffled_2), 1)[2:N[2]]), nrow = N[2], ncol = ncol(m_2))
      ucb_eta_1 <- pmax(matrix(epsilon_1, nrow=N[1],ncol=ncol(m_1)), pmin(ucb_eta_1, m_1 - epsilon_1))
      ucb_eta_2 <- pmax(matrix(epsilon_2, nrow=N[2],ncol=ncol(m_2)), pmin(ucb_eta_2, m_2 - epsilon_2))
      ucb_terms_1 <- (shuffled_1 / m_1) * (ucb_eta_1 - m_1) / (u[1] - m_1) + (u[1] - ucb_eta_1) / (u[1] - m_1)
      ucb_terms_2 <- (shuffled_2 / m_2) * (ucb_eta_2 - m_2) / (u[2] - m_2) + (u[2] - ucb_eta_2) / (u[2] - m_2)
      ucb_terms_1[m_1 < epsilon_1] <- 0
      ucb_terms_2[m_2 < epsilon_2] <- 0
      ucb_terms_1[m_1 > u[1] - epsilon_1] <- Inf
      ucb_terms_2[m_2 > u[2] - epsilon_2] <- Inf
      ucb_mart_1 <- apply(ucb_terms_1, 2, cumprod)
      ucb_mart_2 <- apply(ucb_terms_2, 2, cumprod)
      ucb_mart_1 <- apply(ucb_mart_1, 2, cummax)
      ucb_mart_2 <- apply(ucb_mart_2, 2, cummax)
      
      beta <- .05
      terms_1_stopped <- terms_1
      terms_2_stopped <- terms_2
      terms_1_stopped[ucb_mart_1 > 1/beta] <- 1
      terms_2_stopped[ucb_mart_2 > 1/beta] <- 1
  } else{
    stop("Input valid allocation rule.")
  }
  #stop counting for a given null if it (deterministically) cannot be rejected 
  terms_1_stopped[m_1 > (u[1] - epsilon_1)] <- 1
  terms_2_stopped[m_2 > (u[2] - epsilon_2)] <- 1
  
  mart_1 <- apply(terms_1_stopped, 2, cumprod)
  mart_2 <- apply(terms_2_stopped, 2, cumprod)
  
  if(combine == "product"){
    intersection_martingale <- mart_1 * mart_2
    minimized_martingale <- apply(intersection_martingale , 1, min)
    null_selected <- apply(intersection_martingale, 1, which.min)
    
    #output
    sample_counter_1 <- as.numeric(terms_1_stopped[cbind(1:length(null_selected), null_selected)] != 1)
    sample_counter_2 <- as.numeric(terms_2_stopped[cbind(1:length(null_selected), null_selected)] != 1)
    p_value <- pmin(1, 1 / minimized_martingale)
  } else if(combine == "fisher"){
    pvals_1 <- pmin(matrix(1, nrow=nrow(mart_1), ncol=ncol(mart_1)), 1/mart_1)
    pvals_2 <- pmin(matrix(1, nrow=nrow(mart_2), ncol=ncol(mart_2)), 1/mart_2)
    fisher <- -2 * (log(pvals_1) + log(pvals_2))
    minimized_fisher <- apply(fisher, 1, min)
    null_selected <- apply(fisher, 1, which.min)
    
    #output
    sample_counter_1 <- as.numeric(terms_1_stopped[cbind(1:length(null_selected), null_selected)] != 1)
    sample_counter_2 <- as.numeric(terms_2_stopped[cbind(1:length(null_selected), null_selected)] != 1)
    p_value <- 1 - pchisq(q = minimized_fisher, df = 4)
  }
  
  as.matrix(data.frame(p_value = p_value, counter_1 = sample_counter_1, counter_2 = sample_counter_2))
}




# compute a P-value from two strata using empirical Bernstein with a particular stratumwise sample allocation scheme
#inputs:
  #stratum_1: vector of doubles in [0,u[1]], population values for the first stratum
  #stratum_2: vector of doubles in [0,u[2]], population values for the second stratum
  #replace: boolean, TRUE if sampling is with replacement, FALSE if simple random sampling
  #rule: string, the allocation rule to be used
  #resolution: integer, the resolution of the grid of mu_0, defaults to the equivalent of 1 vote from largest stratum
  #mu_0: double, hypothesized global null, defaults to 1/2 which is the usual global null for RLAs
#output:
#a dataframe with a sequence of P-values and indicators for whether samples were taken from each stratum
get_two_strata_EB <- function(stratum_1, stratum_2, replace, u, rule = "equal", resolution = NULL, mu_0 = 0.5, combine = "product"){
  N <- c(length(stratum_1), length(stratum_2))
  if(is.null(resolution)){resolution <- 2*max(N)}
  w <- N / sum(N)
  
  #we need to rescale everything to [0,1] in order to uses empirical Bernstein
  shuffled_1 <- sample(stratum_1, size = N[1], replace = replace) / u[1]
  shuffled_2 <- sample(stratum_2, size = N[2], replace = replace) / u[2]
  S_1 <- c(0, cumsum(shuffled_1)[-length(shuffled_1)])
  S_2 <- c(0, cumsum(shuffled_2)[-length(shuffled_2)])
  
  
  #epsilon is defined so that eta is always one assorter value (1 invalid vote) above the null mean 
  epsilon_1 <- 1/(2*N[1])
  epsilon_2 <- 1/(2*N[2])
  
  mu_01 <- matrix(seq(epsilon_1, mu_0/w[1] - epsilon_1, length.out = resolution), ncol = resolution, nrow = N[1], byrow = TRUE) / u[1]
  mu_02 <- matrix((mu_0 - w[1] * mu_01[1,]) / w[2], ncol = ncol(mu_01), nrow = N[2], byrow = TRUE) / u[2]
  
  if(replace == FALSE){
    m_1 <- (N[1] * mu_01 - S_1) / (N[1] - 1:N[1] + 1)
    m_2 <- (N[2] * mu_02 - S_2) / (N[2] - 1:N[2] + 1)
  } else{
    m_1 <- mu_01
    m_2 <- mu_02
  }
  

  lambda_1 <- pmin(sqrt(log(2/.05) / (c(1/2, lag(get_sigma_hat_squared(shuffled_1))[2:N[1]]) * 1:N[1] * log(1 + 1:N[1]))), 0.75)
  lambda_2 <- pmin(sqrt(log(2/.05) / (c(1/2, lag(get_sigma_hat_squared(shuffled_2))[2:N[2]]) * 1:N[2] * log(1 + 1:N[2]))), 0.75)
  v_1 <- get_empirical_Vn(shuffled_1)
  v_2 <- get_empirical_Vn(shuffled_2)
  psi_E_1 <- get_psi_E(lambda_1)
  psi_E_2 <- get_psi_E(lambda_2)
  
  log_terms_1 <- lambda_1 * (shuffled_1 - m_1) - v_1 * psi_E_1
  log_terms_2 <- lambda_2 * (shuffled_2 - m_2) - v_2 * psi_E_2
  
  log_terms_1 <- rbind(log_terms_1, matrix(0, ncol = ncol(log_terms_1), nrow = max(0, N[2] - N[1])))
  log_terms_2 <- rbind(log_terms_2, matrix(0, ncol = ncol(log_terms_2), nrow = max(0, N[1] - N[2])))
  
  log_terms_1[m_1 < 0] <- Inf
  log_terms_2[m_2 < 0] <- Inf
  log_terms_1[m_1 > 1] <- -Inf
  log_terms_2[m_2 > 1] <- -Inf
  
  
  #is the running max a valid rule? The running max is a valid rule for a single martingale, but it might not be when the combination is done.
  if(rule == "equal"){
    log_terms_1_stopped <- log_terms_1
    log_terms_2_stopped <- log_terms_2
  } else if(rule == "hard_threshold"){
    allocation <- function(x){
      log_mart <- cumsum(x)
      n <- length(x)
      crossed <- min(which(log_mart < log(.95) & 1:n > 30))
      if(is.finite(crossed)){
        x[crossed:n] <- 0 
      }
      x 
    }
    log_terms_1_stopped <- apply(log_terms_1, 2, allocation)
    log_terms_2_stopped <- apply(log_terms_2, 2, allocation)
  } else if(rule == "gaussian_ucb"){
    allocation <- function(x){
      n <- length(x)
      running_variance <- c(1/2, cumvar(x)[2:n])
      #bound <- sqrt(4 * running_variance * log(n) / 1:n)
      bound <- qnorm(.95) * sqrt(running_variance / 1:n)
      ucb <- cummean(x) + bound
      crossed <- min(which(ucb < 0 & (1:n) > 30))
      if(is.finite(crossed)){
        x[crossed:n] <- 0 
      }
      x 
    } 
    log_terms_1_stopped <- apply(log_terms_1, 2, allocation)
    log_terms_2_stopped <- apply(log_terms_2, 2, allocation)
  } else if(rule == "eb_ucb"){
    #compute terms for lower-sided hypothesis test
    log_ucb_terms_1 <- -lambda_1 * (shuffled_1 - m_1) - v_1 * psi_E_1
    log_ucb_terms_2 <- -lambda_2 * (shuffled_2 - m_2) - v_2 * psi_E_2
    log_ucb_terms_1[m_1 < epsilon_1] <- -Inf
    log_ucb_terms_2[m_2 < epsilon_2] <- -Inf
    log_ucb_terms_1[m_1 > 1 - epsilon_1] <- Inf
    log_ucb_terms_2[m_2 > 1 - epsilon_2] <- Inf
    log_ucb_mart_1 <- apply(log_ucb_terms_1, 2, cumsum)
    log_ucb_mart_2 <- apply(log_ucb_terms_2, 2, cumsum)
    log_ucb_mart_1 <- apply(log_ucb_mart_1, 2, cummax)
    log_ucb_mart_2 <- apply(log_ucb_mart_2, 2, cummax)
    
    beta <- .05
    log_terms_1_stopped <- log_terms_1
    log_terms_2_stopped <- log_terms_2
    log_terms_1_stopped[log_ucb_mart_1 > log(1/beta)] <- 0
    log_terms_2_stopped[log_ucb_terms_2 > log(1/beta)] <- 0
  } else{
    stop("Input valid allocation rule.")
  }
  log_terms_1_stopped[m_1 > (1 - epsilon_1)] <- 0
  log_terms_2_stopped[m_2 > (1 - epsilon_2)] <- 0
  
  log_mart_1 <- apply(log_terms_1_stopped, 2, cumsum)
  log_mart_2 <- apply(log_terms_2_stopped, 2, cumsum)
  
  if(combine == "product"){
    log_intersection_mart <- log_mart_1 + log_mart_2
    log_minimized_martingale <- apply(log_intersection_mart, 1, min)
    null_selected <- apply(log_intersection_mart, 1, which.min)
    p_value <- pmin(1, exp(-log_minimized_martingale))
  } else if(combine == "fisher"){
    pvals_1 <- pmin(matrix(1, nrow=nrow(log_mart_1), ncol=ncol(log_mart_1)), exp(-log_mart_1))
    pvals_2 <- pmin(matrix(1, nrow=nrow(log_mart_2), ncol=ncol(log_mart_2)), exp(-log_mart_2))
    fisher <- -2 * (log(pvals_1) + log(pvals_2))
    minimized_fisher <- apply(fisher, 1, min)
    null_selected <- apply(fisher, 1, which.min)
    p_value <- 1 - pchisq(q = minimized_fisher, df = 4)
  }
  sample_counter_1 <- as.numeric(log_terms_1_stopped[cbind(1:length(null_selected), null_selected)] != 0)
  sample_counter_2 <- as.numeric(log_terms_2_stopped[cbind(1:length(null_selected), null_selected)] != 0)
  
  as.matrix(data.frame(p_value = p_value, counter_1 = sample_counter_1, counter_2 = sample_counter_2))
}





#a function to compute Gaffke bounds or p-values on a stratified population
#inputs:
  #population: a vector of values of a finite population
  #strata: a vector of the same length as population indicating strata membership for each element of population
  #n: a vector of length unique(strata) indicating how many samples to draw from each stratum
  #B: how many monte carlo draws to use when computing Gaffke
  #mu_0: if computing a p-value, the mu_0 to test against
  #alpha: if computing a confidence bound, the level
  #method: the method for dealing with the stratification one of
    #sum: weighted sum across within stratum resamples, weights are the (known) proportions corresponding to each stratums share of the population. Resembles the handling of stratification by the usual bootstrap
    #combine: use the intersection-union strategy, combining P-values across strata and maximizing. Tricky optimization problem, only works for 2 strata right now.
#outputs:
  #a p-value or confidence bound
run_stratified_gaffke <- function(population, strata, n, B = 200, mu_0 = NULL, alpha = .05, method = "sum"){
  strata_names <- unique(strata)
  strata_sizes <- as.numeric(table(strata))
  K <- length(strata_names)
  a <- prop.table(table(strata))
  samples_strata <- list()
  gaffke_strata <- matrix(NA, nrow = B, ncol = K)
  for(k in 1:K){
    samples_strata[[k]] <- sample(population[strata == strata_names[k]], size = n[k], replace = TRUE)
    Z <- matrix(rexp(n = B * n[k]), nrow = B, ncol = n[k]) 
    D <- Z / (rowSums(Z) + rexp(B))
    gaffke_strata[,k] <- D %*% samples_strata[[k]]
  }
  if(method == "sum"){
    gaffke_means <- gaffke_strata %*% as.matrix(a)
    if(!is.null(mu_0)){
      mean(gaffke_means <= mu_0)
    } else{
      quantile(gaffke_means, alpha)
    }
  } else if(method == "combine"){
    if(K > 2){
      stop("Combine only works with 2 strata.")
    }
    if(is.null(mu_0)){
      stop("Combine only works when doing a hypothesis test (supply a value for mu_0).")
    }
    combined_p <- function(mu_01){
      mu_02 <- (mu_0 - a[1] * mu_01) / a[2] 
      p_1 <- mean(c(gaffke_strata[,1] <= mu_01, TRUE))
      p_2 <- mean(c(gaffke_strata[,2] <= mu_02, TRUE))
      combined_test_stat <- -2 * (log(p_1) + log(p_2))
      p_val <- pchisq(q = combined_test_stat, df = 4, lower.tail = FALSE)
      p_val
    }
    #if all the Gaffke samples are equal to 0, the maximum occurs when the null mean for stratum 1 is equal to 0 and the rest is allocated to stratum 2
    if(all(gaffke_strata[,1] == 0)){
      max_p_val <- combined_p(0)
    } else{
      #otherwise, the maximum has to occur within the range of stratum 1 resamples
      max_p_val <- optimize(combined_p, interval = range(gaffke_strata[,1]), maximum = TRUE)$objective
    }
    max_p_val
  }
}


#a function to run the stratified t-test, mirroring the functionality of the ones above
#inputs:
  #population: a vector of values of a finite population
  #strata: a vector of the same length as population indicating strata membership for each element of population
  #n: a vector of length unique(strata) indicating how many samples to draw from each stratum
  #mu_0: if computing a p-value, the mu_0 to test against
#outputs:
  #a p-value or confidence bound
run_stratified_t_test <- function(population, strata, n, mu_0 = 0){
  strata_names <- unique(strata)
  strata_sizes <- as.numeric(table(strata))
  K <- length(strata_names)
  a <- prop.table(table(strata))
  strata_sample_means <- rep(NA, K)
  strata_sample_vars <- rep(NA, K)
  for(k in 1:K){
    stratum_samples <- sample(population[strata == strata_names[k]], size = n[k], replace = TRUE)
    strata_sample_means[k] <- mean(stratum_samples)
    strata_sample_vars[k] <- var(stratum_samples)
  }
  mean_estimate <- sum(a * strata_sample_means)
  std_error <- sqrt(sum(a^2 * strata_sample_vars / n))
  #deals with edge cases where we get 0/0, which returns NaN. 
  #the nature of the null is such that we should treat this as a p value of 1:
  #intuitively the conclusion would be that the mean is equal to the null mean with complete certainty.
  if(mean_estimate == mu_0 & std_error == 0){
    p_value <- 1
  } else{
    p_value <- 1 - pt(q = (mean_estimate - mu_0) / std_error, df = min(n))
  }
  p_value
}
