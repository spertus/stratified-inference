source("stratified_functions.R")

########## ALPHA allocation rules ############
#helper function to get sample sizes from array
get_total_sample_size <- function(x, alpha){
  if(any(x[,1] < alpha)){
    sum(x[1:min(which(x[,1] < alpha)),2:3])
  } else{
    NA
  }
}

run_allocation_simulation <- function(reported_tally, hand_tally, strata, n_sims = 300, alpha = .05){
  v <- tapply(reported_tally, strata, function(x){2 * mean(x) - 1})
  assorter_bound <- 1
  omegas <- reported_tally - hand_tally
  u <-(1 + 1/assorter_bound) / (2 - v/assorter_bound)
  stratum_1 <- (1 - omegas[strata == 1] / assorter_bound) / (2 - v[1] / assorter_bound)
  stratum_2 <- (1 - omegas[strata == 2] / assorter_bound) / (2 - v[2] / assorter_bound)
  
  reported_margin <- 2*mean(reported_tally) - 1
  
  true_margin <- 2*mean(hand_tally) - 1
  understatements <- mean(population == 1)
  overstatements <- mean(population == 0)
  d <- c(20, 20)
  eta_0 <- 1 / (2 - v / assorter_bound)

  results_alpha_equal_product <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, eta_0 = eta_0, u = u, replace = FALSE, combine = "product", rule = "equal"))
  results_alpha_gaussian_product <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, eta_0 = eta_0, u = u, replace = FALSE, combine = "product", rule = "gaussian_ucb"))
  results_alpha_ucb_product <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, eta_0 = eta_0, u = u, replace = FALSE, combine = "product", rule = "alpha_ucb"))
  results_alpha_equal_fisher <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, eta_0 = eta_0, u = u, replace = FALSE, combine = "fisher", rule = "equal"))
  results_alpha_gaussian_fisher <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, eta_0 = eta_0, u = u, replace = FALSE, combine = "fisher", rule = "gaussian_ucb"))
  results_alpha_ucb_fisher <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, eta_0 = eta_0, u = u, replace = FALSE, combine = "fisher", rule = "alpha_ucb"))
  results_eb_equal_product <- replicate(n_sims, get_two_strata_EB(stratum_1, stratum_2, mu_0 = 0.5, u = u, replace = FALSE, combine = "product", rule = "equal"))
  results_eb_gaussian_product <- replicate(n_sims, get_two_strata_EB(stratum_1, stratum_2, mu_0 = 0.5, u = u, replace = FALSE, combine = "product", rule = "gaussian_ucb"))
  results_eb_ucb_product <- replicate(n_sims, get_two_strata_EB(stratum_1, stratum_2, mu_0 = 0.5, u = u, replace = FALSE, combine = "product", rule = "eb_ucb"))
  results_eb_equal_fisher <- replicate(n_sims, get_two_strata_EB(stratum_1, stratum_2, mu_0 = 0.5, u = u, replace = FALSE, combine = "fisher", rule = "equal"))
  results_eb_gaussian_fisher <- replicate(n_sims, get_two_strata_EB(stratum_1, stratum_2, mu_0 = 0.5, u = u, replace = FALSE, combine = "fisher", rule = "gaussian_ucb"))
  results_eb_ucb_fisher <- replicate(n_sims, get_two_strata_EB(stratum_1, stratum_2, mu_0 = 0.5, u = u, replace = FALSE, combine = "fisher", rule = "eb_ucb"))
  
  stopping_times_alpha_equal_product <- apply(results_alpha_equal_product, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_alpha_gaussian_product <- apply(results_alpha_gaussian_product, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_alpha_ucb_product <- apply(results_alpha_ucb_product, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_alpha_equal_fisher <- apply(results_alpha_equal_fisher, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_alpha_gaussian_fisher <- apply(results_alpha_gaussian_fisher, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_alpha_ucb_fisher <- apply(results_alpha_ucb_fisher, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_eb_equal_product <- apply(results_eb_equal_product, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_eb_gaussian_product <- apply(results_eb_gaussian_product, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_eb_ucb_product <- apply(results_eb_ucb_product, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_eb_equal_fisher <- apply(results_eb_equal_fisher, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_eb_gaussian_fisher <- apply(results_eb_gaussian_fisher, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_eb_ucb_fisher <- apply(results_eb_ucb_fisher, 3, function(x){min(which(x[,1] < alpha))})
  
  total_samples_alpha_equal_product <- apply(results_alpha_equal_product, 3, get_total_sample_size, alpha = alpha)
  total_samples_alpha_gaussian_product <- apply(results_alpha_gaussian_product, 3, get_total_sample_size, alpha = alpha)
  total_samples_alpha_ucb_product <- apply(results_alpha_ucb_product, 3, get_total_sample_size, alpha = alpha)
  total_samples_alpha_equal_fisher <- apply(results_alpha_equal_fisher, 3, get_total_sample_size, alpha = alpha)
  total_samples_alpha_gaussian_fisher <- apply(results_alpha_gaussian_fisher, 3, get_total_sample_size, alpha = alpha)
  total_samples_alpha_ucb_fisher <- apply(results_alpha_ucb_fisher, 3, get_total_sample_size, alpha = alpha)
  total_samples_eb_equal_product <- apply(results_eb_equal_product, 3, get_total_sample_size, alpha = alpha)
  total_samples_eb_gaussian_product <- apply(results_eb_gaussian_product, 3, get_total_sample_size, alpha = alpha)
  total_samples_eb_ucb_product <- apply(results_eb_ucb_product, 3, get_total_sample_size, alpha = alpha)
  total_samples_eb_equal_fisher <- apply(results_eb_equal_fisher, 3, get_total_sample_size, alpha = alpha)
  total_samples_eb_gaussian_fisher <- apply(results_eb_gaussian_fisher, 3, get_total_sample_size, alpha = alpha)
  total_samples_eb_ucb_fisher <- apply(results_eb_ucb_fisher, 3, get_total_sample_size, alpha = alpha)
  
  power_frame <- data.frame("stop" = stopping_times_alpha_equal_product, "total_samples" = total_samples_alpha_equal_product, "simulation" = 1:n_sims, "allocation" = "equal", "combined" = "product", "martingale" = "alpha") %>%
    bind_rows(
      data.frame("stop" = stopping_times_alpha_gaussian_product, "total_samples" = total_samples_alpha_gaussian_product, "simulation" = 1:n_sims, "allocation" = "gaussian", "combined" = "product", "martingale" = "alpha"),
      data.frame("stop" = stopping_times_alpha_ucb_product, "total_samples" = total_samples_alpha_ucb_product, "simulation" = 1:n_sims, "allocation" = "ucb", "combined" = "product", "martingale" = "alpha"),
      data.frame("stop" = stopping_times_alpha_equal_fisher, "total_samples" = total_samples_alpha_equal_fisher, "simulation" = 1:n_sims, "allocation" = "equal", "combined" = "fisher", "martingale" = "alpha"),
      data.frame("stop" = stopping_times_alpha_gaussian_fisher, "total_samples" = total_samples_alpha_gaussian_fisher, "simulation" = 1:n_sims, "allocation" = "gaussian", "combined" = "fisher", "martingale" = "alpha"),
      data.frame("stop" = stopping_times_alpha_ucb_fisher, "total_samples" = total_samples_alpha_ucb_fisher, "simulation" = 1:n_sims, "allocation" = "ucb", "combined" = "fisher", "martingale" = "alpha"),
      data.frame("stop" = stopping_times_eb_equal_product, "total_samples" = total_samples_eb_equal_product, "simulation" = 1:n_sims, "allocation" = "equal", "combined" = "product", "martingale" = "eb"),
      data.frame("stop" = stopping_times_eb_gaussian_product, "total_samples" = total_samples_eb_gaussian_product, "simulation" = 1:n_sims, "allocation" = "gaussian", "combined" = "product", "martingale" = "eb"),
      data.frame("stop" = stopping_times_eb_ucb_product, "total_samples" = total_samples_eb_ucb_product, "simulation" = 1:n_sims, "allocation" = "ucb", "combined" = "product", "martingale" = "eb"),
      data.frame("stop" = stopping_times_eb_equal_fisher, "total_samples" = total_samples_eb_equal_fisher, "simulation" = 1:n_sims, "allocation" = "equal", "combined" = "fisher", "martingale" = "eb"),
      data.frame("stop" = stopping_times_eb_gaussian_fisher, "total_samples" = total_samples_eb_gaussian_fisher, "simulation" = 1:n_sims, "allocation" = "gaussian", "combined" = "fisher", "martingale" = "eb"),
      data.frame("stop" = stopping_times_eb_ucb_fisher, "total_samples" = total_samples_eb_ucb_fisher, "simulation" = 1:n_sims, "allocation" = "ucb", "combined" = "fisher", "martingale" = "eb")
    ) %>%
    mutate(overstatements = overstatements, understatements = understatements, reported_margin = round(reported_margin, 3), true_margin = round(true_margin, 3))
  power_frame
}


#allocation simulations
reported_tallies <- rbind(
  c(rep(0, 500), rep(1,500), rep(0, 400), rep(1, 600)),
  c(rep(0, 500), rep(1,500), rep(0, 450), rep(1, 550)),
  c(rep(0, 500), rep(1,500), rep(0, 490), rep(1, 510))
)

#CVR is correct in first instance and understates the margin. CVR is exactly correct in second. CVR overstates in third and 4th while outcome is correct, CVR is wrong in 5th AND the reported outcome is wrong.
#all are 2 vote overstatements (vote for loser recorded as vote for winner)
hand_tallies <- rbind(
  #c(rep(0, 500), rep(1,500), rep(0, 400), rep(1, 600)),
  c(rep(0, 500), rep(1,500), rep(0, 450), rep(1, 550)),
  c(rep(0, 500), rep(1,500), rep(0, 490), rep(1, 510)),
  c(rep(0, 500), rep(1,500), rep(0, 500), rep(1, 500))
)


outer_frames <- list()
for(i in 1:nrow(reported_tallies)){
  inner_frames <- list()
  for(j in 1:nrow(hand_tallies)){
    inner_frames[[j]] <- run_allocation_simulation(
      hand_tally = hand_tallies[j,],
      reported_tally = reported_tallies[i,],
      strata = c(rep(1,1000), rep(2,1000)), 
      n_sims = 300,
      alpha = 0.1
    )
    print(paste("i = ", i, ", j = ",j, sep=""))
  }
  outer_frames[[i]] <- inner_frames %>% reduce(bind_rows)
}

allocation_frame <- outer_frames %>% reduce(bind_rows) %>%
  mutate(finite_stop = ifelse(is.infinite(stop), 1000, stop)) %>%
  mutate(unconditional_total_samples = ifelse(is.na(total_samples), 2000, total_samples)) %>%
  #mutate(reported_margin = paste("Reported Margin =", reported_margin), true_margin = paste("True Margin =", true_margin)) %>%
  #mutate(allocation = recode(allocation, equal = "Proportional", threshold = "Thresholded", ucb = "UCB")) %>%
  as_tibble()


#save(allocation_frame, file = "allocation_power_frame_withnull")
load("allocation_power_frame_withnull")
ggplot(allocation_frame , aes(x = unconditional_total_samples, linetype = combined, color = allocation)) +
  stat_ecdf(size = 1.5) +
  facet_grid(true_margin ~ martingale + reported_margin) +
  xlim(0, 2000) +
  theme_bw() +
  theme(text = element_text(size = 18), axis.text = element_text(size = 14)) +
  labs(x = "Total Samples", y = "Cumulative probability of stopping", color = "Allocation Rule", linetype = "Combination")

allocation_frame %>% 
  group_by(
    martingale, 
    allocation, 
    combined, 
    reported_margin, 
    true_margin) %>% 
  summarize(
    expected_total_samples = mean(unconditional_total_samples), 
    median_total_samples = median(unconditional_total_samples)) %>% 
  filter(reported_margin == 0.1, true_margin == 0.01)

