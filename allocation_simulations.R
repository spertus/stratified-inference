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
  strata_max <-(1 + 1/assorter_bound) / (2 - v/assorter_bound)
  stratum_1 <- (1 - omegas[strata == 1] / assorter_bound) / (2 - v[1] / assorter_bound)
  stratum_2 <- (1 - omegas[strata == 2] / assorter_bound) / (2 - v[2] / assorter_bound)
  
  reported_margin <- 2*mean(reported_tally) - 1
  
  true_margin <- 2*mean(hand_tally) - 1
  understatements <- mean(population == 1)
  overstatements <- mean(population == 0)
  d <- c(20, 20)
  eta_0 <- 1 / (2 - v / assorter_bound)

  
  
  results_equal <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, eta_0 = eta_0, u = strata_max, replace = FALSE, rule = "equal"))
  results_threshold <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, eta_0 = eta_0, u = strata_max, replace = FALSE, rule = "hard_threshold"))
  results_ucb <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, eta_0 = eta_0, u = strata_max, replace = FALSE, rule = "confidence_bound"))
  
  stopping_times_equal <- apply(results_equal, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_threshold <- apply(results_threshold, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_ucb <- apply(results_ucb, 3, function(x){min(which(x[,1] < alpha))})
  
  total_samples_equal <- apply(results_equal, 3, get_total_sample_size, alpha = alpha)
  total_samples_threshold <- apply(results_threshold, 3, get_total_sample_size, alpha = alpha)
  total_samples_ucb <- apply(results_ucb, 3, get_total_sample_size, alpha = alpha)
  
  power_frame <- data.frame("stop" = stopping_times_equal, "total_samples" = total_samples_equal, "simulation" = 1:n_sims, "allocation" = "equal") %>%
    bind_rows(
      data.frame("stop" = stopping_times_threshold, "total_samples" = total_samples_threshold, "simulation" = 1:n_sims, "allocation" = "threshold")
    ) %>%
    bind_rows(
      data.frame("stop" = stopping_times_ucb, "total_samples" = total_samples_ucb, "simulation" = 1:n_sims, "allocation" = "ucb")
    ) %>%
    mutate(overstatements = overstatements, understatements = understatements, reported_margin = round(reported_margin, 3), true_margin = round(true_margin, 3))
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
  c(rep(0, 500), rep(1,500), rep(0, 400), rep(1, 600)),
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
  }
  outer_frames[[i]] <- inner_frames %>% reduce(bind_rows)
}

allocation_frame <- outer_frames %>% reduce(bind_rows) %>%
  mutate(finite_stop = ifelse(is.infinite(stop), 1000, stop)) %>%
  mutate(unconditional_total_samples = ifelse(is.na(total_samples), 2000, total_samples)) %>%
  mutate(reported_margin = paste("Reported Margin =", reported_margin), true_margin = paste("True Margin =", true_margin)) %>%
  as_tibble()


save(allocation_frame, file = "allocation_power_frame")
# load("allocation_power_frame")
# ggplot(allocation_frame, aes(x = total_samples, linetype = allocation, color = allocation)) +
#   stat_ecdf() +
#   facet_grid(true_margin ~ reported_margin)
