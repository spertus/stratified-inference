#script to run a bunch of stratified simulations and output the results
source("stratified_functions.R")


#comparison audits
reported_tallies <- rbind(
  c(rep(0, 5000), rep(1,5000), rep(0, 2500), rep(1, 7500)),
  c(rep(0, 5000), rep(1,5000), rep(0, 4000), rep(1, 6000)),
  c(rep(0, 5000), rep(1,5000), rep(0, 4500), rep(1, 5500)),
  c(rep(0, 5000), rep(1,5000), rep(0, 4900), rep(1, 5100))
  )

#CVR is correct in first instance and understates the margin. CVR is exactly correct in second. CVR overstates in third and 4th while outcome is correct, CVR is wrong in 5th AND the reported outcome is wrong.
#all are 2 vote overstatements (vote for loser recorded as vote for winner)
hand_tallies <- rbind(
  c(rep(0, 5000), rep(1,5000), rep(0, 2500), rep(1, 7500)),
  c(rep(0, 5000), rep(1,5000), rep(0, 3000), rep(1, 7000)),
  c(rep(0, 5000), rep(1,5000), rep(0, 4000), rep(1, 6000)),
  c(rep(0, 5000), rep(1,5000), rep(0, 4500), rep(1, 5500)),
  c(rep(0, 5000), rep(1,5000), rep(0, 4900), rep(1, 5100)),
  c(rep(0, 5000), rep(1,5000), rep(0, 4950), rep(1, 5050)),
  c(rep(0, 5000), rep(1,5000), rep(0, 5000), rep(1, 5000))
)

run_comparison_simulation <- function(reported_tally, hand_tally, strata, sample_sizes, n_sims = 300, alpha = .05){
  v <- 2 * mean(reported_tally) - 1
  u <- 1
  omegas <- reported_tally - hand_tally
  population <- (1 - omegas/u) / (2 - v/u)
  pop_range <- c(0, (1 + 1/u) / (2 - v/u))
  reported_margin <- 2*mean(reported_tally) - 1
  strata_reported_margins <- tapply(reported_tally, strata, function(x){2 * mean(x) - 1})
  true_margin <- 2*mean(hand_tally) - 1
  understatements <- mean(population == 1)
  overstatements <- mean(population == 0)
  eta_0 <- choose_eta_0(audit_type = "comparison", diluted_margin = strata_reported_margins, overstatement_rate = .001)
  pars <- list(d = 100, eta_0 = eta_0)
  
  power_frame <- run_stratified_simulation(
    population = population, 
    strata = strata, 
    sample_sizes = sample_sizes,
    mu_0 = 0.5,
    n_sims = n_sims, 
    alpha = alpha,
    pars = pars,
    bounds = pop_range)
  
  power_frame <- power_frame %>% 
    mutate(overstatements = overstatements, understatements = understatements, reported_margin = round(reported_margin, 3), true_margin = round(true_margin, 3))
}

power_frames <- list()
for(i in 1:nrow(reported_tallies)){
  inner_frames <- list()
  for(j in 1:nrow(hand_tallies)){
    inner_frames[[j]] <- run_comparison_simulation(
      hand_tally = hand_tallies[j,],
      reported_tally = reported_tallies[i,],
      strata = c(rep(1,10000), rep(2,10000)), 
      sample_sizes = matrix(round(10^seq(1, 4, length.out = 30)), nrow = 30, ncol = 2, byrow = F), 
      n_sims = 300,
      alpha = 0.1
    )
  }
  power_frames <- inner_frames %>% reduce(bind_rows)
}

save(power_frames, file = "stratified_comparison_audit_frames")




########## ALPHA allocation rules ############
#helper function to get sample sizes from array
get_total_sample_size <- function(x){
  if(any(x[,1] < alpha)){
    sum(x[1:min(which(x[,1] < alpha)),2:3])
  } else{
    NA
  }
}

run_allocation_simulation <- function(reported_tally, hand_tally, strata, n_sims = 300, alpha = .05){
  v <- 2 * mean(reported_tally) - 1
  u <- 1
  omegas <- reported_tally - hand_tally
  population <- (1 - omegas/u) / (2 - v/u)
  pop_range <- c(0, (1 + 1/u) / (2 - v/u))
  stratum_1 <- population[strata == 1]
  stratum_2 <- population[strata == 2]
  
  reported_margin <- 2*mean(reported_tally) - 1
  strata_reported_margins <- tapply(reported_tally, strata, function(x){2 * mean(x) - 1})
  true_margin <- 2*mean(hand_tally) - 1
  understatements <- mean(population == 1)
  overstatements <- mean(population == 0)
  
  results_equal <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, replace = FALSE, rule = "equal", resolution = length(population)))
  results_threshold <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, replace = FALSE, rule = "hard_threshold", resolution = length(population)))
  
  stopping_times_equal <- apply(results_equal, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_threshold <- apply(results_threshold, 3, function(x){min(which(x[,1] < alpha))})
  
  total_samples_equal <- apply(results_equal, 3, get_total_sample_size)
  total_samples_threshold <- apply(results_threshold, 3, get_total_sample_size)
  
  power_frame <- data.frame("stop" = stopping_times_equal, "total_samples" = total_samples_equal, "simulation" = 1:n_sims, "allocation" = "equal") %>%
    bind_rows(
      data.frame("stop" = stopping_times_threshold, "total_samples" = total_samples_threshold, "simulation" = 1:n_sims, "allocation" = "threshold")
      ) %>%
    mutate(overstatements = overstatements, understatements = understatements, reported_margin = round(reported_margin, 3), true_margin = round(true_margin, 3))
}


#allocation simulations
reported_tallies <- rbind(
  c(rep(0, 500), rep(1,500), rep(0, 300), rep(1, 700)),
  c(rep(0, 500), rep(1,500), rep(0, 400), rep(1, 600)),
  c(rep(0, 500), rep(1,500), rep(0, 490), rep(1, 510))
)

#CVR is correct in first instance and understates the margin. CVR is exactly correct in second. CVR overstates in third and 4th while outcome is correct, CVR is wrong in 5th AND the reported outcome is wrong.
#all are 2 vote overstatements (vote for loser recorded as vote for winner)
hand_tallies <- rbind(
  c(rep(0, 500), rep(1,500), rep(0, 200), rep(1, 800)),
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
        strata = c(rep(1,100), rep(2,100)), 
        n_sims = 1,
        alpha = 0.1
    )
  }
  outer_frames[[i]] <- inner_frames %>% reduce(bind_rows)
}

allocation_frame <- outer_frames %>% reduce(bind_rows) %>%
  mutate(finite_stop = ifelse(is.infinite(stop), 100, stop)) %>%
  mutate(unconditional_total_samples = ifelse(is.na(total_samples), 200, total_samples)) %>%
  mutate(reported_margin = paste("Reported Margin =", reported_margin), true_margin = paste("True Margin =", true_margin)) %>%
  as_tibble()
    
save(allocation_frame, file = "allocation_power_frame")

# ggplot(allocation_frame, aes(x = total_samples, linetype = allocation, color = allocation)) +
#   stat_ecdf() +
#   facet_grid(true_margin ~ reported_margin)

