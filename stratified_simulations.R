#script to run a bunch of stratified simulations and output the results
source("stratified_functions.R")


#comparison audits
reported_tally <- matrix(c(rep(0, 5000), rep(1,5000), rep(0, 2500), rep(1, 7500)), nrow = 1, ncol = 20000, byrow = TRUE)
#CVR is correct in first instance and understates the margin. CVR is exactly correct in second. CVR overstates in third and 4th while outcome is correct, CVR is wrong in 5th AND the reported outcome is wrong.
#all are 2 vote overstatements (vote for loser recorded as vote for winner)
hand_tallies <- rbind(
  c(rep(0, 5000), rep(1,5000), rep(0, 2500), rep(1, 7500)),
  c(rep(0, 5000), rep(1,5000), rep(0, 3000), rep(1, 7000)),
  c(rep(0, 5000), rep(1,5000), rep(0, 4000), rep(1, 6000)),
  c(rep(0, 5000), rep(1,5000), rep(0, 4500), rep(1, 5500)),
  c(rep(0, 5000), rep(1,5000), rep(0, 5000), rep(1, 5000))
)


run_comparison_simulation <- function(hand_tally, reported_tally, strata, sample_sizes, n_sims = 300, alpha = .05){
  v <- 2 * mean(reported_tally) - 1
  u <- 1
  omegas <- reported_tally - hand_tally
  population <- (1 - omegas/u) / (2 - v/u)
  pop_range <- c(0, (1 + 1/u) / (2-v/u))
  mu_0 <- 0.5
  reported_margin <- 2*mean(reported_tally) - 1
  true_margin <- 2*mean(hand_tally) - 1
  understatements <- mean(population == 1)
  overstatements <- mean(population == 0)
  
  power_ppm_unstrat <- matrix(NA, nrow = length(sample_sizes), ncol = 1)
  power_ppm_fisher <- matrix(NA, nrow = length(sample_sizes), ncol = 1)
  power_eb_unstrat <- matrix(NA, nrow = length(sample_sizes), ncol = 1)
  power_eb_fisher <- matrix(NA, nrow = length(sample_sizes), ncol = 1)
  power_eb_martingale <- matrix(NA, nrow = length(sample_sizes), ncol = 1)
  power_hedged_unstrat <- matrix(NA, nrow = length(sample_sizes), ncol = 1)
  power_hedged_fisher <- matrix(NA, nrow = length(sample_sizes), ncol = 1)
  power_hedged_martingale <- matrix(NA, nrow = length(sample_sizes), ncol = 1)

  for(i in 1:length(sample_sizes)){
    replicates_ppm_unstrat <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = rep(1,length(population)), mu_0 = mu_0, n = sample_sizes[i], method = "beta-binomial", bounds = pop_range))
    replicates_ppm_fisher <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = strata, mu_0 = mu_0, n = c(sample_sizes[i]/2,sample_sizes[i]/2), method = "beta-binomial", pool = "fisher", bounds = pop_range))
    replicates_eb_unstrat <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = rep(1,length(population)), mu_0 = mu_0, n = sample_sizes[i], method = "empirical_bernstein", bounds = pop_range))
    replicates_eb_fisher <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = strata, mu_0 = mu_0, n = c(sample_sizes[i]/2,sample_sizes[i]/2), method = "empirical_bernstein", pool = "fisher", bounds = pop_range))
    replicates_eb_martingale <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = strata, mu_0 = mu_0, n = c(sample_sizes[i]/2,sample_sizes[i]/2), method = "empirical_bernstein", pool = "martingale", bounds = pop_range))
    replicates_hedged_unstrat <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = rep(1,length(population)), mu_0 = mu_0, n = sample_sizes[i], method = "hedged", bounds = pop_range))
    replicates_hedged_fisher <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = strata, mu_0 = mu_0, n = c(sample_sizes[i]/2,sample_sizes[i]/2), method = "hedged", pool = "fisher", bounds = pop_range))
    replicates_hedged_martingale <- replicate(n = n_sims, get_stratified_pvalue(population = population, strata = strata, mu_0 = mu_0, n = c(sample_sizes[i]/2,sample_sizes[i]/2), method = "hedged", pool = "martingale", bounds = pop_range))
    
    power_ppm_unstrat[i,1] <- mean(replicates_ppm_unstrat < alpha)
    power_ppm_fisher[i,1] <- mean(replicates_ppm_fisher < alpha)
    power_eb_unstrat[i,1] <- mean(replicates_eb_unstrat < alpha)
    power_eb_fisher[i,1] <- mean(replicates_eb_fisher < alpha)
    power_eb_martingale[i,1] <- mean(replicates_eb_martingale < alpha)
    power_hedged_unstrat[i,1] <- mean(replicates_hedged_unstrat < alpha)
    power_hedged_fisher[i,1] <- mean(replicates_hedged_fisher < alpha)
    power_hedged_martingale[i,1] <- mean(replicates_hedged_martingale < alpha)
  }
  
  colnames(power_ppm_unstrat) <- "power_ppm_unstrat"
  colnames(power_ppm_fisher) <- "power_ppm_strat"
  colnames(power_eb_fisher) <- "power_eb_fisher"
  colnames(power_eb_martingale) <- "power_eb_martingale"
  colnames(power_eb_unstrat) <- "power_eb_unstrat"
  colnames(power_hedged_fisher) <- "power_hedged_fisher"
  colnames(power_hedged_martingale) <- "power_hedged_martingale"
  colnames(power_hedged_unstrat) <- "power_hedged_unstrat"
    

  
  
  
  power_frame <- data.frame(
    power_eb_fisher,
    power_eb_martingale,
    power_eb_unstrat, 
    power_hedged_fisher,
    power_hedged_martingale,
    power_hedged_unstrat,
    power_ppm_unstrat, 
    power_ppm_fisher, 
    sample_sizes) %>%
    pivot_longer(cols = starts_with("power"), names_to = "design", values_to = "power", names_prefix = "power_") %>%
    mutate(overstatements = overstatements, understatements = understatements, reported_margin = reported_margin, true_margin = true_margin)
  power_frame
} 


power_frames <- apply(
  hand_tallies,
  MARGIN = 1,
  FUN = run_comparison_simulation, 
  reported_tally = reported_tally,
  strata = c(rep(1,10000), rep(2,10000)), 
  sample_sizes = round(10^seq(1, 4, length.out = 30)), 
  n_sims = 500,
  alpha = .05
) 
save(power_frames, file = "comparisonaudit_power_frames")
