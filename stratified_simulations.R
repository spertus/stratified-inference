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

power_frames <- apply(
  hand_tallies,
  MARGIN = 1,
  FUN = run_comparison_simulation, 
  reported_tally = reported_tally,
  strata = c(rep(1,10000), rep(2,10000)), 
  sample_sizes = matrix(round(10^seq(1, 4, length.out = 30)), nrow = 30, ncol = 2, byrow = F), 
  n_sims = 500,
  #sample_sizes = matrix(c(100, 100), nrow = 1, ncol = 2, byrow = F),
  #n_sims = 5,
  alpha = .1
) 
save(power_frames, file = "stratified_comparison_audit_frames")
