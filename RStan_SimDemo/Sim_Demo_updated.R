## ------------------------------------------------------------
## Simulating log-RT data and fitting with "hierarchical / non-hierarchical" Stan models
## - set.seed() for reproducibility
## - Options: save each run separately / aggregate into one list
## ------------------------------------------------------------

### Original Code: Udo Boehm
### Adapted Code: Nomiki Koutsoumpari from Yi-Chen (Eden) Hsieh
### CSC 2025 with Jocelyn Halim & Marcell Székely



### Sample means and sds of log-transformed RTs of simulated pps from distributions fitted to real lexical decision data ----

library(rstan)

## Load group-level distribution parameters estimated from real data (contains mu.m, sigma.m, sigma.s) from the Wagenmakers et al., 2008 data
## Ensure the file path is correct
setwd("~/Downloads/RStan_SimDemo")

load("LexDatParams.RData")   # load the objects
ls()                         # list what got loaded
class(mu.m)
str(mu.m)
    # mu.m mean of participant-level means.In hierarchical RT modeling, 
            # each participant has their own average log-RT, and mu.m captures the overall group mean across participants.

    # sigma.m The sd of participant-level means. Tells you how much variability there is between participants’ mean log-RTs.
            # In other words: how spread out the individuals’ average RTs are around mu.m.
  
    # sigma.s The sd of participant-level sds. Each participant not only has a mean RT, but also their own
           # within-person variability (some participants are more consistent, some less).
           # sigma.s quantifies how much these within-person variabilities differ across participants.



## Simulation parameters: parameters are kept in one list (can be modified to whatever you like) ; $ assigns a value here
params <- list()                # creates empty list called params
params$mu1.m   <- mu.m          # group mean for condition 1 (log-RT) (using real data here); as the baseline average log-RT across participants.
params$delta   <- 0             # standardised effect size: condition 2 = mu1.m + delta * sigma.m; Here, delta = 0 → condition 2 is the same as condition 1 (no effect).
                                    # If you set delta = 0.5, condition 2 would be shifted by half a standard deviation.
params$sigma.m <- sigma.m       # group-level SD of subject means (between-subject variability) (using real data here); 
                                    # Not all participants have the same average log-RT → sigma.m models how much participants differ in their means.
params$sigma.s <- sigma.s       # group-level SD of trial-level SDs (between-subject variability in within-subject noise) (using real data here)
params$k       <- 5             # number of trials per subject per condition
params$N       <- 30            # number of subjects
params$NSim    <- 2             # number of repeated simulations

## Saving options (choose one or both) # at the moment they’re just named entries in the params list.
params$save_each_file   <- TRUE   # TRUE: save one .RData  and .rds file per simulation run
params$save_aggregate   <- TRUE   # TRUE: collect all runs in a list (.RData) and save once at the end

## Stan model file names
hier_model_file <- "full model.stan" 
   # file.edit("full model.stan") # open in editor
nonhier_model_file <- "non-hier model.stan"


## Prepare container if aggregate saving is enabled
if (isTRUE(params$save_aggregate)) {
  results_list <- vector("list", params$NSim) # Creates an empty list with params$NSim slots.
}

results_list # at the moment the 2 sims are NULL



### SIMULATION ###

## Main loop: simulate -> fit hierarchical -> fit non-hierarchical -> save
        # Overall purpose: for each simulation (iSim), the script
           # creates group/subject-level parameters,
           # generates trial-level log-RTs for two conditions,
           # packs the data into a 3-D array rt[condition, trial, subject],
           # builds dat_hier (a list) that you pass to the hierarchical Stan model.

for (iSim in 1:params$NSim) { # iSim: the loop index variable that you introduce in the for loop; it is temp
  
  message(sprintf("=== Simulation %d / %d ===", iSim, params$NSim)) # prints which simulation run you’re on
  
  ##  1) Generate simulated data -----
  params$mu2.m <- params$mu1.m + params$delta * params$sigma.m  # group mean for condition 2
                                                                # mu1.m is the group mean (condition 1).
                                                                # delta is a standardized effect size; mu2.m becomes the group mean for condition 2 = mu1.m + delta * sigma.m.
                                                                # If delta = 0, mu2.m == mu1.m (no difference between conditions). We can take it as a control maybe.
                                                                # The mean of condition 2 is the mean of condition 1, shifted by delta SDs.
                                                                # This is a common way in simulations to parameterize group differences
  
  ## Subject-level means (two conditions) - (one mean per subject per condition)
        # M1 and M2 are numeric vectors of length N (number of subjects).
        # They represent each subject’s true mean log-RT in condition 1 and 2 respectively.
        # They are sampled from a normal distribution (rnorm) with group mean mu1.m/mu2.m and between-subject SD sigma.m.
  params$M1 <- rnorm(params$N, mean = params$mu1.m, sd = params$sigma.m)
  params$M2 <- rnorm(params$N, mean = params$mu2.m, sd = params$sigma.m)
  
  ## Subject-level trial SDs (half-normal via abs(N(0, sigma.s)))
        # Each subject has their own trial-to-trial variability (a standard deviation).
        # rnorm(..., mean=0, sd=sigma.s) draws from N(0, sigma.s). Taking abs() makes this a half-normal (positive SD).
        # SIGMA1 and SIGMA2 are length-N vectors giving per-subject within-subject SD for condition 1 and 2.
     
        # Why mean = 0?
        # We’re not simulating data points here, but variability values (the subject-specific SDs).
        # If we just sampled from N(mu, sd), some draws could be negative — which makes no sense for a SD.
        # So the trick is:
        # Draw from N(0, sigma.s) (centered at 0).
        # Then take the absolute value → ensures all SDs are positive.
        # This is equivalent to sampling from a half-normal distribution.
  params$SIGMA1 <- abs(rnorm(params$N, mean = 0, sd = params$sigma.s)) # abs: converts the symmetric bell curve into a distribution that only lives on the positive side.
  params$SIGMA2 <- abs(rnorm(params$N, mean = 0, sd = params$sigma.s))
  
  ## Trial-level data (log-RTs) for each subject, for each condition
  simdat <- NULL # will store the data
  for (isubj in 1:params$N) {
    x <- rnorm(params$k, mean = params$M1[isubj], sd = params$SIGMA1[isubj])
    y <- rnorm(params$k, mean = params$M2[isubj], sd = params$SIGMA2[isubj])
    simdat <- rbind(simdat, data.frame(  # Combines this subject’s data with the previous subjects’ data.
      subj = isubj, x = x, y = y, # Since x and y are vectors of length k, the data frame has k rows per subject.
      m1 = params$M1[isubj], m2 = params$M2[isubj],
      s1 = params$SIGMA1[isubj], s2 = params$SIGMA2[isubj]
    ))
  }
  
  ## Stan hierarchical model requires an array: rt[condition, trial, subject] what stan needs
  rt <- array(NA_real_, dim = c(2, params$k, params$N)) # creates a 3d dimensional array of NAs (numeric type); dim = dimension; so rt[cond, trial, subj] will give the log-RT for a given condition, trial, and subject.
                                                        # This is the format Stan expects for hierarchical models.
  ## use unstack to convert to wide format (ensure subj is consecutive 1..N); long = one trial for one subj; wide (what stan needs) = rows= trial1 columns = Subj1 subj 2 (one matrix per condition)
            # assign to the 3d array; rt[cond, trial, subj] lets you access the log-RT for a specific condition, trial, and subject.
  rt[1, , ] <- as.matrix(unstack(simdat, x ~ subj)) # x ~ subj means: “take the x column (condition 1 trials) and spread it by subj; unstack() groups all rows with the same subj together and makes one column per subject, where rows = trials.
  rt[2, , ] <- as.matrix(unstack(simdat, y ~ subj)) # same as above
  dat_hier <- list(rt = rt, k = params$k, n = params$N)
  
  ##  2) Fit hierarchical Stan model -----
  mod.params.hier <- c("m", "s", "sigma_m", "delta", "mu_m", "scl") # A character vector of parameters to monitor in Stan and it will return posterior samples for them.
  
  # Set initial fitting parameters through a function for Stan's MCMC sampler; 
     # random initialization avoids all chains starting in the exact same place.avoid divergences, pick values that are pausible guven simulated data above.
     # m and s → subject-level parameters
     # mu_m1, delta, sigma_m, scl → group-level / hyperparameters
    # In hierarchical models for two conditions, it’s common to parameterize the second condition relative to the first; that is why no mu_m2 here
  inits_hier <- function() {  
    list(
      mu_m1   = runif(1, 5, 7), # initial guess for the group mean of condition 1; uniform random between 5 and 7.
      delta   = rnorm(1, 0, 1), # Initial guess for effect size (difference between conditions); normal around 0.
      sigma_m = runif(1, 0.001, 4.999), # Initial guess for between-subject SD of means; small positive values allowed.
      scl     = runif(1, 0.001, 9.999), # Initial guess for trial-level SD scaling factor.
      m       = cbind(runif(params$N, 5.5, 6.5), runif(params$N, 5.5, 6.5)), # Initial subject-level means for each condition (matrix N×2).
      s       = matrix(runif(2 * params$N, 0.001, 0.999), ncol = 2, byrow = FALSE) # Initial subject-level SDs for each condition (matrix N×2).
    )
  }
  
  # Fit model using stan model
  fit_hier <- stan(
    file    = hier_model_file,
    data    = dat_hier,
    pars    = mod.params.hier,
    chains  = 3,
    warmup  = 2000, # burn-in iterations per chain (these are discarded)
    thin    = 4,    # Keep every 4th sample (reduces storage & autocorrelation).
    iter    = 20000,
    init    = inits_hier,
    cores   = 3
  )
  
  post_hier <- rstan::extract(fit_hier)
  
  mean(post_hier$delta)            # posterior mean of effect size
  quantile(post_hier$delta, c(0.025, 0.975))  # 95% credible interval
  
  
  ## Posterior summaries
  hier.out <- list(
    mean = list(
      sigma_m = mean(post_hier$sigma_m), # overall posterior mean
      delta   = mean(post_hier$delta), # posterior mean effect
      scl     = mean(post_hier$scl),
      mu_m    = colMeans(post_hier$mu_m), # mean across iterations for each condition
      m       = apply(post_hier$m, c(2, 3), mean), # dimensions 2 & 3 = subjects × conditions; mean over iterations; posterior draws for subject-level means.
      s       = apply(post_hier$s, c(2, 3), mean)  # mean over iterations, dimensions = subjects × conditions; posterior draws for subject-level SDs
    ),
    median = list(
      sigma_m = median(post_hier$sigma_m),
      delta   = median(post_hier$delta),
      scl     = median(post_hier$scl),
      mu_m    = apply(post_hier$mu_m, 2, median),
      m       = apply(post_hier$m, c(2, 3), median),
      s       = apply(post_hier$s, c(2, 3), median)
    )
  )
  
  delta.hier.post <- as.numeric(post_hier$delta)
  
  ## True subject-level means (for comparison)
  od <- cbind( # Combines the two vectors into a matrix
        # tapply() applies a function to subsets of a vector, grouped by some factor
    tapply(simdat$m1, simdat$subj, unique),# Groups simdat$m1 by subj; (simdat = subject-level mean)
                                           # unique returns the unique value for each subject (since all trials for the same subject have the same m1)
                                           # Result: a vector of true means for each subject for condition 1
    tapply(simdat$m2, simdat$subj, unique)
  )
  colnames(od) <- c("true_m1", "true_m2") # Labels the matrix columns
  
  ## ----- 3) Fit non-hierarchical Stan model (using subject-level means only) -----
  m.rt <- array(NA_real_, dim = c(params$N, 2))
    # tapply(..., simdat$subj, mean) = compute average RT per subject for each condition
    # Result: m.rt is a N × 2 matrix, each row = one subject, each column = one condition
    # ✅ This converts the trial-level data to subject means only, which is what the non-hierarchical model expects.
  m.rt[, 1] <- tapply(simdat$x, simdat$subj, mean)# simdat$x = trial-level RTs for condition 1
  m.rt[, 2] <- tapply(simdat$y, simdat$subj, mean)# simdat$y = trial-level RTs for condition 2
    
    # Prepare for Stan; it expects a list of data objects.
  dat_nonhier <- list(m_rt = m.rt, n = params$N) # m_rt = subject means; n=number of subjects
    
    # parameteres to extract
  mod.params.nonhier <- c("mu_m", "delta", "sigma_m")
    
    # initial values for Stan
  inits_nonhier <- function() {
    list(
      mu_m1   = runif(1, 5, 7), # initial guess for condition 1 mean
      delta   = rnorm(1, 0, 1), # initial guess for effect size
      sigma_m = runif(1, 0.001, 4.999) # initial gues for btw-subjects SD
    )
  }
  
  fit_nonhier <- stan( # Calls Stan to run MCMC sampling
    file    = nonhier_model_file,
    data    = dat_nonhier,
    pars    = mod.params.nonhier,
    chains  = 3,
    warmup  = 500,
    thin    = 1,
    iter    = 5000,
    init    = inits_nonhier,
    cores   = 3
  )
  
  post_nonhier <- rstan::extract(fit_nonhier) # list of posterior draws for mu_m, delta, sigma_m
  delta.parthier.post <- as.numeric(post_nonhier$delta) # posterior draws for the effect size
  
  ## 4) Saving -----
  
  ## 4A) Save each run separately (filename includes iSim)
  if (isTRUE(params$save_each_file)) {
    file_i <- sprintf("out/S-demo-%02d.RData", iSim)
    save(list = c("params", "simdat", "od", "hier.out",
                  "delta.hier.post", "delta.parthier.post"),
         file = file_i)
    file_fit_hier <- sprintf("out/fit_hier-%02d.rds", iSim)
    saveRDS(fit_hier, file_fit_hier) # saves individual objects as .rds files (can be read later with readRDS)
    file_fit_nonhier <- sprintf("out/fit_nonhier-%02d.rds", iSim)
    saveRDS(fit_nonhier, file_fit_nonhier)
  }
  
  ## 4B) Save aggregated list (to be written once after loop)
  if (isTRUE(params$save_aggregate)) {
    results_list[[iSim]] <- list( # For each simulation run, we store a list of all important outputs into
      params               = params,
      simdat               = simdat,
      od                   = od,
      hier_summary         = hier.out,
      delta_hier_post      = delta.hier.post,
      delta_nonhier_post   = delta.parthier.post
    )
  }
}

## After loop: if aggregate saving is enabled, save once
if (isTRUE(params$save_aggregate)) {
  save(results_list, file = "out/S-demo-all.RData")
}

message("All simulations completed")



## ------------------------------------------------------------
## Diagnostics for each simulation run
## ------------------------------------------------------------

## Choose which simulation to check
iSim <- 1   # change to 2, 3, ... for other runs

## Load the hierarchical model stan fit for this simulation
ft <- readRDS(sprintf("out/fit_hier-%02d.rds", iSim))

# For the non-hierarchical model, replace with
 ftnh <- readRDS(sprintf("out/fit_nonhier-%02d.rds", iSim))
## and use pars = c("delta","sigma_m","mu_m") (since there’s no scl).
 

## 1) Print parameter summary (posterior means, sds, Rhat, n_eff, etc.)
print(ft, pars = c("delta","sigma_m","scl","mu_m"))
print(ftnh, pars = c("delta", "sigma_m", "mu_m"))


## 2) Trace plots for convergence check (choose parameters you want to plot)
trace_ft <-stan_trace(ft, pars = c("delta","sigma_m"));trace_ft
trace_ftnh <-stan_trace(ftnh, pars = c("delta","sigma_m"));trace_ftnh


## 3) HMC diagnostics (divergent transitions, treedepth, energy, etc.)
check_hmc_diagnostics(ft)
check_hmc_diagnostics(ftnh)


## Posterior distribution of delta

library(ggplot2)

# Extract posterior draws for delta
delta_post <- as.numeric(rstan::extract(ft)$delta)

# Make a data frame for ggplot
df <- data.frame(delta = delta_post)

# Plot posterior
ggplot(df, aes(x = delta)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = mean(delta_post), color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = quantile(delta_post, c(0.025, 0.975)), color = "darkblue", linetype = "dotted", size = 0.8) +
  labs(
    title = "Posterior distribution of delta (effect size)",
    x = expression(delta),
    y = "Density"
  ) +
  theme_minimal()

# Most posterior mass is positive, consistent with a small effect
# Interval includes some negative values → uncertainty, but the posterior leans positive
# Visualizes exactly what the table was showing numerically



# Extract posterior draws for both models(effect size)
delta_hier <- as.numeric(rstan::extract(ft)$delta)
delta_nonhier <- as.numeric(rstan::extract(ftnh)$delta)

# Combine into one data frame
df <- data.frame(
  delta = c(delta_hier, delta_nonhier),
  model = rep(c("Hierarchical", "Non-hierarchical"),
              times = c(length(delta_hier), length(delta_nonhier)))
)

# Plot
ggplot(df, aes(x = delta, fill = model, color = model)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = c(mean(delta_hier), mean(delta_nonhier)),
             linetype = "dashed", size = 0.8, color = c("blue", "red")) +
  labs(
    title = "Posterior distributions of delta (effect size)",
    x = expression(delta),
    y = "Density",
    fill = "Model",
    color = "Model"
  ) +
  theme_minimal()


# Non-hierarchical posterior is narrower (smaller SD = 0.25) → density peak is higher
# Hierarchical posterior is wider (SD = 0.35) → density peak is lower

# ✅ This can happen because:
  # Your hierarchical model includes trial-level variation, which adds uncertainty to the estimate of delta.More realistic
  # Non-hierarchical model uses only subject means, ignoring trial-level noise → delta appears more “precisely” estimated, hence narrower posterior.

