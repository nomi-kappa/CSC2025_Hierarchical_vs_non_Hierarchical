# Simulation Demo (Hierarchical vs Non-Hierarchical RT Models)

## Overview
This folder contains simulation code and results for hierarchical vs non-hierarchical Stan models applied to lexical decision RT data.

The code was modified/adapted by Yi-Chen (Eden) Hsieh and Nomiki Koutsoumpari based on the original OSF materials from Boehm (2018).

---

## Files

- **Sim_Demo.R**  
  Original simulation demo script from Boehm (2018).

- **Sim_Demo_updated.R**  
  Updated version with reproducibility improvements, saving options, and diagnostics.

- **Delta_recovery.R**  
  Utility script to check recovery of the delta parameter (single run and across runs).

- **full model**  
  Stan file: hierarchical RT model.

- **non-hier model**  
  Stan file: non-hierarchical RT model (subject means).

- **LexDatParams.RData**  
  Group-level parameters estimated from Wagenmakers et al. (2008) data.

### Output (`out/` folder)
- **fit_hier-XX.rds** — Stan fit objects for hierarchical model (run XX).  
- **fit_nonhier-XX.rds** — Stan fit objects for non-hierarchical model (run XX).  
- **S-demo-XX.RData** — Saved simulation results for run XX.  
- **S-demo-all.RData** — Aggregated results across runs.  
- Example: `*-02.*` files come from demonstration run with `iSim = 2`.

---

## How to Use

1. **Run simulations with the updated script**  
   Open `Sim_Demo_updated.R` and adjust simulation parameters (e.g., `N`, `k`, `NSim`).  
   
2. **Inspect fits**  
   ```r
   ft <- readRDS("out/fit_hier-02.rds")
   print(ft, pars = c("delta","sigma_m","scl","mu_m"))
   stan_trace(ft, pars = c("delta","sigma_m"))
   check_hmc_diagnostics(ft)
   ```
3. **Check delta recovery**  
   Open `Delta_recovery.R` and run a single run recovery or/and across runs recovery.  
   

# Main differences between the models: 

Main differences
    Feature	                    Non-hierarchical	              Full hierarchical
    Data used	                  subject means	                  trial-level RTs
    Subject-level parameters	  None	                          m[i, cond], s[i, cond]
    Partial pooling	            No	                            Yes (subject means borrow strength from group mean)
    Trial-level variability	    Ignored	                        Explicitly modeled (s[i, cond])
    Complexity	                Simple	                        More complex, but more informative
