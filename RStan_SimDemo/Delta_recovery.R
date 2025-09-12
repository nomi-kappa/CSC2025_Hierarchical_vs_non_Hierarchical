## -------------------------------------------------
## Delta recovery  — from saved files
## -------------------------------------------------

### Original code was generated with ChatGPT, 
### then modified and adapted by Yi-Chen (Eden) Hsieh % Nomiki Koutsoumpari

library(ggplot2)
library (ggridges)

###  Check delta recovery from a single .RData file ----
delta_check_from_file <- function(file_path, true_delta = NULL, title_suffix = NULL) {
  e <- new.env(parent = emptyenv()) # Creates a temporary environment e to load the .RData file.
  load(file_path, envir = e) # loads all objects from the RData file
  if (!all(c("params","delta.hier.post") %in% ls(e))) {
    stop("Th efile should include `params` & `delta.hier.post`")
  }
  post <- as.numeric(e$delta.hier.post) # converts delta.hier.post to a numeric vector of posterior draws.
  dtrue <- if (is.null(true_delta)) e$params$delta else true_delta # determines the “true delta” to check recovery
                                                              # f true_delta argument is given, use that
                                                              # Otherwise, read delta from the simulation parameters in the file (e$params$delta)
  # Compute summary statistics
  mean_post <- mean(post)
  median_post <- median(post)
  sd_post <- sd(post)
  ci <- quantile(post, c(.025,.975), names = FALSE)
  covered <- as.numeric(dtrue >= ci[1] & dtrue <= ci[2]) # is the true delta within the 95% CI of the posterior?
  
  # Summary
  tab <- data.frame(
    true = dtrue,
    post_mean = mean_post,
    post_median = median_post,
    post_sd = sd_post,
    ci_lo = ci[1],
    ci_hi = ci[2],
    covered_95 = covered
  ) 

  # Density
  p <- ggplot(data.frame(post = post), aes(x = post)) +
    geom_density() +
    geom_vline(xintercept = dtrue, linetype = 2) +
    annotate("rect", xmin = ci[1], xmax = ci[2], ymin = -Inf, ymax = Inf, alpha = 0.08) +
    labs(title = paste0("Delta posterior (single)",
                        if (!is.null(title_suffix)) paste0(" — ", title_suffix) else ""),
         x = expression(delta), y = "Density",
         subtitle = sprintf("mean=%.3f, 95%% CI [%.3f, %.3f]; covered=%s",
                            mean_post, ci[1], ci[2],
                            ifelse(covered==1,"YES","NO"))) +
    theme_minimal(base_size = 13)
  
  list(summary_table = tab, plot = p)
}




 
###  Check delta recovery from all runs (results_list) ----
 delta_check_from_results_list <- function(results_list,
                                           title = "Delta posterior across runs",
                                           view = c("overlay", "facets", "ridgeline")) {
   view <- match.arg(view)
   
   # Summary & coverage
       # loop over all runs
   rows <- lapply(seq_along(results_list), function(i){ #results_list a list where each element corresponds to one sim_run
       # extract relevant data
     rr <- results_list[[i]]
     stopifnot(!is.null(rr$delta_hier_post), !is.null(rr$params$delta)) # delta_hier_post = the posterior draws for delta from that run;  and the true delta
       # summarize posterior
     post <- as.numeric(rr$delta_hier_post)
     dtrue <- rr$params$delta
     ci <- quantile(post, c(.025,.975), names = FALSE)
       # Create a summary row
     data.frame(
       sim_id = i,
       true = dtrue,
       post_mean = mean(post),
       post_median = median(post),
       post_sd = sd(post),
       ci_lo = ci[1],
       ci_hi = ci[2],
       covered_95 = as.numeric(dtrue >= ci[1] & dtrue <= ci[2])
     )
   })
   tab <- do.call(rbind, rows) # Takes the list of data frames (one per run) and stacks them into one big data frame. A row for each run
       # Computes the proportion of runs where the true delta was captured by the posterior 95% CI.
       # Example: If 92 out of 100 runs had coverage = 1, then avg_cov = 0.92.
    avg_cov <- mean(tab$covered_95)
   
   # Density
       # seq_along(results_list) → loops over all simulation runs.
        # For each run i, we extract:
        # delta_hier_post → posterior samples of delta for that simulation
        # params$delta → the true delta used in that run
        # data.frame() creates a small table for each run.
        # lapply(..., function(i){...}) returns a list of these tables, one per run.
        # do.call(rbind, ...) stacks all these tables vertically into one big data frame dens_df.
   dens_df <- do.call(rbind, lapply(seq_along(results_list), function(i){ 
     data.frame(sim_id = factor(i),
                post = as.numeric(results_list[[i]]$delta_hier_post),
                true = results_list[[i]]$params$delta)
   }))
   # Every Sim use the same real delta
   dtrue <- unique(dens_df$true)[1]
   
   # plotting
   p <-
     if (view == "overlay") {
       ggplot(dens_df, aes(x = post, colour = sim_id, fill = sim_id)) +
         geom_density(alpha = 0.25) +
         geom_vline(xintercept = dtrue, linetype = 2) +
         guides(fill = "none") +
         labs(title = title,
              subtitle = sprintf("Average 95%%-CI coverage = %.1f%%", 100*avg_cov),
              x = expression(delta), y = "Density") +
         theme_minimal(base_size = 12)
     } else if (view == "facets") {
       ggplot(dens_df, aes(x = post)) +
         geom_density() +
         geom_vline(xintercept = dtrue, linetype = 2) +
         facet_wrap(~ sim_id, scales = "free_y") +
         labs(title = title,
              subtitle = sprintf("Average 95%%-CI coverage = %.1f%%", 100*avg_cov),
              x = expression(delta), y = "Density") +
         theme_minimal(base_size = 12)
     } else { # ridgeline
       if (!requireNamespace("ggridges", quietly = TRUE)) {
         stop("Please install ggridges：install.packages('ggridges')")
       }
       ggplot(dens_df, aes(x = post, y = sim_id, height = ..density.., group = sim_id)) +
         ggridges::geom_ridgeline(stat = "density", scale = 1.2, fill = NA) +
         geom_vline(xintercept = dtrue, linetype = 2) +
         labs(title = title,
              subtitle = sprintf("Average 95%%-CI coverage = %.1f%%", 100*avg_cov),
              x = expression(delta), y = "Simulation run") +
         theme_minimal(base_size = 12)
     }
   
   list(per_run_table = tab,
        avg_coverage = avg_cov,
        plot = p)
 }

 
### Parameter Recovery Check ----
 
 # Example: check 1 saved run
 res1 <- delta_check_from_file("out/S-demo-01.RData", title_suffix = "run 1")
 print(res1$summary_table)
 print(res1$plot)
 
 
 # Parameter recovery check across all runs
 load("out/S-demo-all.RData")  
 
 # 1) Ovelay with diff colors
 res_overlay <- delta_check_from_results_list(results_list, view = "overlay")
 print(res_overlay$per_run_table)
 print(res_overlay$plot)
 
 # 2) Ridgeline（ggridges needed）
 res_ridge <- delta_check_from_results_list(results_list, view = "ridgeline")
 print(res_ridge$plot)
 
 # 3) Facets
 res_facets <- delta_check_from_results_list(results_list, view = "facets")
 print(res_facets$plot)
 
 