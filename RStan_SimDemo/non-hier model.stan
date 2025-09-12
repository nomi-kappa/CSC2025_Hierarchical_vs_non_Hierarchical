# Non-hierarchical = uses subject averages; simpler, no trial-level info.

# Keypoints
  # Data: only subject-level means m_rt[n,2]. No trial-level data.
  # Parameteres: 
     # mu_m1 → group mean for condition 1
     # delta → effect size
    # sigma_m → between-subject SD
  # No subject-level parameters: The model does not explicitly model each subject’s mean and SD.
 # Likelihood: m_rt[isubj, igr] ~ normal(mu_m[igr], sigma_m)
      # Each subject’s mean is drawn from the group-level distribution.
      # This is essentially a simple normal model on subject means, no hierarchical structure at the trial level.

#Group-level
  # mu_m1  delta  sigma_m
  #   |       |
  #   v       v
#Subject-level averages. 
  # m_rt[i, cond]  (i = 1..N)


data{
  int n;
  real m_rt[n, 2];
}
parameters {
  real<lower=0> mu_m1;
  real	delta;
  real<lower=0> sigma_m;
}
transformed parameters {
  real mu_m[2];
  mu_m[1] <- mu_m1;
  mu_m[2] <- mu_m1 + delta*sigma_m;
}
model{
	#put priors on mean and sd of the group-level distributions
	mu_m1 ~ normal(6,0.333)T[0,];
	sigma_m ~ uniform(0,15);

	delta ~ normal(0,1);

	#model the pps
	for (igr in 1:2){
		for (isubj in 1:n){
			m_rt[isubj, igr] ~ normal(mu_m[igr],sigma_m);
		}
	}
}
