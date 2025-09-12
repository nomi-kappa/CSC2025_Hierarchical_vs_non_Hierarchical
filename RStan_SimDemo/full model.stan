data{
  int n;
  int k;
  real rt[2, k, n];

# Full hierarchical = trial-level model with partial pooling → captures both within-subject variability and between-subject variability.
  # Key points
    # Data: trial-level RTs rt[2,k,n]. Each subject × trial × condition.
    # Parameters: 
        # mu_m1, delta, sigma_m → group-level parameters
        # scl → scaling factor for the half-Cauchy prior on subject SDs
        # m[n,2] → subject-level means
        # s[n,2] → subject-level SDs
    # Hierarchical structure: 
      # Group-level: mu_m1, delta, sigma_m, scl
      #  ↓
      # Subject-level: m[i, cond], s[i, cond]
      #  ↓
      # Trial-level: rt[cond, trial, subject] ~ Normal(m[i,cond], s[i,cond])
    # Priors:
      # Subject-level means: m[i,cond] ~ Normal(mu_m[cond], sigma_m) → partial pooling
      # Subject-level SDs: s[i,cond] ~ Half-Cauchy(0, scl) → constrains SDs to be positive
    # Likelihood:
      # Each trial RT is modeled conditionally on the subject’s mean and SD, not just the group.


parameters {
  real<lower=0> mu_m1;
  real	delta;
  real<lower=0> sigma_m;
  real<lower=0> scl;
  real<lower=0> s[n, 2];
  real<lower=0> m[n, 2];
}
transformed parameters {
  real mu_m[2];
  mu_m[1] <- mu_m1;
  mu_m[2] <- mu_m1 + delta*sigma_m;
}
model{	
	#put priors on mean and sd of the group-level normals
	mu_m[1] ~ normal(6,0.333)T[0,];
	sigma_m ~ uniform(0,15)T[0,];
	delta ~ normal(0,1);
	
	scl ~ uniform(0,100)T[0,];
	
	#model the pps
	for(igr in 1:2){
		for (isubj in 1:n){
			s[isubj, igr] ~ cauchy(0, scl)T[0,];# 1/A^2=lambda.s;
			m[isubj, igr] ~ normal(mu_m[igr],sigma_m)T[0,];
			
			for (itrl in 1:k){
				rt[igr, itrl, isubj] ~ normal(m[isubj,igr], s[isubj, igr]);
			}
		}
	}
}
