
// bwm_dirichlet_lognormal.stan — BBWM with Dirichlet center and lognormal ratios
// Data are built from best-to-others and others-to-worst ratios

data {
  int<lower=1> m;               // number of items
  int<lower=1> J;               // number of experts
  matrix[J, m] bos;             // Best-to-Other ratios
  matrix[J, m] ows;             // Other-to-Worst ratios
  int<lower=1,upper=m> best[J]; // index of best item per expert
  int<lower=1,upper=m> worst[J];// index of worst item per expert
  vector<lower=0>[m] m_fuzzy;   // normalized fuzzy-Delphi centers (sum to 1)
  real<lower=0> c;              // Dirichlet concentration parameter
}
parameters {
  simplex[m] alpha;             // item weights
  real<lower=0> sigma_b;        // noise for best-to-other
  real<lower=0> sigma_w;        // noise for other-to-worst
}
model {
  // Prior
  alpha ~ dirichlet(c * m_fuzzy);
  sigma_b ~ normal(0, 1);
  sigma_w ~ normal(0, 1);
  // Likelihood
  for (j in 1:J) {
    for (i in 1:m) {
      bos[j, i] ~ lognormal(log(alpha[best[j]] / alpha[i]), sigma_b);
      ows[j, i] ~ lognormal(log(alpha[i] / alpha[worst[j]]),        sigma_w);
    }
  }
}

