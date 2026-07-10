// =====================================================================
//  cdtmbnma.stan
//  Component dose-response network meta-analysis with a dose-dependent
//  interaction surface.
//
//  Arm-based model. Each arm carries a dose for every component (0 = absent).
//  Continuous outcomes use a normal likelihood with known standard error.
//  Binary outcomes use a binomial-logit likelihood. Study-specific reference
//  levels remain free fixed effects; only component, dose, interaction, and
//  random-effect terms are pooled.
// =====================================================================
data {
  int<lower=1> N;                          // arms
  int<lower=1> S;                          // studies
  int<lower=1> C;                          // components
  int<lower=1> P;                          // interaction pairs; dummy pair when additive

  array[N] int<lower=1, upper=S> study;    // study index per arm
  matrix<lower=0>[N, C] D;                 // component doses per arm
  array[N] int<lower=0, upper=1> is_ref;   // reference arm flag

  array[P] int<lower=1, upper=C> pair_a;
  array[P] int<lower=1, upper=C> pair_b;
  vector<lower=0>[C] dstar;                // bilinear dose normalisers

  int<lower=1, upper=2> outcome;           // 1 = continuous, 2 = binary
  int<lower=0, upper=2> interaction;       // 0 = none, 1 = bilinear, 2 = gpdi

  vector[N] y;
  vector<lower=0>[N] se;
  array[N] int<lower=0> r;
  array[N] int<lower=0> nn;

  real<lower=0> prior_emax_sd;
  real prior_logED50_mean;
  real<lower=0> prior_logED50_sd;
  real<lower=0> prior_int_sd;
  real<lower=0> prior_ref_sd;
  real<lower=0> prior_omega_sd;

  int<lower=1> G;                          // prediction rows; all-zero if none requested
  matrix<lower=0>[G, C] Dpred;
}

parameters {
  vector[S] m;
  vector[C] emax;
  vector[C] logED50;
  vector[P] eta;
  vector[P] INT;
  vector<lower=0>[C] kappa;
  real<lower=0> omega;
  vector[N] z;
}

transformed parameters {
  vector<lower=0>[C] ED50 = exp(logED50);
  vector[N] mu;
  for (i in 1:N) {
    real lp = m[study[i]];
    for (c in 1:C) {
      lp += emax[c] * D[i, c] / (ED50[c] + D[i, c]);
    }
    if (interaction == 1) {
      for (p in 1:P) {
        lp += eta[p]
              * (D[i, pair_a[p]] / dstar[pair_a[p]])
              * (D[i, pair_b[p]] / dstar[pair_b[p]]);
      }
    } else if (interaction == 2) {
      for (p in 1:P) {
        int a = pair_a[p];
        int b = pair_b[p];
        lp += INT[p]
              * (D[i, a] / (kappa[a] + D[i, a]))
              * (D[i, b] / (kappa[b] + D[i, b]));
      }
    }
    if (is_ref[i] == 0) {
      lp += z[i] * omega;
    }
    mu[i] = lp;
  }
}

model {
  m       ~ normal(0, prior_ref_sd);
  emax    ~ normal(0, prior_emax_sd);
  logED50 ~ normal(prior_logED50_mean, prior_logED50_sd);
  eta   ~ normal(0, prior_int_sd);
  INT   ~ normal(0, prior_int_sd);
  kappa ~ lognormal(prior_logED50_mean, prior_logED50_sd);
  omega ~ normal(0, prior_omega_sd);
  z     ~ std_normal();

  if (outcome == 1) {
    y ~ normal(mu, se);
  } else {
    r ~ binomial_logit(nn, mu);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[G] pred;

  for (i in 1:N) {
    if (outcome == 1) {
      log_lik[i] = normal_lpdf(y[i] | mu[i], se[i]);
    } else {
      log_lik[i] = binomial_logit_lpmf(r[i] | nn[i], mu[i]);
    }
  }

  for (g in 1:G) {
    real lp = 0;
    for (c in 1:C) {
      lp += emax[c] * Dpred[g, c] / (ED50[c] + Dpred[g, c]);
    }
    if (interaction == 1) {
      for (p in 1:P) {
        lp += eta[p]
              * (Dpred[g, pair_a[p]] / dstar[pair_a[p]])
              * (Dpred[g, pair_b[p]] / dstar[pair_b[p]]);
      }
    } else if (interaction == 2) {
      for (p in 1:P) {
        int a = pair_a[p];
        int b = pair_b[p];
        lp += INT[p]
              * (Dpred[g, a] / (kappa[a] + Dpred[g, a]))
              * (Dpred[g, b] / (kappa[b] + Dpred[g, b]));
      }
    }
    pred[g] = lp;
  }
}
