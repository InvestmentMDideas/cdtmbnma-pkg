// =====================================================================
//  cdt_stage1.stan
//  Two-component component x dose x time network meta-analysis.
//
//  Exponential time-course; Emax component dose-response on the long-term
//  asymptote; bilinear dose-dependent interaction on the asymptote; additive
//  log-linear dose effects on onset rate; independent random effects on active
//  arms. Continuous aggregate likelihood with known standard errors.
// =====================================================================
data {
  int<lower=1> n_obs;
  int<lower=1> n_arms;
  int<lower=1> n_studies;
  int<lower=1> n_active;

  array[n_obs] int<lower=1, upper=n_arms> row_arm;
  vector<lower=0>[n_obs] time;
  vector[n_obs] ybar;
  vector<lower=0>[n_obs] se;

  array[n_arms] int<lower=1, upper=n_studies> arm_sid;
  vector<lower=0>[n_arms] arm_dA;
  vector<lower=0>[n_arms] arm_dB;
  array[n_active] int<lower=1, upper=n_arms> active_idx;

  real<lower=0> dAstar;
  real<lower=0> dBstar;

  real<lower=0> prior_ref_E_sd;
  real prior_ref_lograte_mean;
  real<lower=0> prior_ref_lograte_sd;
  real<lower=0> prior_emax_sd;
  real prior_logED50_mean;
  real<lower=0> prior_logED50_sd;
  real<lower=0> prior_eta_sd;
  real<lower=0> prior_rate_sd;
  real<lower=0> prior_omega_E_sd;
  real<lower=0> prior_omega_k_sd;

  int<lower=1> n_pred;
  vector<lower=0>[n_pred] pred_dA;
  vector<lower=0>[n_pred] pred_dB;
}

parameters {
  vector[n_studies] mE;
  vector[n_studies] mk;

  real a1_EmaxA;
  real a2_EmaxB;
  real logED50A;
  real logED50B;
  real eta;

  real b1_rateA;
  real b2_rateB;

  real<lower=0> omega_E;
  real<lower=0> omega_k;
  vector[n_active] zE;
  vector[n_active] zk;
}

transformed parameters {
  real<lower=0> ED50A = exp(logED50A);
  real<lower=0> ED50B = exp(logED50B);

  vector[n_arms] uE = rep_vector(0.0, n_arms);
  vector[n_arms] uk = rep_vector(0.0, n_arms);
  vector[n_arms] E_arm;
  vector[n_arms] k_arm;

  for (a in 1:n_active) {
    uE[active_idx[a]] = zE[a] * omega_E;
    uk[active_idx[a]] = zk[a] * omega_k;
  }

  for (j in 1:n_arms) {
    real gA = a1_EmaxA * arm_dA[j] / (ED50A + arm_dA[j]);
    real gB = a2_EmaxB * arm_dB[j] / (ED50B + arm_dB[j]);
    real inter = eta * (arm_dA[j] / dAstar) * (arm_dB[j] / dBstar);
    E_arm[j] = mE[arm_sid[j]] + gA + gB + inter + uE[j];
    k_arm[j] = mk[arm_sid[j]]
               + b1_rateA * log(arm_dA[j] + 1.0)
               + b2_rateB * log(arm_dB[j] + 1.0)
               + uk[j];
  }
}

model {
  mE ~ normal(0, prior_ref_E_sd);
  mk ~ normal(prior_ref_lograte_mean, prior_ref_lograte_sd);
  a1_EmaxA ~ normal(0, prior_emax_sd);
  a2_EmaxB ~ normal(0, prior_emax_sd);
  logED50A ~ normal(prior_logED50_mean, prior_logED50_sd);
  logED50B ~ normal(prior_logED50_mean, prior_logED50_sd);
  eta ~ normal(0, prior_eta_sd);
  b1_rateA ~ normal(0, prior_rate_sd);
  b2_rateB ~ normal(0, prior_rate_sd);
  omega_E ~ normal(0, prior_omega_E_sd);
  omega_k ~ normal(0, prior_omega_k_sd);
  zE ~ std_normal();
  zk ~ std_normal();

  for (i in 1:n_obs) {
    int j = row_arm[i];
    real mu = E_arm[j] * (1.0 - exp(-exp(k_arm[j]) * time[i]));
    ybar[i] ~ normal(mu, se[i]);
  }
}

generated quantities {
  vector[n_obs] log_lik;
  vector[n_pred] pred_effect;

  for (i in 1:n_obs) {
    int j = row_arm[i];
    real mu = E_arm[j] * (1.0 - exp(-exp(k_arm[j]) * time[i]));
    log_lik[i] = normal_lpdf(ybar[i] | mu, se[i]);
  }

  for (g in 1:n_pred) {
    real gA = a1_EmaxA * pred_dA[g] / (ED50A + pred_dA[g]);
    real gB = a2_EmaxB * pred_dB[g] / (ED50B + pred_dB[g]);
    real it = eta * (pred_dA[g] / dAstar) * (pred_dB[g] / dBstar);
    pred_effect[g] = gA + gB + it;
  }
}
