data {
  int<lower=1> N;
  vector[N] theta_o;
  vector[N] theta_r;
  vector<lower=0>[N] var_o;
  vector<lower=0>[N] var_r;
  real<lower=0> a;
  real<lower=0> b;
}

parameters {
  vector<lower=0, upper=1>[N] gamma;
}

transformed parameters {
  vector[N] mu_pred;
  vector<lower=0>[N] var_pred;
  mu_pred = gamma .* theta_o;
  var_pred = var_r + gamma .* var_o;
}

model {
  gamma ~ beta(a, b);
  theta_r ~ normal(mu_pred, sqrt(var_pred));
}

generated quantities {
  vector[N] theta_r_rep;
  vector[N] log_lik;
  for (i in 1:N) {
    theta_r_rep[i] =
      normal_rng(mu_pred[i], sqrt(var_pred[i]));

    log_lik[i] =
      normal_lpdf(theta_r[i] |
                  mu_pred[i],
                  sqrt(var_pred[i]));
  }
}
