functions {
  
  real calc_X(real t,
  real label_end,
  real pb_w,
  real dstar);
  real calc_X(real t,
  real label_end,
  real p,
  real dstar) {
    real X;
    real t_shift;
    
    if(t <= label_end) {
      X = p / dstar * (1 - exp(-dstar * t));
    } else {
      t_shift = t - label_end;
      X = calc_X(label_end, label_end, p, dstar) * exp(-dstar * t_shift);
    }
    return X;
  }
}

// data and fixed parameter values
data {
  int<lower=1> T; // number of sampling times
  real ts[T]; // sampling times
  real L[T]; // fraction of label in granulocytes
  real<lower=0.> label_end; // time at which labelling ends
  real<lower=0.,upper=1> sigma; // standard deviation of observation model
  int evaluate_likelihood; // if 1, sample from posterior, if 0, sample from prior
}

// model parameters
parameters {
  // ranges are ranges of uniform priors
  real<lower = 0., upper = 10.> p;
  real<lower = 0., upper = 10.> dstar;
}

transformed parameters {
  // ranges are ranges of uniform priors
  real Lhat[T];
  
  if(evaluate_likelihood) {
    for (t in 1:T) {
      Lhat[t] = calc_X(ts[t], label_end, p, dstar);
    }
  } else {
    for (t in 1:T) {
      Lhat[t] = 0;
    }
  }
}

// observation model
model {
  if(evaluate_likelihood) {
    for (t in 1:T) {
      L[t] ~ normal(Lhat[t], sigma);
    }
  }
}

generated quantities {
  vector[T] log_lik;
  for (t in 1:T) {
    log_lik[t] = normal_lpdf(L[t] | Lhat[t], sigma);
  }
}
