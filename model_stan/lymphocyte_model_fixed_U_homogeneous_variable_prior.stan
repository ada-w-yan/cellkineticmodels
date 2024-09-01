functions {
  real calc_U(real t,
  real frac,
  real delta,
  real label_end);
  
  real calc_U(real t,
  real frac,
  real delta,
  real label_end) {
    real U;
    real t0;
    real f1;
    if(t == 0) {
      U = 0;
      return U;
    } else if(t <= label_end) {
      f1 = frac;
      t0 = 0;
    } else {
      f1 = 0;
      t0 = label_end;
    }
    U = f1 * (1 - exp(-delta * (t - t0))) + calc_U(t0, frac, delta, label_end) * exp(-delta * (t - t0));
    return U;
  }
  
  real calc_L(real t,
  real frac,
  real delta,
  real label_end,
  real pb_w,
  real dstar,
  int perfect_labelling);
  real calc_L(real t,
  real frac,
  real delta,
  real label_end,
  real pb_w,
  real dstar,
  int perfect_labelling) {
    real L;
    real A;
    real B;
    real t_shift;
      if(perfect_labelling == 1) {
      if(t <= label_end) {
        L = pb_w / dstar - exp(-dstar * t);
      } else {
        t_shift = t - label_end;
        L = calc_L(label_end, frac, delta, label_end, pb_w, dstar, perfect_labelling) * exp(-dstar * t_shift);
      }
      return L;
    }
    if(t <= label_end) {
      L = pb_w * frac / (delta - dstar) * (delta/dstar *(1 - exp(-dstar * t)) - (1 - exp(-delta * t)));
    } else {
      t_shift = t - label_end;
      L = calc_L(label_end, frac, delta, label_end, pb_w, dstar, perfect_labelling) * exp(-dstar * t_shift) + pb_w * exp(-dstar * t_shift) * calc_U(label_end, frac, delta, label_end) / (dstar - delta) * (exp((dstar - delta) * t_shift) - 1);
    }
    return L;
  }

}

// data and fixed parameter values
data {
  int<lower=1> T; // number of sampling times
  real ts[T]; // sampling times
  real L[T]; // fraction of label in granulocytes
  real<lower=0> label_end; // time at which labelling ends
  real<lower=0.,upper=0.1> frac; // 
  real<lower=0.,upper=1.> delta; //
  real<lower=0.,upper=1> sigma; // standard deviation of observation model
  real<lower = 0., upper = 7.> b_w;
  int evaluate_likelihood; // if 1, sample from posterior, if 0, sample from prior
  int<lower = 0, upper = 1> perfect_labelling; // set to 0 if accounting for
  // turnover of label in body water, 1 otherwise
  real<lower = 0.> prior_width;
}

// model parameters
parameters {
  // ranges are ranges of uniform priors
  real<lower = 0., upper = 1.> p_raw;
}

transformed parameters {
  // ranges are ranges of uniform priors]
  real p = p_raw * prior_width;
  real Lhat[T];
  if(evaluate_likelihood) {
    for (t in 1:T) {
      Lhat[t] = calc_L(ts[t], frac, delta, label_end, p * b_w, p, perfect_labelling);
    }
  } else {
    for (t in 1:T) {
      Lhat[t] = 0;
    }
  }
}

// observation model: links viral load predicted by sampled parameter values
// to the observed viral load 
model {
  if(evaluate_likelihood) {
    for (t in 1:T) {
      L[t] ~ normal(Lhat[t], sigma);
    }
  }
}
