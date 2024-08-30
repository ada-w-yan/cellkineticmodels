functions {
  // calculate fraction of label in body water
  real calc_U(real t,
  real f,
  real delta,
  real label_end);
  
  real calc_U(real t,
  real f,
  real delta,
  real label_end) {
    real U;
    real t0;
    real f1;
    if(t == 0) {
      U = 0;
      return U;
    } else if (t <= label_end) {
      f1 = f;
      t0 = 0;
    } else {
      f1 = 0;
      t0 = label_end;
    }
    U = f1 * (1 - exp(-delta * (t - t0))) + calc_U(t0, f, delta, label_end) * exp(-delta * (t - t0));
    return U;
  }
  
  // calculate fraction of label in a single compartment
  real calc_L_single(real t,
  real frac,
  real delta,
  real label_end,
  real pb_w,
  real dstar,
  int perfect_labelling);
  real calc_L_single(real t,
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
    if(t == 0) {
      return 0;
    }
    if(perfect_labelling == 1) {
      if(t <= label_end) {
        L = pb_w / dstar - exp(-dstar * t);
      } else {
        t_shift = t - label_end;
        L = calc_L_single(label_end, frac, delta, label_end, pb_w, dstar, perfect_labelling) * exp(-dstar * t_shift);
      }
      return L;
    }
    if(t <= label_end) {
      L = pb_w * frac / (delta - dstar) * (delta/dstar *(1 - exp(-dstar * t)) - (1 - exp(-delta * t)));
    } else {
      t_shift = t - label_end;
      L = calc_L_single(label_end, frac, delta, label_end, pb_w, dstar, perfect_labelling) * exp(-dstar * t_shift) + pb_w * exp(-dstar * t_shift) * calc_U(label_end, frac, delta, label_end) / (dstar - delta) * (exp((dstar - delta) * t_shift) - 1);
    }
    return L;
  }
  
  // calculate total fraction of label across compartments
  real calc_L(real t, 
  real frac, 
  real delta, 
  real label_end, 
  vector p, 
  real b_w,
  vector frac_pop,
  int n_pop,
  int perfect_labelling) {
    real L;
    L = 0;
    for(i in 1:n_pop) {
      L = L + calc_L_single(t, frac, delta, label_end, p[i] * b_w, p[i], perfect_labelling) * frac_pop[i];
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
  real<lower=0.,upper=0.1> frac; // set to 0 for perfect labelling
  real<lower=0.,upper=1.> delta; // set to 0 for perfect labelling
  real<lower=0.,upper=1> sigma; // standard deviation of observation model
  real<lower = 0., upper = 7.> b_w; // set to 1 for perfect labelling
  int<lower = 1> n_pop; // number of compartments
  int<lower = 0, upper = 1> perfect_labelling; // set to 0 if accounting for
  // turnover of label in body water, 1 otherwise
  real<lower = 0.> prior_width;
  int evaluate_likelihood; // if 1, sample from posterior, if 0, sample from prior
}

// model parameters
parameters {
  // ranges are ranges of uniform priors
  vector<lower = 0., upper = 1.>[n_pop] p_raw; // proliferation rate for each compartment
  positive_ordered[n_pop] frac_pop_raw; // a vector ordered from lowest to highest.
  // dividing each by their sum yields the fraction of cells in each compartment
}

transformed parameters {
  real Lhat[T];
  real bar_p;
  vector[n_pop] p = p_raw .* prior_width;
  simplex[n_pop] frac_pop = frac_pop_raw/sum(frac_pop_raw);
  // dividing frac_pop_raw by its sum yields the fraction of cells in each compartment
  bar_p = sum(p .* frac_pop); // calculate mean proliferation rate as an output
  
  // calculate total fraction of label
  if(evaluate_likelihood) {
    for (t in 1:T) {
      Lhat[t] = calc_L(ts[t], frac, delta, label_end, p, b_w, frac_pop, n_pop, perfect_labelling);
    }
  } else {
    for (t in 1:T) {
      Lhat[t] = 0;
    }
  }
}

// observation model
model {
  
  // this imposes a dirichlet prior with alpha = 1 on frac_pop
  // https://mc-stan.org/docs/2_25/functions-reference/dirichlet-distribution.html
  // which means that the prior is uniform across all simplexes
  // where a simplex is a vector that adds to 1
  frac_pop_raw ~ gamma(1., 1);
  
  // normally distributed noise
  if(evaluate_likelihood) {
    for (t in 1:T) {
      L[t] ~ normal(Lhat[t], sigma);
    }
  }
}

// output log likelihood
generated quantities {
  vector[T] log_lik;
  for (t in 1:T) {
    log_lik[t] = normal_lpdf(L[t] | Lhat[t], sigma);
  }
}
