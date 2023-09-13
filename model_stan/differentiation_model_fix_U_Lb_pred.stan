functions {
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
  
  real calc_C(real t,
  real f,
  real delta,
  real b_w,
  real p_C,
  real dstar_C,
  real Delta,
  real label_end);
  real calc_C(real t,
  real f,
  real delta,
  real b_w,
  real p_C,
  real dstar_C,
  real Delta,
  real label_end) {
    real F_C;
    real t_shift;
    real dstar_C_Delta;
    real B1;
    real C;
    real D;
    real C_end;
    real U_end;
    real t_end;
    real f1;
    dstar_C_Delta = dstar_C + Delta;
    
    if(t == 0) {
      return 0;
    } else if(t <= label_end) {
      t_end = 0;
      f1 = f;
    } else {
      t_end = label_end;
      f1 = 0;
    }
    C_end = calc_C(t_end, f, delta, b_w, p_C, dstar_C, Delta, label_end);
    U_end = calc_U(t_end, f, delta, label_end);
    t_shift = t - t_end;
    B1 = p_C * b_w * f1 / dstar_C_Delta;
    D = p_C * b_w * (U_end - f1) / (dstar_C_Delta - delta);
    C =  C_end - D - B1;
    
    F_C = B1 + C * exp(-dstar_C_Delta * t_shift) + D * exp(-delta * t_shift);
    return F_C;
  }
  
  real calc_E(real t,
  real f,
  real delta,
  real b_w,
  real p_C,
  real p_E,
  real dstar_C,
  real dstar_E,
  real Delta,
  real two_k_minus_one, 
  real C_on_E, 
  real label_end);
  real calc_E(real t,
  real f,
  real delta,
  real b_w,
  real p_C,
  real p_E,
  real dstar_C,
  real dstar_E,
  real Delta,
  real two_k_minus_one, 
  real C_on_E, 
  real label_end) {
    real F_E;
    real t_shift;
    real dstar_C_Delta;
    real A1;
    real B1;
    real C;
    real D;
    real E;
    real H;
    real G;
    real C_end;
    real U_end;
    real E_init;
    real t_end;
    real f1;
    dstar_C_Delta = dstar_C + Delta;
    A1 = (two_k_minus_one * Delta * C_on_E + p_E) * b_w;
    
    if(t == 0) {
      return 0;
    } else if(t <= label_end) {
      t_end = 0;
      f1 = f;
    } else {
      t_end = label_end;
      f1 = 0;
    }
    C_end = calc_C(t_end, f, delta, b_w, p_C, dstar_C, Delta, label_end);
    U_end = calc_U(t_end, f, delta, label_end);
    E_init = calc_E(t_end, f, delta, b_w, p_C, p_E, dstar_C, dstar_E, Delta, two_k_minus_one, C_on_E, label_end);
    t_shift = t - t_end;
    B1 = p_C * b_w * f1 / dstar_C_Delta;
    D = p_C * b_w * (U_end - f1) / (dstar_C_Delta - delta);
    C =  C_end - D - B1;
    E = A1 * f1 + Delta * C_on_E * B1;
    H = Delta * C_on_E * C;
    G = Delta * C_on_E * D + A1 * (U_end - f1);
    
    F_E = E_init * exp(-dstar_E * t_shift) +
    E / dstar_E * (1 - exp(-dstar_E * t_shift)) + 
    H /(dstar_E - dstar_C_Delta) * (exp(-dstar_C_Delta * t_shift) - exp(-dstar_E * t_shift)) + 
    G / (dstar_E - delta) * (exp(-delta * t_shift) - exp(-dstar_E * t_shift));
    return F_E;
  }
}

// data and fixed parameter values
data {
  int<lower=1> N_samples; // number of posterior samples
  int<lower=1> T; // number of sampling times
  real ts[T]; // sampling times
  real<lower=0> label_end; // time at which labelling ends
  real<lower=0.,upper=0.1> f; // 
  real<lower=0.,upper=1.> delta; //
  real<lower = 0.> b_w;
  real<lower = 0.> p_C[N_samples];
  real<lower = 0.> p_E[N_samples];
  real<lower = 0.> dstar_C[N_samples];
  real<lower = 0.> dstar_E[N_samples];
  real<lower = 0.> Delta[N_samples];
  real<lower = 0.> two_k_minus_one[N_samples];
  real<lower = 0.> C_on_E[N_samples];
}

// model parameters
parameters {
}

model {
}

// observation model: links fraction of label predicted by sampled parameter values
// to the observed fraction of label
generated quantities {
  real Chat[T, N_samples];
  real Ehat[T, N_samples];
  real Uhat[T];
  for (t in 1:T) {
    for (n in 1:N_samples) {
      Chat[t,n] = calc_C(ts[t], f, delta, b_w, p_C[n], dstar_C[n], Delta[n], label_end);
      Ehat[t,n] = calc_E(ts[t], f, delta, b_w, p_C[n], p_E[n], dstar_C[n], dstar_E[n], Delta[n], two_k_minus_one[n], C_on_E[n], label_end);
    }
    Uhat[t] = calc_U(ts[t], f, delta, label_end);
  }
}
