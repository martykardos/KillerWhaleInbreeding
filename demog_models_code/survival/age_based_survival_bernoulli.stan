data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for splines
  int Ks;  // number of linear effects
  matrix[N, Ks] Xs;  // design matrix for the linear effects
  // data for spline s(year, k = 7)
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
  // data for spline s(age, k = 3, by = sexF1M2)1
  int nb_2;  // number of bases
  int knots_2[nb_2];  // number of knots
  // basis function matrices
  matrix[N, knots_2[1]] Zs_2_1;
  // data for spline s(age, k = 3, by = sexF1M2)2
  //int nb_3;  // number of bases
  //int knots_3[nb_3];  // number of knots
  // basis function matrices
  //matrix[N, knots_3[1]] Zs_3_1;
  int prior_only;  // should the likelihood be ignored?
}
parameters {
  vector[K] b;  // population-level effects
  //real Intercept;  // temporary intercept for centered predictors
  vector[Ks] bs;  // spline coefficients
  // parameters for spline s(year, k = 7)
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  real<lower=0> sds_1_1;  // standard deviations of spline coefficients
  // parameters for spline s(age, k = 3, by = sexF1M2)1
  // standarized spline coefficients
  vector[knots_2[1]] zs_2_1;
  real<lower=0> sds_2_1;  // standard deviations of spline coefficients
  // parameters for spline s(age, k = 3, by = sexF1M2)2
  // standarized spline coefficients
  //vector[knots_3[1]] zs_3_1;
  //real<lower=0> sds_3_1;  // standard deviations of spline coefficients
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1;
  // actual spline coefficients
  vector[knots_2[1]] s_2_1;
  // actual spline coefficients
  //vector[knots_3[1]] s_3_1;
  // compute actual spline coefficients
  s_1_1 = sds_1_1 * zs_1_1;
  // compute actual spline coefficients
  s_2_1 = sds_2_1 * zs_2_1;
  // compute actual spline coefficients
  //s_3_1 = sds_3_1 * zs_3_1;
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1;// + Zs_3_1 * s_3_1;
    target += bernoulli_logit_glm_lpmf(Y | X, mu, b);
  }
  // priors including constants
  target += student_t_lpdf(sds_1_1 | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
  target += std_normal_lpdf(zs_1_1);
  target += std_normal_lpdf(b);
  target += std_normal_lpdf(bs);
  target += student_t_lpdf(sds_2_1 | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
  target += std_normal_lpdf(zs_2_1);
  //target += student_t_lpdf(sds_3_1 | 3, 0, 1)
  //  - 1 * student_t_lccdf(0 | 3, 0, 1);
  //target += std_normal_lpdf(zs_3_1);
}
generated quantities {
  vector[N] log_lik;
  vector[N] mu;
  vector[N] pred;
  mu = Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1;
  pred = X * b;
  for (n in 1:N) log_lik[n] = bernoulli_logit_lpmf(Y[n] | mu[n]+pred[n]);
}