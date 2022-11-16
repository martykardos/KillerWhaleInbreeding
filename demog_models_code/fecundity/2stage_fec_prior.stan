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
  vector[K] prior_mean;
  matrix[K,K] prior_cov;
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
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1;
  // compute actual spline coefficients
  s_1_1 = sds_1_1 * zs_1_1;
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N) + Xs * bs + Zs_1_1 * s_1_1;
    target += bernoulli_logit_glm_lpmf(Y | X, mu, b);
  }
  // priors including constants
  //target += normal_lpdf(Intercept | -7.811627, 2.885389);
  target += multi_normal_lpdf(b | prior_mean, prior_cov);
  target += std_normal_lpdf(bs);
  target += student_t_lpdf(sds_1_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_1_1);
}
