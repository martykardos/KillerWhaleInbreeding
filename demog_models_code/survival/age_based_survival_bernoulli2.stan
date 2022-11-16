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
  
  vector[knots_1[1]] Zs1_mean;
  matrix[knots_1[1],knots_1[1]] Zs1_cov;
  vector[Ks] bs_mean;
  matrix[Ks,Ks] bs_cov;
  vector[knots_2[1]] Zs2_mean;
  matrix[knots_2[1],knots_2[1]] Zs2_cov;
  vector[K] b_mean;
  matrix[K,K] b_cov;
  //real sd1_mean;
  //real sd1_cov;
  //real sd2_mean;
  //real sd2_cov;  
  int prior_only;  // should the likelihood be ignored?
}
parameters {
  vector[K] b;  // population-level effects
  //real Intercept;  // temporary intercept for centered predictors
  vector[Ks] bs;  // spline coefficients
  // parameters for spline s(year, k = 7)
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  real sds_1_1;  // standard deviations of spline coefficients
  // parameters for spline s(age, k = 3, by = sexF1M2)1
  // standarized spline coefficients
  vector[knots_2[1]] zs_2_1;
  real sds_2_1;  // standard deviations of spline coefficients
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1;
  // actual spline coefficients
  vector[knots_2[1]] s_2_1;
  // compute actual spline coefficients
  s_1_1 = sds_1_1 * zs_1_1;
  // compute actual spline coefficients
  s_2_1 = sds_2_1 * zs_2_1;
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Xs * bs + Zs_1_1 * s_1_1;// + Zs_2_1 * s_2_1;
    target += bernoulli_logit_glm_lpmf(Y | X, mu, b);
  }
  // priors including constants
  target += student_t_lpdf(sds_1_1 | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);

  target += student_t_lpdf(sds_2_1 | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
  //target += normal_lpdf(sds_1_1 | sd1_mean, sd1_cov);
  //target += normal_lpdf(sds_2_1 | sd2_mean, sd2_cov);
  //target += normal_lpdf(zs_2_1 | Zs2_mean, Zs2_cov);
  target += multi_normal_lpdf(b | b_mean, b_cov);
  target += multi_normal_lpdf(zs_1_1 | Zs1_mean, Zs1_cov);
  target += multi_normal_lpdf(zs_2_1 | Zs2_mean, Zs2_cov);
  target += multi_normal_lpdf(bs | bs_mean, bs_cov);

}
