data {
  int<lower=0> p;             // number of variables
  int<lower=0> n;             // number of observations 
  vector[n] y;                // response vector
  matrix[n,p] X;              // design matrix
  positive_ordered[p] lambda; // penalty vector in ascending order
  real<lower=0> a;            // shape parameter for the Gamma prior on sigma
  real<lower=0> g;            // scale parameter for the Gamma prior on sigma
}

parameters {
  vector[p] beta;             // regression coefficients
  real<lower=0> sigma;        // variance 
}

transformed parameters {
    vector[n] eta;
    vector[p] AbsBeta;
    real<lower=0> sigma2;
    sigma2 = square(sigma);
    eta = X * beta;
    for (i in 1:p)
    	AbsBeta[i] = fabs(beta[i]);
}

model {
    sigma2 ~ inv_gamma(a,g); 
    //target += inv_gamma_lpdf(sigma2, a, g)
    target += log(sigma);
    target += -dot_product(lambda,sort_asc(AbsBeta))/sigma - p*log(sigma);   // SLOPE_prior(beta)
    for (j in 1:n)
        target += normal_lpdf(y[j]| eta[j], sigma);
}
