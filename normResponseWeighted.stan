data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Ncol; // total number of columns in model matrix
  matrix[Ntotal, Ncol] X; // model matrix
  real y[Ntotal]; // response variable normal distributed
  real wts[Ntotal]; // weights for each observation
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  vector[Ncol] betas; // regression parameters
  //real populationMean;// constant term
  real<lower=0.01> sigmaPop; // data standard deviation
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  // fitted value
  mu = X * betas;
  //mu = mu + populationMean;
}
model {
  // using diffuse prior
  sigmaPop ~ cauchy(0, 2);
  //populationMean ~ cauchy(0, 1);
  betas ~ cauchy(0, 1); //prior for the betas
  // likelihood function with weights
  for(i in 1:Ntotal){
  target +=  normal_lpdf(y[i] | mu[i], sigmaPop) * wts[i];
  }
}
generated quantities {
  vector[Ntotal] log_lik;
  for (i in 1:Ntotal) log_lik[i] = normal_lpdf(y[i] | mu[i], sigmaPop) * wts[i];
}
