data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> Nclusters; // number of groups for random intercepts
    int<lower=1> NgroupMap[Ntotal]; // mapping variable to map each observation to a group 
    int<lower=1> Ncol; // total number of columns in model matrix
    matrix[Ntotal, Ncol] X; // model matrix
    real y[Ntotal]; // response variable normally distributed
}
// transformed data {
// }
parameters {
  // parameters to estimate in the model
    vector[Ncol] betas; // regression parameters
    real<lower=0> sigmaRan; // random effect standard deviation
    real<lower=0> sigmaPop; // population standard deviation
    vector[Nclusters] rGroupsJitter; // number of random jitters for each cluster member
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  vector[Ntotal] rNewIntercept;
  matrix[Ntotal, (Ncol-1)] mX2;
  // add this random jitter to the population intercept
  vector[Nclusters] rBetaRand;
  rBetaRand = rGroupsJitter + betas[1];
  // new intercept variable equal to number of observations
  // use the mapping variable to expand this vector
  rNewIntercept = rBetaRand[NgroupMap];
  // new model matrix without intercept
  mX2 = X[,2:Ncol];
  // fitted value
  mu = mX2 * betas[2:Ncol]; 
  // add the random intercept
  mu = mu + rNewIntercept;
}
model {
  // using non-informative priors to start with
  //sigmaRan ~ uniform(0, 1e3);
  betas ~ cauchy(0, 10);//prior for the betas
  // random effects sample
  rGroupsJitter ~ normal(0, sigmaRan);
  // likelihood function
  y ~ normal(mu, sigmaPop);
}