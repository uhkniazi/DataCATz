data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> Nclusters; // number of groups for random intercepts
    int<lower=1, upper=Nclusters> NgroupMap[Ntotal]; // mapping variable to map each observation to a group 
    //int<lower=1> Ncol; // total number of columns in model matrix
    vector[Ntotal] X; // model matrix without intercept
    real y[Ntotal]; // response variable normally distributed
}
// transformed data {
// }
parameters {
  // parameters to estimate in the model
    real intercept; // intercept + covariates in model matrix - population parameters
    real slope; 
    real<lower=0> sigmaRan1; // random effect standard deviation
    real<lower=0> sigmaPop; // population standard deviation
    //real<lower=0> sigmaRanSlope1; // random slope standard deviation
    vector[Nclusters] rGroupsJitter; // number of random jitters for each cluster member
    //vector[Nclusters] rGroupsSlope; // number of random slopes for each cluster member
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  // fitted value
  mu = X * slope + intercept + rGroupsJitter[NgroupMap];
  //mu = mu + (betas[1] + rGroupsJitter[NgroupMap]); //(betas[2] + rGroupsSlope[NgroupMap]));
}
model {
  // using non-informative priors to start with
  //sigmaRan ~ uniform(0, 1e3);
  intercept ~ cauchy(0, 2.5);//prior for the betas
  // random effects sample
  rGroupsJitter ~ normal(0, sigmaRan1);
  //rGroupsSlope ~ normal(0, sigmaRanSlope1);
  // likelihood function
  y ~ normal(mu, sigmaPop);
}
