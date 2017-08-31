data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Nvars; // number of variables/vectors in data matrix
  int<lower=1> Neigens; // number of eigen vectors to use
  matrix[Ntotal, Nvars] y; // data matrix where each vector is in column format
}

transformed data {
  // create an identity covariance matrix for rotated components
  // fixed hyperparameter for the covariance matrix 
  matrix[Neigens, Neigens] mIdentity;
  // fixed hyperparameter for mean
  vector[Neigens] mu0;
  mIdentity = diag_matrix(rep_vector(1.0, Neigens));
  mu0 = rep_vector(0.0, Neigens);
}

parameters { // the parameters to track
  // eigen vectors of the covariance matrix of the data matrix
  matrix[Nvars, Neigens] mEigens; // note it is the transformation matrix 
  real<lower=0> sigma; // scale parameter for likelihood function
  vector<lower=0>[Neigens] sigma2; // scale parameter for the eigen vectors
  matrix[Ntotal, Neigens] mComponents; // matrix of components or latent variables
}

model {
  matrix[Nvars, Ntotal] mFitted; // temporary variable to hold fitted values
  // prior for the sigma
  sigma ~ cauchy(0, 2.5); // weak prior
  sigma2 ~ cauchy(0, 2.5); // weak prior
  // sample the eigen vectors
  for (i in 1:Neigens){
    mEigens[,i] ~ normal(0.0, sigma2[i]);  
  }
  
  // sample the latent variables
  to_vector(mComponents) ~ normal(0, 1);
  
  // now take a sample for the data as a function of
  // Eigens X Rotations + Centering
  // i.e. Operations X Inputs
  mFitted = mEigens * mComponents';
  to_vector(y) ~ normal( to_vector(mFitted') , sigma);
}


