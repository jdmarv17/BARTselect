

### Data generation functions ###

sim.xy = function(n, p, nval, rho = 0, s = 5, beta.type = 2, snr = 1, seed = 1) {
  
  set.seed(seed)
  # Generate predictors
  x = matrix(rnorm(n*p),n,p)
  xval = matrix(rnorm(nval*p),nval,p)
  
  # Introduce autocorrelation, if needed
  if (rho != 0) {
    inds = 1:p
    Sigma = rho^abs(outer(inds, inds, "-"))
    obj = svd(Sigma)
    Sigma.half = obj$u %*% (sqrt(diag(obj$d))) %*% t(obj$v)
    x = x %*% Sigma.half
    xval = xval %*% Sigma.half
  }
  else Sigma = diag(1,p)
  
  # Generate underlying coefficients
  s = min(s,p)
  beta = rep(0,p)
  if (beta.type==1) {
    beta[round(seq(1,p,length=s))] = 1
  } else if (beta.type==2) {
    beta[1:s] = 1
  } else if (beta.type==3) {
    beta[1:s] = seq(10,0.5,length=s)
  } else if (beta.type==4) {
    beta[1:6] = c(-10,-6,-2,2,6,10)
  } else {
    beta[1:s] = 1
    beta[(s+1):p] = 0.5^(1:(p-s))
  }
  
  # Set snr based on sample variance on infinitely large test set
  vmu = as.numeric(t(beta) %*% Sigma %*% beta)
  sigma = sqrt(vmu/snr)
  
  
  # Generate responses
  y = as.numeric(x %*% beta + rnorm(n)*sigma)
  yval = as.numeric(xval %*% beta + rnorm(nval)*sigma)
  
  #enlist(x,y,xval,yval,Sigma,beta,sigma)
  
  to.return = list(x, y, xval, yval, Sigma, beta, sigma)
  return(to.return)
}

# THIS FUNCTION MAY NOT BE DONIG SNR CORRECTLY ONCE INTERACTIONS ARE ADDED
sim.xy.interact = function(n = 500, p = 100, nval = 10, rho = 0, s = 5, beta.type = 1, 
                           snr = 1, n.interact = 1, coef.interact = 1, seed = 1) {
  
  ##### This function will take sim.xy arguments and arguments for 
  ##### number of interactions and coefficient for interaction:
  ##### Return 1: a matrix with cols for covariates and interactions
  ##### Return 2: a vector of the outcome y
  ##### Return 3-4: the same as return 1-2 but for test set
  ##### Return 5: true covariance matrix of covariates
  ##### Return 6: true beta vector used to generate outcome from covariates
  ##### Return 7: amount of noise added to outcome y 
  ##### Return 7: pairs of variables selected to interact
  # ARGS:
  # `n` = number of samples
  # `p` = number of covaraites
  # `nval` = number of test samples
  # `rho` = amount of covariate correlation
  # `s` = number of signal variables
  # `beta.type` = beta type from Hastie, Tibshirani, et al
  # `snr` = signal to noise ratio
  # `n.interact` = base value for prior governing tree depth
  # `coef.interact` = base value for prior governing tree depth
  
  
  set.seed(seed)
  # Generate predictors and test set predictors
  x = matrix(rnorm(n * p), n, p)
  xval = matrix(rnorm(nval * p), nval, p)
  
  # Introduce autocorrelation, if needed
  if (rho != 0) {
    inds = 1:p
    Sigma = rho^abs(outer(inds, inds, "-"))
    obj = svd(Sigma)
    Sigma.half = obj$u %*% (sqrt(diag(obj$d))) %*% t(obj$v)
    x = x %*% Sigma.half
    xval = xval %*% Sigma.half
  } else Sigma = diag(1, p)
  
  # Generate underlying coefficients
  s = min(s, p)
  beta = rep(0, p)
  if (beta.type == 1) {
    beta[round(seq(1, p, length = s))] = 1
  } else if (beta.type == 2) {
    beta[1:s] = 1
  } else if (beta.type == 3) {
    beta[1:s] = seq(10, 0.5, length=s)
  } else if (beta.type == 4) {
    beta[1:6] = c(-10, -6, -2, 2, 6, 10)
  } else {
    beta[1:s] = 1
    beta[(s+1):p] = 0.5 ^ (1:(p - s))
  }
  
  
  
  # get indices of non-zero main effects
  possible.ind.logical = c(beta != 0) 
  possible.ind = c(1:p)[possible.ind.logical] # possible col numbers to interact
  possible.pairs = t(combn(possible.ind, 2)) # permutes possible pairs
  
  # randomly select from paired main effect covariates n.interact to interact
  interact.ind = sample(nrow(possible.pairs), size = n.interact)
  pairs = possible.pairs[interact.ind, ]
  
  # save version with indices
  pairs_numeric = pairs
  
  
  # if only one interaction need to put it in matrix in order to loop through its rows
  if (n.interact == 1) {
    pairs = paste0("x", pairs)
    pairs = matrix(nrow = 1, data = c(pairs))
    pairs_numeric = matrix(nrow = 1, c(pairs_numeric))
  } else {
    for (i in 1:nrow(pairs)) {
      pairs[i,] <- paste0("x", pairs[i,])
    }
  }
  
  for (i in 1:nrow(pairs)) {
    
    # grab first pair to interact and get product
    col.ind = pairs_numeric[i,]
    new.interact = x[ ,col.ind[1]] * x[ ,col.ind[2]]
    new.interact.val = xval[ ,col.ind[1]] * xval[ ,col.ind[2]]
    
    # bind new interaction columns to x and xval
    x = cbind(x, unname(new.interact))
    xval = cbind(xval, unname(new.interact.val))
  } # end of i index
  
  
  
  # add some colnames
  names = paste0("x", 1:(p + n.interact))
  colnames(x) <- names
  colnames(xval) <- names
  
  # append interaction coefficients to beta vector and store
  # beta before appending in order to set sigma
  beta.small = beta
  beta = c(beta, rep(coef.interact, each = n.interact))
  
  # Set snr based on sample variance on infinitely large test set
  # DOING THIS WITH SIGMA AND BETA WITHOUT INTERACTIONS 
  # --> THIS MAY BE MAKING vmu/snr/sigma DIFFERENT
  vmu = as.numeric(t(beta.small) %*% Sigma %*% beta.small)
  sigma = sqrt(vmu/snr)
  
  # Generate responses
  y = as.numeric(x %*% beta + rnorm(n)*sigma)
  yval = as.numeric(xval %*% beta + rnorm(nval)*sigma)
  
  # remove cols that are just the interactions
  x = as.data.frame(x[,-c((ncol(x) - (n.interact - 1)):ncol(x))])
  
  # make new col with both vars
  pairs = 
    pairs %>%
    as.data.frame() %>%
    unite(., col = "paired", sep = ":", remove = FALSE)
  
  
  to.return = list(x, y, xval, yval, Sigma, beta, sigma, pairs)
  return(to.return)
}


# similar function to generate binary data from data with an interaction
generate.binary = function(n = 500, p = 100, sigma = 1, seed = 1) {
  
  set.seed(seed)
  # generate predicotrs and label them
  x = matrix(rnorm(n * p), n, p)
  colnames(x) = paste0("x", 1:p)
  
  # form two interactions out of cols 1-4 and bind to 'x'
  interact.col1 = x[,1] * x[,2]
  interact.col2 = x[,3] * x[,4]
  x = cbind(x, interact.col1, interact.col2)
  
  # form linear predictor with X * \beta + \epsilon and draw from binom(1, p) n times
  lin.pred = x %*% c(rep(1, 3), rep(-1, 2), rep(0, (p - 5)), 1, -1) + rnorm(n = n, mean = 0, sd = sigma)
  bern.prob = exp(lin.pred) / (1 + exp(lin.pred))
  y = rbinom(n = n, size = 1, prob = bern.prob)
  x = x[, 1:p]
  
  # return data
  data = data.frame(y = y, x) 
  return(data)
}


# generates friedman function data
friedman_func = function(n = 500, p = 100, sigma = 1, seed = 1){
  
  set.seed(seed)
  # friedman function
  f = function(x){
    10 * sin(pi * x[ ,1] * x[ ,2]) + 20 * (x[ ,3] - 0.5)^2 + 10 * x[ ,4] + 5 * x[ ,5]
  }

  x = matrix(nrow = n, ncol = p, runif(n * p)) # p variables, only first 5 matter
  colnames(x) = paste0("x", 1:p)
  Ey = f(x)
  y = Ey + sigma * rnorm(n)
  
  to_return = list(y, x)
  return(to_return)
}

# data gen function for variable selection, not interaction selection
non.lin.func = function(n = 500, p = 100, sigma = 1, seed = 1) {
  
  set.seed(seed)
  x = matrix(nrow = n, ncol = p, runif(n * p))
  colnames(x) = paste0("x", 1:p)
  
  f.x = log(x[,1] * x[,2]) + x[,3]^2 + 0.5 * x[,4] + (1 / x[,5])
  y = f.x + sigma * rnorm(n)
  
  to.return = list(y, x)
  return(to.return)
}

f1 <- function(x) {
  x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]; x5 <- x[5]
  return(10 * sin(pi * x1 * x2) + 20 * (x3 - 0.5)^2 + 10 * x4 + 5 * x5)
}

f2 <- function(x) {
  x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]; x5 <- x[5]
  return(log(x1 * x2) + x3^2 + 0.5 * x4 + (1 / x5))
}

fried.taylor = function(x) {
  y = -4.03613 + (5 * pi / sqrt(2)) * x[1] + (5 * pi / sqrt(2)) * x[2] + 10 * x[4] + 5 * x[5]
  return(y)
}

nonlin.taylor = function(x) {
  y = 0.363706 + 2 * x[1] + 2 * x[2] + x[3] + 0.5 * x[4] - 4 * x[5]
  return(y)
}
