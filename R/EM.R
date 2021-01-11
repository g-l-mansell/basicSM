#' A function to perform k-means
#' 
#'@param X is a n by d matrix of samples
#'@param k is an integer number of clusters
#'@return list of results including predicted group labels
kmeans <- function(X, k, max_iter=20){
  n <- nrow(X)
  conv <- FALSE
  
  for(iter in 1:max_iter){
    #on the first iteration initialise
    if(iter==1){
      centroids <- X[sample(1:n, k),]
      
      #find distances of all points to centroids
      distances <- matrix(NA, n, k) 
      for(i in 1:n){
        for(j in 1:k){
          distances[i,j] <- sum((X[i,] - centroids[j,])^2) 
        }
      }
      groups <- apply(distances, 1, which.min)
    } else{
      #find centroids
      for(j in 1:k){
        centroids[j,] <- colMeans(X[groups==j,])
      }
      
      #find distances of all points to centroids
      for(i in 1:n){
        for(j in 1:k){
          distances[i,j] <- sum((X[i,] - centroids[j,])^2) 
        }
      }
      
      #store the last grouping
      prev_groups <- groups
      
      #allocate points to the closet group
      groups <- apply(distances, 1, which.min)
      
      #check for convergence
      if(all(groups == prev_groups)){
        conv <- TRUE
        break
      }
    }
  }
  
  cost <- sum(apply(distances, 1, min))
  return(list("Groups"=groups, "Cost"=cost, "Iterations"=iter, "Converged"=conv))
}

#' A function to calculate p(x | z) = N(x | mu, Sigma) for gaussian_EM
pdf <- function(x, mu, Sigma){
  D <- length(x)
  Sigma <- matrix(Sigma, D, D)
  dist <- as.numeric(x)-as.numeric(mu)
  exponent <- -0.5 * crossprod(dist, solve(Sigma, dist))
  p <- exp(exponent) / ((2*pi)^(D/2) * det(Sigma)^(1/2))
  return(as.numeric(p))
}

#' A function to preform EM for Gaussian Mixtures
#' 
#'@param X is a n by d matrix of samples
#'@param k is the number of component gaussains
#'@param max_iter is an integer limiting the number of EM iterations
#'@param init is a list providing initial parameters, e.g from k-means
#'@param tol algorithm is considered to converge when the change 
#in log-likelihood change is lower than this tolerance 

#init$means is a k by d matrix of component means
#init$vars is a k length list of d by d component covariance matrices
#init$pis is a k length vector of mixing coefficients/group priors
#if init=NULL the means are selected as random samples, 
#vars as the identity matrix and pies as 1/k

gaussian_EM <- function(X, k, max_iter=40, init=NULL, tol=1e-3){
  
  # Initialise parameters 
  n <- nrow(X)
  d <- ncol(X)
  p <- matrix(NA, nrow=n, ncol=k)
  conv <- FALSE
  
  if(is.null(init)){
    means <- X[sample(1:n, k),]
    vars <- rep(list(diag(d)), k) 
    pis <- rep(1/k, k)
  } else{
    means <- init$means
    vars <- init$vars
    pis <- init$pis
  }
  
  # Begin EM iterations
  for(iter in 1:max_iter){
    # E-step
    #find p(x,z) for each sample and each group, store in matrix p
    for(j in 1:k){
      for(s in 1:n){
        p[s,j] <- pdf(X[s,], means[j,], vars[[j]]) * pis[j]
      }
    }
    #calculate responsibilities and store in matrix r (each row sums to 1)
    r <- p / rowSums(p)
    
    #M-step
    #re-estimate parameters using current responsibilities
    Nk <- colSums(r)
    pis[] <- Nk / nsamples
    for(j in 1:k){
      means[j,] <- colSums(r[,j] * X) / Nk[j]
      sigma <- matrix(0, d, d)
      for(s in 1:n){
        sigma <- sigma + r[s,j] * tcrossprod(X[s,] - means[j,])
      }
      vars[[j]] <- sigma / Nk[j]
    }
    
    #Store last loglikelihood
    prev_ll <- ll
    
    #Calculate log likelihood of current expectations
    ll <- sum(log(rowSums(p)))    
    
    #Check for convergence (after 2 rounds)
    if(iter > 1){
      ll_change <- ll - prev_ll
      if(abs(ll_change) < tol){
        conv <- TRUE
        break
      }
    }
  }
  
  #returning a list of group allocations
  groups <- apply(r, 1, which.max)
  return(list("Groups"=groups, "LogLikelihood"=ll, "Iterations"=iter, 
              "Converged"=conv, "Means"=means, "CovarianceMatrices"=vars,
              "MixingCoeffs"=pis))
}