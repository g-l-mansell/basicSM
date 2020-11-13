#' Perform ordinary least squares
#'
#'@param X is a design matrix with samples as rows and variables as columns
#'@param y is a vector containing the target outcomes
#'
#'@return a vector of best fit regression coefficients
#'
#'@examples
#'ord_least_squares(X=matrix(c(1, 1.5, 3, 2.5, 3, 5), ncol=2), y=c(10, 20, 30))
ord_least_squares<- function(X, y){
  if(nrow(X) != length(y)) stop("X and y have different numbers of samples")
  w <-  solve(t(X) %*% X) %*% (t(X) %*% y)
  return(w)
}

#' Perform regularised least squares
#'
#'@param X is a design matrix with samples as rows and variables as columns
#'@param y is a vector containing the target outcomes
#'@param lambda is the regularsation constant
#'
#'@return a vector of best fit regression coefficients
#'
#'@examples
#'reg_least_squares(X=matrix(c(1, 1.5, 3, 2.5, 3, 5), ncol=2), y=c(10, 20, 30), lambda=2)
reg_least_squares<- function(X, y, lambda){
  if(nrow(X) != length(y)) stop("X and y have different numbers of samples")
  if(length(lambda) != 1) stop("Regularisation coefficient must be a scalar")
  w <- solve((lambda * diag(ncol(X))) + (t(X) %*% X)) %*% t(X) %*% y
  return(w)
}

#' Calculate regression error
#'
#'@param X is a design matrix with samples as rows and variables as columns
#'@param y is a vector containing the target outcomes
#'@param w is a vector of best fit regression coefficients
#'
#'@return the residual sum of squares
#'
#'@examples
#'error_function(X=matrix(c(1, 1.5, 3, 2.5, 3, 5), ncol=2), y=c(10, 20, 30), w=c(2, 4))
error_function <- function(X, y, w){
  diff <- y - (X %*% w)
  E <- t(diff) %*% diff
  return(as.numeric(E))
}

#' Perform simple polynomial transform
#'
#'@param x is a vector of inputs
#'@param b is the polyomial degree
#'
#'@return a desgin matrix with polynomial variables as columns
#'
#'@examples
#'poly_transform(x=c(3, 5, 1, 2.5), b=3)
poly_transform <- function(x, b){
  X <- x
  for(i in 2:b){
    X <- cbind(X, x^i)
  }
  return(X)
}

#' Perform cross validation
#'
#'@param X is a design matrix with samples as rows and variables as columns
#'@param y is a vector containing the target outcomes
#'@param k is the number of groups
#'@param FUN is the function to use for the regression, by default ord_least_squares
#'@... Additional parameter will be passed onto FUN
#'
#'@return the average testing error
#'
#'@examples
#'cross_validation(X=matrix(c(1, 1.5, 3, 7, 3, 1, 2.5, 3, 5, 10, 6.7, 3), ncol=2), 
#'                 y=c(10, 20, 30, 60, 45, 15), k=2)
#'cross_validation(X=matrix(c(1, 1.5, 3, 7, 3, 1, 2.5, 3, 5, 10, 6.7, 3), ncol=2), 
#'                 y=c(10, 20, 30, 60, 45, 15), k=2, FUN=reg_least_squares, lambda=2)

cross_validation <- function(X, y, k, FUN=ord_least_squares, ...){
  if(nrow(X) != length(y)) stop("X and y have different numbers of samples")
  
  #create the vector groups to use to allocate samples to k groups
  size <- floor(nrow(X)/k)
  remainder <- nrow(X) %% k
  groups <- c(rep(seq(1:k), size), seq(1, remainder, length.out = remainder)) 
  groups <- groups[sample(groups, length(groups))]
  
  #for each group, it the model on all the other data, then calculate the testing error
  testing_error <- rep(NA, k)
  for(i in 1:k){
    training_data_x <- X[groups!=i,]
    training_data_y <- y[groups!=i]
    coefficients <- FUN(training_data_x, training_data_y, ...)
    
    testing_data_x <- X[groups==i,]
    testing_data_y <- y[groups==i]
    testing_error[i] <- error_function(testing_data_x, testing_data_y, coefficients)
  }
  return(mean(testing_error))
}
