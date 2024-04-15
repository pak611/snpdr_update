#' Cox Regression
#'
#' This script defines an objective function for a Cox regression model and uses
#' the `optim` function to estimate the model's beta coefficients.
#'
#' @author Your Name
#' @seealso \code{\link{optim}}
#' @keywords internal
#'
#' @examples
#' # Define your data matrix X here
#' # X <- ...
#' 
#' # Initial guess for beta
#' beta_init <- rep(0, ncol(X))
#' 
#' # Use optim to find the beta that minimizes -L (maximizes L)
#' optim_result <- optim(par = beta_init, fn = objective_function, X = X, method = "BFGS")
#' 
#' # Extract the estimated beta coefficients
#' beta_est <- optim_result$par
#' 
#' @export
#' # Define the objective function L(beta)
objectiveFunction <- function(beta, y, X) {
  #browser()

  y <- as.numeric(as.character(y))
  X <- as.matrix(sapply(X, as.numeric))
  # Compute the linear term beta^T * X for each i
  linear.terms <- y%*%X * beta
  
  # Compute the log-sum-exp term
  log.sum.exp <- log(sum(exp(-(y-1)%*%X * beta)))
  
  # Sum over all i (assuming X is a matrix where each row is an observation)
  L <- sum(linear.terms) - log.sum.exp
  
  # Return negative because optim() minimizes by default
  return(-L)
}

cox_mle <- function(y, X) {
  # Use optim to find the beta that minimizes -L (maximizes L)

  # Initial guess for beta
 
  optimResult <- optim(par = 0, fn = objectiveFunction, y = y, X = X, method = "BFGS")
  
  # Extract the estimated beta coefficients
  betaEst <- optimResult$par
  
  return(betaEst)
}

