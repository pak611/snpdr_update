#' Generate Survival Data
#'
#' This function generates synthetic survival data based on a combination of main effects,
#' interaction effects, and random noise. It also visualizes the generated Kaplan-Meier survival curves.
#'
#' @param n Number of observations. Default is 1000.
#' @param p Number of features. Default is 5.
#' @param max_observation_time Maximum observation time for censoring. Default is 100.
#' @param intEff Include interactions
#'
#' @return A dataframe containing the generated survival data, with columns for the features, survival times, censoring status, and the survival object.
#' @export
#'
#' @examples
#' \dontrun{
#' data <- generate_survival_data(n=1000, p=5, max_observation_time=100)
#' head(data)
#' }


generate_survival_data <- function(n = 1000, p = 5, max_observation_time = 100, intEff = FALSE) {
  
  # Load required libraries
  library(survival)
  
  # Generate random features
  set.seed(42)
  data <- as.data.frame(matrix(rnorm(n * p), ncol=p))
  
  # Coefficients for main effects for each feature
  beta <- runif(p, -2, 2) # Random coefficients between -2 and 2 for simplicity
  
  # Get the order for beta values in descending order
  order_beta_desc <- order(beta, decreasing = TRUE)
  
  # Get corresponding attribute names based on the order
  attributes_ordered <- colnames(data)[order_beta_desc]
  
  # Print ordered beta values and their corresponding attributes
  cat("Ordered attributes:", attributes_ordered, "\n")
  cat("Corresponding betas:", beta[order_beta_desc], "\n")
  

  # Calculate main effects
  main_effects <- as.matrix(data[,1:p]) %*% unlist(beta)

  
  # Simulate epistasis effect
  if (intEff == TRUE) {
    interaction_effect <- ifelse(data$V1 > 1 & data$V2 > 1, -15, 0)
  } else {
    interaction_effect <- 0.0
  }
  

  
  # Generate survival times (with main effects, interaction effect, and some randomness)
  data$survival_time <- rexp(n, rate=0.05) + main_effects + interaction_effect + rnorm(n, mean=0, sd=5)
  
  data$survival_time[data$survival_time < 0] <- 0.01
  
  # Introduce censoring
  data$censored <- ifelse(data$survival_time > max_observation_time, 1, 0)
  data$survival_time[data$censored == 1] <- max_observation_time
  
  # Create a Surv object
  data$surv_obj <- with(data, Surv(survival_time, censored == 0))
  
  # Visualize the survival curves for different groups
  km_fit_all <- survfit(data$surv_obj ~ 1, data=data)
  km_fit_interaction <- survfit(data$surv_obj ~ (V1 > 1 & V2 > 1), data=data)
  
  plot(km_fit_all, col="black", lty=1, main="Kaplan-Meier survival curves", xlab="Time", ylab="Survival probability")
  lines(km_fit_interaction, col=c("red", "blue"), lty=1:2)
  legend("bottomleft", legend=c("All", "Interaction: FALSE", "Interaction: TRUE"), col=c("black", "red", "blue"), lty=1:3)
  
  return(data)
}
