#------------------------------
# Preliminaries, load packages
#------------------------------

# Install the required packages
install.packages("MASS")
install.packages("flexsurv")
install.packages("rstpm2")
install.packages("coxme")
install.packages("simsurv")
install.packages("ggplot2")
install.packages("ggpubr")

library("MASS")
library("flexsurv")
library("rstpm2")
library("coxme")
library("simsurv")
library("ggplot2")
library("ggpubr")


#--------------------------------------------
# Plot four example baseline hazards for the
# two-component mixture models
#--------------------------------------------

haz <- function(t, lambdas, gammas, pi) {
  pi1 <- pi
  pi2 <- 1 - pi
  numer_1 <- (pi1 * gammas[1] * lambdas[1] * (t ^ (gammas[1] - 1)) * exp(-lambdas[1] * (t ^ gammas[1])))
  numer_2 <- (pi2 * gammas[2] * lambdas[2] * (t ^ (gammas[2] - 1)) * exp(-lambdas[2] * (t ^ gammas[2])))
  denom_1 <- (pi1 * exp(-lambdas[1] * (t ^ gammas[1])))
  denom_2 <- (pi2 * exp(-lambdas[2] * (t ^ gammas[2])))
  (numer_1 + numer_2) / (denom_1 + denom_2)
}
gg <-
  ggplot2::ggplot(data.frame(t = c(0.03, 10)), ggplot2::aes(t)) +
  ggplot2::theme_set(ggplot2::theme_bw()) +
  ggplot2::xlab("Time") +
  ggplot2::ylab("Hazard rate") +
  ggplot2::ylim(0, 1.5)
pars1 <- list(lambdas = c(1.0, 1.0), gammas = c(1.5, 0.5), pi = 0.5)
pars2 <- list(lambdas = c(0.1, 0.1), gammas = c(3.0, 1.6), pi = 0.8)
pars3 <- list(lambdas = c(1.4, 0.1), gammas = c(1.3, 0.5), pi = 0.9)
pars4 <- list(lambdas = c(1.5, 0.5), gammas = c(0.2, 0.1), pi = 0.1)
ggs <- lapply(list(pars1, pars2, pars3, pars4), function(x) {
  gg + ggplot2::stat_function(fun = haz, args = x, n = 1000) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
})
ggpubr::ggarrange(plotlist = ggs, ncol = 2, nrow = 2,
                  label.x = .53, label.y = 0.95,
                  font.label = list(size = 12, face = "plain"),
                  labels = LETTERS[1:4])


#---------------------------------------------
# Example: Weibull proportional hazards model
#---------------------------------------------

# Define function for simulation study
sim_run <- function() {
  cov <- data.frame(id = 1:200, 
                    trt = rbinom(200, 1, 0.5))
  dat <- simsurv(lambdas = 0.1, 
                 gammas = 1.5, 
                 betas = c(trt = -0.5), 
                 x = cov, 
                 maxt = 5)
  dat <- merge(cov, dat)
  
  mod <- flexsurvspline(Surv(eventtime, status) ~ trt, data = dat)
  
  est <- mod$coefficients[["trt"]]
  ses <- sqrt(diag(mod$cov))[["trt"]]
  cil <- est + qnorm(.025) * ses
  ciu <- est + qnorm(.975) * ses
  
  c(bias = est - (-0.5), 
    coverage = ((-0.5 > cil) && (-0.5 < ciu)))
}

# Set seed for reproducibility
set.seed(908070)

# Run 100 replicates in simulation study and report
# the estimated mean bias and coverage probability
rowMeans(replicate(100, sim_run()))


#----------------------------------------------------
# Example: Clustered event times (IPD meta-analysis)
#----------------------------------------------------

# Define dimensions for simulations
num_studies  <- 50
num_patients <- 200
tot_patients <- num_studies * num_patients

# Define covariate data
cov <- data.frame(
  id = 1:tot_patients,                              # patient IDs
  study = rep(1:num_studies, each = num_patients),  # study IDs
  treat = rbinom(tot_patients, 1, 0.5))             # treatment covariate

# Define 'true' parameters
pop_treat_effect <- -0.5
study_treat_effect <- pop_treat_effect + rnorm(num_studies, 0, 0.5)
pars <- data.frame(treat = rep(study_treat_effect, each = num_patients))

# Set seed for reproducibility
set.seed(908070)

# call simsurv to simulate event times
dat <- simsurv(dist = "gompertz", 
               lambdas = 0.1, 
               gammas = 0.05, 
               x = cov, 
               betas = pars)

# Merge event times onto covariate data
dat <- merge(cov, dat) 

# Display first few rows of simulated data
head(dat)

# Fit analysis model with random effect for treatment
mod <- coxme(Surv(eventtime, status) ~ treat + (treat | study), data = dat)

# Summarize analysis model (should show SD of random effect ~= 0.5)
summary(mod)


#----------------------------------------------------------------------
# Example: Flexible parametric model (ie. user-defined log cum hazard)
#----------------------------------------------------------------------

# Load and examine German Breast Cancer dataset
data("brcancer", package = "simsurv")
head(brcancer)

# Fit a Weibull model to the data
mod_weib <- flexsurvspline(Surv(rectime, censrec) ~ hormon, 
                           data = brcancer, k = 0)

# Fit a FPM model to the data
mod_flex <- flexsurvspline(Surv(rectime, censrec) ~ hormon, 
                           data = brcancer, k = 3)

# Plot the survival curve for each model, overlaid with Kaplan-Meier
par(mfrow = c(1, 2), cex = 0.85) # graphics parameters
plot(mod_weib, 
     main = "Weibull model", 
     ylab = "Survival probability", 
     xlab = "Time") 
plot(mod_flex, 
     main = "Flexible parametric model", 
     ylab = "Survival probability", 
     xlab = "Time")

# Define a function that returns the log cumulative hazard for an FPM model
logcumhaz <- function(t, x, betas, knots) {
  basis <- flexsurv::basis(knots, log(t))
  res <- 
    betas[["gamma0"]] * basis[[1]] + 
    betas[["gamma1"]] * basis[[2]] +
    betas[["gamma2"]] * basis[[3]] +
    betas[["gamma3"]] * basis[[4]] +
    betas[["gamma4"]] * basis[[5]] +
    betas[["hormon"]] * x[["hormon"]]
  res
}

# Fit the FPM model to the data, to obtain the values that will be 
# used as 'true' parameter values when simulating the event times
true_mod <- flexsurvspline(Surv(rectime, censrec) ~ hormon, 
                           data = brcancer, k = 3)

# Define a function for conducting the simulation study
sim_run <- function(true_mod) {
  
  # Covariate data
  cov <- data.frame(id = 1:200, hormon = rbinom(200, 1, 0.5))
  
  # Simulate event times using simsurv
  dat <- simsurv(betas = true_mod$coefficients, 
                 x = cov, 
                 knots = true_mod$knots, 
                 logcumhazard = logcumhaz, 
                 maxt = NULL, 
                 interval = c(1E-8, 100000))
  
  # Merge covariate data and event times
  dat <- merge(cov, dat)

  # Fit the two analysis models to the simulated dataset
  weib_mod <- flexsurvspline(Surv(eventtime, status) ~ hormon, 
                             data = dat, k = 0)

  flex_mod <- flexsurvspline(Surv(eventtime, status) ~ hormon, 
                             data = dat, k = 3)
  
  # Return the log HR for each model
  true_loghr <- true_mod$coefficients[["hormon"]]
  weib_loghr <- weib_mod$coefficients[["hormon"]]
  flex_loghr <- flex_mod$coefficients[["hormon"]]
 
  # Return the bias in the log HR under each model
  c(true_loghr = true_loghr, 
    weib_bias  = weib_loghr - true_loghr, 
    flex_bias  = flex_loghr - true_loghr)
}

# Set seed for reproducibility
set.seed(543543)

# Run 100 replicates in simulation study and report the 
# estimated mean bias under the Weibull and FPM models
rowMeans(replicate(100, sim_run(true_mod = true_mod)))


#--------------------------------------------------------------
# Example: Piecewise hazard function (ie. user-defined hazard)
#--------------------------------------------------------------

# Set seed for reproducibility
set.seed(1729)

# Specify numer of cutpoints for the piecewise hazard
ncuts <- 19

# Simulate cutpoints for the piecewise hazard
cuts  <- sort(rexp(ncuts, rate = 0.1))

# Create vector with lower limits for each time interval in piecewise hazard
pw_times <- c(0, cuts)

# Simulate a hazard rate for each time interval (increasing piecewise hazard)
N <- length(pw_times)
pw_haz <- sort(abs(rnorm(N)))

# Transform the increasing piecewise hazard to a bathtub-shaped piecewise hazard
pw_haz <- abs(pw_haz - median(pw_haz))

# Plot the bathtub-shaped piecewise hazard
dd <- data.frame(x = c(pw_times, 20), 
                 y = c(pw_haz, tail(pw_haz, 1)))
ggplot2::ggplot() + 
  ggplot2::theme_set(ggplot2::theme_bw()) +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::geom_step(ggplot2::aes(y = y, x = x), dd) +
  ggplot2::xlab("Time") +
  ggplot2::ylab("Hazard rate") +
  ggplot2::ylim(0, 2)

# Define an R function returning the true hazard for any time t
haz <- function(t, x, betas, lb_interval, haz_interval, ...){
  haz_interval[findInterval(t, lb_interval)]
}

# Simulate event times based on the bathtub-shaped piecewise hazard 
cov <- data.frame(id = 1:5000)
dat <- simsurv(x = cov, 
               hazard = haz, 
               lb_interval  = pw_times, 
               haz_interval = pw_haz)
summary(dat$eventtime)


#----------------------------------------------------------------------------
# Example: Weibull model with time-dependent effects (ie. non-prop. hazards)
#----------------------------------------------------------------------------

# Define covariate data
cov <- data.frame(id = 1:10000, 
                  trt = rbinom(10000, 1, 0.5))

# Simulate event times using simsurv
dat <- simsurv(dist = "weibull",      # baseline haz distribution
               lambdas = 0.1,         # scale parameter
               gammas = 1.5,          # shape parameter
               betas = c(trt = -0.5), # 'true' par for time-fixed coefficient
               x = cov,               # covariate data
               tde = c(trt = 0.15),   # 'true' par for interaction with log(t)
               tdefunction = "log",   # function of time to use in interaction term
               maxt = 5)              # max follow up time

# Merge covariate data and event times
dat <- merge(cov, dat)

# Display first few rows of simulated data
head(dat)

# Fit a Weibull model (with time-dependent effect) to the simulated data
mod_tvc <- stpm2(Surv(eventtime, status) ~ trt, 
                 data = dat, 
                 tvc = list(trt = 1))

# Fit a Weibull model (without time-dependent effect) to the simulated data
mod_ph <- stpm2(Surv(eventtime, status) ~ trt, 
                data = dat)

# Plot the time-fixed and time-varying hazard ratios, based on the
# two analysis models that were estimated in the previous steps.
# Also plot the 'true' time-varying hazard ratio.
plot(mod_tvc, newdata = data.frame(trt = 0), type = "hr", 
     var = "trt", ylim = c(0, 1), ci = TRUE, rug = FALSE, 
     main = "Time-dependent hazard ratio", 
     ylab = "Hazard ratio", xlab = "Time")
plot(mod_ph,  newdata = data.frame(trt = 0), type = "hr", 
     var = "trt", ylim = c(0, 1), add = TRUE, ci = FALSE, lty = 5)
curve(exp(-0.5 + 0.15 * log(x)), 0.1, to = 5, add = TRUE, lty = 3, lwd = 2)


#-------------------------------------------------------------------------
# Example: Joint longitudinal-survival model (ie. time-varying covariate)
#-------------------------------------------------------------------------

# Define a function that returns the hazard (at time t) under
# the joint longitudinal-survival model
haz <- function(t, x, betas, ...) {
  betas[["delta"]] * (t ^ (betas[["delta"]] - 1)) * exp(
    betas[["gamma_0"]] +
    betas[["gamma_1"]] * x[["x1"]] +
    betas[["gamma_2"]] * x[["x2"]] +
    betas[["alpha"]] * (
      betas[["beta_0i"]] +
      betas[["beta_1i"]] * t +
      betas[["beta_2"]]  * x[["x1"]] +
      betas[["beta_3"]]  * x[["x2"]]
    )
  )
}

# Set seed for reproducibility
set.seed(5454)

# Set number of patients
N <- 200

# Define 'true' parameters for the longitudinal and survival submodels
betas <- data.frame(
  delta   = rep(2,    N), # weibull shape parameter
  gamma_0 = rep(-11.9, N), # pop-level coef for intercept            in surv. submodel
  gamma_1 = rep(0.6,  N), # pop-level coef for binary covariate     in surv. submodel
  gamma_2 = rep(0.08, N), # pop-level coef for continuous covariate in surv. submodel
  alpha   = rep(0.03, N), # association parameter in survival submodel
  beta_0  = rep(90,   N), # pop-level coef for intercept            in long. submodel
  beta_1  = rep(2.5,  N), # pop-level coef for slope                in long. submodel
  beta_2  = rep(-1.5, N), # pop-level coef for binary covariate     in long. submodel
  beta_3  = rep(1,    N)  # pop-level coef for continuous covariate in long. submodel
)
b_corrmat <- matrix(c(1, 0.5, 0.5, 1), 2, 2) # corr matrix for patient-level ranefs
b_sds     <- c(20, 3)                        # sds         for patient-level ranefs
b_means   <- rep(0, 2)                       # means       for patient-level ranefs
b_z       <- MASS::mvrnorm(n = N, mu = b_means, Sigma = b_corrmat)
b         <- sapply(1:length(b_sds), 
                    FUN = function(x) b_sds[x] * b_z[, x]) # patient-level ranefs
betas$beta_0i <- betas$beta_0 + b[, 1] # combine pop-level and patient-level intercept
betas$beta_1i <- betas$beta_1 + b[, 2] # combine pop-level and patient-level slope

# Define covariate data
covdat <- data.frame(
  x1 = stats::rbinom(N, 1, 0.45), 
  x2 = stats::rnorm(N, 44, 8.5)
)

# Simulate event times using simsurv, by specifying the 
# user-defined hazard function (with a maximum follow up of 10 years)
times <- simsurv(hazard = haz, x = covdat, betas = betas, maxt = 10)

# Display first few simulated event times
head(times)



