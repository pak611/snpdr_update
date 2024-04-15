install.packages('coxed')
library(coxed)



num.covariates <- 1000
# Lets try user defined beta values (coefficients)
factor <- 10

coefficients <- rep(1,num.covariates)

coefficients[c(3,5)] <- coefficients[c(3,5)] * factor

# Make T bigger and decrease the variance of the X variables
# N = 200: number of observations
# T = 5: maximum time
# xvars = 19: number of covariates
# censor = 0.2: censoring rate (proportion of observations that are censored)
# num.data.frames = 1: number of datasets to simulate
simdata <- sim.survdata(N=200, T=20, xvars = num.covariates, censor = 0.2, num.data.frames=1, beta=coefficients/2)

# View the first few rows of the simulated data
str(simdata)

# round simsuve to the nearest integer


# output simdata to a csv
write.csv(simdata$data, "C:/Users/patri/Dropbox/Ph.D/Research/Doctoral_Research/Bioinformatics/sNPDR/data/simulatedData/simdata.csv")


