
#//////////////////////////////////////////////////////////////////////////////////////////////////#
#------------------------------------------ SIMULATE SURVIVAL DATA -------------------------------#
#//////////////////////////////////////////////////////////////////////////////////////////////////#


# for windows
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/simSurvData/cox_epistasis")

# for mac
#devtools::load_all("/Users/patrickkampmeyer/Desktop/sNPDR/simSurvData/cox_epistasis")


#------------------------------------------ SINGLE SET ----------------------------------------------#

num.covariates <- 100
num.inst <- 100
# Lets try user defined beta values (coefficients)
factor <- 30

coefficients <- rep(0.1,num.covariates)
coefficients[3:4] <- coefficients[3:4] * factor



# specify the interactions matrix as a [1:num.inst, 1:num.covariates] matrix
# Initialize the matrix with zeros
inter.mat <- matrix(0, nrow=num.covariates, ncol=num.covariates)

# Set the elements at [1,2] and [2,1] to 1
inter.mat[1,2] <- 1
inter.mat[2,1] <- 1

# Make T bigger and decrease the variance of the X variables
# N = 200: number of observations
# T = 5: maximum time
# xvars = 19: number of covariates
# censor = 0.2: censoring rate (proportion of observations that are censored)
# num.data.frames = 1: number of datasets to simulate

#For interaction effects
#simdata <- coxed::sim.survdata(N=num.inst, T=40, xvars = num.covariates, censor = 0.1, num.data.frames=1, beta=coefficients/2, interactions=TRUE, inter.mat=inter.mat, mu=0, sd=0.5)

#For no interaction effects
simdata <- coxed::sim.survdata(N=num.inst, T=20, xvars = num.covariates, censor = 0.1, num.data.frames=1, beta=coefficients/2, interactions=FALSE, inter.mat=inter.mat, mu=0, sd=0.5)

# output simdata to a csv
write.csv(simdata$data, "/Users/patrickkampmeyer/Desktop/sNPDR/data/simulatedData/simdata_new.csv")



#------------------------------------------ MULTIPLE SETS ----------------------------------------------#


num.covariates <- 100
num.inst <- 400
num.sig <- 2
# Lets try user defined beta values (coefficients)
factor <- 30
coefficients <- rep(0.1,num.covariates)
coefficients[1:num.sig] <- coefficients[1:num.sig] * factor



# specify the interactions matrix as a [1:num.inst, 1:num.covariates] matrix
# Initialize the matrix with zeros
inter.mat <- matrix(0, nrow=num.covariates, ncol=num.covariates)

# Set the elements at [1,2] and [2,1] to 1
inter.mat[1,2] <- 1
inter.mat[2,1] <- 1

# Generate data in a loop

for (i in seq(1,10)) {
        
    #For no interaction effects
    simdata <- sim.survdata(N=num.inst, T=100, xvars = num.covariates, censor = 0.1, num.data.frames=1, beta=coefficients/2, interactions=FALSE, inter.mat=inter.mat, mu=0, sd=0.5)

    # output simdata to a csv

    # for windows
    write.csv(simdata$data, paste0("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData/data_",as.character(i),".csv"))
}





