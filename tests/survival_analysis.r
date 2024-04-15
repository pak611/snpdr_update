#----------------------------------------------------- SURVIVAL ANALYSIS -----------------------------------------------------------

# 1. Filter + Binomial Regression: Use filtering to create neighborhood datasets with time to event
#   respective to Ri. Results in design matrix with pheno = (time, status) and geno = (Ri, Ni). Use binomial 
#   regression on the snpdr design matrix to identify important features.


# 2. Filter + Binomial Regression: Use filtering to create neighborhood datasets with time to event
# respective to Ri. Results in design matrix with pheno = (time, status) and SNP = (Ri, Ni). Use binomial
# regression on the snpdr design matrix to identify important features. 


# 3. GLMNET (Cluster = True) Cox Regression: Use glmnet with cluster = True to perform cox regression on the design matrix

# 4. GLMNET Cox Regression

# 5. Accelerated Failure Time Model (AFT): Use the survreg function to fit an accelerated failure time model to the data

# 6. Model Evaluation: Calculate the True Positive Rate (TPR) and False Positive Rate (FPR) for the top 10% of genes selected by sNPDR

#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#-------------------------------------- FILTER + BINOMIAL REGRESSION ON GENOTYPE SURVIVAL DATA  ----------------------------------
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

library(glmnet)
library(survival)
#install.packages('devtools')
# for windows
#devtools::load_all("C:/Users/patri/Dropbox/Ph.D/Research/Doctoral_Research/Bioinformatics/sNPDR/sNPDR")

# for mac
devtools::load_all("/Users/patrickkampmeyer/Desktop/sNPDR/sNPDR")

# for windows
#file_path <- "C:/Users/patri/Dropbox/Ph.D/Research/Doctoral_Research/Bioinformatics/sNPDR/data/simulatedData/simdata.csv"

# for mac
file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/simulatedData/simdata_new.csv"


# Read in gene expression data (The important features are (3,5) by a factor of 10)
survival.data <- read.table(file_path, sep = ",", header = T)
attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
outcome <- survival.data[, ncol(survival.data)]
time <- survival.data[, ncol(survival.data)-1]

print(str(attr.mat)) ## 200 individuals, 1000 features, Tstudy = 20 (years), 20% censoring rate

# call npdr
npdr_res <- sNPDR::npdr(time, outcome, attr.mat,
                        regression.type = "binomial-surv", attr.diff.type = "numeric-abs",
                        nbd.method = "relieff", nbd.metric = "euclidean",
                        knn = 200, msurf.sd.frac = 0.5,
                        covars = "none", covar.diff.type = "numeric-abs",
                        padj.method = "bonferroni", verbose = FALSE,
                        use.glmnet = FALSE, glmnet.alpha = 0.1, glmnet.lower = 0,
                        glmnet.lam = "lambda.1se",
                        rm.attr.from.dist = c(), neighbor.sampling = "none",
                        separate.hitmiss.nbds = FALSE,
                        corr.attr.names = NULL,
                        fast.reg = FALSE, fast.dist = FALSE,
                        dopar.nn = FALSE, dopar.reg = FALSE,
                        unique.dof = FALSE, external.dist=NULL, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=0.75)


# sort by abs(beta.Z.att) in descending order, keeping the original row names
ranked.vars <- npdr_res[order(-abs(npdr_res$beta.Z.att)),]
print(head(ranked.vars, 10))


# > print(head(ranked.vars, 10))
#     att      pval.adj      pval.att beta.raw.att beta.Z.att     beta.0 pval.0
# 1    X4  0.000000e+00  0.000000e+00   -0.6366129  -38.94911 -0.7127937      0
# 2  X570 4.861792e-168 4.861792e-171   -0.4503735  -27.92600 -0.8083471      0
# 3  X938 2.651602e-151 2.651602e-154   -0.4471788  -26.50266 -0.8208364      0
# 4  X439 3.821811e-150 3.821811e-153   -0.4346320  -26.40132 -0.8214569      0
# 5  X776 1.600876e-146 1.600876e-149   -0.4172878  -26.08204 -0.8266977      0
# 6    X3 1.204743e-140 1.204743e-143   -0.4113745  -25.55569 -0.8326882      0
# 7  X227 1.496462e-132 1.496462e-135   -0.3846034  -24.81271 -0.8390880      0
# 8  X260 8.157002e-132 8.157002e-135   -0.3925143  -24.74402 -0.8398064      0
# 9  X769 8.864762e-132 8.864762e-135   -0.3956962  -24.74064 -0.8402213      0
# 10 X800 2.959393e-129 2.959393e-132   -0.3977277  -24.50379 -0.8411331      0


#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#-------------------------------------- FILTER + BINOMIAL REGRESSION ON SNP SURVIVAL DATA  ----------------------------------
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# library(glmnet)
# library(survival)
# #install.packages('devtools')
# # for windows
# #devtools::load_all("C:/Users/patri/Dropbox/Ph.D/Research/Doctoral_Research/Bioinformatics/sNPDR/sNPDR")

# # for mac
# devtools::load_all("/Users/patrickkampmeyer/Desktop/sNPDR/sNPDR")


# # read in SNP data
# file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/sReliefF_data/Episet2.csv"
# survival.data <- read.table(file_path, sep = ",", header = T)
# attr.mat <- survival.data[, -c(1, 2)]
# outcome <- survival.data[, ncol(survival.data)]
# time <- survival.data[, ncol(survival.data)-1]

# # call npdr
# npdr_res <- sNPDR::npdr(time, outcome, attr.mat,
#                         regression.type = "binomial-surv", attr.diff.type = "numeric-abs",
#                         nbd.method = "relieff", nbd.metric = "euclidean",
#                         knn = 200, msurf.sd.frac = 0.5,
#                         covars = "none", covar.diff.type = "numeric-abs",
#                         padj.method = "bonferroni", verbose = FALSE,
#                         use.glmnet = FALSE, glmnet.alpha = 0.1, glmnet.lower = 0,
#                         glmnet.lam = "lambda.1se",
#                         rm.attr.from.dist = c(), neighbor.sampling = "none",
#                         separate.hitmiss.nbds = FALSE,
#                         corr.attr.names = NULL,
#                         fast.reg = FALSE, fast.dist = FALSE,
#                         dopar.nn = FALSE, dopar.reg = FALSE,
#                         unique.dof = FALSE, external.dist=NULL, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=0.75)


# # sort by abs(beta.Z.att) in descending order, keeping the original row names
# ranked.vars <- npdr_res[order(-abs(npdr_res$beta.Z.att)),]
# print(head(ranked.vars, 5))

# # # Successfull in identifying 14 and 15 as the top two features
# # #> print(head(ranked.vars, 5))
# # #       pval.att beta.raw.att  beta.Z.att     beta.0        pval.0
# # #15 0.000000e+00  -2.06990975 -114.489792 -0.0261107  1.304219e-02
# # #14 0.000000e+00  -1.20126717  -80.018725 -0.3718214 1.336942e-285
# # #13 5.091485e-08   0.07724988    5.448593 -1.0371470  0.000000e+00
# # #6  5.176381e-07   0.06979277    5.020056 -1.0334348  0.000000e+00
# # #11 5.784996e-06   0.06006434    4.534396 -1.0349130  0.000000e+00


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#-------------------------------------- GLMNET (CLUSTER = TRUE) COX REGRESSION ON DESIGN MATRIX  ------------------
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# library(glmnet)
# library(survival)
# #install.packages('devtools')
# # for windows
# #devtools::load_all("C:/Users/patri/Dropbox/Ph.D/Research/Doctoral_Research/Bioinformatics/sNPDR/sNPDR")

# # for mac
# devtools::load_all("/Users/patrickkampmeyer/Desktop/sNPDR/sNPDR")

# # for windows
# #file_path <- "C:/Users/patri/Dropbox/Ph.D/Research/Doctoral_Research/Bioinformatics/sNPDR/data/simulatedData/simdata.csv"

# # for mac
# file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/simulatedData/simdata_new.csv"


# # Read in gene expression data
# survival.data <- read.table(file_path, sep = ",", header = T)
# attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
# outcome <- survival.data[, ncol(survival.data)]
# time <- survival.data[, ncol(survival.data)-1]

# print(str(attr.mat))

# # call npdr
# npdr_res <- sNPDR::npdr(time, outcome, attr.mat,
#                         regression.type = "cox", attr.diff.type = "numeric-abs",
#                         nbd.method = "relieff", nbd.metric = "euclidean",
#                         knn = 50, msurf.sd.frac = 0.5,
#                         covars = "none", covar.diff.type = "numeric-abs",
#                         padj.method = "bonferroni", verbose = FALSE,
#                         use.glmnet = TRUE, glmnet.alpha = "cluster", glmnet.lower = 0,
#                         glmnet.lam = "lambda.1se",
#                         rm.attr.from.dist = c(), neighbor.sampling = "none",
#                         separate.hitmiss.nbds = FALSE,
#                         corr.attr.names = NULL,
#                         fast.reg = FALSE, fast.dist = FALSE,
#                         dopar.nn = FALSE, dopar.reg = FALSE,
#                         unique.dof = FALSE, external.dist=NULL, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=0.75)


# y <- npdr_res[, -((ncol(npdr_res)-2:ncol(npdr_res)))]
# colnames(y) <- c("time", "status")
# x <- npdr_res[, ((ncol(npdr_res)-2:ncol(npdr_res)))]

# mod <- glmnet(as.matrix(x), as.matrix(y), family = "cox")

# cv.fit <- cv.glmnet(as.matrix(x), as.matrix(y), family = "cox")

# plot(cv.fit)

# coef <- coef(mod, s = -5.5)

# # Get the coefficients and their names
# coef.values <- coef@x
# coef.names <- coef@Dimnames[[1]][coef@i + 1]

# # Order the coefficients by their absolute values in decreasing order
# indices <- order(abs(coef.values), decreasing = TRUE)

# # Select the top 5 coefficients and their names
# top.coefs <- coef.values[indices[1:10]]
# top.names <- coef.names[indices[1:10]]

# # Combine the names and coefficients into a named vector
# top.coefs.named <- setNames(top.coefs, top.names)

# print(top.coefs.named)


#////////////////////////////////////////////////////////////////////////////////////////////////////////////
#---------------------------------------------- GLMNET COX MODEL --------------------------------------------
#////////////////////////////////////////////////////////////////////////////////////////////////////////////

# # Load the necessary libraries
# library(glmnet)
# library(survival)
# library(dplyr)
# # for mac
# file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/simulatedData/simdata_new.csv"

# # Read in gene expression data
# survival.data <- read.table(file_path, sep = ",", header = T)
# attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
# outcome <- survival.data[, ncol(survival.data)]
# time <- survival.data[, ncol(survival.data)-1]

# cv.fit <- cv.glmnet(x = attr.mat , y = Surv(time, outcome), family = "cox")

# best.lambda <- cv.fit$lambda.min

# fit <- glmnet(x = attr.mat, y = Surv(time, outcome), family = "cox", lambda = best.lambda)

# coefficients <- coef(fit, s = 0.5)

# # Convert to a data frame (ensure it's a data frame)
# coefficients_df <- as.data.frame(as.matrix(coefficients))

# # Add row names as a column
# coefficients_df$Feature <- row.names(coefficients_df)

# # Sort the features by the absolute values of their coefficients in descending order
# sorted_coefficients <- coefficients_df %>%
#   mutate(AbsoluteCoefficient = abs(coefficients_df[,1])) %>%
#   arrange(desc(AbsoluteCoefficient))

# # Display top ranking features
# print(head(sorted_coefficients,10))

# # > print(head(sorted_coefficients,10))
# #              1 Feature AbsoluteCoefficient
# # X4   10.428360      X4           10.428360
# # X3    8.009181      X3            8.009181
# # X717  2.299106    X717            2.299106
# # X202  2.160551    X202            2.160551
# # X82   2.105866     X82            2.105866
# # X945  1.981550    X945            1.981550
# # X708  1.909794    X708            1.909794
# # X608  1.864089    X608            1.864089
# # X604  1.650644    X604            1.650644
# # X380  1.559240    X380            1.559240


#////////////////////////////////////////////////////////////////////////////////////////////////////////////
#-------------------------------------- ACCELERATED FAILURE TIME MODEL (AFT) --------------------------------
#////////////////////////////////////////////////////////////////////////////////////////////////////////////

# # for mac
# file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/simulatedData/simdata_new.csv"


# survival.data <- read.table(file_path, sep = ",", header = T)
# attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
# outcome <- survival.data[, ncol(survival.data)]
# time <- survival.data[, ncol(survival.data)-1]

# aft.model <- survreg(Surv(time,outcome) ~ ., data = attr.mat, dist = "weibull")

# model.summary <- summary(aft.model)
# print(model.summary)

# # Get the coefficients from the model summary
# coef <- model.summary$coefficients

# # Convert the coefficients to a data frame
# coef_df <- as.data.frame(coef) %>% 
#   tibble::rownames_to_column(var = "feature")

# # Add a column with the absolute values of the coefficients
# coef_df <- coef_df %>% mutate(abs_coef = abs(coef))

# # Rank the features based on the absolute values of the coefficients
# coef_df <- coef_df %>% arrange(desc(abs_coef))

# # Print the ranked features
# print(head(coef_df, 10))

# # > print(head(coef_df, 10))
# #    feature      coef abs_coef
# # 1     X125 -81.38469 81.38469
# # 2      X78  61.08043 61.08043
# # 3     X107  54.32749 54.32749
# # 4     X168 -46.26383 46.26383
# # 5     X153 -41.39945 41.39945
# # 6     X146  39.53909 39.53909
# # 7       X8 -39.06747 39.06747
# # 8      X82 -38.71182 38.71182
# # 9      X16  37.27626 37.27626
# # 10     X58 -35.98975 35.98975

#//////////////////////////////////////////////////////////////////////////////////////////////////////
#----------------------------------------- MODEL EVALUATION ------------------------------------------
#//////////////////////////////////////////////////////////////////////////////////////////////////////
# Calculate the True Positive Rate (TPR) and False Positive Rate (FPR) for the top 10% of genes selected by sNPDR

# # Get a list of all CSV files in the simulatedData directory
# file_list <- list.files(path = "C:/Users/patri/Dropbox/Ph.D/Research/Doctoral_Research/Bioinformatics/sNPDR/data/simulatedData", pattern = "*.csv", full.names = TRUE)

# # Initialize an empty list to store the results
# npdr_results <- list()

# # Loop over the files
# for(i in seq_along(file_list)) {


# survival.data <- read.table(file_list[i], sep = ",", header = T)
# attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
# outcome <- survival.data[, ncol(survival.data)]
# time <- survival.data[, ncol(survival.data)-1]

#   # Call npdr
#   npdr_res <- sNPDR::npdr(time, outcome, attr.mat,
#                           regression.type = "binomial", attr.diff.type = "numeric-abs",
#                           nbd.method = "relieff", nbd.metric = "euclidean",
#                           knn = 100, msurf.sd.frac = 0.5,
#                           covars = "none", covar.diff.type = "numeric-abs",
#                           padj.method = "bonferroni", verbose = FALSE,
#                           use.glmnet = FALSE, glmnet.alpha = 1, glmnet.lower = 0,
#                           glmnet.lam = "lambda.1se",
#                           rm.attr.from.dist = c(), neighbor.sampling = "none",
#                           separate.hitmiss.nbds = FALSE,
#                           corr.attr.names = NULL,
#                           fast.reg = FALSE, fast.dist = FALSE,
#                           dopar.nn = FALSE, dopar.reg = FALSE,
#                           unique.dof = FALSE, external.dist=NULL, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=0.75)

#   # Sort npdr_res by beta.Z.att in descending order
#   npdr_res <- npdr_res[order(-abs(npdr_res$beta.Z.att)),]


#   # Append the result to the list
#   npdr_results[[i]] <- head(npdr_res$att,220)
# }




# tp_rate <- function(df){
#     feature.numbers  <- as.numeric(gsub("X", "", rownames(df)))
#     conv.numbers <- ifelse(feature.numbers <= 2200, 0, 1)
#     trp <- sum(conv.numbers == 1) / lenght(conv.numbers)
#     return(tpr)
# }

# tpr <- lapply(npdr_results, tp_rate)
# # Use the concordance index C-index

# # Calculate the predicted risk scores


# # Calculate the concordance index (C-index)



# Calculate the Integrated Brier Score (IBS)


