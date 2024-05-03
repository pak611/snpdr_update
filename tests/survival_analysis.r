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
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/sNPDR")

# for mac
#devtools::load_all("/Users/patrickkampmeyer/Desktop/sNPDR/sNPDR")


#-------------------------------------------Genotype Data---------------------------------------------

file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/sReliefF_data/EpiSet2.txt"
# Read in the genotype data

survival.data <- read.table(file_path, sep = "\t", header = FALSE, skip = 1, quote = "")

# Read the first line to get the column names as a single string
headers <- readLines(file_path, n = 1)

# Split the header string into individual column names using strsplit
column_names <- strsplit(headers, "\t", fixed = TRUE)[[1]]

colnames(survival.data) <- column_names

attr.mat <- survival.data[, -c(1,2)]
outcome <- as.numeric(survival.data[, 1])
time <- survival.data[, 2]


#-------------------------------------------Gene-Expression Data---------------------------------------------

# #for mac
# #file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/simulatedData/simdata_new.csv"


# # for windows
# file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData/data_2.csv"


# # Read in gene expression data (The important features are (3,5) by a factor of 10)
# survival.data <- read.table(file_path, sep = ",", header = T)
# attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
# outcome <- survival.data[, ncol(survival.data)]
# time <- survival.data[, ncol(survival.data)-1]

# call npdr
npdr_res <- sNPDR::npdr(time, outcome, attr.mat,
                        regression.type = "binomial-surv", attr.diff.type = "numeric-abs",
                        nbd.method = "relieff", nbd.metric = "euclidean",
                        knn = 50, msurf.sd.frac = 0.5,
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
#        pval.att beta.raw.att beta.Z.att     beta.0 pval.0
# 17 4.454345e-11  -0.20768294  -6.591273 -0.8082828      0
# 12 1.941346e-06  -0.15179405  -4.760650 -0.8228472      0
# 3  1.408919e-05  -0.11247348  -4.343396 -0.8047119      0
# 9  1.820300e-03  -0.08615807  -3.118446 -0.8286300      0
# 15 7.708959e-03  -0.06966812  -2.664791 -0.8286299      0
# 7  9.209130e-03  -0.06753177  -2.604410 -0.8277701      0
# 2  2.005770e-02  -0.06099279  -2.325427 -0.8368737      0
# 5  5.471251e-02  -0.05700426  -1.921249 -0.8471894      0
# 13 6.891331e-02   0.05753530   1.819072 -0.8756835      0
# 6  1.237605e-01   0.04607367   1.539235 -0.8754902      0

# Rank	Variable		Score
# 1	Var12		0.0759
# 2	Var19		0.043
# 3	Var9		0.0407
# 4	Var4		0.0343
# 5	Var11		0.0173
# 6	Var5		-0.0055
# 7	Var16		-0.0142
# 8	Var15		-0.0503
# 9	Var0		-0.0553
# 10	Var8		-0.062
# 11	Var1		-0.1321
# 12	Var13		-0.1412
# 13	Var2		-0.1534
# 14	Var14		-0.1538
# 15	Var10		-0.166
# 16	Var3		-0.1686
# 17	Var7		-0.1855
# 18	Var6		-0.1865
# 19	Var18		-0.1882
# 20	Var17		-0.2117



#-------------------------------------------Gene-Expression Data---------------------------------------------
# > print(head(ranked.vars, 10))
#        pval.att beta.raw.att beta.Z.att    beta.0 pval.0
# 74 2.330272e-10    0.3216577   6.341351 -2.411518      0
# 56 2.758856e-10    0.3314397   6.315237 -2.415650      0
# 43 4.224716e-08   -0.3278016  -5.483489 -2.043583      0
# 73 1.074449e-07   -0.2883822  -5.315758 -2.052037      0
# 48 2.074789e-07   -0.2846058  -5.194469 -2.054083      0
# 46 2.492382e-07   -0.2933707  -5.160196 -2.055562      0
# 78 5.084550e-07   -0.2819113  -5.024873 -2.058599      0
# 82 1.503624e-06   -0.2836475  -4.812321 -2.062809      0
# 91 1.719346e-06   -0.2690923  -4.785433 -2.069016      0
# 75 4.726613e-06   -0.2445048  -4.577926 -2.072387      0
#-------------------------------------------Gene-Expression Data (revised) ---------------------------------------------
# > print(head(ranked.vars, 10))
#        pval.att beta.raw.att beta.Z.att    beta.0 pval.0
# 87 9.091136e-17    0.3099250   8.323973 -1.568922      0
# 48 7.556039e-14   -0.3006592  -7.483577 -1.201872      0
# 19 9.523886e-12    0.2435935   6.817884 -1.528022      0
# 49 9.429869e-11   -0.2702944  -6.479574 -1.221949      0
# 96 1.141282e-10   -0.2597223  -6.450649 -1.227314      0
# 17 6.646412e-10   -0.2436142  -6.177545 -1.236406      0
# 61 1.540816e-09    0.2364560   6.043091 -1.517799      0
# 56 2.867202e-09    0.2394930   5.941935 -1.517156      0
# 89 1.081049e-08    0.2194685   5.720097 -1.511783      0
# 64 2.242962e-08   -0.2222094  -5.594570 -1.245819      0



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

# # srelieff: 
# 1	V7		0.0417
# 2	V4		0.0356
# 3	V9		0.0285
# 4	V1		0.0229
# 5	V2		0.0177
# 6	V8		0.0167
# 7	V13		0.0113
# 8	V6		0.0053
# 9	V10		0.0022
# 10	V11		-0.0174
# 11	V3		-0.0183
# 12	V12		-0.0348
# 13	V14		-0.0387
# 14	V15		-0.0591
# 15	V5		-0.0653


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#-------------------------------------- GLMNET (CLUSTER = TRUE) COX REGRESSION ON DESIGN MATRIX  ------------------
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


library(glmnet)
library(survival)
#install.packages('devtools')
# for windows
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/sNPDR")

# for mac
#devtools::load_all("/Users/patrickkampmeyer/Desktop/sNPDR/sNPDR")

# for windows
file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData/data_1.csv"

# for mac
#file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/simulatedData/simdata_new.csv"


file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/sReliefF_data/EpiSet2.txt"
# Read in the genotype data

file_path <- "C:/Users/patri/OneDrive/Desktop/Episet_Datasets/dataset_0.txt"

survival.data <- read.table(file_path, sep = "\t", header = FALSE, skip = 1, quote = "")

# Read the first line to get the column names as a single string
headers <- readLines(file_path, n = 1)

# Split the header string into individual column names using strsplit
column_names <- strsplit(headers, "\t", fixed = TRUE)[[1]]

colnames(survival.data) <- column_names

attr.mat <- survival.data[, -c(1,2)]
outcome <- as.numeric(survival.data[, 1])
time <- survival.data[, 2]

# # Read in gene expression data
# survival.data <- read.table(file_path, sep = ",", header = T)
# attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
# outcome <- survival.data[, ncol(survival.data)]
# time <- survival.data[, ncol(survival.data)-1]

# print(str(attr.mat))

# call npdr
npdr_res <- sNPDR::npdr(time, outcome, attr.mat,
                        regression.type = "cox", attr.diff.type = "numeric-abs",
                        nbd.method = "relieff", nbd.metric = "euclidean",
                        knn = 50, msurf.sd.frac = 0.5,
                        covars = "none", covar.diff.type = "numeric-abs",
                        padj.method = "bonferroni", verbose = FALSE,
                        use.glmnet = TRUE, glmnet.alpha = "cluster", glmnet.lower = 0,
                        glmnet.lam = "lambda.1se",
                        rm.attr.from.dist = c(), neighbor.sampling = "none",
                        separate.hitmiss.nbds = FALSE,
                        corr.attr.names = NULL,
                        fast.reg = FALSE, fast.dist = FALSE,
                        dopar.nn = FALSE, dopar.reg = FALSE,
                        unique.dof = FALSE, external.dist=NULL, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=0.75)


y <- npdr_res[, -((ncol(npdr_res)-2:ncol(npdr_res)))]
colnames(y) <- c("time", "status")
x <- npdr_res[, ((ncol(npdr_res)-2:ncol(npdr_res)))]

mod <- glmnet(as.matrix(x), as.matrix(y), family = "cox")

cv.fit <- cv.glmnet(as.matrix(x), as.matrix(y), family = "cox")

plot(cv.fit)

coef <- coef(mod, s = -6.0)

# Get the coefficients and their names
coef.values <- coef@x
coef.names <- coef@Dimnames[[1]][coef@i + 1]

# Order the coefficients by their absolute values in decreasing order
indices <- order(abs(coef.values), decreasing = TRUE)

# Select the top 5 coefficients and their names
top.coefs <- coef.values[indices[1:10]]
top.names <- coef.names[indices[1:10]]

# Combine the names and coefficients into a named vector
top.coefs.named <- setNames(top.coefs, top.names)

print(top.coefs.named)


#////////////////////////////////////////////////////////////////////////////////////////////////////////////
#---------------------------------------------- GLMNET COX MODEL --------------------------------------------
#////////////////////////////////////////////////////////////////////////////////////////////////////////////

# # Load the necessary libraries
# library(glmnet)
# library(survival)
# library(dplyr)
# # for mac
# file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData/data_1.csv"

# # Read in gene expression data
# survival.data <- read.table(file_path, sep = ",", header = T)
# attr.mat <- as.matrix(survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))])
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

# library(survival)
# # # for mac
# # file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/simulatedData/simdata_new.csv"

# # for windows
# file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData/simdata_1.csv"


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


# coef_df$feature <- sub("X", "", coef_df$feature)

# rownames(coef_df) <- coef_df$feature

# coef_df$feature <- NULL

# # Remove the row with the name "(Intercept)"
# coef_df <- subset(coef_df, rownames(coef_df) != "(Intercept)")



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


# # #------------------------------------------------- SNPDR ------------------------------------------------
# #Calculate the True Positive Rate (TPR) and False Positive Rate (FPR) for the top 10% of genes selected by sNPDR

# #Get a list of all CSV files in the simulatedData directory

# #load the sNPDR package

# #for mac
# #devtools::load_all("/Users/patrickkampmeyer/Desktop/sNPDR/sNPDR")

# # for windows
# devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/sNPDR")

# # for mac
# #file_list <- list.files(path = "C:/Users/patri/Dropbox/Ph.D/Research/Doctoral_Research/Bioinformatics/sNPDR/data/simulatedData", pattern = "*.csv", full.names = TRUE)

# # for windows
# file_list <- list.files(path = "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData", pattern = "*.csv", full.names = TRUE)


# # number of significant features
# num.sig <- 10

# # Initialize an empty list to store the results
# npdr_results <- list()

# # Loop over the files
# for(i in seq_along(file_list)) {


#     survival.data <- read.table(file_list[i], sep = ",", header = T)
#     attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
#     outcome <- survival.data[, ncol(survival.data)]
#     time <- survival.data[, ncol(survival.data)-1]


#     npdr_res <- sNPDR::npdr_revised(time, outcome, attr.mat,
#                             regression.type = "binomial-surv", attr.diff.type = "numeric-abs",
#                             nbd.method = "relieff", nbd.metric = "euclidean",
#                             knn = 50, msurf.sd.frac = 0.5,
#                             covars = "none", covar.diff.type = "numeric-abs",
#                             padj.method = "bonferroni", verbose = FALSE,
#                             use.glmnet = FALSE, glmnet.alpha = 0.1, glmnet.lower = 0,
#                             glmnet.lam = "lambda.1se",
#                             rm.attr.from.dist = c(), neighbor.sampling = "none",
#                             separate.hitmiss.nbds = FALSE,
#                             corr.attr.names = NULL,
#                             fast.reg = FALSE, fast.dist = FALSE,
#                             dopar.nn = FALSE, dopar.reg = FALSE,
#                             unique.dof = FALSE, external.dist=NULL, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=0.75)

#     # Sort npdr_res by beta.Z.att in descending order
#     npdr_res <- npdr_res[order(-abs(npdr_res$beta.Z.att)),]


#     # Append the result to the list
#     npdr_results[[i]] <- head(npdr_res,num.sig)
# }


# # source in model_val function
# source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")

# results <- model_eval(feature.ranking = npdr_results, important.features = seq(1,num.sig))

# print(results)


#*****************************************      DATASET 1      ***************************************************

#--- 100 covariates, 100 individuals, effect size = 30, 10 significant features, KM.weight=FALSE)
# snpdr
# > print(results)
#   precision recall
# 1      0.18   0.18
# npdr_revised
# > print(results)
#   precision recall
# 1      0.11   0.11

#--- 100 covariates, 100 individuals, effect size = 30, 10 significant features, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=0.75)
# snpdr
# > print(results)
#   precision recall
# 1       0.2    0.2
# npdr_revised
# > print(results)
#   precision recall
# 1      0.12   0.12

#--- 100 covariates, 100 individuals, effect size = 30, 10 significant features, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=0.25)
# > print(results)
#   precision recall
# 1       0.1    0.1

#--- 100 covariates, 100 individuals, effect size = 30, 10 significant features, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=2.0)
# > print(results)
#   precision recall
# 1      0.19   0.19

#*****************************************      DATASET 2     ***************************************************

# REDUCE THE NUMBER OF SIGNIFICANT FEATURES

#--- 100 covariates, 100 individuals, effect size = 30, 5 significant features, KM.weight=FALSE)
# > print(results)
#   precision recall
# 1      0.16   0.16

#--- 100 covariates, 100 individuals, effect size = 30, 5 significant features, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=0.75)
# > print(results)
#   precision recall
# 1       0.2    0.2

#*****************************************      DATASET 3     ***************************************************

# INCREASE THE NUMBER OF INSTANCES AND FEATURES

#--- 1000 covariates, 500 individuals, effect size = 30, 10 significant features, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=0.75)
# > print(results)
#   precision recall
# 1      0.06   0.06

#--------------------------------------------------------- ACCELERATED FAILURE TIME MODEL (AFT) ------------------------------------------------

# library(survival)
# library(dplyr)

# source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")

# # for windows
# file_list <- list.files(path = "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData", pattern = "*.csv", full.names = TRUE)

# # number of significant features
# num.sig <- 10

# # Initialize an empty list to store the results
# aft_results <- list()


# for (i in seq_along(file_list)) {
    
#     print(i)

#     survival.data <- read.table(file_list[i], sep = ",", header = T)
#     attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
#     outcome <- survival.data[, ncol(survival.data)]
#     time <- survival.data[, ncol(survival.data)-1]

#     aft.model <- survreg(Surv(time,outcome) ~ ., data = attr.mat, dist = "weibull")

#     model.summary <- summary(aft.model)

#     # Get the coefficients from the model summary
#     coef <- model.summary$coefficients

#     # Convert the coefficients to a data frame
#     coef_df <- as.data.frame(coef) %>% 
#     tibble::rownames_to_column(var = "feature")

#     # Add a column with the absolute values of the coefficients
#     coef_df <- coef_df %>% mutate(abs_coef = abs(coef))

#     # Rank the features based on the absolute values of the coefficients
#     coef_df <- coef_df %>% arrange(desc(abs_coef))

#     coef_df$feature <- sub("X", "", coef_df$feature)

#     rownames(coef_df) <- coef_df$feature

#     # Remove the row with the name "(Intercept)"
#     coef_df <- subset(coef_df, rownames(coef_df) != "(Intercept)")

#     coef_df$feature <- NULL

#     # Append the result to the list
#     aft_results[[i]] <- head(coef_df,num.sig)
# }


# # source in model_val function

# source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")

# results <- model_eval(feature.ranking = aft_results, important.features = seq(1,num.sig))

# print(results)

#*****************************************      DATASET 1     ***************************************************
#--- 100 covariates, 100 individuals, effect size = 30, 10 significant features)

# > print(results)
#   precision recall
# 1      0.62   0.62
#*****************************************      DATASET 2     ***************************************************




#*****************************************      DATASET 3     ***************************************************



# #--------------------------------------------------------- GLMNET COX MODEL ------------------------------------------------

# # Load the necessary libraries
# library(glmnet)
# library(survival)
# library(dplyr)



# # for windows
# file_list <- list.files(path = "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData", pattern = "*.csv", full.names = TRUE)


# # number of significant features
# num.sig <- 10

# # Initialize an empty list to store the results
# glm_results <- list()



# for (i in seq_along(file_list)){

#     # Read in gene expression data
#     survival.data <- read.table(file_list[i], sep = ",", header = T)
#     attr.mat <- as.matrix(survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))])
#     outcome <- survival.data[, ncol(survival.data)]
#     time <- survival.data[, ncol(survival.data)-1]

#     cv.fit <- cv.glmnet(x = attr.mat , y = Surv(time, outcome), family = "cox")

#     best.lambda <- cv.fit$lambda.min

#     fit <- glmnet(x = attr.mat, y = Surv(time, outcome), family = "cox", lambda = best.lambda)

#     coefficients <- coef(fit, s = 0.5)

#     # Convert to a data frame (ensure it's a data frame)
#     coefficients_df <- as.data.frame(as.matrix(coefficients))

#     # Add row names as a column
#     coefficients_df$Feature <- row.names(coefficients_df)

#     # Sort the features by the absolute values of their coefficients in descending order
#     sorted_coefficients <- coefficients_df %>%
#     mutate(AbsoluteCoefficient = abs(coefficients_df[,1])) %>%
#     arrange(desc(AbsoluteCoefficient))


#     rownames(sorted_coefficients) <- sub("X", "", rownames(sorted_coefficients))
#     # Display top ranking features
#     sorted_coefficients$Feature <- NULL

#     glm_results[[i]] <- head(sorted_coefficients,num.sig)

# }

# # source in model_val function
# source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")

# results <- model_eval(feature.ranking = glm_results, important.features = c(1,num.sig))

# print(results)

#*****************************************      DATASET 1     ***************************************************

#--- 100 covariates, 100 individuals, effect size = 30, 10 significant features, s=0.5)
# > print(results)
#   precision recall
# 1       0.8    0.8

#*****************************************      DATASET 2     ***************************************************




#*****************************************      DATASET 3     ***************************************************


#--------------------------------------------COX REGRESSION ON NPDR DESIGN MATRIX (ALPHA = CLUSTER) -------------------------------------------------------


# library(glmnet)
# library(survival)
# #install.packages('devtools')
# # for windows
# devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/sNPDR")

# # for windows
# file_list <- list.files(path = "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData", pattern = "*.csv", full.names = TRUE)


# # number of significant features
# num.sig <- 10

# # Initialize an empty list to store the results
# glm_results <- list()


# # for mac
# #file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/simulatedData/simdata_new.csv"


# for (i in seq_along(file_list)) {

#         # Read in gene expression data
#     survival.data <- read.table(file_list[i], sep = ",", header = T)
#     attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
#     outcome <- survival.data[, ncol(survival.data)]
#     time <- survival.data[, ncol(survival.data)-1]


#     # call npdr
#     npdr_res <- sNPDR::npdr(time, outcome, attr.mat,
#                             regression.type = "cox", attr.diff.type = "numeric-abs",
#                             nbd.method = "relieff", nbd.metric = "euclidean",
#                             knn = 50, msurf.sd.frac = 0.5,
#                             covars = "none", covar.diff.type = "numeric-abs",
#                             padj.method = "bonferroni", verbose = FALSE,
#                             use.glmnet = TRUE, glmnet.alpha = "cluster", glmnet.lower = 0,
#                             glmnet.lam = "lambda.1se",
#                             rm.attr.from.dist = c(), neighbor.sampling = "none",
#                             separate.hitmiss.nbds = FALSE,
#                             corr.attr.names = NULL,
#                             fast.reg = FALSE, fast.dist = FALSE,
#                             dopar.nn = FALSE, dopar.reg = FALSE,
#                             unique.dof = FALSE, external.dist=NULL, KM.weight=FALSE, KM.kernal.type="gaussian", KM.kernal.sigma=0.75)


#     y <- npdr_res[, -((ncol(npdr_res)-2:ncol(npdr_res)))]
#     colnames(y) <- c("time", "status")
#     x <- npdr_res[, ((ncol(npdr_res)-2:ncol(npdr_res)))]

#     mod <- glmnet(as.matrix(x), as.matrix(y), family = "cox")

#     cv.fit <- cv.glmnet(as.matrix(x), as.matrix(y), family = "cox")

#     coef <- coef(mod, s = -6.0)

#     # Get the coefficients and their names
#     coef.values <- coef@x

#     # Order the coefficients by their absolute values in decreasing order
#     indices <- order(abs(coef.values), decreasing = TRUE)

#     # Create a dataframe
#     coef_df <- data.frame(coef.values)
#     rownames(coef_df) <- indices


#     glm_results[[i]] <- head(coef_df,num.sig)

# }

# # source in model_val function

# source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")

# results <- model_eval(feature.ranking = glm_results, important.features = seq(1,num.sig))

# print(results)



