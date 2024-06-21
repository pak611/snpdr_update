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
#-------------------------------------- LINEAR REGRESSION ON DESIGN MATRIX  ----------------------------------
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

library(glmnet)
library(survival)
#install.packages('devtools')
# for windows
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/sNPDR")

# for windows
file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData/data_1.csv"


# Read in gene expression data (The important features are (3,5) by a factor of 10)
survival.data <- read.table(file_path, sep = ",", header = T)
attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
outcome <- survival.data[, ncol(survival.data)]
time <- survival.data[, ncol(survival.data)-1]

library(survivalsvm)



# call npdr
npdr_res <- sNPDR::npdr(time, outcome, attr.mat,
                        regression.type = "cox", attr.diff.type = "numeric-abs",
                        nbd.method = "relieff", nbd.metric = "euclidean",
                        knn = 200, msurf.sd.frac = 0.5,
                        covars = "none", covar.diff.type = "numeric-abs",
                        padj.method = "bonferroni", verbose = FALSE,
                        use.glmnet = TRUE, glmnet.alpha = 'cluster', glmnet.lower = 0,
                        glmnet.lam = "lambda.1se",
                        rm.attr.from.dist = c(), neighbor.sampling = "none",
                        separate.hitmiss.nbds = FALSE,
                        corr.attr.names = NULL,
                        fast.reg = FALSE, fast.dist = FALSE,
                        dopar.nn = FALSE, dopar.reg = FALSE,
                        unique.dof = FALSE, external.dist=NULL, KM.weight=FALSE, KM.kernal.type="gaussian", KM.kernal.sigma=2.0)


y <- npdr_res[, -((ncol(npdr_res)-2:ncol(npdr_res)))]
colnames(y) <- c("time", "status")
x <- npdr_res[, ((ncol(npdr_res)-2:ncol(npdr_res)))]


# Fit a linear regression model to the data
mod <- lm(y$time ~ ., data = x)

# Get the coefficients and their names
coef <- coef(mod)

# Rank the coefficients by their absolute values
sorted_coef <- sort(abs(coef), decreasing = TRUE)

# Rank the coefficients by their value
sorted_coef2 <- sort(coef, decreasing = TRUE)

print(sorted_coef)

#> print(sorted_coef)
#          X1 (Intercept)         X10          X2          X6         X41
# 92.36142679 48.60347066 14.23731124 11.62141987  9.85661924  9.05748589
#         X49          X3         X18         X48         X24         X31
#  8.34568368  7.74924526  6.17118503  6.02187568  5.54452101  5.46812186
#         X42         X15         X27          X5         X23         X21
#  5.30712581  5.18464182  5.14767216  4.91018759  4.60932351  4.52297756
#         X47         X37         X45          X7         X32         X44
#  4.47342563  4.21447786  3.82694923  3.68802626  3.60247361  3.38933340
#         X36         X20         X29         X39         X40         X14
#  2.97182086  2.80419014  2.74605348  2.70671170  2.40040093  2.33324065
#         X16         X25         X38         X11         X13         X12
#  2.21925763  2.12200446  2.02199507  1.87383319  1.76655285  1.64946383
#         X17         X26          X8         X34         X22         X35
#  1.31322198  1.25992994  1.22448916  1.11320388  0.97203553  0.86995288
#         X19         X28         X43         X33         X30          X9
#  0.77334799  0.53041301  0.48631958  0.40580693  0.34756103  0.22636219
#         X50         X46          X4
#  0.19644087  0.09627586  0.08120409


print(sorted_coef2)
# now try on just the data itself
mod <- lm(time ~ ., data = attr.mat)


# Get the coefficients and their names
coef <- coef(mod)

# Rank the coefficients by their absolute values
sorted_coef <- sort(abs(coef), decreasing = TRUE)

print(sorted_coef)

#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#-----------------------------------------------------------------  SURVIVAL SVM  ----------------------------------
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# #for mac
# #file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/simulatedData/simdata_new.csv"


# # for windows
# file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData/data_1.csv"


# # Read in gene expression data (The important features are (3,5) by a factor of 10)
# survival.data <- read.table(file_path, sep = ",", header = T)
# attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
# outcome <- survival.data[, ncol(survival.data)]
# time <- survival.data[, ncol(survival.data)-1]

# library(survivalsvm)

# fit <- survivalsvm(
#   formula = Surv(time, outcome) ~ .,
#   data = attr.mat,
#   type = "regression",
#   gamma.mu = 0.1,
#   opt.meth = "ipop",
#   kernel = "lin_kernel"
# )

# str(fit)


# support_vectors <- fit$model.fit$SV  
# coefficients <- fit$model.fit$Beta


# # Compute the weight vector (w)
# weights <- t(support_vectors) %*% coefficients

# # Assign names to the weights
# names(weights) <- fit$var.names

# # Rank the coefficients by their absolute values 
# sorted_weights <- sort(abs(weights), decreasing = TRUE)

# print(sorted_weights)



# #/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# #-----------------------------------------------------------------  NPDR + SURVIVAL SVM  ----------------------------------
# #/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

library(glmnet)
library(survival)
#install.packages('devtools')
# for windows
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/sNPDR")

# for windows
file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData/data_1.csv"


# Read in gene expression data (The important features are (3,5) by a factor of 10)
survival.data <- read.table(file_path, sep = ",", header = T)
attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
outcome <- survival.data[, ncol(survival.data)]
time <- survival.data[, ncol(survival.data)-1]

library(survivalsvm)



# call npdr
npdr_res <- sNPDR::npdr(time, outcome, attr.mat,
                        regression.type = "cox", attr.diff.type = "numeric-abs",
                        nbd.method = "relieff", nbd.metric = "euclidean",
                        knn = 50, msurf.sd.frac = 0.5,
                        covars = "none", covar.diff.type = "numeric-abs",
                        padj.method = "bonferroni", verbose = FALSE,
                        use.glmnet = TRUE, glmnet.alpha = 'cluster', glmnet.lower = 0,
                        glmnet.lam = "lambda.1se",
                        rm.attr.from.dist = c(), neighbor.sampling = "none",
                        separate.hitmiss.nbds = FALSE,
                        corr.attr.names = NULL,
                        fast.reg = FALSE, fast.dist = FALSE,
                        dopar.nn = FALSE, dopar.reg = FALSE,
                        unique.dof = FALSE, external.dist=NULL, KM.weight=FALSE, KM.kernal.type="gaussian", KM.kernal.sigma=2.0)


y <- npdr_res[, -((ncol(npdr_res)-2:ncol(npdr_res)))]
colnames(y) <- c("time", "status")
x <- npdr_res[, ((ncol(npdr_res)-2:ncol(npdr_res)))]



fit <- survivalsvm(
  formula = Surv(y$time, y$status) ~ .,
  data = x,
  type = "regression",
  gamma.mu = 0.1,
  opt.meth = "ipop",
  kernel = "lin_kernel"
)

str(fit)


support_vectors <- fit$model.fit$SV  
coefficients <- fit$model.fit$Beta


# Compute the weight vector (w)
weights <- t(support_vectors) %*% coefficients

# Assign names to the weights
names(weights) <- fit$var.names

# Rank the coefficients by their absolute values 
sorted_weights <- sort(abs(weights), decreasing = TRUE)

print(sorted_weights)





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


# #-------------------------------------------Genotype Data---------------------------------------------

# file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/sReliefF_data/EpiSet2.txt"
# # Read in the genotype data

# survival.data <- read.table(file_path, sep = "\t", header = FALSE, skip = 1, quote = "")

# # Read the first line to get the column names as a single string
# headers <- readLines(file_path, n = 1)

# # Split the header string into individual column names using strsplit
# column_names <- strsplit(headers, "\t", fixed = TRUE)[[1]]

# colnames(survival.data) <- column_names

# attr.mat <- survival.data[, -c(1,2)]
# outcome <- as.numeric(survival.data[, 1])
# time <- survival.data[, 2]


# #-------------------------------------------Gene-Expression Data---------------------------------------------

# #for mac
# #file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/simulatedData/simdata_new.csv"


# # for windows
# file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData/data_1.csv"


# # Read in gene expression data (The important features are (3,5) by a factor of 10)
# survival.data <- read.table(file_path, sep = ",", header = T)
# attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
# outcome <- survival.data[, ncol(survival.data)]
# time <- survival.data[, ncol(survival.data)-1]


#------------------------------------------- Real Data --------------------------------------------

library(Survival)


# Extract the expression data 
attr.mat <- as.matrix(exprs(eset))


# Extract the phenotype data
pheno_data <- pData(eset)



time <- as.numeric(trimws(pheno_data[["Decease.delay.after.surgery..months."]]))
max_time <- max(time, na.rm = TRUE)
time[is.na(time)] <- max_time
outcome <- ifelse(pheno_data[["State.of.health"]] == "deceased", 0, 1)
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



# #-------------------------------------------Gene-Expression Data---------------------------------------------
# # > print(head(ranked.vars, 10))
# #        pval.att beta.raw.att beta.Z.att    beta.0 pval.0
# # 74 2.330272e-10    0.3216577   6.341351 -2.411518      0
# # 56 2.758856e-10    0.3314397   6.315237 -2.415650      0
# # 43 4.224716e-08   -0.3278016  -5.483489 -2.043583      0
# # 73 1.074449e-07   -0.2883822  -5.315758 -2.052037      0
# # 48 2.074789e-07   -0.2846058  -5.194469 -2.054083      0
# # 46 2.492382e-07   -0.2933707  -5.160196 -2.055562      0
# # 78 5.084550e-07   -0.2819113  -5.024873 -2.058599      0
# # 82 1.503624e-06   -0.2836475  -4.812321 -2.062809      0
# # 91 1.719346e-06   -0.2690923  -4.785433 -2.069016      0
# # 75 4.726613e-06   -0.2445048  -4.577926 -2.072387      0
# #-------------------------------------------Gene-Expression Data (revised) ---------------------------------------------
# # > print(head(ranked.vars, 10))
# #        pval.att beta.raw.att beta.Z.att    beta.0 pval.0
# # 87 9.091136e-17    0.3099250   8.323973 -1.568922      0
# # 48 7.556039e-14   -0.3006592  -7.483577 -1.201872      0
# # 19 9.523886e-12    0.2435935   6.817884 -1.528022      0
# # 49 9.429869e-11   -0.2702944  -6.479574 -1.221949      0
# # 96 1.141282e-10   -0.2597223  -6.450649 -1.227314      0
# # 17 6.646412e-10   -0.2436142  -6.177545 -1.236406      0
# # 61 1.540816e-09    0.2364560   6.043091 -1.517799      0
# # 56 2.867202e-09    0.2394930   5.941935 -1.517156      0
# # 89 1.081049e-08    0.2194685   5.720097 -1.511783      0
# # 64 2.242962e-08   -0.2222094  -5.594570 -1.245819      0



#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#-------------------------------------- FILTER + BINOMIAL REGRESSION ON SNP SURVIVAL DATA  ----------------------------------
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# library(glmnet)
# library(survival)
# #install.packages('devtools')
# # for windows
#devtools::load_all("C:/Users/patri/Dropbox/Ph.D/Research/Doctoral_Research/Bioinformatics/sNPDR/sNPDR")

# # for mac
# devtools::load_all("/Users/patrickkampmeyer/Desktop/sNPDR/sNPDR")


# # read in SNP data
# # for mac
# #file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/sReliefF_data/Episet2.csv"
# # for windows
# file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/sReliefF_data/EpiSet2.txt"

# # Another dataset
# file_path <- "C:/Users/patri/OneDrive/Desktop/Episet_Datasets/Episet_Dataset_0.txt"
# survival.data <- read.table(file_path, sep = "\t", header = FALSE, skip = 1, quote = "")
# #survival.data <- read.table(file_path, sep = ",", header = T)
# attr.mat <- survival.data[, -c(1, 2)]
# outcome <- survival.data[, ncol(survival.data)]
# time <- survival.data[, ncol(survival.data)-1]

# # call npdr
# npdr_res <- sNPDR::npdr(time, outcome, attr.mat,
#                         regression.type = "binomial-surv", attr.diff.type = "numeric-abs",
#                         nbd.method = "relieff", nbd.metric = "euclidean",
#                         knn = 50, msurf.sd.frac = 0.5,
#                         covars = "none", covar.diff.type = "numeric-abs",
#                         padj.method = "bonferroni", verbose = FALSE,
#                         use.glmnet = FALSE, glmnet.alpha = 0.1, glmnet.lower = 0,
#                         glmnet.lam = "lambda.1se",
#                         rm.attr.from.dist = c(), neighbor.sampling = "none",
#                         separate.hitmiss.nbds = FALSE,
#                         corr.attr.names = NULL,
#                         fast.reg = FALSE, fast.dist = FALSE,
#                         dopar.nn = FALSE, dopar.reg = FALSE,
#                         unique.dof = FALSE, external.dist=NULL, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=2.0)


# # sort by abs(beta.Z.att) in descending order, keeping the original row names
# ranked.vars <- npdr_res[order(-abs(npdr_res$beta.Z.att)),]
# print(head(ranked.vars, 5))

# # Successfull in identifying 14 and 15 as the top two features
# #> print(head(ranked.vars, 5))
# #       pval.att beta.raw.att  beta.Z.att     beta.0        pval.0
# #15 0.000000e+00  -2.06990975 -114.489792 -0.0261107  1.304219e-02
# #14 0.000000e+00  -1.20126717  -80.018725 -0.3718214 1.336942e-285
# #13 5.091485e-08   0.07724988    5.448593 -1.0371470  0.000000e+00
# #6  5.176381e-07   0.06979277    5.020056 -1.0334348  0.000000e+00
# #11 5.784996e-06   0.06006434    4.534396 -1.0349130  0.000000e+00

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




# Now for the generated episet

# snpdr:

# > print(head(ranked.vars, 5))
#        pval.att beta.raw.att  beta.Z.att    beta.0 pval.0
# 20 0.000000e+00   -8.6685057 -123.662188 5.0856604      0
# 8  2.007303e-65    0.3159400   17.098564 0.9882829      0
# 16 4.095403e-23    0.1642019    9.904801 1.0298208      0
# 14 1.671757e-21    0.1393316    9.526488 1.0032593      0
# 18 8.005679e-21    0.1319419    9.362257 1.0049787      0


# sReliefF:

# Rank	Variable		Score
# 1	Var14		0.0715
# 2	Var19		0.0528
# 3	Var18		0.0394
# 4	Var5		0.0285
# 5	Var3		0.0242
# 6	Var15		0.0084
# 7	Var7		-0.008
# 8	Var2		-0.0552
# 9	Var8		-0.0725
# 10	Var1		-0.0738
# 11	Var0		-0.0889
# 12	Var13		-0.1248
# 13	Var17		-0.159
# 14	Var4		-0.1599
# 15	Var9		-0.168
# 16	Var16		-0.1733
# 17	Var6		-0.1857
# 18	Var12		-0.1906
# 19	Var11		-0.1916
# 20	Var10		-0.2199







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
file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData/data_2.csv"

# for mac
#file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/simulatedData/simdata_new.csv"


#file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/sReliefF_data/EpiSet2.txt"
# Read in the genotype data

# Another episet dataset
#file_path <- "C:/Users/patri/OneDrive/Desktop/Episet_Datasets/Episet_Dataset_0.txt"

# file_path <- "C:/Users/patri/OneDrive/Desktop/Episet_Datasets/dataset_2.txt"

# survival.data <- read.table(file_path, sep = "\t", header = FALSE, skip = 1, quote = "")

# # Read the first line to get the column names as a single string
# headers <- readLines(file_path, n = 1)

# # Split the header string into individual column names using strsplit
# column_names <- strsplit(headers, "\t", fixed = TRUE)[[1]]

# colnames(survival.data) <- column_names

# attr.mat <- survival.data[, -c(1,2)]
# outcome <- as.numeric(survival.data[, 1])
# time <- survival.data[, 2]

# # Read in gene expression data
survival.data <- read.table(file_path, sep = ",", header = T)
attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
outcome <- survival.data[, ncol(survival.data)]
time <- survival.data[, ncol(survival.data)-1]

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
                        unique.dof = FALSE, external.dist=NULL, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=5.0)


y <- npdr_res[, -((ncol(npdr_res)-2:ncol(npdr_res)))]
colnames(y) <- c("time", "status")
x <- npdr_res[, ((ncol(npdr_res)-2:ncol(npdr_res)))]

fraction_0 <- rep(1 - sum(y$status == 0) / nrow(y), sum(y$status == 0)) * 1.0
fraction_1 <- rep(1 - sum(y$status == 1) / nrow(y), sum(y$status == 1))

weights <- numeric(nrow(y))

weights[y$status == 0] <- fraction_0
weights[y$status == 1] <- fraction_1
 

fit <- survreg(Surv(y$time, y$status) ~ ., data = x, dist = "weibull", weights = weights)

fit_summary <- summary(fit)

coef_df <- as.data.frame(fit_summary$coefficients)

coef_df$Feature <- rownames(coef_df)

names(coef_df) <- c("z", "Feature")

coef_df <- coef_df[-1, ]

coef_df_sorted <- coef_df[order(-abs(coef_df$z)), ]

coef_df_sorted$Feature <- sub("X", "", coef_df_sorted$Feature)

rownames(coef_df_sorted) <- coef_df_sorted$Feature

# Display the sorted data frame
print(coef_df_sorted)


predict(cv.fit, s = cv.fit$lambda[length(cv.fit$lambda)], type = "coefficients", lower.limits = 0.0)

#cv.fit <- cv.glmnet(as.matrix(x), as.matrix(y), family = "cox")

#plot(cv.fit)

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


# # survreg(formula = Surv(y$time, y$status) ~ ., data = x, dist = "weibull")
# #                Value Std. Error     z       p
# # (Intercept)  3.08365    0.11654 26.46 < 2e-16
# # Var19       -0.00826    0.05429 -0.15 0.87913
# # Var18        0.24412    0.05493  4.44 8.8e-06
# # Var17       -0.06519    0.04436 -1.47 0.14164
# # Var16       -0.10575    0.04518 -2.34 0.01926
# # Var15       -0.19491    0.04896 -3.98 6.9e-05
# # Var14        0.12918    0.05856  2.21 0.02740
# # Var13        0.02216    0.04488  0.49 0.62151
# # Var12       -0.05732    0.04490 -1.28 0.20168
# # Var11       -0.01999    0.04408 -0.45 0.65025
# # Var10       -0.09671    0.04433 -2.18 0.02916
# # Var9        -0.04204    0.04420 -0.95 0.34158
# # Var8        -0.17143    0.04643 -3.69 0.00022
# # Var7        -0.08589    0.05318 -1.62 0.10631
# # Var6        -0.08105    0.04469 -1.81 0.06975
# # Var5        -0.08000    0.05305 -1.51 0.13158
# # Var4        -0.04861    0.04497 -1.08 0.27973
# # Var3        -0.20293    0.05119 -3.96 7.4e-05
# # Var2        -0.06814    0.04701 -1.45 0.14727
# # Var1        -0.01367    0.04602 -0.30 0.76634
# # Var0        -0.04378    0.04552 -0.96 0.33618
# # Log(scale)   0.71345    0.01028 69.38 < 2e-16


# # NPDR no time reclassification
# #                   z Feature
# # Var3   0.0443860896    Var3
# # Var19 -0.0293031376   Var19
# # Var2  -0.0251185272    Var2
# # Var13  0.0234811157   Var13
# # Var7   0.0230309151    Var7
# # Var0   0.0223468739    Var0
# # Var8   0.0211175533    Var8
# # Var6   0.0187145412    Var6
# # Var15  0.0173721612   Var15
# # Var4   0.0148467616    Var4
# # Var18 -0.0116934200   Var18
# # Var14  0.0104353251   Var14
# # Var11  0.0097397908   Var11
# # Var9  -0.0096306075    Var9
# # Var17  0.0094155964   Var17
# # Var5   0.0088854346    Var5
# # Var1   0.0060362063    Var1
# # Var12 -0.0038051500   Var12
# # Var10 -0.0023267694   Var10
# # Var16 -0.0003794072   Var16


# # NPDR with time reclassification
# #                 z Feature
# # Var2  -0.17801854    Var2
# # Var9   0.16499869    Var9
# # Var1   0.08240912    Var1
# # Var11 -0.07989843   Var11
# # Var13 -0.07895185   Var13
# # Var19 -0.07416082   Var19
# # Var10 -0.07164858   Var10
# # Var16 -0.07068879   Var16
# # Var17  0.06973748   Var17
# # Var0   0.06694037    Var0
# # Var6  -0.06549161    Var6
# # Var8  -0.06175918    Var8
# # Var14  0.05499545   Var14
# # Var18 -0.04949380   Var18
# # Var15 -0.04241567   Var15
# # Var12 -0.04029568   Var12
# # Var4  -0.03414764    Var4
# # Var3   0.03280781    Var3
# # Var5  -0.02769813    Var5
# # Var7   0.01989020    Var7

#////////////////////////////////////////////////////////////////////////////////////////////////////////////
#---------------------------------------------- GLMNET COX MODEL --------------------------------------------
#////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Load the necessary libraries
library(glmnet)
library(survival)
library(dplyr)

# # episet data
# file_path <- "C:/Users/patri/OneDrive/Desktop/Episet_Datasets/dataset_2.txt"

# survival.data <- read.table(file_path, sep = "\t", header = FALSE, skip = 1, quote = "")

# # Read the first line to get the column names as a single string
# headers <- readLines(file_path, n = 1)

# # Split the header string into individual column names using strsplit
# column_names <- strsplit(headers, "\t", fixed = TRUE)[[1]]

# colnames(survival.data) <- column_names

# attr.mat <- survival.data[, -c(1,2)]
# outcome <- as.numeric(survival.data[, 1])
# time <- survival.data[, 2]


# for mac
file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData/data_1.csv"

# Read in gene expression data
survival.data <- read.table(file_path, sep = ",", header = T)
attr.mat <- as.matrix(survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))])
outcome <- survival.data[, ncol(survival.data)]
time <- survival.data[, ncol(survival.data)-1]

cv.fit <- cv.glmnet(x = as.matrix(attr.mat) , y = Surv(time, outcome), family = "cox")

best.lambda <- cv.fit$lambda.min

fit <- glmnet(x = attr.mat, y = Surv(time, outcome), family = "cox", lambda = best.lambda)

coefficients <- coef(fit, s = 0.5)

# Convert to a data frame (ensure it's a data frame)
coefficients_df <- as.data.frame(as.matrix(coefficients))

# Add row names as a column
coefficients_df$Feature <- row.names(coefficients_df)

# Sort the features by the absolute values of their coefficients in descending order
sorted_coefficients <- coefficients_df %>%
  mutate(AbsoluteCoefficient = abs(coefficients_df[,1])) %>%
  arrange(desc(AbsoluteCoefficient))

# Display top ranking features
print(head(sorted_coefficients,20))

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



# # # Episet Dataset
# # > print(head(sorted_coefficients,10))
# #                  1 Feature AbsoluteCoefficient
# # Var19 -0.241236112   Var19         0.241236112
# # Var9   0.188016672    Var9         0.188016672
# # Var16 -0.149724180   Var16         0.149724180
# # Var2  -0.108448621    Var2         0.108448621
# # Var18 -0.093409392   Var18         0.093409392
# # Var7   0.089584818    Var7         0.089584818
# # Var0   0.079081356    Var0         0.079081356
# # Var17  0.052090170   Var17         0.052090170
# # Var5   0.020770258    Var5         0.020770258
# # Var11  0.005377796   Var11         0.005377796



#////////////////////////////////////////////////////////////////////////////////////////////////////////////
#-------------------------------------- ACCELERATED FAILURE TIME MODEL (AFT) --------------------------------
#////////////////////////////////////////////////////////////////////////////////////////////////////////////

library(survival)
# # for mac
# file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/simulatedData/simdata_new.csv"

# for windows
file_path <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData/data_1.csv"


survival.data <- read.table(file_path, sep = ",", header = T)
attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
outcome <- survival.data[, ncol(survival.data)]
time <- survival.data[, ncol(survival.data)-1]

aft.model <- survreg(Surv(time,outcome) ~ ., data = attr.mat, dist = "weibull")

model.summary <- summary(aft.model)
print(model.summary)

# Get the coefficients from the model summary
coef <- model.summary$coefficients

# Convert the coefficients to a data frame
coef_df <- as.data.frame(coef) %>% 
  tibble::rownames_to_column(var = "feature")

# Add a column with the absolute values of the coefficients
coef_df <- coef_df %>% mutate(abs_coef = abs(coef))

# Rank the features based on the absolute values of the coefficients
coef_df <- coef_df %>% arrange(desc(abs_coef))


coef_df$feature <- sub("X", "", coef_df$feature)

rownames(coef_df) <- coef_df$feature

coef_df$feature <- NULL

# Remove the row with the name "(Intercept)"
coef_df <- subset(coef_df, rownames(coef_df) != "(Intercept)")



# Print the ranked features
print(head(coef_df, 10))

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

# # # Episet Dataset
# # > print(head(coef_df, 10))
# #              coef   abs_coef
# # Var19  0.28029085 0.28029085
# # Var9  -0.25083413 0.25083413
# # Var2   0.19068399 0.19068399
# # Var16  0.18143042 0.18143042
# # Var18  0.16313533 0.16313533
# # Var7  -0.15229937 0.15229937
# # Var17 -0.12522466 0.12522466
# # Var0  -0.11956210 0.11956210
# # Var11 -0.09941089 0.09941089
# # Var6  -0.07931023 0.07931023


#////////////////////////////////////////////////////////////////////////////////////////////////////////////

#//////////////////////////////////////////////////////////////////////////////////////////////////////
#----------------------------------------- MODEL EVALUATION ------------------------------------------
#//////////////////////////////////////////////////////////////////////////////////////////////////////


#------------------------------------------------- SNPDR ------------------------------------------------
#Calculate the True Positive Rate (TPR) and False Positive Rate (FPR) for the top 10% of genes selected by sNPDR

#Get a list of all CSV files in the simulatedData directory

#load the sNPDR package

library(ggplot2)
library(PRROC)
library(precrec)
library(gridExtra)
library(cowplot)

#for mac
#devtools::load_all("/Users/patrickkampmeyer/Desktop/sNPDR/sNPDR")

# for windows
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/sNPDR")

# for mac
#file_list <- list.files(path = "C:/Users/patri/Dropbox/Ph.D/Research/Doctoral_Research/Bioinformatics/sNPDR/data/simulatedData", pattern = "*.csv", full.names = TRUE)

# for windows
file_list <- list.files(path = "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData", pattern = "*.csv", full.names = TRUE)


# number of significant features
num.sig <- 20

# Initialize an empty list to store the results
npdr_results <- list()

# Loop over the files
for(i in seq_along(file_list)) {


    survival.data <- read.table(file_list[i], sep = ",", header = T)
    attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
    outcome <- survival.data[, ncol(survival.data)]
    time <- survival.data[, ncol(survival.data)-1]


    npdr_res <- sNPDR::npdr_revised(time, outcome, attr.mat,
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

    # Sort npdr_res by beta.Z.att in descending order
    npdr_res <- npdr_res[order(-abs(npdr_res$beta.Z.att)),]


    # Append the result to the list
    #npdr_results[[i]] <- head(npdr_res,num.sig)
    npdr_results[[i]] <- npdr_res
}

print(head(npdr_results))


# source in model_val function
source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")
source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/auRC_data.r")


pr_df <- calculate_auRC(npdr_results, seq(3.0, 5.0, 0.1))

scores <- pr_df$score
labels <- pr_df$label

# Compute PR curve using PRROC
pr1 <- pr.curve(scores.class0 = scores[labels == 1], scores.class1 = scores[labels == 0], curve = TRUE)

# label columns of pr$curve
colnames(pr1$curve) <- c("V1", "V2", "V3")

# Plot PR curve
pr_plot <- ggplot(data.frame(pr1$curve), aes(x = V1, y = V2)) +
  geom_line(size = 1.5) +
    theme(
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold"),
    legend.key.size = unit(1.5, "lines"),
    panel.grid.major = element_line(size = 0.8),
    panel.grid.minor = element_line(size = 0.5)
  ) +
  labs(x = "Recall", y = "Precision") +
  ggtitle("sNPDR: Precision-Recall Curve")

print(pr_plot)


results_df <- npdr_results
threshold_range <- seq(3.0, 6.0, 0.1)

npdr_results_list <- lapply(threshold_range, function(threshold_value) {
    results <- lapply(results_df, function(df) {
        #browser()
        df$predicted_class <- ifelse(abs(df$beta.Z.att) > threshold_value, 1, 0)
        df$actual_class <- ifelse(rownames(df) %in% seq(1, num.sig), 1, 0)
        return(df)
    })

    results_list <- results
    return(results_list)
})



# Apply model_eval2 to each dataframe in npdr_results_class
results_list <- lapply(npdr_results_list, function(df) {
  model_eval2(feature.ranking = df, important.features = seq(1,num.sig), num.sig = num.sig)
})

# Combine the results into a single dataframe
results <- do.call(rbind, results_list)

# Combine the results into a single dataframe
results <- do.call(rbind, results_list)

print(results)

# > print(results)
#  precision     recall        mcc
# 0.17000000 0.53941755 0.07773311

# *****************************************      DATASET 1      ***************************************************

# --- 100 covariates, 100 individuals, effect size = 30, 10 significant features, KM.weight=FALSE)
# snpdr
# > print(results)
#   precision recall
# 1      0.18   0.18
# npdr_revised
# > print(results)
#   precision recall
# 1      0.11   0.11

# --- 100 covariates, 100 individuals, effect size = 30, 10 significant features, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=0.75)
# snpdr
# > print(results)
#   precision recall
# 1       0.2    0.2
# npdr_revised
# > print(results)
#   precision recall
# 1      0.12   0.12

# --- 100 covariates, 100 individuals, effect size = 30, 10 significant features, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=0.25)
# > print(results)
#   precision recall
# 1       0.1    0.1

# --- 100 covariates, 100 individuals, effect size = 30, 10 significant features, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=2.0)
# > print(results)
#   precision recall
# 1      0.19   0.19

# *****************************************      DATASET 2     ***************************************************

# REDUCE THE NUMBER OF SIGNIFICANT FEATURES

# --- 100 covariates, 100 individuals, effect size = 30, 5 significant features, KM.weight=FALSE)
# > print(results)
#   precision recall
# 1      0.16   0.16

# --- 100 covariates, 100 individuals, effect size = 30, 5 significant features, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=0.75)
# > print(results)
#   precision recall
# 1       0.2    0.2

# *****************************************      DATASET 3     ***************************************************

# INCREASE THE NUMBER OF INSTANCES AND FEATURES

# --- 1000 covariates, 500 individuals, effect size = 30, 10 significant features, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=0.75)
# > print(results)
#   precision recall
# 1      0.06   0.06

#--------------------------------------------------------- ACCELERATED FAILURE TIME MODEL (AFT) ------------------------------------------------

library(survival)
library(dplyr)

source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")

# for windows
file_list <- list.files(path = "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData", pattern = "*.csv", full.names = TRUE)

# number of significant features
num.sig <- 20

# Initialize an empty list to store the results
aft_results <- list()


for (i in seq_along(file_list)) {
    
    print(i)

    survival.data <- read.table(file_list[i], sep = ",", header = T)
    attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
    outcome <- survival.data[, ncol(survival.data)]
    time <- survival.data[, ncol(survival.data)-1]

    model <- survreg(Surv(time, outcome) ~ ., data = attr.mat, dist = "weibull", control = list(maxiter = 1000, rel.tolerance = 1e-2))
    coef_df <- as.data.frame(summary(model)$coefficients)
    coef_df$beta.Z.att <- coef_df[,1]
    coef_df <- coef_df[order(-coef_df$beta.Z.att),]
    rownames(coef_df) <- sub("X", "", rownames(coef_df))
    coef_df <- subset(coef_df, rownames(coef_df) != "(Intercept)")
    aft_results[[i]] <- coef_df
}


# source in model_val function

source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")

results <- model_eval(feature.ranking = aft_results, important.features = seq(1,num.sig), num.sig = num.sig)

print(results)

# > print(results)
# precision    recall       mcc
# 1.0000000 0.9890110 0.9482093


pr_df <- calculate_auRC(aft_results, seq(0, 0.5, 0.01))

scores <- pr_df$score
labels <- pr_df$label

# Compute PR curve using PRROC
pr2 <- pr.curve(scores.class0 = scores[labels == 1], scores.class1 = scores[labels == 0], curve = TRUE)

# label columns of pr$curve
colnames(pr2$curve) <- c("V1", "V2", "V3")

# plot the PR curve
pr_plot <- ggplot(data.frame(pr2$curve), aes(x = V1, y = V2)) +
  geom_line(size = 1.5) +
    theme(
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold"),
    legend.key.size = unit(1.5, "lines"),
    panel.grid.major = element_line(size = 0.8),
    panel.grid.minor = element_line(size = 0.5)
  ) +
  labs(x = "Recall", y = "Precision") +
  ggtitle("AFT: Precision-Recall Curve")

print(pr_plot)




results <- model_eval2(feature.ranking = npdr_results_class, important.features = seq(1,num.sig), num.sig = num.sig)

print(results)

#*****************************************      DATASET 1     ***************************************************
#--- 100 covariates, 100 individuals, effect size = 30, 10 significant features)

# > print(results)
#   precision recall
# 1      0.62   0.62
#*****************************************      DATASET 2     ***************************************************




#*****************************************      DATASET 3     ***************************************************



# # #--------------------------------------------------------- GLMNET COX MODEL ------------------------------------------------

# Load the necessary libraries
library(glmnet)
library(survival)
library(dplyr)
library(ggplot2)
library(PRROC)
library(precrec)
library(gridExtra)
library(cowplot)

# for windows
file_list <- list.files(path = "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData", pattern = "*.csv", full.names = TRUE)


# number of significant features
num.sig <- 20

# Initialize an empty list to store the results
glm_results <- list()



for (i in seq_along(file_list)){

    # Read in gene expression data
    survival.data <- read.table(file_list[i], sep = ",", header = T)
    attr.mat <- as.matrix(survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))])
    outcome <- survival.data[, ncol(survival.data)]
    time <- survival.data[, ncol(survival.data)-1]

    

    cv.fit <- cv.glmnet(x = attr.mat , y = Surv(time, outcome), family = "cox")

    best.lambda <- cv.fit$lambda.min

    fit <- glmnet(x = attr.mat, y = Surv(time, outcome), family = "cox", lambda = best.lambda)

    coefficients <- coef(fit, s = 0.5)

    # Convert to a data frame (ensure it's a data frame)
    coefficients_df <- as.data.frame(as.matrix(coefficients))

    # Add row names as a column
    coefficients_df$Feature <- row.names(coefficients_df)

    # Sort the features by the absolute values of their coefficients in descending order
    sorted_coefficients <- coefficients_df %>%
    mutate(AbsoluteCoefficient = abs(coefficients_df[,1])) %>%
    arrange(desc(AbsoluteCoefficient))


    rownames(sorted_coefficients) <- sub("X", "", rownames(sorted_coefficients))
    # Display top ranking features
    sorted_coefficients$Feature <- NULL

    glm_results[[i]] <- sorted_coefficients

}

# change column name from AbsoluteCoefficient to beta.Z.att
glm_results <- lapply(glm_results, function(df) {
  names(df)[names(df) == "AbsoluteCoefficient"] <- "beta.Z.att"
  return(df)
})

# source in the calculate_auRC function
source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/auRC_data.r")
source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")



pr_df <- calculate_auRC(glm_results, seq(0, 0.1, 0.01))

scores <- as.numeric(pr_df$score)
labels <- as.numeric(pr_df$label)


# Compute PR curve using PRROC
pr3 <- pr.curve(scores.class0 = scores[labels == 1], scores.class1 = scores[labels == 0], curve = TRUE)

# label columns of pr$curve
colnames(pr3$curve) <- c("V1", "V2", "V3")

# plot the PR curve
pr_plot <- ggplot(data.frame(pr3$curve), aes(x = V1, y = V2)) +
  geom_line(size = 1.5) +
    theme(
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold"),
    legend.key.size = unit(1.5, "lines"),
    panel.grid.major = element_line(size = 0.8),
    panel.grid.minor = element_line(size = 0.5)
  ) +
  labs(x = "Recall", y = "Precision") +
  ggtitle("COX: Precision-Recall Curve")

print(pr_plot)



results <- model_eval(feature.ranking = glm_results, important.features = seq(1,num.sig), num.sig = num.sig)

print(results)

# > print(results)
# precision    recall       mcc
# 1.0000000 0.9890110 0.9482093


#*****************************************      DATASET 1     ***************************************************

#--- 100 covariates, 100 individuals, effect size = 30, 10 significant features, s=0.5)
# > print(results)
#   precision recall
# 1       0.8    0.8

#*****************************************      DATASET 2     ***************************************************




#*****************************************      DATASET 3     ***************************************************


#--------------------------------------------COX REGRESSION ON NPDR DESIGN MATRIX (ALPHA = CLUSTER) -------------------------------------------------------


library(glmnet)
library(survival)
#install.packages('devtools')
# for windows
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/sNPDR")

# for windows
file_list <- list.files(path = "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData", pattern = "*.csv", full.names = TRUE)



# number of significant features
num.sig <- 20

# Initialize an empty list to store the results
glm_results <- list()


# for mac
#file_path <- "/Users/patrickkampmeyer/Desktop/sNPDR/data/simulatedData/simdata_new.csv"


for (i in seq_along(file_list)) {

        # Read in gene expression data
    survival.data <- read.table(file_list[i], sep = ",", header = T)
    attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
    outcome <- survival.data[, ncol(survival.data)]
    time <- survival.data[, ncol(survival.data)-1]


    # call npdr
    npdr_res <- sNPDR::npdr_revised(time, outcome, attr.mat,
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
                            unique.dof = FALSE, external.dist=NULL, KM.weight=TRUE, KM.kernal.type="gaussian", KM.kernal.sigma=2.0)


    y <- npdr_res[, -((ncol(npdr_res)-2:ncol(npdr_res)))]
    colnames(y) <- c("time", "status")
    x <- npdr_res[, ((ncol(npdr_res)-2:ncol(npdr_res)))]


    fraction_0 <- (rep(1 - sum(y$status == 0) / nrow(y), sum(y$status == 0)) * 0.001) + 0.01
    fraction_1 <- rep(1 - sum(y$status == 1) / nrow(y), sum(y$status == 1))

    weights <- numeric(nrow(y))

    weights[y$status == 0] <- fraction_0
    weights[y$status == 1] <- fraction_1


    fit <- survreg(Surv(y$time, y$status) ~ ., data = x, dist = "weibull", weights = weights)

    fit_summary <- summary(fit)

    coef_df <- as.data.frame(fit_summary$coefficients)

    coef_df$Feature <- rownames(coef_df)

    names(coef_df) <- c("z", "Feature")

    coef_df <- coef_df[-1, ]

    coef_df_sorted <- coef_df[order(-abs(coef_df$z)), ]

    coef_df_sorted$Feature <- sub("X", "", coef_df_sorted$Feature)

    rownames(coef_df_sorted) <- coef_df_sorted$Feature

    glm_results[[i]] <- coef_df_sorted

}

# change column name from AbsoluteCoefficient to beta.Z.att
glm_results <- lapply(glm_results, function(df) {
  names(df)[names(df) == "z"] <- "beta.Z.att"
  return(df)
})

# source in the calculate_auRC function
source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/auRC_data.r")
source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")



pr_df <- calculate_auRC(glm_results, seq(0, 1.0, 0.01))

scores <- as.numeric(pr_df$score)
labels <- as.numeric(pr_df$label)


# Compute PR curve using PRROC
pr4 <- pr.curve(scores.class0 = scores[labels == 1], scores.class1 = scores[labels == 0], curve = TRUE)

# label columns of pr$curve
colnames(pr4$curve) <- c("V1", "V2", "V3")

# plot the PR curve
pr_plot <- ggplot(data.frame(pr4$curve), aes(x = V1, y = V2)) +
  geom_line(size = 1.5) +
    theme(
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold"),
    legend.key.size = unit(1.5, "lines"),
    panel.grid.major = element_line(size = 0.8),
    panel.grid.minor = element_line(size = 0.5)
  ) +
  labs(x = "Recall", y = "Precision") +
  ggtitle("COX-NPDR: Precision-Recall Curve")

print(pr_plot)

# source in model_val function

source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")

results <- model_eval(feature.ranking = glm_results, important.features = seq(1,num.sig), num.sig = num.sig)

print(results)

# npdr Kernal sigma = 2.0
# > print(results) 
# precision    recall       mcc
# 0.2700000 0.6771057 0.1849167

# npdr no kernel
# > print(results)
# precision    recall       mcc
# 0.2400000 0.6904654 0.1535939

# npdr_revised Kernal sigma = 2.0
# > print(results)
# precision    recall       mcc
# 0.3100000 0.7728411 0.2287121

# npdr_revised no kernel
# > print(results) 
# precision    recall       mcc
# 0.2600000 0.6947746 0.1762496


# original censor weight

# > print(results)
# precision    recall       mcc
# 0.2300000 0.6506244 0.1390051

# cencor weight * 0.1
# > print(results)
# precision    recall       mcc
# 0.3600000 0.8199871 0.2848068

# cencor weight * 0.01
# > print(results)
# precision    recall       mcc
# 0.3800000 0.8239834 0.3036571

# cencor weight * 0.001
# > print(results)
# precision    recall       mcc
# 0.3800000 0.8211997 0.3027413

# Assuming 'pr1', 'pr2', 'pr3', 'pr4' are your PR curve data frames
pr_plot <- ggplot() +
  geom_line(data = data.frame(pr1$curve), aes(x = V1, y = V2, colour = "pr1")) +
  geom_line(data = data.frame(pr2$curve), aes(x = V1, y = V2, colour = "pr2")) +
  geom_line(data = data.frame(pr3$curve), aes(x = V1, y = V2, colour = "pr3")) +
  geom_line(data = data.frame(pr4$curve), aes(x = V1, y = V2, colour = "pr4")) +
  scale_colour_manual(values = c("pr1" = "red", "pr2" = "blue", "pr3" = "green", "pr4" = "purple"),
                      name = "Legend Title",
                      labels = c("SNPDR", "AFT", "COX", "COX-NPDR")) +
  theme(
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold"),
    legend.key.size = unit(1.5, "lines"),
    panel.grid.major = element_line(size = 0.8),
    panel.grid.minor = element_line(size = 0.5)
  ) +
  labs(x = "Recall", y = "Precision") +
  ggtitle("sNPDR")  # This is the title you want

print(pr_plot)
