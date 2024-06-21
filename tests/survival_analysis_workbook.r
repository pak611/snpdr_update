#//////////////////////////////////////////////////////////////////////////////////////////////////////
#----------------------------------------- MODEL EVALUATION ------------------------------------------
#//////////////////////////////////////////////////////////////////////////////////////////////////////


# source in model_val function
source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")
source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/auRC_data.r")

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
num.sig <- 10

# Initialize an empty list to store the results
npdr_results <- list()

# Loop over the files
for(i in seq_along(file_list)) {


    survival.data <- read.table(file_list[i], sep = ",", header = T)
    attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
    outcome <- survival.data[, ncol(survival.data)]
    time <- survival.data[, ncol(survival.data)-1]


    npdr_res <- sNPDR::npdr(time, outcome, attr.mat,
                            regression.type = "keep-same", attr.diff.type = "numeric-abs",
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
    npdr_results[[i]] <- npdr_res
}


print(head(npdr_results))


pr1_df <- calculate_auRC(npdr_results, seq(0, 9, 0.01))

scores <- pr1_df$score
labels <- pr1_df$label

# Compute PR curve using PRROC
pr1 <- pr.curve(scores.class0 = scores[labels == 1], scores.class1 = scores[labels == 0], curve = TRUE)

# label columns of pr$curve
colnames(pr1$curve) <- c("V1", "V2", "V3")

# Plot PR curve
pr_plot <- ggplot(data.frame(pr1$curve), aes(x = V1, y = V2)) +
  geom_line(size = 1.5) +  # Increase line thickness
  labs(x = "Recall", y = "Precision") +
  ggtitle("Precision-Recall Curve") +
  theme(axis.text = element_text(size = 12),  # Increase axis label size
        axis.title = element_text(size = 14),  # Increase axis title size
        axis.ticks = element_line(size = 1.5))  # Increase axis tick size

print(pr_plot)


results_df <- npdr_results
threshold_range <- seq(0, 8, 0.01)

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


print(results_list)

# > print(results_list)
# [[1]]
#     precision    recall         mcc
# 1  0.33333333 0.5000000  0.13678242
# 2         NaN 0.0000000         NaN
# 3  0.05882353 0.5000000 -0.06211740
# 4  0.00000000 0.0000000 -0.07647191
# 5  0.00000000 0.0000000 -0.03350126
# 6  0.00000000 0.0000000 -0.04761905
# 7         NaN 0.0000000         NaN
# 8         NaN 0.0000000         NaN
# 9  0.00000000 0.0000000 -0.05862104
# 10 0.12765957 0.9310345  0.08682309


#--------------------------------------------------------- ACCELERATED FAILURE TIME MODEL (AFT) ------------------------------------------------


library(survival)
library(dplyr)
library(tibble)



# For Windows
file_list <- list.files(path = "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData", pattern = "*.csv", full.names = TRUE)

# Number of significant features
num.sig <- 10

# Initialize an empty list to store the results
aft_results <- list()

for (i in seq_along(file_list)) {
    #browser()
    print(i)
    
    survival.data <- read.table(file_list[i], sep = ",", header = TRUE)
    attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
    outcome <- survival.data[, ncol(survival.data)]
    time <- survival.data[, ncol(survival.data)-1]
    
    tryCatch({
        #browser()
        aft.model <- survreg(Surv(time, outcome) ~ ., data = attr.mat, dist = "weibull", control = list(maxiter = 1000, rel.tolerance = 1e-2))
        model.summary <- summary(aft.model)
        
        # Get the coefficients from the model summary
        coef <- model.summary$coefficients
        
        # Convert the coefficients to a data frame
        coef_df <- as.data.frame(coef) %>% 
            rownames_to_column(var = "feature") %>%
            mutate(abs_coef = abs(coef))
        
        # Rank the features based on the absolute values of the coefficients
        coef_df <- coef_df %>% arrange(desc(abs_coef))
        
        coef_df$feature <- sub("X", "", coef_df$feature)
        rownames(coef_df) <- coef_df$feature
        
        # Remove the row with the name "(Intercept)"
        coef_df <- subset(coef_df, rownames(coef_df) != "(Intercept)")
         # change column name of coef_df$z to beta.Z.att
        names(coef_df)[2] <- "beta.Z.att"
        coef_df$feature <- NULL
        
        # Append the result to the list
        aft_results[[i]] <- coef_df
    }, error = function(e) {
        message("Error fitting model for file ", file_list[i], ": ", e)
    })
}

# Evaluate the results
results <- model_eval(feature.ranking = aft_results, important.features = seq(1, num.sig), num.sig = num.sig)
print(results)

# > print(results)
# precision    recall       mcc
# 1.0000000 0.9890110 0.9482093


pr2_df <- calculate_auRC(aft_results, seq(0, 9, 0.1))

scores <- pr2_df$score
labels <- pr2_df$label

# Compute PR curve using PRROC
pr2 <- pr.curve(scores.class0 = scores[labels == 1], scores.class1 = scores[labels == 0], curve = TRUE)

# label columns of pr$curve
colnames(pr2$curve) <- c("V1", "V2", "V3")

# Plot PR curve
pr_plot <- ggplot(data.frame(pr2$curve), aes(x = V1, y = V2)) +
  geom_line() +
  labs(x = "Recall", y = "Precision") +
  ggtitle("Precision-Recall Curve")

print(pr_plot)



results <- model_eval2(feature.ranking = npdr_results_class, important.features = seq(1,num.sig), num.sig = num.sig)

print(results)


results_df <- aft_results
threshold_range <- seq(0.5, 0.5, 1)

aft_results_list <- lapply(threshold_range, function(threshold_value) {
    results <- lapply(results_df, function(df) {
        #browser()
        df$predicted_class <- ifelse(abs(df$abs_coef) > threshold_value, 1, 0)
        df$actual_class <- ifelse(rownames(df) %in% seq(1, num.sig), 1, 0)
        return(df)
    })

    results_list <- results
    return(results_list)
})



# Apply model_eval2 to each dataframe in npdr_results_class
results_list <- lapply(aft_results_list, function(df) {
  model_eval2(feature.ranking = df, important.features = seq(1,num.sig), num.sig = num.sig)
})

print(results_list)

# > print(results_list)
# [[1]]
#    precision    recall       mcc
# 1          1 1.0000000 1.0000000
# 2        NaN 0.0000000       NaN
# 3          1 1.0000000 1.0000000
# 4          1 1.0000000 1.0000000
# 5          1 0.9545455 0.8230549
# 6          1 1.0000000 1.0000000
# 7          1 0.9878049 0.9434564
# 8          1 0.8571429 0.6123724
# 9          1 0.9310345 0.7579367
# 10       NaN 0.0000000       NaN

#--------------------------------------------------------- COX REGRESSION ------------------------------------------------
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
num.sig <- 10

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



pr3_df <- calculate_auRC(glm_results, seq(0, 9, 0.1))

scores <- as.numeric(pr3_df$score)
labels <- as.numeric(pr3_df$label)


# Compute PR curve using PRROC
pr3 <- pr.curve(scores.class0 = scores[labels == 1], scores.class1 = scores[labels == 0], curve = TRUE)

# label columns of pr$curve
colnames(pr3$curve) <- c("V1", "V2", "V3")

# Plot PR curve
pr_plot <- ggplot(data.frame(pr3$curve), aes(x = V1, y = V2)) +
  geom_line() +
  labs(x = "Recall", y = "Precision") +
  ggtitle("Precision-Recall Curve")

print(pr_plot)



results_df <- glm_results
threshold_range <- seq(0.8, 1.0, 0.01)

glm_results_list <- lapply(threshold_range, function(threshold_value) {
    results <- lapply(results_df, function(df) {
        #browser()
        df$predicted_class <- ifelse(abs(df$beta.Z.att) > threshold_value, 1, 0)
        df$actual_class <- ifelse(rownames(df) %in% seq(1, num.sig), 1, 0)
        return(df)
    })

    results_list <- results
    return(results_list)
})

# Apply model_eval2 to each dataframe in glm_results_list
results_list <- lapply(glm_results_list, function(df) {
  model_eval2(feature.ranking = df, important.features = seq(1,num.sig), num.sig = num.sig)
})


print(results_list)

#    precision    recall       mcc
# 1        NaN 0.0000000       NaN
# 2          1 0.9545455 0.8230549
# 3        NaN 0.0000000       NaN
# 4        NaN 0.0000000       NaN
# 5          1 0.6923077 0.4285714
# 6          1 0.9310345 0.7579367
# 7          1 1.0000000 1.0000000
# 8          1 0.6923077 0.4285714
# 9          1 0.7941176 0.5275893
# 10         1 0.5000000 0.3015113


#--------------------------------------------COX REGRESSION ON NPDR DESIGN MATRIX (ALPHA = CLUSTER) -------------------------------------------------------


library(glmnet)
library(survival)
#install.packages('devtools')
# for windows
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/sNPDR")

# for windows
file_list <- list.files(path = "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData", pattern = "*.csv", full.names = TRUE)



# number of significant features
num.sig <- 10

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


    npdr_res <- sNPDR::npdr(time, outcome, attr.mat, regression.type = "cox", attr.diff.type = "numeric-abs",
                                    nbd.method = "relieff", nbd.metric = "euclidean", knn = 50, msurf.sd.frac = 0.5,
                                    covars = "none", covar.diff.type = "numeric-abs", padj.method = "bonferroni", verbose = FALSE,
                                    use.glmnet = TRUE, glmnet.alpha = "cluster", glmnet.lower = 0, glmnet.lam = "lambda.1se",
                                    rm.attr.from.dist = c(), neighbor.sampling = "none", separate.hitmiss.nbds = FALSE,
                                    corr.attr.names = NULL, fast.reg = FALSE, fast.dist = FALSE, dopar.nn = FALSE,
                                    dopar.reg = FALSE, unique.dof = FALSE, external.dist=NULL, KM.weight=TRUE,
                                    KM.kernal.type="gaussian", KM.kernal.sigma=2.0)
    
    y <- npdr_res[, ((ncol(npdr_res)-1):ncol(npdr_res))]
    colnames(y) <- c("time", "status")
    x <- npdr_res[, 1:(ncol(npdr_res)-2)]
    
    fraction_0 <- rep(1 - sum(y$status == 0) / nrow(y), sum(y$status == 0)) * 0.001
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
    coef_df <- coef_df[order(-abs(coef_df$z)), ]
    coef_df$Feature <- sub("X", "", coef_df$Feature)
    rownames(coef_df) <- coef_df$Feature
    coef_df <- subset(coef_df, rownames(coef_df) != "(Intercept)")
    # change column name of coef_df$z to beta.Z.att
    names(coef_df)[1] <- "beta.Z.att"

}

# source in model_val function

source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")

results <- model_eval(feature.ranking = glm_results, important.features = seq(1,num.sig), num.sig = num.sig)

print(results)


results_df <- glm_results
threshold_range <- seq(0.3, 0.3, 0.01)

glm_results_list <- lapply(threshold_range, function(threshold_value) {
    results <- lapply(results_df, function(df) {
        #browser()
        df$predicted_class <- ifelse(abs(df$z) > threshold_value, 1, 0)
        df$actual_class <- ifelse(rownames(df) %in% seq(1, num.sig), 1, 0)
        return(df)
    })

    results_list <- results
    return(results_list)
})

# Apply model_eval2 to each dataframe in glm_results_list
results_list <- lapply(glm_results_list, function(df) {
  model_eval2(feature.ranking = df, important.features = seq(1,num.sig), num.sig = num.sig)
})


print(results_list)

# > print(results_list)
# [[1]]
#    precision    recall         mcc
# 1  0.5714286 0.8571429  0.43112399
# 2  1.0000000 0.5000000  0.30151134
# 3  1.0000000 0.5000000  0.30151134
# 4  0.2500000 0.9000000  0.25000000
# 5  0.0000000 0.0000000 -0.03350126
# 6  0.4000000 0.6923077  0.22941573
# 7        NaN 0.0000000         NaN
# 8  0.4000000 0.6923077  0.22941573
# 9  0.6666667 0.6923077  0.33218588
# 10 0.5555556 0.9000000  0.47755198

# source in the calculate_auRC function
source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/auRC_data.r")
source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")



pr3_df <- calculate_auRC(glm_results, seq(0, 9, 0.1))

scores <- as.numeric(pr_df3$score)
labels <- as.numeric(pr_df3$label)


# Compute PR curve using PRROC
pr3 <- pr.curve(scores.class0 = scores[labels == 1], scores.class1 = scores[labels == 0], curve = TRUE)

# label columns of pr$curve
colnames(pr$curve) <- c("V1", "V2", "V3")

# Plot PR curve
pr_plot <- ggplot(data.frame(pr$curve), aes(x = V1, y = V2)) +
  geom_line() +
  labs(x = "Recall", y = "Precision") +
  ggtitle("Precision-Recall Curve")

print(pr_plot)



results_df <- glm_results
threshold_range <- seq(0.8, 1.0, 0.01)

glm_results_list <- lapply(threshold_range, function(threshold_value) {
    results <- lapply(results_df, function(df) {
        #browser()
        df$predicted_class <- ifelse(abs(df$beta.Z.att) > threshold_value, 1, 0)
        df$actual_class <- ifelse(rownames(df) %in% seq(1, num.sig), 1, 0)
        return(df)
    })

    results_list <- results
    return(results_list)
})

# Apply model_eval2 to each dataframe in glm_results_list
results_list <- lapply(glm_results_list, function(df) {
  model_eval2(feature.ranking = df, important.features = seq(1,num.sig), num.sig = num.sig)
})


print(results_list)



#-------------------------------------------------- CUMULATIVE ------------------------------------------------

# Assume pr1, pr2, pr3 are your PR curve data
# Combine all PR curves into one data frame with an additional 'curve' column
all_pr_data <- rbind(
  data.frame(Recall = pr1$curve[,1], Precision = pr1$curve[,2], Curve = 'Curve 1'),
  data.frame(Recall = pr2$curve[,1], Precision = pr2$curve[,2], Curve = 'Curve 2')
  #data.frame(Recall = pr3$curve[,1], Precision = pr3$curve[,2], Curve = 'Curve 3')
)

# Plot PR curves
pr_plot <- ggplot(all_pr_data, aes(x = Recall, y = Precision, color = Curve, group = Curve)) +
  geom_line() +
  labs(x = "Recall", y = "Precision") +
  ggtitle("Precision-Recall Curves") +
  scale_color_discrete(name = "Curves")

print(pr_plot)
