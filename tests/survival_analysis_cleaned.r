# Load libraries
library(ggplot2)
library(PRROC)
library(precrec)
library(gridExtra)
library(cowplot)
library(survival)
library(dplyr)
library(tibble)
library(glmnet)

# Load custom functions
source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/model_eval.r")
source("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/tests/auRC_data.r")
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/sNPDR")

# Load data files
file_list <- list.files(path = "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/simulatedData", pattern = "*.csv", full.names = TRUE)
num.sig <- 20



evaluate_model <- function(file, model_type) {
  survival.data <- read.table(file, sep = ",", header = TRUE)
  attr.mat <- survival.data[, -c(1, ncol(survival.data)-1, ncol(survival.data))]
  outcome <- survival.data[, ncol(survival.data)]
  time <- survival.data[, ncol(survival.data)-1]
  
    if (model_type == "sNPDR") {
        #browser()
        model <- sNPDR::npdr(time, outcome, attr.mat, regression.type = "binomial-surv", attr.diff.type = "numeric-abs",
                            nbd.method = "relieff", nbd.metric = "euclidean", knn = 50, msurf.sd.frac = 0.5, 
                            covars = "none", covar.diff.type = "numeric-abs", padj.method = "bonferroni", verbose = FALSE)
        model <- model[order(-abs(model$beta.Z.att)),]
        model$feature <- rownames(model)
        coef_df <- model[, c("feature", "beta.Z.att")]
  } else if (model_type == "AFT") {
    #browser()
    model <- survreg(Surv(time, outcome) ~ ., data = attr.mat, dist = "weibull", control = list(maxiter = 1000, rel.tolerance = 1e-2))
    coef_df <- as.data.frame(summary(model)$coefficients)
    coef_df$beta.Z.att <- coef_df[,1]
    coef_df <- coef_df[order(-coef_df$beta.Z.att),]
    rownames(coef_df) <- sub("X", "", rownames(coef_df))
    coef_df <- subset(coef_df, rownames(coef_df) != "(Intercept)")

  } else if (model_type == "Cox") {
    cv.fit <- cv.glmnet(x = as.matrix(attr.mat), y = Surv(time, outcome), family = "cox")
    fit <- glmnet(x = as.matrix(attr.mat), y = Surv(time, outcome), family = "cox", lambda = cv.fit$lambda.min)
    coef_df <- as.data.frame(as.matrix(coef(fit)))
    coef_df$beta.Z.att <- abs(coef_df[,1])
    coef_df <- coef_df[order(-coef_df$beta.Z.att),]
    rownames(coef_df) <- sub("X", "", rownames(coef_df))
    coef_df <- subset(coef_df, rownames(coef_df) != "(Intercept)")
  } else if (model_type == "cox_npdr") {
    #browser()
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
  
  return(coef_df)
}

calculate_metrics <- function(results, threshold_range) {
  #browser()
  results_list <- lapply(threshold_range, function(threshold_value) {
    lapply(results, function(df) {
      #browser()
      df$predicted_class <- ifelse(abs(df$beta.Z.att) > threshold_value, 1, 0)
      df$actual_class <- ifelse(rownames(df) %in% seq(1, num.sig), 1, 0)
      return(df)
    })
  })
  return(results_list)
}


#model_types <- c("Cox")
model_types <- c("cox_npdr", "Cox")
all_results <- list()
pr_curves <- list()
#pr_dfs <- list()

for (model_type in model_types) {
  #browser()
  #browser()
  model_results <- lapply(file_list, evaluate_model, model_type = model_type)
  all_results[[model_type]] <- model_results # list of dataframes containing feature scores
  
  # determine the threshold range based on range of the beta values
  beta.range.ul <- max(sapply(model_results, function(df) max(abs(df$beta.Z.att))))
  beta.range.ll <- min(sapply(model_results, function(df) min(abs(df$beta.Z.att))))
  # choose the threshold range top 25% of the range
  threshold.range.ll <- beta.range.ll + (beta.range.ul - beta.range.ll) * 0.10
  threshold.range.ul <- beta.range.ul * 0.80
  print(c("Threshold range: ", threshold.range.ll, threshold.range.ul))

  threshold_range <- seq(threshold.range.ll, threshold.range.ul, 0.01)


  results_list <- calculate_metrics(model_results, threshold_range)

  evaluation <- lapply(results_list, function(df) {
  model_eval2(feature.ranking = df, important.features = seq(1, num.sig), num.sig = num.sig)
   })

  print(paste("Results for", model_type, ":"))
  print(evaluation)
  
  # Plot PR curve
  pr_df <- calculate_auRC(model_results, threshold_range)
  scores <- pr_df$score
  labels <- pr_df$label
  pr_curve <- pr.curve(scores.class0 = scores[labels == 1], scores.class1 = scores[labels == 0], curve = TRUE)


  #pr_dfs[[model_type]] <- pr_df
  pr_curves[[model_type]] <- pr_curve
  
  
}


# Initialize an empty ggplot object
plot <- ggplot() +
  labs(x = "Recall", y = "Precision") +
  ggtitle("Precision-Recall Curves") +
  scale_color_discrete(name = "Curves")

# Create an empty data frame to store all curves
all_curves <- data.frame()

# Iterate over pr_curves and add each one to the data frame
for (model_type in names(pr_curves)) {
  pr_curve <- pr_curves[[model_type]]
  data <- data.frame(Recall = pr_curve$curve[,1], Precision = pr_curve$curve[,2], Model = model_type)
  all_curves <- rbind(all_curves, data)
}
# Add the curves to the plot
plot <- plot + geom_line(data = all_curves, aes(x = Recall, y = Precision, color = Model), size = 1.5) +
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
  guides(color = guide_legend(override.aes = list(size = 1.5)))

# Print the plot
print(plot)


# plot the PR curve from the evaluation list

# append all the PR curves to a single dataframe
#eval.df <- do.call(rbind, evaluation)
eval.df <- data

plot <- ggplot() +
  geom_line(data = eval.df, aes(x = Recall, y = Precision), size = 1.5) +
  labs(x = "Recall", y = "Precision") +
  ggtitle("Precision-Recall Curves") +
  theme(
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold"),
    legend.key.size = unit(1.5, "lines"),
    panel.grid.major = element_line(size = 0.8),
    panel.grid.minor = element_line(size = 0.5)
  )

print(plot)
