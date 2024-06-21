## load libraries required for analysis
library(MLeval)
library(caret)

## simulate data
im <- twoClassSim(2000, intercept = -25, linearVars = 20)
table(im$Class)

## run caret
fitControl <- trainControl(method = "cv",summaryFunction=prSummary,
classProbs=T,savePredictions = T,verboseIter = F)
im_fit <- train(Class ~ ., data = im,method = "ranger",metric = "AUC",
trControl = fitControl)
im_fit2 <- train(Class ~ ., data = im,method = "xgbTree",metric = "AUC",
trControl = fitControl)

## run MLeval
x <- evalm(list(im_fit,im_fit2))

## curves and metrics are in the 'x' list object

# Load the ggplot2 package
library(ggplot2)

# Extract the PR curves from the 'x' list object
pr_curves <- x$PRC

# Initialize an empty ggplot object
plot <- ggplot() +
  labs(x = "Recall", y = "Precision") +
  ggtitle("Precision-Recall Curves") +
  scale_color_discrete(name = "Models")

# Iterate over pr_curves and add each one to the plot
for (model_name in names(pr_curves)) {
  pr_curve <- pr_curves[[model_name]]
  plot <- plot + geom_line(data = pr_curve, aes(x = Recall, y = Precision, color = model_name))
}

# Print the plot
print(plot)