# Install necessary packages if not already installed
install.packages(c("ggplot2", "PRROC", "precrec", "gridExtra", "cowplot"))

# Load packages
library(ggplot2)
library(PRROC)
library(precrec)
library(gridExtra)
library(cowplot)

# Generate some example data
set.seed(123)
scores <- runif(100)
labels <- sample(0:1, 100, replace = TRUE)

# Compute PR curve using PRROC
pr <- pr.curve(scores.class0 = scores[labels == 1], scores.class1 = scores[labels == 0], curve = TRUE)

# label columns of pr$curve
colnames(pr$curve) <- c("V1", "V2", "V3")

# Plot PR curve
pr_plot <- ggplot(data.frame(pr$curve), aes(x = V1, y = V2)) +
  geom_line() +
  labs(x = "Recall", y = "Precision") +
  ggtitle("Precision-Recall Curve")

# Generate example boxplot data
data <- data.frame(
  Method = rep(c("Random Forest", "Relief", "NPDR"), each = 100),
  auPRC = c(runif(100, 0.4, 0.6), runif(100, 0.5, 0.7), runif(100, 0.6, 0.8))
)

# Plot boxplot
box_plot <- ggplot(data, aes(x = Method, y = auPRC, color = Method)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(y = "Area under the PR Curve") +
  ggtitle("auPRC for Different Methods")

# Combine plots
combined_plot <- plot_grid(pr_plot, box_plot, labels = c("A", "B"))

# Display combined plot
print(combined_plot)
