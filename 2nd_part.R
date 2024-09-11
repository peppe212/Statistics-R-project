# Author: Giuseppe Muschetta
# File: 2nd_part.R


# INTRODUCTION TO PART 2
# Modified Modulus Transformation:
# A modified modulus transformation is introduced with different parameters for 
# positive and negative values of x. It is defined as:
#  ((x + 1)^λ+ -1)/λ+    if x >= 0 and λ+ != 0
#  log(x + 1)            if x >= 0 and λ+ = 0,
# -(-(x + 1)^λ- -1)/λ-   if x < 0  and λ- != 0
# -log(-x + 1)           if x < 0  and λ- = 0
# 
# If we assume that the second derivative is continuous at x=0, which we will 
# remember as the point where the change of sign occurs, and carry out the 
# calculations, we obtain the following relationship: λ- + λ+ = 2, 
# ergo  λ- = 2 - λ+
# In this way we can rewrite the transformation function so that it depends 
# on a single parameter.
#
# The ultimate YEO-JOHNSON POWER TRANSFORMATION IS:
#  ((x + 1)^λ -1)/λ                 if x >= 0 and λ != 0
#  log(x + 1)                       if x >= 0 and λ = 0,
# -(-(x + 1)^(2 - λ) -1)/(2 - λ)    if x < 0  and λ != 2
# -log(-x + 1)                      if x < 0  and λ = 2
#
# Application and Impact:
# These transformations can potentially handle dataset with both positive and 
# negative entries more effectively, accommodating the inherent asymmetry in such data. 
# The adaptability of the transformation parameters based on the sign of x 
# values is particularly advantageous for real-world data, which often exhibits
# mixed-sign skewness.



# First, we need the file script "functions.R" containing all the functions:
# box_cox(x, λ), yeo_johnson(x, λ), yeo_johnson_derivative(x, λ)
# test_continuity(...), plot_continuity(...), test_increasing(...)
file_path <- "functions.R"
if (file.exists(file_path)) {
  source(file_path)
} else {
  stop("File has not been found:", file_path)
}

# Libraries needed in this file script:
library(ggplot2)
library(gridExtra)
library(stats)

# Set the seed for reproducibility
set.seed(123)


#########################
# CONCAVITY AND CONVEXITY
# This part of the script demonstrates that the transformation is 
# concave for λ < 1 and convex for λ > 1.
concavity_convexity <- function() {
  data_pos <- 0.7 * rnorm(500, mean = 2, sd = 1)  
  # Generate negatively skewed normal data
  data_neg <- 0.3 * rnorm(500, mean = -2, sd = 1) 
  # Combine the positive and negative data into a single vector
  data <- c(data_pos, data_neg)  
  
  # Define lambda values for transformation demonstration
  lambda_values <- c(0, 0.5, 1, 1.5, 2)
  
  # Apply the Yeo-Johnson transformation across all data points
  yeojohnson_data <- sapply(lambda_values, function(l) yeo_johnson(data, l))
  
  # Create a dataframe for plotting
  yeojohnson_plot_data <- expand.grid(X = data, Lambda = lambda_values)
  yeojohnson_plot_data$Transformed <- as.vector(yeojohnson_data)
  
  # Set line type adjustments as a named vector for plot aesthetics
  linetypes <- setNames(rep(c("solid", "dashed", "dotted", "twodash", "longdash"), length.out = length(lambda_values)), lambda_values)
  
  # Create the plot with appropriate styles using ggplot2
  yeojohnson_plot <- ggplot(yeojohnson_plot_data, aes(x = X, y = Transformed, color = factor(Lambda), linetype = factor(Lambda))) +
    geom_line() +
    scale_color_manual(values = rainbow(length(lambda_values))) +
    scale_linetype_manual(values = linetypes) +
    labs(title = "Yeo-Johnson Transformations on Mixed Normal Data", x = "Original Values", y = "Transformed Values") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  
  # Display the plot
  print(yeojohnson_plot)
}

### Execute concavity_convexity function
concavity_convexity()



######################################
# VIEW ONLY YEO-JOHNSON TRANSFORMATION 
view_yeo_johnson <- function() {
  
  # Define lambda values for transformation demonstration (matching those in the figure)
  lambda_values <- c(0, 0.5, 1, 1.5, 2)
  
  # Generate data in a range suitable for the plot, between -4 and 4
  # We're contemplating a symmetrical interval
  data <- seq(-4, 4, length.out = 1000)  
  
  # Apply the Yeo-Johnson transformation across all data points for each lambda
  yeojohnson_data <- sapply(lambda_values, function(l) yeo_johnson(data, l))
  
  # Create a dataframe for plotting
  yeojohnson_plot_data <- expand.grid(X = data, Lambda = lambda_values)
  yeojohnson_plot_data$Transformed <- as.vector(yeojohnson_data)
  
  # Filter the data to ensure all transformed values lie within the range [-4, 4]
  yeojohnson_plot_data <- subset(yeojohnson_plot_data, Transformed >= -4 & Transformed <= 4)
  
  # Set line type adjustments as a named vector for plot aesthetics
  # Set lambda = 1 as solid, and all other values as dashed or dotted
  linetypes <- setNames(c("dotted", "dashed", "solid", "dashed", "dotted"), lambda_values)
  
  # Set color adjustments: all lines black but lambda = 1 stands out as a solid black line
  yeojohnson_plot <- ggplot(yeojohnson_plot_data, aes(x = X, y = Transformed, color = factor(Lambda), linetype = factor(Lambda))) +
    geom_line(linewidth = 1) +  # Set all lines to the same width
    scale_color_manual(values = rep("black", length(lambda_values))) +  # All lines in black
    scale_linetype_manual(values = linetypes) +  # Set the correct linetypes
    labs(title = "Yeo-Johnson Transformations", x = "Original values", y = "") +
    xlim(-4, 4) + ylim(-4, 4) +  # Set axis limits to match the figure
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right", 
          legend.title = element_blank())  # Remove the legend title for consistency with the figure
  
  # Display the plot
  print(yeojohnson_plot)
}

### Execute view_yeo_johnson function
view_yeo_johnson()



###################################
# REPRODUCING FIGURE 2 OF THE PAPER
# BOX-COX and YEO_JOHNSON transformations
# This part of the script demonstrates a comparative analysis of the 
# Box-Cox vs Yeo-Johnson power transformations.
# The Box-Cox transformation is applied only to positive data, whereas the 
# Yeo-Johnson transformation is applicable to both positive and negative values, 
# making it more versatile for real-world data.
# The transformations are explored over a range of x values from -4 to 6, 
# with specified lambda values of 0, 0.5, 1, 1.5, and 2. 
# The lambda parameter dictates the nature of the transformation:
# - A concave transformation (lambda < 1) addresses right skewness.
# - A convex transformation (lambda > 1) addresses left skewness.
box_cox_vs_yeo_johnson <- function() {
  # Set the seed for reproducibility
  set.seed(123)
  
  # Define the x values range
  x_vals <- seq(-4, 6, by = 0.1)
  # Positive values for Box-Cox
  x_vals_positive <- x_vals[x_vals > 0]  
  
  # Define lambda values for transformation
  lambda_values <- c(0, 0.5, 1, 1.5, 2)
  
  # Calculate Box-Cox transformation only on positive values
  boxcox_data <- sapply(lambda_values, function(l) box_cox(x_vals_positive, l))
  
  # Calculate Yeo-Johnson transformation for all x_vals
  yeojohnson_data <- sapply(lambda_values, function(l) yeo_johnson(x_vals, l))
  
  # Create a dataframe for ggplot for Box-Cox
  boxcox_plot_data <- expand.grid(X = x_vals_positive, Lambda = lambda_values)
  boxcox_plot_data$Transformed <- as.vector(boxcox_data)
  
  # Create a dataframe for Yeo-Johnson
  yeojohnson_plot_data <- expand.grid(X = x_vals, Lambda = lambda_values)
  yeojohnson_plot_data$Transformed <- as.vector(yeojohnson_data)
  
  # Set linetype adjustments as named vector
  linetypes <- setNames(c("solid", "longdash", "dashed", "dotted", "twodash"), lambda_values)
  
  # Plotting both transformations with appropriate styles
  boxcox_plot <- ggplot(boxcox_plot_data, aes(x = X, y = Transformed, linetype = factor(Lambda))) +
    geom_line(color = "black") +
    labs(title = "Box-Cox Transformations", x = "Original Values", y = "Transformed Values") +
    scale_linetype_manual(values = linetypes) +
    ylim(-4, 4) + 
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  
  yeojohnson_plot <- ggplot(yeojohnson_plot_data, aes(x = X, y = Transformed, linetype = factor(Lambda))) +
    geom_line(color = "black") +
    labs(title = "Yeo-Johnson Transformations", x = "Original Values", y = "Transformed Values") +
    scale_linetype_manual(values = linetypes) +
    ylim(-4, 4) + # set the same limit as the paper
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  
  grid.arrange(boxcox_plot, yeojohnson_plot, ncol = 2) 
  #plot(boxcox_plot)
  #plot(yeojohnson_plot)
}

### Executes box_cox_vs_yeo_johnson
box_cox_vs_yeo_johnson()



###################################
# REPRODUCING FIGURE 3 of the paper
# APPLYING YEO-JOHNSON TRANSFORMATION to the mixture density
# Compute the Kullback-Leibler divergence between the transformed and normal densities
compute_kl_divergence <- function(lambda, sampled_x) {
  transformed_x <- yeo_johnson(sampled_x, lambda)
  
  # Calculate normal density for the transformed data
  mean_transformed <- mean(transformed_x)
  sd_transformed <- sd(transformed_x)
  normal_density <- dnorm(transformed_x, mean = mean_transformed, sd = sd_transformed)
  
  # Calculate mixture density after transformation
  transformed_density <- density(transformed_x, n = length(transformed_x))$y
  
  # Ensure densities are positive and avoid NaN values
  valid_indices <- which(transformed_density > 0 & normal_density > 0)
  
  # Compute KL divergence on valid indices
  kl_divergence <- sum(transformed_density[valid_indices] * log(transformed_density[valid_indices] / normal_density[valid_indices]))
  
  return(kl_divergence)
}

# Function to apply Yeo-Johnson transformation to the mixture density and minimize KL divergence
apply_yeo_johnson_transformation <- function() {
  # Generate x values and calculate the mixture density
  x_vals <- seq(-5, 10, length.out = 3000)
  mixture_vals <- mixture_density(x_vals)
  
  # Sample from the mixture density
  sampled_x <- sample(x_vals, size = 10000, replace = TRUE, prob = mixture_vals)
  
  # Optimize lambda to minimize KL divergence using optim
  optim_result <- optim(par = 0.3, fn = function(lambda) compute_kl_divergence(lambda, sampled_x), 
                        method = "L-BFGS-B", lower = 0.01, upper = 2.0)
  opt_lambda <- optim_result$par
  
  # Recalculate the optimal mean and standard deviation for the transformed distribution
  transformed_x_optimal <- yeo_johnson(sampled_x, opt_lambda)
  opt_mean <- mean(transformed_x_optimal)
  opt_sd <- sd(transformed_x_optimal)
  
  cat(sprintf("Best lambda = %f\n", opt_lambda))
  cat(sprintf("Mean = %f, Standard Deviation = %f, Variance = %f\n", opt_mean, opt_sd, opt_sd^2))
  
  # Plot original mixture, transformed, and normal densities
  transformed_density <- density(transformed_x_optimal, from = -5, to = 10)
  
  plot <- ggplot() +
    geom_line(aes(x = x_vals, y = mixture_vals), color = "black", linetype = "solid", linewidth = 1) +
    geom_line(aes(x = transformed_density$x, y = transformed_density$y), color = "black", linetype = "dotted", linewidth = 1) +
    geom_line(aes(x = x_vals, y = dnorm(x_vals, mean = opt_mean, sd = opt_sd)), color = "black", linetype = "dashed", linewidth = 1) +
    xlim(-5, 10) + ylim(0, 0.4) +
    labs(x = "X", y = "Density") +
    theme_minimal()
  
  print(plot)
}

# Run the function to apply the transformation
apply_yeo_johnson_transformation()


# CONCLUSION OF THE SECOND PART
# The graphs illustrate the transformation effects clearly, showing how the 
# Yeo-Johnson transformation provides a more normal-like distribution compared 
# to the bimodal outcome of the signed power transformation. 
# This is a significant finding as it demonstrates the Yeo-Johnson 
# transformation's capability to handle skewed data effectively, 
# thus supporting the hypothesis that it could be better suited for 
# distributions with pronounced skewness.
