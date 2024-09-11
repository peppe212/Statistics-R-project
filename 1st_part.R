# Author: Giuseppe Muschetta
# File: 1st_part.R

# INTRODUCTION TO PART 1:
# In this section, we explore the limitations of existing power transformations,
# specifically the Box-Cox and Signed-Power transformations. 
# The goal is to transform a random variable X such that the transformed
# distribution approximates normality. 
# We will minimize the Kullback-Leibler information number to select the optimal 
# transformation parameter λ for the signed power transformation. 
# Additionally, we will demonstrate that despite optimization, the transformed
# density may still exhibit bimodality, highlighting the need for new 
# transformation families.

# We will use the Signed Power transformation applied to the mixture density 
# f(x) = 0.3Φ(x) + 0.7γ(x) 
# The Box-Cox transformation is limited to positive values, while the signed 
# power transformation can handle the entire real line.

# Loading file script "functions.R" containing all the functions we need:
# box_cox(x, λ), signed_power(x, λ) and mixture_density(x)
file_path <- "functions.R"
if (file.exists(file_path)) {
  source(file_path)
} else {
  stop("File has not been found...", file_path)
}

# Load required libraries
library(ggplot2)
library(stats)

# Set seed for reproducibility
set.seed(123)

run_part_one <- function() {
  
  # Generate x values and calculate the corresponding mixture density
  x_vals <- seq(-10, 10, length.out = 3000)
  mixture_vals <- mixture_density(x_vals)
  
  # Sample 10,000 values from x_vals with replacement, using probabilities 
  # from mixture_vals. Values with higher mixture density are more likely 
  # to be sampled, simulating data from the mixture distribution.
  sampled_x <- sample(x_vals, size = 10000, replace = TRUE, prob = mixture_vals)
  
  # Estimate the density of the sampled data (sampled_x) over the range from -10 to 10, 
  # producing a smooth approximation of the data distribution.
  sampled_density <- density(sampled_x, from = -10, to = 10)
  
  
  # Function to compute the KL divergence
  compute_kl_divergence <- function(lambda) {
    
    # Apply the Signed Power transformation to the sampled data
    transformed_x <- signed_power(sampled_x, lambda)
    
    # Calculate the sample mean and empirical standard deviation 
    # for the transformed data
    mean_ <- mean(transformed_x)
    standard_deviation_ <- sd(transformed_x)
    
    # Generate the normal density based on the parameters of the transformed data
    normal_density <- dnorm(sampled_x, mean = mean_, sd = standard_deviation_)
    
    # Filter out non-positive values that might cause NaN in log calculation
    valid_indices <- which(transformed_x > 0 & normal_density > 0)
    
    # Calculate the Kullback-Leibler divergence on valid indices in order to avoid warnings
    kl_divergence <- sum(log(transformed_x[valid_indices] / normal_density[valid_indices]))
    
    # Return the KL divergence, mean, and standard deviation
    return(list(kl_divergence = kl_divergence, mean = mean_, sd = standard_deviation_))
  }
  
  
  # Optimize lambda to minimize KL divergence using optim
  # "L-BFGS-B" is an optimization algorithm that handles bounded parameter 
  # constraints, suitable for optimization over limited intervals.
  starting_lambda_value = 0.5
  optim_result <- optim(
                        par = starting_lambda_value, 
                        fn = function(lambda) compute_kl_divergence(lambda)$kl_divergence, 
                        method = "L-BFGS-B", lower = 0.01, upper = 2.0
                        )
  opt_lambda <- optim_result$par
  
  # Recalculate the optimal mean and standard deviation 
  # for the transformed distribution
  optimal_params <- compute_kl_divergence(opt_lambda)
  opt_mean <- optimal_params$mean
  opt_sd <- optimal_params$sd
  
  # Print the optimal lambda value and the corresponding mean, 
  # standard deviation, and variance for the transformed 
  # Signed Power density distribution.
  cat(sprintf("Best lambda = %f\n", opt_lambda))
  cat("For the transformed Signed Power density we have:\n")
  cat(sprintf("mean = %f\nstandard_deviation = %f\nvariance = %f\n", 
              opt_mean, opt_sd, opt_sd * opt_sd))
  

  # Recalculate the transformed density for plotting
  transformed_x_optimal <- signed_power(sampled_x, opt_lambda)
  opt_transformed_density <- density(transformed_x_optimal, from = -10, to = 10)
  
  # Generate the plots to visually compare the original, transformed, and normal densities
  plot <- ggplot() +
    geom_line(aes(x = x_vals, y = mixture_vals), color = "black", linetype = "solid", linewidth = 1.1) +
    geom_line(aes(x = opt_transformed_density$x, y = opt_transformed_density$y), color = "black", linetype = "dotted", linewidth = 1) +
    geom_line(aes(x = x_vals, y = dnorm(x_vals, mean = opt_mean, sd = opt_sd)), color = "black", linetype = "dashed", linewidth = 1) +
    xlim(-10, 10) + ylim(0, 0.3) +
    labs(title = "Optimized Density Transformations", x = "X", y = "Density") +
    theme_minimal()
  
  print(plot)
}

# Run the complete Part One process
run_part_one()


# CONCLUSION:
# This outcome underscores the inherent limitations of the signed power 
# transformation in handling data that are asymmetrically distributed across 
# the entire real line. 
# The observed bimodality, particularly pronounced under conditions of 
# extreme skewness, clearly illustrates the need for the development of a new 
# family of power transformations.
