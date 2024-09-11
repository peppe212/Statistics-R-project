# Author: Giuseppe Muschetta
# File: 3rd_part.R


# INTRODUCTION TO PART 3:
# This script performs a comprehensive statistical analysis based on 
# Charles Darwin's (1876) study on cross-fertilization versus self-fertilization 
# in plant growth. Darwin examined 15 pairs of plants, with one plant cross-fertilized 
# and the other self-fertilized. The goal was to determine whether cross-fertilized 
# plants exhibit superior growth.

# This script carries out several key tasks:
# 
# 1. Initial Statistical Analysis: A one-tailed paired t-test evaluates if the 
#    mean height difference between cross-fertilized and self-fertilized plants 
#    is significantly greater than zero, confirming Darwin's hypothesis.
# 
# 2. Normality Check: A Q-Q plot, Shapiro-Wilk test, and skewness calculation 
#    are used to assess the normality of the data.
# 
# 3. Yeo-Johnson Transformation: To correct for non-normality, the script applies 
#    the Yeo-Johnson transformation, optimizing the transformation parameter λ 
#    by maximizing the log-likelihood function.
# 
# 4. Parameter Estimation: The transformed data's mean, variance, and covariance 
#    matrix are estimated, aligning the results with Darwin's original study.
# 
# 5. Likelihood Ratio Test: A likelihood ratio test compares the model with and 
#    without the transformation to evaluate its effect.
# 
# 6. Post-Transformation Analysis: The script re-evaluates the normality of the 
#    transformed data and recalculates skewness and Shapiro-Wilk statistics.
# 
# 7. Quantile Estimation: The script estimates the 0.01-quantile of the mean 
#    height difference in the transformed data to assess the impact of cross-fertilization.


# LOAD DEPENDENCIES:
# Load external functions contained in the script file "functions.R" 
# used for Yeo-Johnson transformation, likelihood calculation, and skewness estimation.
file_path <- "functions.R"
if (!file.exists(file_path)) {
  stop("The file does not exist in the specified path:", file_path)
}
source(file_path)

# Data: Height differences between cross-fertilized and self-fertilized plants
heights_diff <- c(-8.4, -6.0, 0.7, 1.0, 1.8, 2.0, 2.9, 3.0, 3.5, 3.6, 5.1, 6.1, 7.0, 7.5, 9.3)

# Function for initial statistical analysis
# Performs a one-tailed paired t-test to see if the mean height difference is greater than zero,
# followed by a skewness check and Shapiro-Wilk test for normality.
initial_stat_analysis <- function(data) {
  # Paired t-test (one-tailed to match the paper)
  paired_t_test <- t.test(data, mu = 0, alternative = "greater")
  
  # Skewness and Shapiro-Wilk test for normality
  skewness_stat <- calculate_skewness(data)
  shapiro_test <- shapiro.test(data)
  
  # Density plot for original data
  plot(density(data), main = "Original Sample Density")
  
  # Print initial statistics
  cat("Paired t-statistic:", sprintf("%.5f", paired_t_test$statistic), "\n")
  cat(sprintf("One-tailed p-value (df = %d): %.5f\n", paired_t_test$parameter, paired_t_test$p.value))
  cat(sprintf("Skewness statistic (original data): %.5f\n", skewness_stat))
  cat(sprintf("Shapiro-Wilk test statistic: %.5f\n", shapiro_test$statistic))
  cat(sprintf("Shapiro-Wilk test p-value: %.5f\n", shapiro_test$p.value))
}

### Perform initial statistical analysis on Darwin's data
initial_stat_analysis(heights_diff)


# Function to maximize log-likelihood for Yeo-Johnson transformation
# This optimizes the λ parameter by maximizing the likelihood function.
objective_function <- function(lambda, data) {
  transformed_data <- yeo_johnson(data, lambda)
  mu_init <- mean(transformed_data)
  sigma2_init <- var(transformed_data)
  log_likelihood(lambda, mu_init, sigma2_init, data)
}

# Function to optimize the λ parameter using log-likelihood maximization
optimize_lambda_loglikelihood <- function(data) {
  optim_result <- optim(par = 0, fn = objective_function, data = data, method = "L-BFGS-B", lower = -2, upper = 2)
  return(optim_result$par)
}

# Function to estimate parameters and covariance matrix for transformed data
# This optimizes λ, μ, and σ², and calculates the covariance matrix using the Hessian matrix.
estimate_parameters_and_covariance <- function(data) {
  best_lambda <- optimize_lambda_loglikelihood(data)
  transformed_data <- yeo_johnson(data, best_lambda)
  
  # Estimate mean and variance for the transformed data
  mu_hat <- mean(transformed_data)
  sigma2_hat <- var(transformed_data)
  
  # Perform optimization to refine λ, μ, and σ²
  optim_result <- optim(par = c(best_lambda, mu_hat, sigma2_hat), fn = function(params) {
    lambda <- params[1]
    mu <- params[2]
    sigma2 <- params[3]
    # Debug information on the correct convergence of values
    # cat("Current lambda:", lambda, "\n")
    # cat("Current mu:", mu, "\n")
    # cat("Current sigma^2:", sigma2, "\n")
    log_likelihood(lambda, mu, sigma2, data)
  }, hessian = TRUE, method = "L-BFGS-B")
  
  # Extract optimized parameters
  best_params <- optim_result$par
  best_lambda <- best_params[1]
  mu_hat <- best_params[2]  
  sigma2_hat <- best_params[3]  
  
  # Compute covariance matrix from Hessian
  hessian_matrix <- optim_result$hessian
  covariance_matrix <- solve(hessian_matrix) * 15  # Scaling factor
  
  return(list(
    lambda = best_lambda,
    mu = mu_hat,
    sigma2 = sigma2_hat,
    covariance_matrix = round(covariance_matrix, 3)
  ))
}

# Function to perform the Yeo-Johnson transformation and subsequent analysis
transformation_analysis <- function(data) {
  param_estimates <- estimate_parameters_and_covariance(data)
  best_lambda <- param_estimates$lambda
  mu_hat <- param_estimates$mu
  sigma2_hat <- param_estimates$sigma2
  covariance_matrix <- param_estimates$covariance_matrix
  
  # Apply Yeo-Johnson transformation and perform normality tests on transformed data
  transformed_data <- yeo_johnson(data, best_lambda)
  shapiro_test <- shapiro.test(transformed_data)
  transformed_skewness <- calculate_skewness(transformed_data)
  
  # Perform paired t-test on transformed data
  paired_t_test_transformed <- t.test(transformed_data, mu = 0, alternative = "greater")
  
  # Likelihood ratio test comparing λ = 1 (no transformation) vs optimal λ
  log_likelihood_original <- -log_likelihood(1, mean(data), var(data), data)
  log_likelihood_transformed <- -log_likelihood(best_lambda, mu_hat, sigma2_hat, data)
  log_likelihood_ratio <- 2 * (log_likelihood_transformed - log_likelihood_original)
  p_value_likelihood_ratio <- 1 - pchisq(log_likelihood_ratio, df = 1)
  
  return(list(
    lambda = best_lambda,
    transformed_data = transformed_data,
    mu = mu_hat,
    sigma2 = sigma2_hat,
    covariance_matrix = covariance_matrix,
    transformed_skewness = transformed_skewness,
    shapiro_test = shapiro_test,
    t_statistic_transformed = paired_t_test_transformed$statistic,
    p_value_transformed = paired_t_test_transformed$p.value,
    log_likelihood_ratio = log_likelihood_ratio,
    p_value_likelihood_ratio = p_value_likelihood_ratio
  ))
}

# Function to print the final results after Yeo-Johnson transformation
print_final_statistics <- function(result) {
  cat("Best lambda:", result$lambda, "\n")
  cat("Estimated mu:", result$mu, "\n")
  cat("Estimated sigma^2:", result$sigma2, "\n")
  
  # Plot the transformed data's density
  plot(density(yeo_johnson(heights_diff, result$lambda)), 
       main = paste("Density after Yeo-Johnson Transformation (λ =", round(result$lambda, 3), ")"))
  
  # Print the final statistics and test results
  cat("Paired t-statistic (transformed data):", result$t_statistic_transformed, "\n")
  cat("p-value (transformed data):", result$p_value_transformed, "\n")
  cat("Shapiro-Wilk test (transformed data): W =", result$shapiro_test$statistic, ", p-value =", result$shapiro_test$p.value, "\n")
  cat("Skewness statistic (transformed data):", result$transformed_skewness, "\n")
  cat("Covariance matrix:\n")
  print(result$covariance_matrix)
  
  # Likelihood ratio test result
  cat("For λ = 1 \n")
  cat("Likelihood ratio chi-squared statistic:", result$log_likelihood_ratio, "\n")
  cat("p-value of the likelihood ratio test:", result$p_value_likelihood_ratio, "\n")
  
  # Estimate and print the quantile for the mean difference
  estimated_quantile <- estimate_quantile_mean_difference(result$lambda, 
                                                          result$mu, 
                                                          result$sigma2, 
                                                          length(heights_diff))
  cat("Quantile 0.01 (mu:", result$mu, "sigma^2:", result$sigma2, "):", estimated_quantile, "\n")
}

### Performing the transformation analysis on Darwin's data
result <- transformation_analysis(heights_diff)

### Printing the final results
print_final_statistics(result)

# PROJECT ENDS HERE.
