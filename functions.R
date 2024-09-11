# Author: Giuseppe Muschetta
# File: functions.R

################################# BOX-COX 
#' Apply the original Box-Cox transformation to stabilize variance and make data 
#' more normally distributed.
#'
#' @param x Numeric vector, strictly positive, the input data for which to apply the Box-Cox transformation.
#' @param lambda Numeric scalar, the transformation parameter. If lambda equals 0, the transformation
#'        defaults to the natural logarithm of x.
#'
#' @return Numeric vector containing the Box-Cox transformation for each element in x.
#'
#' @examples
#' # Apply Box-Cox transformation to a vector of positive numbers with lambda = 0.5
#' box_cox(c(1, 2, 3, 4, 5), 0.5)
#'
#' @export
box_cox <- function(x, lambda) {
  # Defensive programming approach:
  
  # Making sure x is a numeric vector with positive values.
  if (!is.numeric(x) || !is.vector(x) || any(x <= 0)) {
    stop("Input 'x' must be a numeric vector with all elements positive.")
  }
  
  # Making sure lambda is a numeric scalar
  if (!is.numeric(lambda) || length(lambda) != 1) {
    stop("Parameter 'lambda' must be a numeric scalar.")
  }
  
  # Initialize a vector for results of the same length as x
  y <- numeric(length(x))
  
  # Apply the Box-Cox transformation
  pos_values <- x > 0  # Only positive values are valid for Box-Cox
  
  y[pos_values] <- if (lambda == 0) {
    log(x[pos_values])
  } else {
    (x[pos_values] ^ lambda - 1) / lambda
  }
  
  # Returns a numeric vector with transformed data
  return(y)
}


######################## SIGNED POWER TRANSFORMATION
#' Apply the Signed Power Transformation to data.
#'
#' This function applies the Signed Power Transformation to a numeric vector. The transformation
#' adjusts data values using a power function, taking into account the sign of each value, 
#' which makes it suitable for data that includes both positive and negative values.
#'
#' @param x Numeric vector, the input data for which to apply the signed power transformation.
#'          It should be a vector containing any real numbers.
#' @param lambda Numeric scalar, the transformation parameter, which should be strictly greater than zero.
#'          This parameter controls the power to which the absolute value of each data point is raised.
#'
#' @return Numeric vector containing the signed power transformation for each element in x. The transformation
#'         is performed separately for positive and negative values, maintaining their signs post-transformation.
#'
#' @examples
#' # Apply the signed power transformation to a vector of mixed sign numbers with lambda = 0.5
#' signed_power_transformation(c(-2, -1, 0, 1, 2), 0.5)
#'
#' @export
signed_power <- function(x, lambda) {
  # Defensive programming approach:
  
  # Making sure x is a numeric vector
  if (!is.numeric(x) || !is.vector(x)) {
    stop("Input 'x' must be a numeric vector.")
  }
  
  # Making sure lambda is a numeric scalar and greater than 0
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0) {
    stop("Parameter 'lambda' must be a numeric scalar greater than 0.")
  }
  
  # Initialize a vector for results of the same length as x
  y <- numeric(length(x))
  
  # Vectors to identify positive and negative values
  pos_values <- x >= 0
  neg_values <- x < 0
  
  # Apply the transformation for positive elements.
  # For positive values, a straightforward power transformation is applied.
  y[pos_values] <- (x[pos_values]^lambda - 1) / lambda
  
  # Apply the transformation for negative elements.
  # For negative values, the Signed Power transformation preserves the sign 
  # and applies a power transformation based on the absolute value.
  y[neg_values] <- (sign(x[neg_values]) * abs(x[neg_values])^(lambda) - 1) / lambda
  
  # Return the transformed data
  return(y)
}

######################## MIXTURE DENSITY FUNCTION
#' Calculate the mixture density function consisting of a standard normal and 
#' a shifted gamma distribution.
#'
#' @param x Numeric vector representing the values at which to evaluate the mixed density.
#' @return Numeric vector containing the density values for each element of x.
#' @description
#' The function computes a mixture density that is 30% standard normal and 70% a specific shifted gamma distribution.
#' The gamma distribution component is adjusted to be defined for x > -2, using a specified transformation.
#' 
#' @examples
#' # Calculate the mixture density at a range of x values
#' x_vals <- seq(-10, 10, length.out = 1000)
#' densities <- mixture_density(x_vals)
#' plot(x_vals, densities, type = 'l')
#'
#' @export
mixture_density <- function(x) {
  # Calculate standard normal density
  phi <- dnorm(x)
  
  # Pre-allocate gamma density vector
  gamma_shifted <- numeric(length(x))
  
  # Compute gamma density only for x > -2
  valid_values <- x > -2
  gamma_shifted[valid_values] <- (1/6) * (x[valid_values] + 2)^3 * exp(-(x[valid_values] + 2))
  
  # Calculate the mixture density
  f_x <- 0.3 * phi + 0.7 * gamma_shifted
  
  return(f_x)
}



########################## YEO JOHNSON POWER TRANSFOMTATION
#' Yeo-Johnson Power Transformation function R x R --> R
#' This function applies the Yeo-Johnson transformation to a numeric vector. 
#' It is designed to transform non-normal dependent variables into a normal shape. 
#' Yeo-Johnson can handle both positive and negative data by adjusting 
#' the lambda parameter.
#'
#' @param x A numeric vector that contains the data to be transformed.
#' @param lambda A numeric scalar that controls the nature of the transformation.
#'  Specific cases:
#'  if (x >= 0):
#'    if(lambda == 0)
#'      "Transforms data using a logarithmic scale: log(x + 1)"
#'    else
#'      "Performs a Box-Cox transformation with a shifting -1 constant"
#'  if (x < 0):
#'    if(lambda == 2)
#'      "Transforms data using -log(-x + 1) to be less sensitive at the lower bound"
#'    else
#'      "Power transformation where original lambda has become lambda - 2
#'  End_specific_cases.  
#' @return A numeric vector with the transformed data.
#' @examples of usage:
#' set.seed(123)
#' data <- rlnorm(100, meanlog = 0, sdlog = 1) - 50
#' transformed_data <- yeo_johnson(data, lambda = 0.5)
#'
#' @importFrom stats rlnorm
#' @export
yeo_johnson <- function(x, lambda) {
  # Defensive programming approach:
  
  # Making sure x is a numeric vector
  if (!is.numeric(x) || !is.vector(x)) {
    stop("Input 'x' must be a numeric vector.")
  }
  
  # Making sure lambda is a numeric scalar
  if (!is.numeric(lambda) || length(lambda) != 1) {
    stop("Parameter 'lambda' must be numeric scalar.")
  }
  
  # Initialize a vector for results of the same length as x
  y <- numeric(length(x))
  pos_values <- x >= 0
  neg_values <- x < 0
  
  # Apply the transformation for negative elements.
  # The Yeo-Johnson transformation handles negative values by applying a modified 
  # power transformation, allowing for continuous transformation across zero. 
  # The sign is preserved during the process.
  y[neg_values] <- if (lambda == 2) {
    -log(-x[neg_values] + 1)
  } else {
    -((-x[neg_values] + 1)^(2 - lambda) - 1) / (2 - lambda)
  }
  
  # Apply the transformation for positive elements.
  # For positive values, the transformation is similar to the Box-Cox transformation, 
  # ensuring that the data can be transformed to approximate normality while 
  # preserving the continuity at zero.
  y[pos_values] <- if (lambda == 0) {
    log(x[pos_values] + 1)
  } else {
    ((x[pos_values] + 1)^lambda - 1) / lambda
  }
  
  # Returns a numeric vector with transformed data
  return(y)
}



########################## INVERSE YEO-JOHNSON POWER TRANSFORMATION
#' Inverse Yeo-Johnson Power Transformation function R x R --> R
#' This function applies the inverse Yeo-Johnson transformation to a numeric vector.
#' It is designed to revert data that has been previously transformed using the
#' Yeo-Johnson transformation back to its original scale.
#' @param y A numeric vector that contains the data to be inversely transformed.
#' @param lambda A numeric scalar that controls the nature of the transformation.
#' @return A numeric vector with the inversely transformed data.
#' @examples of usage:
#' set.seed(123)
#' transformed_data <- c(-0.5, 0, 0.5, 1)
#' original_data <- inv_yeo_johnson(transformed_data, lambda = 0.5)
#'
#' @export
inv_yeo_johnson <- function(y, lambda) {
  # Defensive programming approach:
  
  # Making sure y is a numeric vector
  if (!is.numeric(y) || !is.vector(y)) {
    stop("Input 'y' must be a numeric vector.")
  }
  
  # Making sure lambda is a numeric scalar
  if (!is.numeric(lambda) || length(lambda) != 1) {
    stop("Parameter 'lambda' must be numeric scalar.")
  }
  
  # Initialize a vector for results of the same length as y
  x <- numeric(length(y))
  pos_values <- y >= 0
  neg_values <- y < 0
  
  # Deals with positive elements.
  x[pos_values] <- if (lambda == 0) {
    exp(y[pos_values]) - 1
  } else {
    (y[pos_values] * lambda + 1)^(1 / lambda) - 1
  }
  
  # Deals with negative elements.
  x[neg_values] <- if (lambda == 2) {
    -exp(-y[neg_values]) + 1
  } else {
    -((1 - lambda * y[neg_values])^(1 / (2 - lambda)) - 1)
  }
  
  # Returns a numeric vector with inversely transformed data
  return(x)
}



############################ YEO JOHNSON DERIVATIVES
#  Compute the k-th derivative of the Yeo-Johnson transformation
#' @param x Numeric vector, the input data for which to compute the transformation.
#' @param lambda Numeric scalar, the transformation parameter, not equal to 0 or 2 for recursive calls.
#' @param k Integer scalar, the order of the derivative to compute, must be >= 1.
#' @return Numeric vector containing the k-th derivative of the Yeo-Johnson transformation for each element in x.
#' @examples
#' # Calculates the 2nd derivative for x = 1.5 and lambda = 0.5
#' yeo_johnson_derivative(1.5, 0.5, 2)
#' @export
yeo_johnson_derivative <- function(x, lambda, k) {
  # Validate inputs
  if (!is.numeric(x) || !is.numeric(lambda) || !is.numeric(k) || length(lambda) != 1 || length(k) != 1 || k < 1) {
    stop("Inputs must be numeric and lambda, k must be numeric scalars with k >= 1.")
  }
  
  # Initialize a vector for the derivative results
  y_derivative <- numeric(length(x))
  pos_values <- x >= 0
  neg_values <- x < 0
  
  # Base case: k = 1 uses the original Yeo-Johnson function for k-1
  if (k == 1) {
    yeo_johnson_base <- yeo_johnson(x, lambda)
    
    # Calculate the first derivative for positive and negative x values
    y_derivative[pos_values] <- if (lambda != 0) {
      (x[pos_values] + 1)^lambda * log(x[pos_values] + 1) - yeo_johnson_base[pos_values]
    } else {
      log(x[pos_values] + 1)^2 / 2  # Special case when lambda = 0
    }
    
    y_derivative[neg_values] <- if (lambda != 2) {
      -((-x[neg_values] + 1)^(2 - lambda) * log(-x[neg_values] + 1) - yeo_johnson_base[neg_values]) / (2 - lambda)
    } else {
      -log(-x[neg_values] + 1)^2 / 2  # Special case when lambda = 2
    }
  } else {
    # Recursive call: Compute the (k-1)th derivative using the same function
    previous_derivative <- yeo_johnson_derivative(x, lambda, k - 1)
    
    # Calculate the k-th derivative for positive and negative x values
    y_derivative[pos_values] <- if (lambda != 0) {
      (x[pos_values] + 1)^lambda * log(x[pos_values] + 1)^k - k * previous_derivative[pos_values] / lambda
    } else {
      log(x[pos_values] + 1)^(k+1) / (k + 1)  # When lambda = 0
    }
    
    y_derivative[neg_values] <- if (lambda != 2) {
      -((-x[neg_values] + 1)^(2 - lambda) * log(-x[neg_values] + 1)^k - k * previous_derivative[neg_values]) / (2 - lambda)
    } else {
      -log(-x[neg_values] + 1)^(k+1) / (k + 1)  # When lambda = 2
    }
  }
  
  return(y_derivative)
}



############################ TEST CONTINUITY FUNCTION
#' Test Continuity of a Function at a Specified Point
#'
#' This function tests the continuity of either the Yeo-Johnson transformation or its derivatives
#' at a specified point. It evaluates the function within a defined range around the point,
#' and checks if all values are close to the function value at the point itself within a specified tolerance.
#'
#' @param func Function to test, either the Yeo-Johnson transformation or its derivative.
#' @param x0 Numeric, the x-coordinate at which to test continuity.
#' @param lambda0 Numeric, the lambda parameter at which to test continuity.
#' @param k Integer, optional; specifies the order of the derivative to test. If NULL, tests the base function.
#' @param delta Numeric, the radius around (x0, lambda0) to test for continuity; defaults to 0.1.
#' @param n_points Integer, the number of points in each dimension (x and lambda) to generate in the grid; defaults to 300.
#' @param tolerance Numeric, the maximum allowed difference for continuity; defaults to 0.1.
#' @param convergence_zone_factor Numeric, the factor to define a smaller zone around the test point for a precise continuity check; defaults to 10.
#'
#' @return A list with elements:
#'   - `Continuous`: Logical, TRUE if the function is continuous at the point within the tolerance, FALSE otherwise.
#'   - `Expected`: Numeric, the expected function value at the point (x0, lambda0).
#'
#' @examples
#' test_continuity(yeo_johnson, 1.0, 0.5)  # Test continuity of the Yeo-Johnson function
#' test_continuity(yeo_johnson_derivative, 1.0, 0.5, k = 2)  # Test continuity of the second derivative
#'
#' @export
test_continuity <- function(func, x0, lambda0, k = NULL, delta = 0.01, n_points = 300, tolerance = 0.1, convergence_zone_factor = 10) {
  
  # Generate sequences of x and lambda values around the test point.
  # This creates a grid of values centered around (x0, lambda0) to evaluate the 
  # function's behavior in a small neighborhood, checking for potential discontinuities.
  x_seq <- seq(x0 - delta, x0 + delta, length.out = n_points)
  lambda_seq <- seq(lambda0 - delta, lambda0 + delta, length.out = n_points)
  
  # Calculate the function values for a grid of (x, lambda) values
  results <- outer(x_seq, lambda_seq, Vectorize(function(x, lambda) {
    if (is.null(k)) {
      func(x, lambda)  # Testing the base function
    } else {
      func(x, lambda, k)  # Testing the derivative function
    }
  }))
  
  # Determine the expected value at the test point
  expected_value <- if (is.null(k)) {
    func(x0, lambda0)
  } else {
    func(x0, lambda0, k)
  }
  
  # Select values from the results that are within a defined smaller zone for more precise continuity check
  converging_values <- results[abs(x_seq - x0) < delta / convergence_zone_factor & abs(lambda_seq - lambda0) < delta / convergence_zone_factor]
  
  # Check if all selected values are within the tolerance of the expected value
  continuity_check <- all(abs(converging_values - expected_value) < tolerance)
  
  # Ulterior explanations on the use of the "convergence_zone_factor" parameter:
  # A high value such as 100 is used to focus closely around the test point,
  # allowing for precise assessment of continuity in that very restricted region.
  # This is particularly useful for functions with complex behaviors near test points.
  return(list(Expected = expected_value, Check = continuity_check))
}



########################## PLOT CONTINUITY FUNCTION
#' Visualize the Continuity of Yeo-Johnson or Its Derivatives
#'
#' This function plots the Yeo-Johnson transformation or one of its k-th order derivatives
#' over a specified range of x values for a given lambda value,
#' facilitating the graphical visualization of continuity or discontinuities.
#'
#' @param func The function to plot; can be either yeo_johnson or yeo_johnson_derivative.
#' @param x_range Numeric vector with two elements specifying the minimum and maximum x values to plot.
#' @param lambda The value of the lambda parameter for which to plot the function.
#' @param k The degree of the derivative to plot; if NULL, plots the base function.
#'
#' @return Does not return a value but produces a plot.
#' @examples
#' plot_continuity(yeo_johnson, c(-2, 2), 0.5)
#' plot_continuity(yeo_johnson_derivative, c(-2, 2), 0.5, k = 2)
#'
#' @importFrom graphics plot lines
#' @export
plot_continuity <- function(func, x_range, lambda, k = NULL) {
  x_values <- seq(from = x_range[1], to = x_range[2], length.out = 100)
  y_values <- if (is.null(k)) {
    sapply(x_values, function(x) func(x, lambda))
  } else {
    sapply(x_values, function(x) func(x, lambda, k))
  }
  
  plot(x_values, y_values, type = 'l', main = paste("Continuity Plot for k =", ifelse(is.null(k), "Function", k), 
                                                    "and lambda =", lambda),
       xlab = "x", ylab = "Function Value", col = "blue")
}



##################### TEST INCREASING FUNCTION
# Function to test if the transformation function is increasing with respect to lambda and x
#' Test Increasing Behavior of a Function with Respect to Lambda and X
#'
#' This function evaluates whether a given transformation function is monotonically increasing
#' with respect to lambda and x by comparing the function's output at a base point and 
#' after a small increment in lambda or x.
#'
#' @param func Function to be tested, should take x and lambda as inputs.
#' @param x Numeric, the x-coordinate at which the function's increasing behavior is to be tested.
#' @param lambda Numeric, the lambda value at which the function's increasing behavior is to be tested.
#' @param delta_lambda Numeric, the small increment to lambda to test for increase, default is 0.1.
#' @param delta_x Numeric, the small increment to x to test for increase, default is 0.1.
#'
#' @return A list containing two logical values:
#'   - `IncreasingLambda`: TRUE if the function value increases with an increase in lambda,
#'     FALSE otherwise.
#'   - `IncreasingX`: TRUE if the function value increases with an increase in x, FALSE otherwise.
#'
#' @examples
#' test_increasing(yeo_johnson, x = 1, lambda = 0.5)
#' test_increasing(yeo_johnson_derivative, x = 1, lambda = 0.5, k = 1) # If testing derivatives
#'
#' @export
test_increasing <- function(func, x, lambda, delta_x = 0.1, delta_lambda = 0.1) {
  # Evaluate the function at the original lambda and x
  initial_value <- func(x, lambda)
  
  # Check if the function is zero when x = 0, which is a special case in which 
  # the function remains constant at 0, so it cannot increase!
  if (x == 0 && initial_value == 0) {
    print("At x = 0, the function is constant at zero.")
    return(list(
      IncreasingLambda = NA,
      IncreasingX = NA
    ))
  }
  
  # Evaluate the function at incremented lambda and x
  increased_value_lambda <- func(x, lambda + delta_lambda)
  increased_value_x <- func(x + delta_x, lambda)
  
  # Determine if the function values are increasing
  is_increasing_lambda <- increased_value_lambda > initial_value
  is_increasing_x <- increased_value_x > initial_value
  
  return(list(
    IncreasingLambda = is_increasing_lambda,
    IncreasingX = is_increasing_x
  ))
}




##################### TEST CONVEXITY CONCAVITY FUNCTION
# Function to test convexity/concavity of the Yeo-Johnson transformation with respect to lambda
#' Test Convexity or Concavity in Lambda for Yeo-Johnson Transformation
#'
#' This function calculates the second derivative of the Yeo-Johnson transformation
#' with respect to lambda and determines if the function is convex (for x > 0)
#' or concave (for x < 0).
#'
#' @param x Numeric, the value of x to test.
#' @param lambda Numeric, the value of lambda to test.
#' @return A list containing the second derivative and the convexity/concavity property.
#' @examples
#' x_values = seq(-2, 2, by = 0.5)
#' lambda_test = 0.5  # A typical value of lambda for testing
#' for (x in x_values) {
#'   result = test_convexity_concavity(x, lambda_test)
#'   cat(sprintf("For x = %f: The function is %s, second derivative = %f\n", x, result$Property, result$SecondDerivative))
#' }
#' @export
test_convexity_concavity <- function(x, lambda) {
  # Calculate the second derivative with respect to lambda
  second_derivative <- yeo_johnson_derivative(x, lambda, 2)
  
  # Determine convexity for x > 0 (second derivative >= 0)
  if (x > 0) {
    is_convex = (second_derivative >= 0)
    property = ifelse(is_convex, "convex", "not convex")
  }
  
  # Determine concavity for x < 0 (second derivative <= 0)
  if (x < 0) {
    is_concave = (second_derivative <= 0)
    property = ifelse(is_concave, "concave", "not concave")
  }
  
  # Special consideration for x = 0.
  # At x = 0, determining convexity or concavity is difficult because the 
  # second derivative may not provide clear information, and other methods may 
  # be required to assess the function's properties.
  if (x == 0) {
    property = "check manually"
  }
  
  return(list(
    SecondDerivative = second_derivative,
    Property = property
  ))
}



##################### LOG-LIKELIHOOD FUNCTION TO MAXIMIZE (by minimizing the minimum)
#' Calculate the log-likelihood for the Yeo-Johnson transformation.
#'
#' @param lambda A numeric scalar, the transformation parameter lambda.
#' @param mu A numeric scalar, the mean of the transformed data.
#' @param sigma2 A numeric scalar, the variance of the transformed data.
#' @param x A numeric vector, the original data to be transformed.
#' 
#' @return A numeric scalar, the negative value of the log-likelihood function.
#' 
#' @examples
#' set.seed(123)
#' data <- rnorm(100)
#' log_likelihood(lambda = 0.5, mu = 0, sigma2 = 1, x = data)
#' 
#' @export
log_likelihood <- function(lambda, mu, sigma2, x) {
  # Defensive programming: Check that x is a numeric vector.
  if (!is.numeric(x) || !is.vector(x)) {
    stop("Input 'x' must be a numeric vector.")
  }
  
  # Check that lambda, mu, and sigma2 are numeric scalars.
  if (!is.numeric(lambda) || length(lambda) != 1) {
    stop("Parameter 'lambda' must be a numeric scalar.")
  }
  if (!is.numeric(mu) || length(mu) != 1) {
    stop("Parameter 'mu' must be a numeric scalar.")
  }
  if (!is.numeric(sigma2) || length(sigma2) != 1) {
    stop("Parameter 'sigma2' must be a numeric scalar.")
  }
  
  # Apply the Yeo-Johnson transformation to the data.
  transformed_x <- yeo_johnson(x, lambda)
  
  # Calculate the number of data points.
  n <- length(x)
  
  # Calculate the log-likelihood.
  log_likelihood_value <- -(n / 2) * log(2 * pi) - (n / 2) * log(sigma2) - (1 / (2 * sigma2)) * sum((transformed_x - mu)^2) + 
    (lambda - 1) * sum(sign(x) * log(abs(x) + 1))
  
  # Return the negative log-likelihood for minimization purposes.
  return(-log_likelihood_value)
}




############################# CALCULATE SKEWNESS 
#' This function calculates the skewness using Minitab formula: 
#' g1 = (sum((x - mean(x))^3) / n) / (sum((x - mean(x))^2) / n)^1.5
#' b1 = g1 * ((1 - 1 / n))^1.5
#' @param x A numeric vector for which to calculate the skewness.
#' @return A numeric value representing the skewness of the input vector.
#' @examples
#' data <- c(1, 2, 3, 4, 5)
#' calculate_skewness(data)
#' @export
calculate_skewness <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  sd_x <- sd(x)
  g1 <- (sum((x - mean_x)^3) / n) / (sd_x^3)
  b1 <- g1 * ((n - 1) / n)^1.5
  return(b1)
}




##################### ESTIMATE QUANTILE OF MEAN DIFFERENCE
#' 
#' This function estimates the \eqn{\gamma}-th quantile of the mean difference in the transformed data.
#' 
#' @param lambda The transformation parameter obtained from the Yeo-Johnson transformation.
#' @param mu The estimated mean of the transformed data.
#' @param sigma2 The estimated variance of the transformed data.
#' @param n The sample size of the original data.
#' @param gamma The quantile to estimate, with a default value of 0.01.
#' 
#' @return The estimated \eqn{\gamma}-th quantile of the mean difference in the transformed data.
#' 
#' @examples
#' lambda <- 1.305  # Example lambda value from transformation
#' mu <- 4.570      # Example estimated mean of transformed data
#' sigma2 <- 29.786 # Example estimated variance of transformed data
#' n <- 15          # Sample size
#' gamma <- 0.01    # Quantile to estimate
#' 
#' estimated_quantile <- estimate_quantile(lambda, mu, sigma2, n, gamma)
#' print(estimated_quantile)
#' 
#' @export
estimate_quantile_mean_difference <- function(lambda, mu, sigma2, n, gamma = 0.01) {
  # Calculate the t quantile for the given gamma and degrees of freedom
  t_quantile <- qt(gamma, df = n - 1)
  
  # Estimate the quantile of the mean difference in the transformed data
  estimated_quantile <- inv_yeo_johnson(mu + t_quantile * (sqrt(sigma2) / sqrt(n)), lambda)
  
  # Return the estimated quantile
  return(estimated_quantile)
}



