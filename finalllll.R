###############################################################################
# Section 1: Data Loading and Preprocessing
###############################################################################

# Load necessary library
library(fda)
options(warn=-1)

# Set working directory to where the data files are located
setwd("~/Desktop/New Folder With Items")

# Load data from CSV files
grw = read.csv("GDP g.csv")  # GDP growth data
ipc = read.csv("GDP y.csv")  # GDP per capita data
fdi = read.csv("FDI r.csv")  # Foreign Direct Investment data

# Convert data frames to matrices and remove the last row if necessary
grw = as.matrix(grw)
ipc = as.matrix(ipc[-nrow(ipc), ])
fdi = as.matrix(fdi[-nrow(fdi), ])

# Define the time sequence from 1980 to 2020
Time = 1980:2020

###############################################################################
# Section 2: Defining Basis Functions for Functional Data Smoothing
###############################################################################

# Define the range of the data
rangeval <- range(Time)

# Set the number of basis functions to use
nbasis <- 35  # Adjust this number as needed

# Create Fourier basis functions for each data set
grw_basis <- create.fourier.basis(rangeval, nbasis)  # For GDP growth
ipc_basis <- create.fourier.basis(rangeval, nbasis)  # For GDP per capita
fdi_basis <- create.fourier.basis(rangeval, nbasis)  # For FDI

# Define functional parameter objects with smoothing parameter lambda
grw_fdPar <- fdPar(grw_basis, Lfdobj = 2, lambda = 0.01)  # For GDP growth
ipc_fdPar <- fdPar(ipc_basis, Lfdobj = 2, lambda = 0.01)  # For GDP per capita
fdi_fdPar <- fdPar(fdi_basis, Lfdobj = 2, lambda = 0.01)  # For FDI

###############################################################################
# Section 3: Smoothing the Functional Data
###############################################################################

# Smooth the data using the defined basis and functional parameters
grw_fd <- smooth.basis(Time, grw, grw_fdPar)$fd  # Smoothed GDP growth functions
ipc_fd <- smooth.basis(Time, ipc, ipc_fdPar)$fd  # Smoothed GDP per capita functions
fdi_fd <- smooth.basis(Time, fdi, fdi_fdPar)$fd  # Smoothed FDI functions

# Obtain the y to c mapping matrices (used later for inference)
grw_y2c <- smooth.basis(Time, grw, grw_basis)$y2cMap
ipc_y2c <- smooth.basis(Time, ipc, ipc_basis)$y2cMap
fdi_y2c <- smooth.basis(Time, fdi, fdi_basis)$y2cMap

# Get the country names from the cleaned data
countries <- colnames(grw)

# Load regional data and align with country data
regional <- read.csv("Country_Income_Region_Info_Completed_excel.csv", stringsAsFactors = FALSE)
regional <- regional[regional$Country %in% countries, ]  # Filter to include only present countries
regional <- regional[match(countries, regional$Country), ]  # Ensure the order matches

###############################################################################
# Section 4: Registration of Functional Data
###############################################################################

# Define a more flexible basis for the warping functions
Wnbasis <- 7  # Number of basis functions for warping; adjust as needed
Wbasis <- create.bspline.basis(rangeval = range(Time), nbasis = Wnbasis)
Wfd0 <- fd(matrix(0, Wnbasis, 1), Wbasis)

# Set up the functional parameter object for warping functions
WfdParobj <- fdPar(Wfd0, Lfdobj = 2, lambda = 0.01)

#---------------------------- Registration of FDI Data ----------------------------#

# Remove columns with NA values from FDI data
fdi <- fdi[, colSums(is.na(fdi)) == 0]

# Remove columns with zero variance (constant functions)
fdi_variances <- apply(fdi, 2, var)
fdi <- fdi[, fdi_variances > 0]

# Define a sequence of lambda values to try for smoothing
lambdas <- exp(seq(log(0.1), log(0.001), length.out = 10))

# Function to perform registration given a lambda value
perform_registration <- function(lambda) {
  fdi_fdPar <- fdPar(fdi_basis, Lfdobj = 2, lambda = lambda)
  fdi_fd <- smooth.basis(Time, fdi, fdi_fdPar)$fd
  
  # Find the curves with the minimum and maximum area under the curve (AUC)
  auc <- apply(fdi, 2, sum)
  min_auc_curve <- which.min(auc)
  max_auc_curve <- which.max(auc)
  
  # Create a target function that combines the min and max curves
  target_fd <- fdi_fd[, c(min_auc_curve, max_auc_curve)]
  target_fd <- mean.fd(target_fd)  # Average of min and max curves
  
  # Perform registration using the min-max target function
  fdi_fd_registered <- register.fd(yfd = fdi_fd, y0fd = target_fd, WfdParobj = WfdParobj)
  
  return(list(original = fdi_fd, registered = fdi_fd_registered))
}

# Loop through lambda values to find the best registration
for (lambda in lambdas) {
  cat("Trying lambda:", lambda, "\n")
  tryCatch({
    result <- perform_registration(lambda)
    cat("Registration successful with lambda:", lambda, "\n")
    
    # Plot the results
    par(mfrow = c(2, 2))
    plot(result$original, main = "Original FDI Functions")
    plot(result$registered$regfd, main = "Registered FDI Functions")
    plot(mean.fd(result$registered$regfd), main = "Mean Function After Registration")
    plot(result$registered$warpfd, main = "Warping Functions")
    
    # Exit the loop if successful
    break
    
  }, error = function(e) {
    cat("Error occurred. Trying next lambda.\n")
  })
}

# Save the registered result
fdi_fd_registered <- result

# Check if registration was successful
if (!exists("result")) {
  cat("Registration failed for all lambda values.\n")
} else {
  cat("Final lambda used:", lambda, "\n")
}

#---------------------------- Registration of IPC Data ----------------------------#

# Remove columns with NA values from IPC data
ipc <- ipc[, colSums(is.na(ipc)) == 0]

# Remove columns with zero variance (constant functions)
ipc_variances <- apply(ipc, 2, var)
ipc <- ipc[, ipc_variances > 0]

# Define a sequence of lambda values to try for smoothing
lambdas <- exp(seq(log(0.1), log(0.001), length.out = 10))

# Function to perform registration given a lambda value
perform_registration <- function(lambda) {
  ipc_fdPar <- fdPar(ipc_basis, Lfdobj = 2, lambda = lambda)
  ipc_fd <- smooth.basis(Time, ipc, ipc_fdPar)$fd
  
  # Find the curves with the minimum and maximum area under the curve (AUC)
  auc <- apply(ipc, 2, sum)
  min_auc_curve <- which.min(auc)
  max_auc_curve <- which.max(auc)
  
  # Create a target function that combines the min and max curves
  target_fd <- ipc_fd[, c(min_auc_curve, max_auc_curve)]
  target_fd <- mean.fd(target_fd)  # Average of min and max curves
  
  # Perform registration using the min-max target function
  ipc_fd_registered <- register.fd(yfd = ipc_fd, y0fd = target_fd, WfdParobj = WfdParobj)
  
  return(list(original = ipc_fd, registered = ipc_fd_registered))
}

# Loop through lambda values to find the best registration
for (lambda in lambdas) {
  cat("Trying lambda:", lambda, "\n")
  tryCatch({
    result <- perform_registration(lambda)
    cat("Registration successful with lambda:", lambda, "\n")
    
    # Plot the results
    par(mfrow = c(2, 2))
    plot(result$original, main = "Original IPC Functions")
    plot(result$registered$regfd, main = "Registered IPC Functions")
    plot(mean.fd(result$registered$regfd), main = "Mean Function After Registration")
    plot(result$registered$warpfd, main = "Warping Functions")
    
    # Exit the loop if successful
    break
    
  }, error = function(e) {
    cat("Error occurred. Trying next lambda.\n")
  })
}

ipc_fd_registered <- result

# Check if registration was successful
if (!exists("result")) {
  cat("Registration failed for all lambda values.\n")
} else {
  cat("Final lambda used:", lambda, "\n")
}

#---------------------------- Registration of GRW Data ----------------------------#

# Remove columns with NA values from GRW data
grw <- grw[, colSums(is.na(grw)) == 0]

# Remove columns with zero variance (constant functions)
grw_variances <- apply(grw, 2, var)
grw <- grw[, grw_variances > 0]

# Define a sequence of lambda values to try for smoothing
lambdas <- exp(seq(log(0.1), log(0.001), length.out = 10))

# Function to perform registration given a lambda value
perform_registration <- function(lambda) {
  grw_fdPar <- fdPar(grw_basis, Lfdobj = 2, lambda = lambda)
  grw_fd <- smooth.basis(Time, grw, grw_fdPar)$fd
  
  # Find the curves with the minimum and maximum area under the curve (AUC)
  auc <- apply(grw, 2, sum)
  min_auc_curve <- which.min(auc)
  max_auc_curve <- which.max(auc)
  
  # Create a target function that combines the min and max curves
  if (min_auc_curve != max_auc_curve) {
    target_fd <- grw_fd[, c(min_auc_curve, max_auc_curve)]
    target_fd <- mean.fd(target_fd)  # Average of min and max curves
  } else {
    target_fd <- grw_fd[, min_auc_curve]  # Use the single curve if min and max are the same
  }
  
  # Perform registration using the min-max target function
  grw_fd_registered <- register.fd(yfd = grw_fd, y0fd = target_fd, WfdParobj = WfdParobj)
  
  return(list(original = grw_fd, registered = grw_fd_registered))
}

# Loop through lambda values to find the best registration
for (lambda in lambdas) {
  cat("Trying lambda:", lambda, "\n")
  tryCatch({
    result <- perform_registration(lambda)
    cat("Registration successful with lambda:", lambda, "\n")
    
    # Plot the results
    par(mfrow = c(2, 2))
    plot(result$original, main = "Original GRW Functions")
    plot(result$registered$regfd, main = "Registered GRW Functions")
    plot(mean.fd(result$registered$regfd), main = "Mean Function After Registration")
    if (!is.null(result$registered$warpfd)) {
      plot(result$registered$warpfd, main = "Warping Functions")
    } else {
      cat("Warping function is not available for lambda:", lambda, "\n")
    }
    
    # Exit the loop if successful
    break
    
  }, error = function(e) {
    cat("Error occurred. Trying next lambda.\n")
  })
}

grw_fd_registered <- result

# Check if registration was successful
if (!exists("result")) {
  cat("Registration failed for all lambda values.\n")
} else {
  cat("Final lambda used:", lambda, "\n")
}

###############################################################################
# Section 5: Functional Regression Analysis
###############################################################################

# Re-define the number of basis functions if needed
nbasis <- 35  # Adjust as needed; 21 was correct in previous runs

# Create the Fourier basis for each functional object
grw_basis <- create.fourier.basis(rangeval, nbasis)
ipc_basis <- create.fourier.basis(rangeval, nbasis)
fdi_basis <- create.fourier.basis(rangeval, nbasis)

# Use the registered functional data
grw_fd <- grw_fd_registered$registered$regfd  # Registered GDP growth functions
ipc_fd <- ipc_fd_registered$registered$regfd  # Registered GDP per capita functions
fdi_fd <- fdi_fd_registered$registered$regfd  # Registered FDI functions

# Check dimensions and ensure all functional objects are aligned
print(dim(eval.fd(Time, grw_fd)))
print(dim(eval.fd(Time, ipc_fd)))
print(dim(eval.fd(Time, fdi_fd)))

# Create the interaction term pointwise
ipc_eval <- eval.fd(Time, ipc_fd)
fdi_eval <- eval.fd(Time, fdi_fd)
ipc_fdi_interaction <- ipc_eval * fdi_eval  # Interaction term ipc * fdi

# Convert the interaction back to a functional data object
ipc_fdi_fd <- smooth.basis(Time, ipc_fdi_interaction, ipc_basis)$fd

# Define the intercept as a functional data object
const_basis <- create.constant.basis(rangeval)
const_fd <- fd(matrix(1, 1, length(countries)-1), const_basis)

# Define the list of predictors (including the interaction term)
X_list <- list(
  const = const_fd, 
  ipc_fd = ipc_fd,               # IPC predictor
  fdi_fd = fdi_fd,               # FDI predictor
  ipc_fdi_fd = ipc_fdi_fd        # IPC * FDI interaction predictor
)

# Define basis for beta functions
beta_basis <- create.fourier.basis(rangeval, nbasis)
betalist <- list(
  const = fdPar(const_basis),   # Beta for intercept
  ipc_fd = fdPar(beta_basis),      # Beta for IPC
  fdi_fd = fdPar(beta_basis),      # Beta for FDI
  ipc_fdi_fd = fdPar(beta_basis)   # Beta for IPC * FDI interaction
)

# Perform the functional regression
fRegress_result <- fRegress(grw_fd, X_list, betalist)

# Extract and plot the beta functions
beta_estimates <- fRegress_result$betaestlist

# Plot beta functions
plot(beta_estimates[[1]], main = "Beta 0 (Intercept)")
plot(beta_estimates[[2]], main = "Beta 1 (IPC)")
plot(beta_estimates[[3]], main = "Beta 2 (FDI)")
plot(beta_estimates[[4]], main = "Beta 3 (IPC * FDI Interaction)")

###############################################################################
# Section 6: Functional Regression using 'refund' Package
###############################################################################

# Load necessary libraries
library(lme4)
library(refund)
library(plotly)
library(fda)
library(mgcv)

# Define a fine grid over the time domain
finegrid <- seq(min(Time), max(Time), length.out = 100)

# Evaluate the functional data objects over the fine grid
Y_mat <- eval.fd(finegrid, grw_fd)    # Response matrix (t x n)
ipc_mat <- eval.fd(finegrid, ipc_fd)
fdi_mat <- eval.fd(finegrid, fdi_fd)
ipc_fdi_mat <- ipc_mat * fdi_mat      # Interaction term

# Transpose the matrices to get (n x t)
Y_mat <- t(Y_mat)
ipc_mat <- t(ipc_mat)
fdi_mat <- t(fdi_mat)
ipc_fdi_mat <- t(ipc_fdi_mat)

# Prepare the data list for 'pffr' function
data_pffr <- list(
  Y = Y_mat,          # Response matrix (n x t)
  ipc = ipc_mat,      # Predictor matrices (n x t)
  fdi = fdi_mat,
  ipc_fdi = ipc_fdi_mat
)

# Fit the model using 'pffr' function
model_pffr <- pffr(
  formula = Y ~
    ff(ipc, xind = finegrid ) +
    ff(fdi, xind = finegrid ) +
    ff(ipc_fdi, xind = finegrid),
  yind = finegrid,
  data = data_pffr,
  algorithm = "gam",   # Use 'bam' for large datasets
  method = "REML"     # Use 'fREML' method suitable for 'bam'
)

# Summarize the model
summary(model_pffr)

coef(
  model_pffr,
  raw = FALSE,
  se = TRUE,
  freq = FALSE,
  sandwich = FALSE,
  seWithMean = TRUE)

plot(model_pffr)

# Flatten the response and predictors for 'gam' model
Y_vector <- as.vector(t(Y_mat))
ipc_vector <- as.vector(t(ipc_mat))
fdi_vector <- as.vector(t(fdi_mat))
ipc_fdi_vector <- as.vector(t(ipc_fdi_mat))

# Create a data frame with repeated time points
n <- nrow(Y_mat)
t_points <- ncol(Y_mat)

data_long <- data.frame(
  Y = Y_vector,
  t = rep(finegrid, times = n),
  ipc = ipc_vector,
  fdi = fdi_vector,
  ipc_fdi = ipc_fdi_vector
)

# Fit the function-on-function regression model using tensor product smooths
gam_model <- gam(
  Y ~ te(ipc, t, bs = c("ps", "ps"), k = c(10, 10)) +
    te(fdi, t, bs = c("ps", "ps"), k = c(10, 10)) +
    te(ipc_fdi, t, bs = c("ps", "ps"), k = c(10, 10)),
  data = data_long,
  method = "REML"
)

# Summarize the model
summary(gam_model)

plot(gam_model)

# Function to plot the estimated coefficient surfaces from 'gam' model
plot_beta_surface_gam <- function(gam_model, term_label, predictor_name) {
  # Generate a grid for predictor and time
  predictor_vals <- seq(min(data_long[[predictor_name]]), max(data_long[[predictor_name]]), length.out = 100)
  t_vals <- finegrid
  
  pred_grid <- expand.grid(
    predictor = predictor_vals,
    t = t_vals
  )
  
  # Set up the new data for prediction
  newdata <- data.frame(
    ipc = 0,
    fdi = 0,
    ipc_fdi = 0,
    t = pred_grid$t
  )
  
  newdata[[predictor_name]] <- pred_grid$predictor
  
  # Predict the effect of the smooth term
  pred_effect <- predict(gam_model, newdata = newdata, type = "terms", terms = term_label)
  
  # Reshape to matrix
  effect_matrix <- matrix(pred_effect, nrow = length(predictor_vals), ncol = length(t_vals))
  
  # Plot the surface using plotly
  library(plotly)
  
  plot_ly(
    x = predictor_vals,
    y = t_vals,
    z = effect_matrix,
    type = "surface"
  ) %>%
    layout(
      title = paste("Estimated Beta(s, t) Surface for", predictor_name),
      scene = list(
        xaxis = list(title = paste(predictor_name, "(s)")),
        yaxis = list(title = "Time (t)"),
        zaxis = list(title = "Beta(s, t)")
      )
    )
}

# Plot for each predictor
plot_beta_surface_gam(gam_model, "te(ipc,t)", "ipc")
plot_beta_surface_gam(gam_model, "te(fdi,t)", "fdi")
plot_beta_surface_gam(gam_model, "te(ipc_fdi,t)", "ipc_fdi")

###############################################################################
# Section 7: Statistical Analysis and Inference
###############################################################################

# Compute residuals from the functional regression
residuals_fd <- grw_fd - fRegress_result$yhatfdobj
residuals_matrix <- eval.fd(Time, residuals_fd)

# Estimate the variance of the residuals
sigma_squared <- var(as.vector(residuals_matrix))

# Create the variance-covariance matrix
n <- length(Time)
SigmaE <- sigma_squared * diag(n)

# Compute standard errors of the estimated beta functions
stderr_list <- fRegress.stderr(fRegress_result, y2cMap = grw_y2c, SigmaE = SigmaE)

# For each coefficient, compute t-statistics and p-values
p_values_list <- list()  # To store p-values for each coefficient

for (i in 1:length(fRegress_result$betaestlist)) {
  beta_fd <- fRegress_result$betaestlist[[i]]$fd
  beta_stderr_fd <- stderr_list$betastderrlist[[i]]
  
  # Evaluate functions at time points
  beta_values <- eval.fd(Time, beta_fd)
  beta_stderr_values <- eval.fd(Time, beta_stderr_fd)
  
  # Compute t-statistics
  t_values <- beta_values / beta_stderr_values
  
  # Compute p-values
  p_values <- 2 * (1 - pnorm(abs(t_values)))
  
  # Store p-values
  p_values_list[[i]] <- p_values
  
  # Plot p-values over time
  plot(Time, p_values, type = 'l', main = paste("P-values for Beta", i - 1), ylab = "P-value")
  # Add red line where p-value is 0.1
  abline(h = 0.1, col = "red", lty = 2)
}

# Load regional data
regional <- read.csv("Country_Income_Region_Info_Completed_excel.csv", stringsAsFactors = FALSE)

# After data cleaning, get the list of countries
countries <- colnames(grw)

# Update regional data to include only the countries present
regional <- regional[regional$Country %in% countries, ]

# Ensure the order matches
regional <- regional[match(countries, regional$Country), ]

# Extract regions and income levels
regions <- regional$Region
Income.Range <- regional$Income.level

# Perform depth analysis
library(fda.usc)
ops.fda.usc()

depth_result <- depth.FM(grw_fd)

# Get unique regions
unique_regions <- unique(regions)

# Perform Wilcoxon tests for depth differences between regions
wilcox_results_region <- list()

for (region in unique_regions) {
  indices_region <- which(regions == region)
  indices_other <- which(regions != region)
  
  depth_region <- depth_result$dep[indices_region]
  depth_other <- depth_result$dep[indices_other]
  
  # Check for sufficient data
  if (length(depth_region) >= 1 && length(depth_other) >= 1) {
    test_result <- wilcox.test(depth_region, depth_other)
    wilcox_results_region[[region]] <- test_result
  } else {
    wilcox_results_region[[region]] <- NA
  }
}

# View Wilcoxon test results for regions
print(wilcox_results_region)

# Perform Wilcoxon tests for depth differences between income levels
unique_Income.Range <- unique(regional$Income.Range)
Income = regional$Income.Range
wilcox_results_income <- list()

for (level in unique_Income.Range) {
  indices_level <- which(Income == level)
  indices_other <- which(Income != level)
  
  depth_level <- depth_result$dep[indices_level]
  depth_other <- depth_result$dep[indices_other]
  
  # Check for sufficient data
  if (length(depth_level) >= 1 && length(depth_other) >= 1) {
    test_result <- wilcox.test(depth_level, depth_other)
    wilcox_results_income[[level]] <- test_result
  } else {
    wilcox_results_income[[level]] <- NA
  }
}

# View Wilcoxon test results for income levels
print(wilcox_results_income)

# Plot coefficients with confidence bands
for (i in 1:length(fRegress_result$betaestlist)) {
  beta_fd <- fRegress_result$betaestlist[[i]]$fd
  beta_stderr_fd <- stderr_list$betastderrlist[[i]]
  
  beta_values <- eval.fd(Time, beta_fd)
  stderr_values <- eval.fd(Time, beta_stderr_fd)
  
  # Compute confidence bands
  upper_band <- beta_values + 2 * stderr_values
  lower_band <- beta_values - 2 * stderr_values
  
  # Plot beta function and confidence bands
  plot(Time, beta_values, type = 'l', lwd = 2, main = paste("Beta", i - 1, "with 2-SE Confidence Band"), ylab = "Coefficient")
  lines(Time, upper_band, col = 'blue', lty = 2)
  lines(Time, lower_band, col = 'blue', lty = 2)
  legend("topright", legend = c("Estimate", "Confidence Band"), col = c("black", "blue"), lty = c(1, 2))
}

# Boxplot visualization of depth values by region and income level
regional <- read.csv("Country_Income_Region_Info_Completed_excel.csv", stringsAsFactors = FALSE)

regional <- regional[regional$Country %in% countries, ]
regional <- regional[match(countries, regional$Country), ]

# Extract region and income level vectors
region <- regional$Region
income_level <- regional$Income.Range

b1 <- boxplot(grw_fd, method = "MBD")
DM <- b1$depth

# Create data frame for plotting
depth_data <- data.frame(
  Depth = DM,
  Region = factor(region),
  IncomeLevel = factor(income_level)
)

# Load ggplot2 for plotting
library(ggplot2)

# Boxplot of depth by region
ggplot(depth_data, aes(x = Region, y = Depth)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "Depth Values by Region", x = "Region", y = "Depth") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot of depth by income level
ggplot(depth_data, aes(x = IncomeLevel, y = Depth)) +
  geom_boxplot(fill = "lightgreen") +
  labs(title = "Depth Values by Income Level", x = "Income Level", y = "Depth") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Use ggpubr for adding statistical comparisons
library(ggpubr)

# For regions
ggboxplot(depth_data, x = "Region", y = "Depth", fill = "Region") +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  labs(title = "Depth Values by Region", x = "Region", y = "Depth") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# For income levels
ggboxplot(depth_data, x = "IncomeLevel", y = "Depth", fill = "IncomeLevel") +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  labs(title = "Depth Values by Income Level", x = "Income Level", y = "Depth") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Compute and plot acceleration (second derivative) for GRW data
accel_fd <- deriv.fd(grw_fd, deriv = 2)
accel_values <- eval.fd(Time, accel_fd)

# Plot acceleration curves for all countries
matplot(Time, accel_values, type = 'l', lty = 1, col = 'grey',
        xlab = 'Years', ylab = 'Acceleration', main = 'Acceleration Curves for All Countries (GRW)')

# Compute and plot the mean acceleration curve
mean_accel <- rowMeans(accel_values)
lines(Time, mean_accel, col = 'red', lwd = 2)
legend("topright", legend = c("Mean Acceleration"), col = c("red"), lty = 1, lwd = 2)

# Compute and plot acceleration for FDI data
accel_fd <- deriv.fd(fdi_fd, deriv = 2)
accel_values <- eval.fd(Time, accel_fd)
matplot(Time, accel_values, type = 'l', lty = 1, col = 'grey',
        xlab = 'Years', ylab = 'Acceleration', main = 'Acceleration Curves for All Countries (FDI)')
mean_accel <- rowMeans(accel_values)
lines(Time, mean_accel, col = 'red', lwd = 2)
legend("topright", legend = c("Mean Acceleration"), col = c("red"), lty = 1, lwd = 2)

# Compute and plot acceleration for IPC data
accel_fd <- deriv.fd(ipc_fd, deriv = 2)
accel_values <- eval.fd(Time, accel_fd)
matplot(Time, accel_values, type = 'l', lty = 1, col = 'grey',
        xlab = 'Years', ylab = 'Acceleration', main = 'Acceleration Curves for All Countries (IPC)')
mean_accel <- rowMeans(accel_values)
lines(Time, mean_accel, col = 'red', lwd = 2)
legend("topright", legend = c("Mean Acceleration"), col = 'red', lty = 1, lwd = 2)

# Compute mean acceleration for each country
mean_accel_countries <- colMeans(accel_values)

# Create data frame with mean acceleration and grouping variables
accel_data <- data.frame(
  Country = countries,
  MeanAcceleration = mean_accel_countries,
  Region = factor(region),
  IncomeLevel = factor(income_level)
)

# Plot mean acceleration by region
ggplot(accel_data, aes(x = Region, y = MeanAcceleration)) +
  geom_boxplot(fill = "orange") +
  labs(title = "Mean Acceleration by Region", x = "Region", y = "Mean Acceleration") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot mean acceleration by income level
ggplot(accel_data, aes(x = IncomeLevel, y = MeanAcceleration)) +
  geom_boxplot(fill = "purple") +
  labs(title = "Mean Acceleration by Income Level", x = "Income Level", y = "Mean Acceleration") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Perform ANOVA for region
anova_accel_region <- aov(MeanAcceleration ~ Region, data = accel_data)
summary(anova_accel_region)

# Perform ANOVA for income level
anova_accel_income <- aov(MeanAcceleration ~ IncomeLevel, data = accel_data)
summary(anova_accel_income)

###############################################################################
# Section 8: Functional Regression by Regions and Income Levels
###############################################################################

# Function to perform regression for a subset of countries
perform_regression <- function(country_indices, group_name) {
  # Subset the functional data objects
  grw_fd_subset <- grw_fd[country_indices]
  ipc_fd_subset <- ipc_fd[country_indices]
  fdi_fd_subset <- fdi_fd[country_indices]
  
  # Recompute the interaction term
  ipc_eval_subset <- eval.fd(Time, ipc_fd_subset)
  fdi_eval_subset <- eval.fd(Time, fdi_fd_subset)
  ipc_fdi_interaction_subset <- ipc_eval_subset * fdi_eval_subset
  ipc_fdi_fd_subset <- smooth.basis(Time, ipc_fdi_interaction_subset, ipc_basis)$fd
  
  # Define the intercept
  const_basis <- create.constant.basis(rangeval)
  const_fd <- fd(matrix(1, 1, length(country_indices)), const_basis)
  
  # Define the list of predictors
  X_list_subset <- list(
    const = const_fd,
    ipc_fd = ipc_fd_subset,
    fdi_fd = fdi_fd_subset,
    ipc_fdi_fd = ipc_fdi_fd_subset
  )
  
  # Define basis for beta functions
  beta_basis <- create.fourier.basis(rangeval, nbasis)
  betalist <- list(
    beta_basis_const = fdPar(const_basis),
    ipc_fd = fdPar(beta_basis),
    fdi_fd = fdPar(beta_basis),
    ipc_fdi_fd = fdPar(beta_basis)
  )
  
  # Perform the functional regression
  fRegress_result <- fRegress(grw_fd_subset, X_list_subset, betalist)
  
  # Extract beta functions
  beta_estimates <- fRegress_result$betaestlist
  
  # Return results
  list(
    fRegress_result = fRegress_result,
    beta_estimates = beta_estimates,
    group_name = group_name
  )
}

# Get unique regions
unique_regions <- levels(regions)

# Initialize a list to store regression results
regression_results_regions <- list()

for (region in unique_regions) {
  # Get indices of countries in this region
  country_indices <- which(regions == region)
  
  # Check if there are enough countries
  if (length(country_indices) >= 2) {
    result <- perform_regression(country_indices, region)
    regression_results_regions[[region]] <- result
  } else {
    cat("Not enough countries in region:", region, "\n")
  }
}

# Get unique income levels
unique_income_levels <- levels(income_levels)

# Initialize a list to store regression results
regression_results_income <- list()

for (level in unique_income_levels) {
  # Get indices of countries in this income level
  country_indices <- which(income_levels == level)
  
  # Check if there are enough countries
  if (length(country_indices) >= 2) {
    result <- perform_regression(country_indices, level)
    regression_results_income[[level]] <- result
  } else {
    cat("Not enough countries in income level:", level, "\n")
  }
}

# Plot interaction beta functions for regions
par(mfrow = c(2, 2))  # Adjust layout as needed

for (region in names(regression_results_regions)) {
  beta_estimates <- regression_results_regions[[region]]$beta_estimates
  beta3_fd <- beta_estimates[[4]]$fd  # Interaction term
  
  plot(beta3_fd, main = paste("Beta 3 (Interaction) for", region))
}

par(mfrow = c(1, 1))  # Reset layout

# Plot interaction beta functions for income levels
par(mfrow = c(2, 2))  # Adjust layout as needed

for (level in names(regression_results_income)) {
  beta_estimates <- regression_results_income[[level]]$beta_estimates
  beta3_fd <- beta_estimates[[4]]$fd  # Interaction term
  
  plot(beta3_fd, main = paste("Beta 3 (Interaction) for", level))
}

par(mfrow = c(1, 1))  # Reset layout

beta3_fd <- fRegress_result$betaestlist[[4]]$fd  # Interaction term beta3
beta3_values <- eval.fd(Time, beta3_fd)

###############################################################################
# Section 9: Principal Component Analysis (PCA)
###############################################################################

# Perform PCA on the registered GRW data
pca_grw <- pca.fd(grw_fd_registered$registered$regfd, nharm = 3)

# Plot the principal component functions (harmonics)
plot(pca_grw$harmonics, main = "Principal Components of Registered GRW Data")

# Display the proportion of variance explained by each component
print(pca_grw$varprop)

# Perform PCA on the registered IPC data
pca_ipc <- pca.fd(ipc_fd_registered$registered$regfd, nharm = 3)
plot(pca_ipc$harmonics, main = "Principal Components of Registered IPC Data")
print(pca_ipc$varprop)

# Perform PCA on the registered FDI data
pca_fdi <- pca.fd(fdi_fd_registered$registered$regfd, nharm = 3)
plot(pca_fdi$harmonics, main = "Principal Components of Registered FDI Data")
print(pca_fdi$varprop)

###############################################################################
# Section 10: Clustering Analysis
###############################################################################

# Analyze the IPC * FDI interaction term from 'gam' model

# Extract the smooth terms from the model
summary(gam_model)$s.table

# Get the estimated coefficients for the interaction term
interaction_term_label <- "te(ipc_fdi,t)"
interaction_smooth <- gam_model$smooth[[which(sapply(gam_model$smooth, function(x) x$label) == interaction_term_label)]]

# Generate a grid over s (ipc_fdi) and t (time)
s_vals <- seq(min(data_long$ipc_fdi), max(data_long$ipc_fdi), length.out = 100)
t_vals <- finegrid  # Time grid

# Create a grid for prediction
pred_grid <- expand.grid(ipc_fdi = s_vals, t = t_vals)
pred_grid$ipc <- 0
pred_grid$fdi <- 0

# Predict the effect of the interaction term
interaction_effect <- predict(
  gam_model,
  newdata = pred_grid,
  type = "terms",
  terms = interaction_term_label
)

# Reshape the predicted effects into a matrix
interaction_matrix <- matrix(interaction_effect, nrow = length(s_vals), ncol = length(t_vals))

# Identify where the interaction effect is negative
negative_indices <- interaction_matrix < 0

# Calculate the proportion of negative values
proportion_negative <- sum(negative_indices) / length(interaction_matrix)
print(paste("Proportion of negative beta(s, t) values:", round(proportion_negative, 4)))

# Plot the regions where beta(s, t) is negative
library(plotly)

plot_ly(
  x = s_vals,
  y = t_vals,
  z = interaction_matrix,
  type = "surface",
  surfacecolor = interaction_matrix < 0,
  colorscale = list(c(0, 1), c("red", "blue")),
  showscale = FALSE
) %>%
  layout(
    title = "Beta(s, t) Surface for IPC * FDI Interaction (Blue regions are negative)",
    scene = list(
      xaxis = list(title = "IPC * FDI (s)"),
      yaxis = list(title = "Time (t)"),
      zaxis = list(title = "Beta(s, t)")
    )
  )

# Number of countries and time points
n_countries <- nrow(Y_mat)
n_time <- ncol(Y_mat)

# Initialize a matrix to store the interaction effects for each country over time
country_interaction_effects <- matrix(0, nrow = n_countries, ncol = n_time)

for (i in 1:n_countries) {
  # Extract the predictor function for the interaction term for country i
  ipc_fdi_i <- ipc_fdi_mat[i, ]
  
  # For each time point t, compute the effect
  for (t_idx in 1:n_time) {
    # For the predictor value ipc_fdi_i at s and time t
    s_val <- ipc_fdi_i
    t_val <- finegrid[t_idx]
    
    # Create a data frame for prediction
    pred_data <- data.frame(
      ipc_fdi = s_val,
      t = rep(t_val, length(s_val)),
      ipc = rep(0, length(s_val)),
      fdi = rep(0, length(s_val))
    )
    
    # Predict the effect
    effect <- predict(
      gam_model,
      newdata = pred_data,
      type = "terms",
      terms = interaction_term_label
    )
    
    # Integrate over s (assuming equal spacing)
    delta_s <- s_vals[2] - s_vals[1]
    country_interaction_effects[i, t_idx] <- sum(effect) * delta_s
  }
}

# Compute summary statistics for each country
country_negative_proportions <- apply(country_interaction_effects < 0, 1, mean)
countries <- countries[countries != "United.States"]

country_summary <- data.frame(
  Country = countries,
  NegativeProportion = country_negative_proportions
)

# Perform PCA on the interaction effects
library(refund)
library(fda)

# Create a functional data object from the country interaction effects
interaction_fd <- Data2fd(argvals = finegrid, y = t(country_interaction_effects))

# Perform FPCA
pca_results <- pca.fd(interaction_fd, nharm = 3)

# Plot the principal components
plot(pca_results$harmonics, main = "Principal Components of Interaction Effects")

# Proportion of variance explained
print(pca_results$varprop)

# Get the scores for each country
country_scores <- pca_results$scores

# Combine scores with country summaries
country_data <- data.frame(
  Country = countries,
  NegativeProportion = country_negative_proportions,
  PC1 = country_scores[, 1],
  PC2 = country_scores[, 2],
  PC3 = country_scores[, 3]
)

# Determine the number of clusters (e.g., 3)
set.seed(123)
num_clusters <- 3

# Perform K-means clustering
kmeans_result <- kmeans(country_data[, c("NegativeProportion", "PC1", "PC2", "PC3")], centers = num_clusters)

# Add cluster assignments to the data frame
country_data$Cluster <- factor(kmeans_result$cluster)

library(ggplot2)

# Plot the clusters using the first two principal components
ggplot(country_data, aes(x = PC1, y = PC2, color = Cluster, label = Country)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, size = 3) +
  labs(title = "Country Clusters Based on Interaction Effects") +
  theme_minimal()

# Load region and income level data
regional <- read.csv("Country_Income_Region_Info_Completed_excel.csv", stringsAsFactors = FALSE)

# Ensure alignment
regional <- regional[match(country_data$Country, regional$Country), ]

# Add region and income level to country_data
country_data$Region <- factor(regional$Region)
country_data$IncomeLevel <- factor(regional$Income.Range)  # Adjust the column name as needed

# Cross-tabulation of clusters by region
table_clusters_region <- table(country_data$Cluster, country_data$Region)
print(table_clusters_region)

# Cross-tabulation of clusters by income level
table_clusters_income <- table(country_data$Cluster, country_data$IncomeLevel)
print(table_clusters_income)

# Plot clusters colored by region
ggplot(country_data, aes(x = PC1, y = PC2, color = Region, shape = Cluster)) +
  geom_point(size = 3) +
  labs(title = "Country Clusters Colored by Region") +
  theme_minimal()

# Plot clusters colored by income level
ggplot(country_data, aes(x = PC1, y = PC2, color = IncomeLevel, shape = Cluster)) +
  geom_point(size = 3) +
  labs(title = "Country Clusters Colored by Income Level") +
  theme_minimal()
