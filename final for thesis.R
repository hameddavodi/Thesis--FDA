# Start redirecting all output to "output.txt"
sink("output.txt")

# All output below will be saved to "output.txt"
cat("This is a cat() output.\n")
print("This is a print() output.")

###############################################################################
# Section 1: Data Loading and Preprocessing
###############################################################################

# Load necessary library
library(fda)
library(lme4)
library(refund)
library(plotly)
library(mgcv)

options(warn=-1)

# Set working directory to where the data files are located
setwd("C:/Users/MEDIA MARK/Desktop/thesis/")

# Load data from CSV files
grw = read.csv("GDP g.csv")  # Diff GDP growth data vs USA
ipc = read.csv("GDP y.csv")  # Diff GDP per capita data vs USA
fdi = read.csv("FDI.csv")  # Diff Financial Development index data vs USA

# Convert data frames to matrices and remove the last row if necessary
grw = as.matrix(grw)
ipc = as.matrix(ipc[-nrow(ipc), ])
fdi = as.matrix(fdi[-nrow(fdi), ])

# Define the time sequence from 1980 to 2020
Time = 1980:2020

# Define the range of the data
rangeval <- range(Time)

# Load raaw data for comparison
ipc_raw = read.csv("GDP.csv")  # GDP per capita data
fdi_raw = read.csv("FDI.csv")  # financial development data
ipc_raw = as.matrix(ipc_raw[-nrow(ipc_raw), ])
fdi_raw = as.matrix(fdi_raw[-nrow(fdi_raw), ])

dim(grw)
dim(ipc)
dim(fdi)
dim(ipc_raw)
dim(fdi_raw)
###############################################################################
# Section 2: Using GCV to Select Optimal Number of Basis Functions
###############################################################################

# Define a range of possible nbasis values to evaluate
nbasis_range <- 3:40  # You can adjust this range as needed

# Initialize vectors to store GCV values
grw_gcv <- numeric(length(nbasis_range))
ipc_gcv <- numeric(length(nbasis_range))
fdi_gcv <- numeric(length(nbasis_range))
ipc_raw_gcv <- numeric(length(nbasis_range))
fdi_raw_gcv <- numeric(length(nbasis_range))

# Loop over nbasis values to compute GCV for each data set
for (i in seq_along(nbasis_range)) {
  # Create Fourier basis for each nbasis value
  grw_basis <- create.fourier.basis(rangeval, nbasis_range[i])
  ipc_basis <- create.fourier.basis(rangeval, nbasis_range[i])
  fdi_basis <- create.fourier.basis(rangeval, nbasis_range[i])
  ipc_raw_basis <- create.fourier.basis(rangeval, nbasis_range[i])
  fdi_raw_basis <- create.fourier.basis(rangeval, nbasis_range[i])
  
  # Smooth the data using the defined basis and functional parameters
  grw_gcv[i] <- smooth.basis(Time, grw, fdPar(grw_basis, Lfdobj = 2, lambda = 0.0001))$gcv
  ipc_gcv[i] <- smooth.basis(Time, ipc, fdPar(ipc_basis, Lfdobj = 2, lambda = 0.0001))$gcv
  fdi_gcv[i] <- smooth.basis(Time, fdi, fdPar(fdi_basis, Lfdobj = 2, lambda = 0.0001))$gcv
  ipc_raw_gcv[i] <- smooth.basis(Time, ipc_raw, fdPar(ipc_raw_basis, Lfdobj = 2, lambda = 0.0001))$gcv
  fdi_raw_gcv[i] <- smooth.basis(Time, fdi_raw, fdPar(fdi_raw_basis, Lfdobj = 2, lambda = 0.0001))$gcv
}

# Plot GCV values to visualize the optimal nbasis
par(mfrow=c(5,1))  # Set up a 5x1 plotting layout

plot(nbasis_range, grw_gcv, type = "b", main = "GCV for GDP Growth", xlab = "Number of Basis Functions", ylab = "GCV")
plot(nbasis_range, ipc_gcv, type = "b", main = "GCV for GDP per Capita", xlab = "Number of Basis Functions", ylab = "GCV")
plot(nbasis_range, fdi_gcv, type = "b", main = "GCV for FDI", xlab = "Number of Basis Functions", ylab = "GCV")
plot(nbasis_range, ipc_raw_gcv, type = "b", main = "GCV for Raw GDP per Capita", xlab = "Number of Basis Functions", ylab = "GCV")
plot(nbasis_range, fdi_raw_gcv, type = "b", main = "GCV for Raw FDI", xlab = "Number of Basis Functions", ylab = "GCV")

par(mfrow=c(1,1))  # Reset layout

# Find the nbasis that minimizes GCV for each data set
optimal_grw_nbasis <- nbasis_range[which.min(grw_gcv)]
optimal_ipc_nbasis <- nbasis_range[which.min(ipc_gcv)]
optimal_fdi_nbasis <- nbasis_range[which.min(fdi_gcv)]
optimal_ipc_raw_nbasis <- nbasis_range[which.min(ipc_raw_gcv)]
optimal_fdi_raw_nbasis <- nbasis_range[which.min(fdi_raw_gcv)]

# Create Fourier basis functions with the optimal number of basis
grw_basis_opt <- create.fourier.basis(rangeval, optimal_grw_nbasis)
ipc_basis_opt <- create.fourier.basis(rangeval, optimal_ipc_nbasis)
fdi_basis_opt <- create.fourier.basis(rangeval, optimal_fdi_nbasis)
ipc_raw_basis_opt <- create.fourier.basis(rangeval, optimal_ipc_raw_nbasis)
fdi_raw_basis_opt <- create.fourier.basis(rangeval, optimal_fdi_raw_nbasis)

# Define functional parameter objects
grw_fdPar_opt <- fdPar(grw_basis_opt, Lfdobj = 2, lambda = 0.0001)
ipc_fdPar_opt <- fdPar(ipc_basis_opt, Lfdobj = 2, lambda = 0.0001)
fdi_fdPar_opt <- fdPar(fdi_basis_opt, Lfdobj = 2, lambda = 0.0001)
ipc_raw_fdPar_opt <- fdPar(ipc_raw_basis_opt, Lfdobj = 2, lambda = 0.0001)
fdi_raw_fdPar_opt <- fdPar(fdi_raw_basis_opt, Lfdobj = 2, lambda = 0.0001)

# Smooth the data using the optimal number of basis functions
grw_fd_opt <- smooth.basis(Time, grw, grw_fdPar_opt)$fd
ipc_fd_opt <- smooth.basis(Time, ipc, ipc_fdPar_opt)$fd
fdi_fd_opt <- smooth.basis(Time, fdi, fdi_fdPar_opt)$fd
ipc_raw_fd_opt <- smooth.basis(Time, ipc_raw, ipc_raw_fdPar_opt)$fd
fdi_raw_fd_opt <- smooth.basis(Time, fdi_raw, fdi_raw_fdPar_opt)$fd

# Calculate the mean of each optimized functional data object
grw_mean_fd_opt <- mean.fd(grw_fd_opt)
ipc_mean_fd_opt <- mean.fd(ipc_fd_opt)
fdi_mean_fd_opt <- mean.fd(fdi_fd_opt)
ipc_raw_mean_fd_opt <- mean.fd(ipc_raw_fd_opt)
fdi_raw_mean_fd_opt <- mean.fd(fdi_raw_fd_opt)

# Plot the mean functions
par(mfrow=c(5,1))  # Set up a 5x1 plotting layout

plot(grw_mean_fd_opt, main = "Mean of Optimized GDP Growth Functions", xlab = "Year", ylab = "GDP Growth")
plot(ipc_mean_fd_opt, main = "Mean of Optimized GDP per Capita Functions", xlab = "Year", ylab = "GDP per Capita")
plot(fdi_mean_fd_opt, main = "Mean of Optimized FDI Functions", xlab = "Year", ylab = "FDI")
plot(ipc_raw_mean_fd_opt, main = "Mean of Optimized Raw GDP per Capita Functions", xlab = "Year", ylab = "Raw GDP per Capita")
plot(fdi_raw_mean_fd_opt, main = "Mean of Optimized Raw FDI Functions", xlab = "Year", ylab = "Raw FDI")



###############################################################################
# Section 3: Summary Statistics for GRW, FDI, and IPC
###############################################################################

# Function to compute summary statistics
compute_summary_stats <- function(fd_obj, var_name) {
  # Evaluate the functional data over the time grid
  fd_values <- eval.fd(Time, fd_obj)
  
  # Compute statistics at each time point
  mean_values <- rowMeans(fd_values)
  median_values <- apply(fd_values, 1, median)
  var_values <- apply(fd_values, 1, var)
  sd_values <- sqrt(var_values)
  
  # Overall statistics (averaged over time)
  overall_mean <- mean(mean_values)
  overall_median <- median(median_values)
  overall_variance <- mean(var_values)
  overall_sd <- sqrt(overall_variance)
  
  # Compile into a data frame
  summary_df <- data.frame(
    Statistic = c("Overall Mean", "Overall Median", "Overall Variance", "Overall Standard Deviation"),
    Value = c(overall_mean, overall_median, overall_variance, overall_sd)
  )
  
  # Print the summary table
  cat("\nSummary Statistics for", var_name, ":\n")
  print(summary_df)
  
  # Plot the mean function with confidence bands
  # Compute 95% confidence intervals
  n_countries <- ncol(fd_values)
  se_values <- sd_values / sqrt(n_countries)
  ci_upper <- mean_values + 1.96 * se_values
  ci_lower <- mean_values - 1.96 * se_values
  
  # Create a data frame for plotting
  plot_df <- data.frame(
    Time = Time,
    Mean = mean_values,
    UpperCI = ci_upper,
    LowerCI = ci_lower
  )
  
  # Load ggplot2
  library(ggplot2)
  
  # Plot the mean function with confidence intervals
  ggplot(plot_df, aes(x = Time, y = Mean)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), alpha = 0.2, fill = "blue") +
    labs(title = paste("Mean Function with 95% CI for", var_name),
         x = "Time", y = var_name) +
    theme_minimal()
}

# Compute and display summary statistics for GRW
compute_summary_stats(grw_fd_opt, "GDP Growth (GRW)")

# Compute and display summary statistics for IPC
compute_summary_stats(ipc_fd_opt, "GDP per Capita (IPC)")

# Compute and display summary statistics for FDI
compute_summary_stats(fdi_fd_opt, "Financial Development Index (FDI)")

# Compute and display summary statistics for IPC
compute_summary_stats(ipc_raw_fd_opt, "GDP per Capita (IPC)")

# Compute and display summary statistics for FDI
compute_summary_stats(fdi_raw_fd_opt, "Financial Development Index (FDI)")

###############################################################################
# OR Defining Basis Functions for Functional Data Smoothing
###############################################################################

# Set the number of basis functions to use
nbasis <- 33  # Adjust this number as needed

# Create Fourier basis functions for each data set
grw_basis <- create.fourier.basis(rangeval, nbasis)  # For GDP growth
ipc_basis <- create.fourier.basis(rangeval, nbasis)  # For GDP per capita
fdi_basis <- create.fourier.basis(rangeval, nbasis)  # For FDI

# Define functional parameter objects with smoothing parameter lambda
grw_fdPar <- fdPar(grw_basis, Lfdobj = 2, lambda = 0.0001)  # For GDP growth
ipc_fdPar <- fdPar(ipc_basis, Lfdobj = 2, lambda = 0.0001)  # For GDP per capita
fdi_fdPar <- fdPar(fdi_basis, Lfdobj = 2, lambda = 0.0001)  # For FDI

###############################################################################
# Section 4: Smoothing the Functional Data
###############################################################################

# Smooth the data using the defined basis and functional parameters
grw_fd <- smooth.basis(Time, grw, grw_fdPar)$fd  # Smoothed GDP growth functions
ipc_fd <- smooth.basis(Time, ipc, ipc_fdPar)$fd  # Smoothed GDP per capita functions
fdi_fd <- smooth.basis(Time, fdi, fdi_fdPar)$fd  # Smoothed FDI functions

# Obtain the y to c mapping matrices (used later for inference)
grw_y2c <- smooth.basis(Time, grw, grw_basis)$y2cMap
ipc_y2c <- smooth.basis(Time, ipc, ipc_basis)$y2cMap
fdi_y2c <- smooth.basis(Time, fdi, fdi_basis)$y2cMap


# Load regional data and align with country data
regional <- read.csv("Category.csv")

# Replace periods (.) with spaces in the grw column names
colnames(grw) <- gsub("\\.", " ", colnames(grw))

# Now, match the country names with the regional data
countries <- colnames(grw)
regional_cleaned <- regional[regional$Country %in% countries, ]  # Filter to include only present countries
regional_cleaned <- regional_cleaned[match(countries, regional_cleaned$Country), ]  # Ensure the order matches

# Check if there are any NA values after the fix
sum(is.na(regional_cleaned$Country))  # This should be 0 if all matches are correct

###############################################################################
# Section 5: Registration of Functional Data
###############################################################################

# Define a more flexible basis for the warping functions
Wnbasis <- 7  # Number of basis functions for warping;
Wbasis <- create.bspline.basis(rangeval = range(Time), nbasis = Wnbasis)
Wfd0 <- fd(matrix(0, Wnbasis, 1), Wbasis)

# Set up the functional parameter object for warping functions
WfdParobj <- fdPar(Wfd0, Lfdobj = 2, lambda = 0.0001)

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
# Section 6: Functional Regression Analysis
###############################################################################

# Re-define the number of basis functions if needed
nbasis <- 33

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

# Extract the predicted values (yhat)
yhat_fd <- fRegress_result$yhatfdobj

# yhat_fd is a functional data object (fd class), and you can evaluate it over the time grid
Time <- seq(1980, 2020, length.out = 100)  # Adjust the time grid as per your data

# Evaluate the predicted values over the time grid
yhat_values <- eval.fd(Time, yhat_fd)

# Now you can work with the predicted values (yhat_values)
print(head(yhat_values))  # Print a preview of the predicted values

# Optionally, plot the predicted values
# matplot(Time, yhat_values, type = "l", main = "Predicted Yhat (GRW)", xlab = "Time", ylab = "Predicted Value")
# plot(grw_fd_registered$original, main = "Original GRW Functions")


# Plot beta functions
plot(beta_estimates[[1]], main = "Beta 0 (Intercept)")
plot(beta_estimates[[2]], main = "Beta 1 (IPC)")
plot(beta_estimates[[3]], main = "Beta 2 (FDI)")
plot(beta_estimates[[4]], main = "Beta 3 (IPC * FDI Interaction)")

###############################################################################
# Section 7: Functional Regression using 'refund' Package
###############################################################################

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

# plot(model_pffr)

qq.pffr(
  model_pffr,
  rep = 0,
  level = 0.9,
  s.rep = 10,
  type = c("deviance", "pearson", "response"),
  pch = ".",
  rl.col = 2,
  rep.col = "gray80")

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

qq.gam(
  gam_model,
  rep = 0,
  level = 0.9,
  s.rep = 10,
  type = c("deviance", "pearson", "response"),
  pch = ".",
  rl.col = 2,
  rep.col = "gray80")

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
# Section 8: Statistical Analysis and Inference
###############################################################################

# Load regional data
regional <- read.csv("Category.csv")

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
  }
}

# View Wilcoxon test results for income levels
print(wilcox_results_income)

# Boxplot visualization of depth values by region and income level
regional <- read.csv("Category.csv")

regional <- regional[regional$Country %in% countries, ]
regional <- regional[match(countries, regional$Country), ]

# Extract region and income level vectors
region <- regional$Region
income_level <- regional$Income.Range

b1 <- boxplot(grw_fd, method = "Both", xlab = "Time", ylab = "Grw")
b2 <- boxplot(ipc_fd, method = "Both", xlab = "Time", ylab = "Ipc")
b3 <- boxplot(fdi_fd, method = "Both", xlab = "Time", ylab = "Fdi")

# DM <- b1$depth


# Compute and plot acceleration (second derivative) for GRW data
accel_fd <- deriv.fd(ipc_fd, deriv = 2)
accel_values <- eval.fd(Time, accel_fd)

# Plot acceleration curves for all countries
matplot(Time, accel_values, type = 'l', lty = 1, col = 'grey',
        xlab = 'Years', ylab = 'Acceleration', main = 'Acceleration Curves for All Countries (GRW)')


# Compute and plot the mean acceleration curve
mean_accel <- rowMeans(accel_values)
lines(Time, mean_accel, col = 'red', lwd = 2)
legend("topright", legend = c("Mean Acceleration"), col = c("red"), lty = 1, lwd = 2)
par(mfrow = c(1, 1))  # Reset layout

# Compute and plot acceleration for FDI data
accel_fd <- deriv.fd(fdi_fd, deriv = 2)
accel_values <- eval.fd(Time, accel_fd)
matplot(Time, accel_values, type = 'l', lty = 1, col = 'grey',
        xlab = 'Years', ylab = 'Acceleration', main = 'Acceleration Curves for All Countries (FDI)')
mean_accel <- rowMeans(accel_values)
lines(Time, mean_accel, col = 'red', lwd = 2)
legend("topright", legend = c("Mean Acceleration"), col = c("red"), lty = 1, lwd = 2)
par(mfrow = c(1, 1))  # Reset layout

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
par(mfrow = c(1, 1))  # Reset layout

# Plot mean acceleration by income level
ggplot(accel_data, aes(x = IncomeLevel, y = MeanAcceleration)) +
  geom_boxplot(fill = "purple") +
  labs(title = "Mean Acceleration by Income Level", x = "Income Level", y = "Mean Acceleration") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
par(mfrow = c(1, 1))  # Reset layout


###############################################################################
# Section 9: Functional Regression by Regions and Income Levels
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

regional <- read.csv("Category.csv")

regional <- regional[regional$Country %in% countries, ]
regional <- regional[match(countries, regional$Country), ]

# Extract region and income level vectors
region <- regional$Region
income_level <- regional$Income.Range

# Get unique regions
unique_regions <- unique(region)

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
unique_income_levels <- unique(income_level)

# Initialize a list to store regression results
regression_results_income <- list()

for (level in unique_income_levels) {
  # Get indices of countries in this income level
  country_indices <- which(income_level == level)
  
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
# Section 10: Principal Component Analysis (PCA)
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
# Section 11: Clustering Analysis
###############################################################################

# Perform PCA on the interaction effects

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

# Plot the clusters using the first two principal components
ggplot(country_data, aes(x = PC1, y = PC2, color = Cluster, label = Country)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, size = 3) +
  labs(title = "Country Clusters Based on Interaction Effects") +
  theme_minimal()

# Load region and income level data
regional <- read.csv("Category.csv")

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


###############################################################################
# Section 12: Heatmaps of Beta Surfaces
###############################################################################

# Function to plot heatmap of beta surface
plot_beta_heatmap <- function(gam_model, term_label, predictor_name) {
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
  
  # Create a data frame for plotting
  beta_df <- data.frame(
    Predictor = rep(predictor_vals, times = length(t_vals)),
    Time = rep(t_vals, each = length(predictor_vals)),
    BetaValue = as.vector(effect_matrix)
  )
  
  # Plot heatmap
  ggplot(beta_df, aes(x = Time, y =Predictor , fill = BetaValue)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    labs(title = paste("Heatmap of Beta(s, t) for", predictor_name), x =  "Time (t)" , y =paste(predictor_name, "(s)")) +
    theme_minimal()
}

# Plot heatmaps for each predictor
plot_beta_heatmap(gam_model, "te(ipc,t)", "ipc")
plot_beta_heatmap(gam_model, "te(fdi,t)", "fdi")
plot_beta_heatmap(gam_model, "te(ipc_fdi,t)", "ipc_fdi")



























































































library(plotly)
library(ggplot2)
library(mgcv)
library(reshape2)

# Extract beta functions from fRegress_result
beta_estimates <- fRegress_result$betaestlist

# Beta functions
beta0_fd <- beta_estimates$const$fd       # Intercept
beta_ipc_fd <- beta_estimates$ipc_fd$fd   # Beta for IPC
beta_fdi_fd <- beta_estimates$fdi_fd$fd   # Beta for FDI
beta_inter_fd <- beta_estimates$ipc_fdi_fd$fd  # Beta for IPC * FDI interaction

# Define the time grid
Time <- seq(1980, 2020, length.out = 100)

# Evaluate beta functions over the time grid
beta_ipc_values <- eval.fd(Time, beta_ipc_fd)
beta_fdi_values <- eval.fd(Time, beta_fdi_fd)
beta_inter_values <- eval.fd(Time, beta_inter_fd)


# Load necessary library
library(plotly)

# Create data frame for plotting
plot_data <- data.frame(
  Time = Time,
  Beta_FDI = beta_fdi_values,
  Beta_Interaction = beta_inter_values
)


# Create scatter plot with horizontal line at y = 0
ggplot(heatmap_data, aes(x = Time, y = Beta_FDI, color = Beta_Interaction)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, color = "red", linetype = "solid", size = 0.8) +  # Adds red line at y = 0
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "Beta Interaction vs Beta FDI over Time",
    x = "Time",
    y = "Beta FDI",
    color = "Beta Interaction"
  ) +
  theme_minimal()


# Assuming beta_ipc_values, beta_inter_values, and Time are already defined in your environment
# Assuming beta_ipc_values, beta_inter_values, and Time are already defined in your environment

# Create data frame for plotting
plot_data <- data.frame(
  Time = Time,
  Beta_c = - (beta_ipc_values / beta_inter_values)
)

# Limit Beta_c to the boundaries of ±2
plot_data$Beta_c <- pmax(pmin(plot_data$Beta_c, 2), -2)

# Create a column to indicate positive or negative values
plot_data$Sign <- ifelse(plot_data$Beta_c >= 0, "Positive", "Negative")

# Load ggplot2 library
library(ggplot2)

# Plot Beta_c over Time with color based on positive or negative values (Scatter Plot)
ggplot(plot_data, aes(x = Time, y = Beta_c, color = Sign)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Positive" = "blue", "Negative" = "red")) +
  ylim(-2, 2) +
  labs(
    title = expression("Scatter Plot of " ~ F[c] ~ " over Time"),
    x = "Time",
    y = expression(F[c] == -beta[IPC]/beta[Interaction])
  ) +
  theme_minimal()
