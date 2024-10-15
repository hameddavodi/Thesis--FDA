library(fda)

# Provide directory of the RData
setwd("~/Desktop/New Folder With Items")

# Load data
grw = read.csv("GDP g.csv")
ipc = read.csv("GDP y.csv")
fdi = read.csv("FDI.csv")

# Convert to matrix and remove headers/rows if necessary
grw = as.matrix(grw)
ipc = as.matrix(ipc[-1, ])
fdi = as.matrix(fdi[-1, ])

# Time sequence
Time = 1980:2020

# Define a Fourier basis for functional data smoothing
rangeval <- range(Time)
nbasis <- 35  # Adjust as needed 21 was correct

# Create the basis for each functional object
grw_basis <- create.fourier.basis(rangeval, nbasis)
ipc_basis <- create.fourier.basis(rangeval, nbasis)
fdi_basis <- create.fourier.basis(rangeval, nbasis)

# Convert data to functional objects using smooth.basis
grw_fd <- smooth.basis(Time, grw, grw_basis)$fd
ipc_fd <- smooth.basis(Time, ipc, ipc_basis)$fd
fdi_fd <- smooth.basis(Time, fdi, fdi_basis)$fd

grw_y2c <- smooth.basis(Time, grw, grw_basis)$y2cMap
ipc_y2c <- smooth.basis(Time, ipc, ipc_basis)$y2cMap
fdi_y2c <- smooth.basis(Time, fdi, fdi_basis)$y2cMap

# Check dimensions and ensure all functional objects are aligned
print(dim(eval.fd(Time, grw_fd)))
print(dim(eval.fd(Time, ipc_fd)))
print(dim(eval.fd(Time, fdi_fd)))

# Create the interaction term pointwise
ipc_eval <- eval.fd(Time, ipc_fd)
fdi_eval <- eval.fd(Time, fdi_fd)
ipc_fdi_interaction <- ipc_eval * fdi_eval  # Pointwise multiplication

# Convert the interaction back to a functional object
ipc_fdi_fd <- smooth.basis(Time, ipc_fdi_interaction, ipc_basis)$fd

# Define the intercept as a functional data object
const_basis <- create.constant.basis(rangeval)
# Define the intercept as a functional data object
const_fd <- fd(matrix(1, 1, 94), const_basis)


# Define the list of predictors (including the interaction term)
X_list <- list(
  const = const_fd, 
  ipc_fd = ipc_fd,               # ipc predictor
  fdi_fd = fdi_fd,               # fdi predictor
  ipc_fdi_fd = ipc_fdi_fd        # ipc*fdi interaction predictor
)

# Define basis for beta functions
beta_basis <- create.fourier.basis(rangeval, nbasis)
betalist <- list(
  beta_basis_const <- const_basis,
  ipc_fd = fdPar(beta_basis),  # Beta for ipc
  fdi_fd = fdPar(beta_basis),  # Beta for fdi
  ipc_fdi_fd = fdPar(beta_basis)  # Beta for ipc*fdi interaction
)

# Perform the functional regression
fRegress_result <- fRegress(grw_fd, X_list, betalist)

# Extract and plot the beta functions
beta_estimates <- fRegress_result$betaestlist

# Plot beta functions
plot(beta_estimates[[1]], main = "Beta 0 (Intercept)")
plot(beta_estimates[[2]], main = "Beta 1 (ipc)")
plot(beta_estimates[[3]], main = "Beta 2 (fdi)")
plot(beta_estimates[[4]], main = "Beta 3 (ipc*fdi interaction)")



##################  


# Load country information
regional <- read.csv("Country_Income_Region_Info_Completed_excel.csv", stringsAsFactors = FALSE)

# Extract country names, income levels, and regions
countries <- regional$Country
income_levels <- regional$Income.level
regions <- regional$Region

# Calculate the second derivative of the growth rate functional data
grw_fd_acceleration <- deriv.fd(grw_fd, deriv=2)
# Plot the acceleration data
plot(grw_fd_acceleration, main = "Second Derivative of Growth Rate (Acceleration)")

# Calculate the second derivative of the growth rate functional data
ipc_fd_acceleration <- deriv.fd(ipc_fd, deriv=2)
# Plot the acceleration data
plot(ipc_fd_acceleration, main = "Second Derivative of Income Per Capita (Acceleration)")

# Calculate the second derivative of the growth rate functional data
fdi_fd_acceleration <- deriv.fd(fdi_fd, deriv=2)
# Plot the acceleration data
plot(fdi_fd_acceleration, main = "Second Derivative of Growth Rate (Acceleration)")


# Define a more flexible basis for the warping functions
Wnbasis <- 7  # You can adjust this number
Wbasis <- create.bspline.basis(rangeval = range(Time), nbasis = Wnbasis)
Wfd0 <- fd(matrix(0, Wnbasis, 1), Wbasis)

# Set up the functional parameter object for warping functions
WfdParobj <- fdPar(Wfd0, Lfdobj = 2, lambda = 0.01)

# Define a functional parameter object with adjusted lambda
fdi_fdPar <- fdPar(fdi_basis, Lfdobj = 2, lambda = 0.01)  # Adjust lambda as needed
fdi_fd <- smooth.basis(Time, fdi, fdi_fdPar)$fd

# Remove columns with NA values
fdi <- fdi[, colSums(is.na(fdi)) == 0]

# Remove constant functions (if any)
fdi_variances <- apply(fdi, 2, var)
fdi <- fdi[, fdi_variances > 0]

# First registration
fdi_fd_registered1 <- register.fd(yfd = fdi_fd, y0fd = mean.fd(fdi_fd), WfdParobj = WfdParobj)

# Update the target function to the mean of registered curves
fdi_mean_registered <- mean.fd(fdi_fd_registered1$regfd)

# Second registration using the updated mean
fdi_fd_registered2 <- register.fd(yfd = fdi_fd, y0fd = fdi_mean_registered, WfdParobj = WfdParobj)


# Plot registered FDI functions
plot(fdi_fd_registered2$regfd, main = "Registered FDI Functions")

# Plot the new mean function after registration
plot(mean.fd(fdi_fd_registered2$regfd), main = "Mean Function After Registration")

# Plot warping functions
plot(fdi_fd_registered2$warpfd, main = "Warping Functions")

#####
# Define a functional parameter object with adjusted lambda
ipc_fdPar <- fdPar(ipc_basis, Lfdobj = 2, lambda = 0.001)  # Adjust lambda as needed
ipc_fd <- smooth.basis(Time, ipc, ipc_fdPar)$fd

# Remove columns with NA values
ipc <- ipc[, colSums(is.na(ipc)) == 0]

# Remove constant functions (if any)
ipc_variances <- apply(ipc, 2, var)
ipc <- ipc[, ipc_variances > 0]

# First registration
ipc_fd_registered1 <- register.fd(yfd = ipc_fd, y0fd = mean.fd(ipc_fd), WfdParobj = WfdParobj)

# Update the target function to the mean of registered curves
ipc_mean_registered <- mean.fd(ipc_fd_registered1$regfd)

# Second registration using the updated mean
ipc_fd_registered2 <- register.fd(yfd = ipc_fd, y0fd = ipc_mean_registered, WfdParobj = WfdParobj)


# Plot registered ipc functions
plot(ipc_fd_registered2$regfd, main = "Registered ipc Functions")

# Plot the new mean function after registration
plot(mean.fd(ipc_fd_registered2$regfd), main = "Mean Function After Registration")

# Plot warping functions
plot(ipc_fd_registered2$warpfd, main = "Warping Functions")
#####
# Define a functional parameter object with adjusted lambda
grw_fdPar <- fdPar(grw_basis, Lfdobj = 2, lambda = 0.01)  # Adjust lambda as needed
grw_fd <- smooth.basis(Time, grw, grw_fdPar)$fd

# Remove columns with NA values
grw <- grw[, colSums(is.na(grw)) == 0]

# Remove constant functions (if any)
grw_variances <- apply(grw, 2, var)
grw <- grw[, grw_variances > 0]

# First registration
grw_fd_registered1 <- register.fd(yfd = grw_fd, y0fd = mean.fd(grw_fd), WfdParobj = WfdParobj)

# Update the target function to the mean of registered curves
grw_mean_registered <- mean.fd(grw_fd_registered1$regfd)

# Second registration using the updated mean
grw_fd_registered2 <- register.fd(yfd = grw_fd, y0fd = grw_mean_registered, WfdParobj = WfdParobj)


# Plot registered grw functions
plot(grw_fd_registered2$regfd, main = "Registered grw Functions")

# Plot the new mean function after registration
plot(mean.fd(grw_fd_registered2$regfd), main = "Mean Function After Registration")

# Plot warping functions
plot(grw_fd_registered2$warpfd, main = "Warping Functions")
grw_fd_registered = grw_fd_registered2
#####


# Perform PCA on the registered data
pca_grw <- pca.fd(grw_fd_registered$regfd, nharm = 6)
# Plot the principal component functions (harmonics)
plot(pca_grw$harmonics, main = "Principal Components of Registered Growth Rate")
# Display the proportion of variance explained by each component
print(pca_grw$varprop)



# Calculate residuals
# Compute residuals from the regression
residuals_fd <- grw_fd - fRegress_result$yhatfdobj
residuals_matrix <- eval.fd(Time, residuals_fd)

# Estimate the variance of the residuals
sigma_squared <- var(as.vector(residuals_matrix))

# Create the variance-covariance matrix (assuming independence)
n <- length(Time)
SigmaE <- sigma_squared * diag(n)


# Compute Mean Squared Error (MSE)
mse <- mean((eval.fd(Time, residuals_fd))^2)
print(paste("Mean Squared Error:", mse))

# Compute R-squared value
ss_total <- sum((eval.fd(Time, grw_fd) - mean(eval.fd(Time, grw_fd)))^2)
ss_residual <- sum((eval.fd(Time, residuals_fd))^2)
r_squared <- 1 - (ss_residual / ss_total)
print(paste("R-squared:", r_squared))


# Compute standard errors of regression coefficients
stderr_list <- fRegress.stderr(fRegress_result, y2cMap = grw_y2c, SigmaE=SigmaE)


# For each coefficient, compute t-statistics and p-values
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
  
  # Plot p-values over time
  plot(Time, p_values, type = 'l', main = paste("P-values for Beta", i - 1), ylab = "P-value")
}

#############


library(fda)
library(fda.usc)
# Perform depth analysis using Modified Band Depth (MBD)
depth_result <- depth.FM(grw_fd)
# Plot depth values
plot(depth_result$dep, type = 'h', main = "Functional Depth of Growth Rate Functions", xlab = "Country Index", ylab = "Depth")
# Identify the median curve (highest depth)
median_index <- which.max(depth_result$dep)
median_curve <- grw_fd[median_index]

# Plot all curves with the median curve highlighted
plot(grw_fd, col = "grey", main = "Growth Rate Functions with Median Curve")
lines(median_curve, col = "red", lwd = 2)

# Identify potential outliers (lowest depth)
outlier_threshold <- quantile(depth_result$dep, 0.1)  # Lower 10% as outliers
outlier_indices <- which(depth_result$dep <= outlier_threshold)

# Highlight outlier curves
for (i in outlier_indices) {
  lines(grw_fd[i], col = "blue", lty = 2)
}

legend("topright", legend = c("Median Curve", "Outliers"), col = c("red", "blue"), lty = c(1, 2), cex = 0.8)



# Use your growth rate functional data object
temp_fd <- grw_fd

# Plot the functional data
plot(temp_fd, main = "Growth Rate Functional Data")

# Perform functional boxplot using MBD
boxplot_result <- boxplot(temp_fd, method = "MBD", main = "Functional Boxplot of Growth Rates")
boxplot_result <- boxplot(temp_fd, method = "BD2", main = "Functional Boxplot of Growth Rates")
boxplot_result <- boxplot(temp_fd, method = "Both", main = "Functional Boxplot of Growth Rates")




# Get unique regions
unique_regions <- unique(regions)

# Perform Wilcoxon tests
wilcox_results_region <- list()

for (region in unique_regions) {
  indices_region <- which(regions == region)
  indices_other <- which(regions != region)
  
  depth_region <- depth_result$dep[indices_region]
  depth_other <- depth_result$dep[indices_other]
  
  # Check for sufficient data
  if (length(depth_region) >= 2 && length(depth_other) >= 2) {
    test_result <- wilcox.test(depth_region, depth_other)
    wilcox_results_region[[region]] <- test_result
  } else {
    wilcox_results_region[[region]] <- NA
  }
}

# View Wilcoxon test results for regions
print(wilcox_results_region)


# Get unique income levels
unique_income_levels <- unique(income_levels)

# Perform Wilcoxon tests
wilcox_results_income <- list()

for (level in unique_income_levels) {
  indices_level <- which(income_levels == level)
  indices_other <- which(income_levels != level)
  
  depth_level <- depth_grw$dep[indices_level]
  depth_other <- depth_grw$dep[indices_other]
  
  # Check for sufficient data
  if (length(depth_level) >= 2 && length(depth_other) >= 2) {
    test_result <- wilcox.test(depth_level, depth_other)
    wilcox_results_income[[level]] <- test_result
  } else {
    wilcox_results_income[[level]] <- NA
  }
}

# View Wilcoxon test results for income levels
print(wilcox_results_income)



# FANOVA for regions
group_region <- as.factor(regions)
fanova_region <- fanova.onefactor(grw_fd, group_region)
print(fanova_region)

# FANOVA for income levels
group_income <- as.factor(income_levels)
fanova_income <- fanova.onefactor(grw_fd, group_income)
print(fanova_income)



# Compute standard errors
stderr_list <- fRegress.stderr(fRegress_result)

# Plot coefficients with confidence bands
for (i in 1:length(fRegress_result$betaestlist)) {
  beta_fd <- fRegress_result$betaestlist[[i]]$fd
  beta_stderr_fd <- stderr_list$betastderrlist[[i]]$fd
  
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


# Evaluate functional data at time points
grw_values <- eval.fd(Time, grw_fd)
ipc_values <- eval.fd(Time, ipc_fd)
fdi_values <- eval.fd(Time, fdi_fd)

# Example: Use average values over time for each country
mean_grw <- colMeans(grw_values)
mean_ipc <- colMeans(ipc_values)
mean_fdi <- colMeans(fdi_values)

# Create a data frame
df <- data.frame(IPC = mean_ipc, FDI = mean_fdi, GRW = mean_grw)

# 3D Scatter plot
library(scatterplot3d)
scatterplot3d(df$IPC, df$FDI, df$GRW, xlab = "IPC", ylab = "FDI", zlab = "GRW", main = "3D Scatter Plot of IPC, FDI, and GRW")


# Evaluate beta_4 function
beta_4_fd <- fRegress_result$betaestlist[[4]]$fd
beta_4_values <- eval.fd(Time, beta_4_fd)

# Compute interaction term for each country
ipc_values <- eval.fd(Time, ipc_fd)
fdi_values <- eval.fd(Time, fdi_fd)
interaction_term <- ipc_values * fdi_values

# Compute interaction effect
interaction_effect <- sweep(interaction_term, 1, beta_4_values, '*')

# Sum interaction effects over time for each country
total_effect <- colSums(interaction_effect)

# Classify countries based on the sign of the total effect
cluster_labels <- ifelse(total_effect > 0, "Positive", ifelse(total_effect < 0, "Negative", "Neutral"))

# Create a data frame with countries and their cluster labels
clustering_results <- data.frame(Country = countries, Cluster = cluster_labels)

# View clustering results
print(clustering_results)
