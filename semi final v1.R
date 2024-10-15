library(fda)

# Provide directory of the RData
setwd("~/Desktop/New Folder With Items")

# Load data
grw = read.csv("GDP g.csv")
ipc = read.csv("GDP y.csv")
fdi = read.csv("FDI r.csv")



# Convert to matrix and remove headers/rows if necessary
grw = as.matrix(grw)

ipc = as.matrix(ipc[-nrow(ipc), ])

fdi = as.matrix(fdi[-nrow(fdi), ])



# Time sequence
Time = 1980:2020


# Define Fourier basis for functional data smoothing
rangeval <- range(Time)
nbasis <- 35  # Adjust as needed

# Create the basis for each functional object
grw_basis <- create.fourier.basis(rangeval, nbasis)
ipc_basis <- create.fourier.basis(rangeval, nbasis)
fdi_basis <- create.fourier.basis(rangeval, nbasis)

# Define functional parameter objects with adjusted lambda
grw_fdPar <- fdPar(grw_basis, Lfdobj = 2, lambda = 0.01)
ipc_fdPar <- fdPar(ipc_basis, Lfdobj = 2, lambda = 0.01)
fdi_fdPar <- fdPar(fdi_basis, Lfdobj = 2, lambda = 0.01)

# Smooth the cleaned data
grw_fd <- smooth.basis(Time, grw, grw_fdPar)$fd
ipc_fd <- smooth.basis(Time, ipc, ipc_fdPar)$fd
fdi_fd <- smooth.basis(Time, fdi, fdi_fdPar)$fd

# Get the country names from the cleaned data
countries <- colnames(grw)

regional <- read.csv("Country_Income_Region_Info_Completed_excel.csv", stringsAsFactors = FALSE)

# Filter the regional data to include only the countries present in the cleaned data
regional <- regional[regional$Country %in% countries, ]

# Ensure the order of countries matches
regional <- regional[match(countries, regional$Country), ]

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
grw_fdPar <- fdPar(grw_basis, Lfdobj = 2, lambda = 0.001)  # Adjust lambda as needed
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
####

nbasis <- 35  # Adjust as needed 21 was correct
#
#
# get back usa to grw data
# Load data
#


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

ipc_fdi_interaction <- ipc_eval * fdi_eval

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



# Smooth the grw data to get y2cMap
grw_smooth <- smooth.basis(Time, grw, grw_basis)
grw_fd <- grw_smooth$fd
grw_y2c <- grw_smooth$y2cMap

# Compute residuals
residuals_fd <- grw_fd - fRegress_result$yhatfdobj
residuals_matrix <- eval.fd(Time, residuals_fd)

# Estimate the variance of the residuals
sigma_squared <- var(as.vector(residuals_matrix))

# Create the variance-covariance matrix
n <- length(Time)
SigmaE <- sigma_squared * diag(n)
# Compute standard errors
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
}




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
depth_result <- depth.FM(grw_fd)

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
  if (length(depth_region) >= 1 && length(depth_other) >= 1) {
    test_result <- wilcox.test(depth_region, depth_other)
    wilcox_results_region[[region]] <- test_result
  } else {
    wilcox_results_region[[region]] <- NA
  }
}

# View Wilcoxon test results for regions
print(wilcox_results_region)

regional <- read.csv("Country_Income_Region_Info_Completed_excel.csv", stringsAsFactors = FALSE)


# After data cleaning, get the list of countries
countries <- colnames(grw)

# Update regional data to include only the countries present
regional <- regional[regional$Country %in% countries, ]

# Ensure the order matches
regional <- regional[match(countries, regional$Country), ]

# Get unique income levels
unique_Income.Range <- unique(regional$Income.Range)
Income = regional$Income.Range
# Perform Wilcoxon tests
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




# Compute standard errors (if not already computed)
# stderr_list <- fRegress.stderr(fRegress_result, y2cMap = grw_y2c, SigmaE = SigmaE)

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
ggplot(depth_data, aes(x = income_level, y = Depth)) +
  geom_boxplot(fill = "lightgreen") +
  labs(title = "Depth Values by Income Level", x = "Income Level", y = "Depth") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

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


# Assuming temp_fd is your functional data object for temperature
# Compute the second derivative (acceleration) for all countries
accel_fd <- deriv.fd(grw_fd, deriv = 2)

# Evaluate the acceleration at the time points
accel_values <- eval.fd(Time, accel_fd)

# Plot acceleration curves for all countries
matplot(Time, accel_values, type = 'l', lty = 1, col = 'grey',
        xlab = 'Years', ylab = 'Acceleration', main = 'Acceleration Curves for All Countries')

# Compute and plot the mean acceleration curve
mean_accel <- rowMeans(accel_values)
lines(Time, mean_accel, col = 'red', lwd = 2)
legend("topright", legend = c("Mean Acceleration"), col = c("red"), lty = 1, lwd = 2)



accel_fd <- deriv.fd(fdi_fd, deriv = 2)

# Evaluate the acceleration at the time points
accel_values <- eval.fd(Time, accel_fd)

# Plot acceleration curves for all countries
matplot(Time, accel_values, type = 'l', lty = 1, col = 'grey',
        xlab = 'Years', ylab = 'Acceleration', main = 'Acceleration Curves for All Countries')

# Compute and plot the mean acceleration curve
mean_accel <- rowMeans(accel_values)
lines(Time, mean_accel, col = 'red', lwd = 2)
legend("topright", legend = c("Mean Acceleration"), col = c("red"), lty = 1, lwd = 2)



accel_fd <- deriv.fd(ipc_fd, deriv = 2)

# Evaluate the acceleration at the time points
accel_values <- eval.fd(Time, accel_fd)

# Plot acceleration curves for all countries
matplot(Time, accel_values, type = 'l', lty = 1, col = 'grey',
        xlab = 'Years', ylab = 'Acceleration', main = 'Acceleration Curves for All Countries')

# Compute and plot the mean acceleration curve
mean_accel <- rowMeans(accel_values)
lines(Time, mean_accel, col = 'red', lwd = 2)
legend("topright", legend = c("Mean Acceleration"), col = c("red"), lty = 1, lwd = 2)




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


# Perform ANOVA
anova_accel_region <- aov(MeanAcceleration ~ Region, data = accel_data)
summary(anova_accel_region)


# Perform ANOVA
anova_accel_income <- aov(MeanAcceleration ~ IncomeLevel, data = accel_data)
summary(anova_accel_income)



# Assuming 'regional' data frame is already loaded and aligned with 'countries'
countries <- colnames(grw)

# Ensure 'regional' data frame is aligned
regional <- regional[regional$Country %in% countries, ]
regional <- regional[match(countries, regional$Country), ]

# Extract regions and income levels
regions <- factor(regional$Region)
income_levels <- factor(regional$Income.level)

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
    beta_basis_const = const_basis,
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

par(mfrow = c(2, 2))  # Adjust layout as needed

for (region in names(regression_results_regions)) {
  beta_estimates <- regression_results_regions[[region]]$beta_estimates
  beta3_fd <- beta_estimates[[4]]$fd  # Interaction term
  
  plot(beta3_fd, main = paste("Beta 3 (Interaction) for", region))
}

par(mfrow = c(1, 1))  # Reset layout


par(mfrow = c(2, 2))  # Adjust layout as needed

for (level in names(regression_results_income)) {
  beta_estimates <- regression_results_income[[level]]$beta_estimates
  beta3_fd <- beta_estimates[[4]]$fd  # Interaction term
  
  plot(beta3_fd, main = paste("Beta 3 (Interaction) for", level))
}

par(mfrow = c(1, 1))  # Reset layout

beta3_fd <- fRegress_result$betaestlist[[4]]$fd  # Interaction term beta3
beta3_values <- eval.fd(Time, beta3_fd)


# Create a data frame with Time and beta3_values
beta3_data <- data.frame(Time = Time, Beta3 = beta3_values)

# Initialize plot
plot(Time, beta3_values, type = 'n', xlab = 'Time', ylab = 'Beta3 (Interaction)', main = 'Beta3 Interaction Term')

# Identify positive and negative sections
positive_indices <- which(beta3_values >= 0)
negative_indices <- which(beta3_values < 0)

# Plot positive sections in red
plot(Time[positive_indices], beta3_values[positive_indices], col = 'red', lwd = 2)

# Plot negative sections in green
plot(Time[negative_indices], beta3_values[negative_indices], col = 'green', lwd = 2)

# Add a horizontal line at zero
abline(h = 0, lty = 2)

# Add legend
legend('topright', legend = c('Positive', 'Negative'), col = c('red', 'green'), lty = 1, lwd = 2)

# Perform PCA on the registered data
pca_grw <- pca.fd(grw_fd_registered1$regfd, nharm = 3)
# Plot the principal component functions (harmonics)
plot(pca_grw$harmonics, main = "Principal Components of Registered Growth Rate")
# Display the proportion of variance explained by each component
print(pca_grw$varprop)


# Perform PCA on the registered data
pca_grw <- pca.fd(ipc_fd_registered2$regfd, nharm = 3)
# Plot the principal component functions (harmonics)
plot(pca_grw$harmonics, main = "Principal Components of Registered Growth Rate")
# Display the proportion of variance explained by each component
print(pca_grw$varprop)

# Perform PCA on the registered data
pca_grw <- pca.fd(fdi_fd_registered2$regfd, nharm = 3)
# Plot the principal component functions (harmonics)
plot(pca_grw$harmonics, main = "Principal Components of Registered Growth Rate")
# Display the proportion of variance explained by each component
print(pca_grw$varprop)



