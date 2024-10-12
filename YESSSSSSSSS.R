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
nbasis <- 25  # Adjust as needed 21 was correct

# Create the basis for each functional object
grw_basis <- create.fourier.basis(rangeval, nbasis)
ipc_basis <- create.fourier.basis(rangeval, nbasis)
fdi_basis <- create.fourier.basis(rangeval, nbasis)

# Convert data to functional objects using smooth.basis
grw_fd <- smooth.basis(Time, grw, grw_basis)$fd
ipc_fd <- smooth.basis(Time, ipc, ipc_basis)$fd
fdi_fd <- smooth.basis(Time, fdi, fdi_basis)$fd

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
