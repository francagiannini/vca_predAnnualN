library(copula) # For copula fitting and sampling
library(dplyr)   # For data manipulation (e.g., %>% and select)
library(stats)   # For empirical CDF/quantile functions

# --- 0. Pre-computation / Setup ---
# Your original data source (assuming it's already loaded)
# sasoutput_fil_NLES5 <- read.csv("path/to/your/sasoutput_fil_NLES5.csv")

# Identify the columns you want to vary in your sensitivity analysis
# Make sure to exclude any fixed parameters if you've decided to keep them constant.
input_names_nles5 <- c(
  "Y", "NT", "MNCS", "MNCA", "MNudb", "M1", "M2", "F0", "F1", "F2",
  "G0", "G1", "G2", "WC", "M", "W", "MP", "WP", "CU", "jbnr",
  "AAa", "AAb", "APb"
)
# Ensure these columns are numeric (remove NAs if any)
original_inputs_df <- sasoutput_fil_NLES5 %>%
  select(all_of(input_names_nles5)) %>%
  na.omit() # Remove rows with NA in any selected column

# Number of simulations needed for Sobol' (e.g., 2000 * (D + 2) for D inputs)
# Ensure n_simulations is the same as what you used in soboljansen
# n_simulations <- 2000 # Example, use your actual n_simulations

# --- 1. Estimate Empirical Marginals (ECDFs) ---
# For each input, create its empirical cumulative distribution function (ECDF)
# and its quantile function (inverse ECDF).
# These will be used to transform uniform samples back to the original marginals.

ecdfs <- list()
quants <- list()
for (col_name in colnames(original_inputs_df)) {
  ecdfs[[col_name]] <- ecdf(original_inputs_df[[col_name]])
  quants[[col_name]] <- quantile(original_inputs_df[[col_name]], probs = seq(0, 1, length.out = length(original_inputs_df[[col_name]])), type = 1)
  # A simpler approximation for quantile: approxfun(ecdfs[[col_name]])$f(p) which is what `q_dist` does in `copula` for empirical.
}

# --- 2. Fit a Copula Model to Your Data ---
# Convert data to pseudo-observations (values between 0 and 1, representing ranks)
# This is essentially applying the ECDF to each data point.
pseudo_obs <- as.data.frame(lapply(original_inputs_df, function(x) {
  rank(x, ties.method = "average") / (length(x) + 1)
}))

# Choose a copula family. Common choices:
# - NormalCopula (Gaussian copula): models linear correlation
# - tCopula (Student-t copula): models linear and tail dependency (useful for extreme events)
# - For highly complex/non-linear dependencies, you might explore others or an empirical copula
#   However, for this application, a Gaussian or t-copula is a good starting point.

# For simplicity, let's start with a Gaussian copula.
# It uses Pearson's correlation coefficient as its parameter.
# dim = number of input variables
d <- ncol(original_inputs_df)
gaussian_cop <- normalCopula(dim = d, dispstr = "un") # "un" for unstructured correlation matrix

# Fit the copula to the pseudo-observations
# This estimates the correlation matrix (for Gaussian/t-copulas) or other parameters.
# method = "ml" (maximum likelihood) is generally robust.
# Use start = P2kendall(as.matrix(pseudo_obs)) for a better starting point for ML
# Or method = "itau" for inversion of Kendall's tau (simpler, but less efficient for Gaussian/t)

fit_cop <- tryCatch({
  # Use empirical Kendall's tau for a starting point
  P_kendall_mat <- cor(pseudo_obs, method = "kendall")
  # Convert Kendall's tau to Pearson's rho for Gaussian copula starting value
  # rho = sin(tau * pi / 2)
  start_rho <- sin(P_kendall_mat * pi / 2)
  start_rho[is.na(start_rho)] <- 0 # Handle potential NA for constant columns

  # Ensure the starting matrix is positive definite
  start_rho_chol <- try(chol(start_rho), silent = TRUE)
  if (inherits(start_rho_chol, "try-error")) {
    warning("Starting correlation matrix not positive definite. Using default start.")
    start_rho <- diag(d) # Start with identity matrix if problematic
  }
  start_param_vec <- start_rho[upper.tri(start_rho)] # Extract upper triangle for parameter vector

  fitCopula(copula = gaussian_cop, data = as.matrix(pseudo_obs), method = "ml",
            start = start_param_vec)
}, error = function(e) {
  warning(paste("Error fitting copula, attempting with simpler approach or default start:", e$message))
  # Fallback to simple method or no start
  fitCopula(copula = gaussian_cop, data = as.matrix(pseudo_obs), method = "ml")
})


# Check if the fitting was successful
if (inherits(fit_cop, "try-error")) {
  stop("Failed to fit the copula model. Please check your input data or try a different copula family.")
} else {
  cat("\nCopula model fitted successfully.\n")
  print(summary(fit_cop))
}

# --- 3. Generate Correlated Samples from the Fitted Copula ---
# These samples will be uniform (U[0,1]) but with the specified dependency structure.
# You need 2 * n_simulations for soboljansen (X1 and X2 are stacked)
num_samples_needed <- 2 * n_simulations

# Set seed for reproducibility
set.seed(123)
uniform_correlated_samples <- rCopula(n = num_samples_needed, copula = fit_cop@copula)

# --- 4. Transform Uniform Samples Back to Original Marginals ---
# Apply the inverse of the ECDF (quantile function) for each variable.
# This ensures the generated samples have the original marginal distributions
# while retaining the correlation structure from the copula.

X_correlated_data <- as.data.frame(matrix(NA, nrow = num_samples_needed, ncol = d))
colnames(X_correlated_data) <- colnames(original_inputs_df)

for (j in 1:d) {
  col_name <- colnames(original_inputs_df)[j]
  # Use the quantile function for the empirical distribution
  # For data with many unique values, `quants[[col_name]]` can be large.
  # A faster way for the empirical quantile function:
  X_correlated_data[[col_name]] <- quantile(original_inputs_df[[col_name]], probs = uniform_correlated_samples[, j], type = 1)
}

# --- 5. Split into X1_data and X2_data for Sobol' ---
X1_data_correlated <- X_correlated_data[1:n_simulations, ]
X2_data_correlated <- X_correlated_data[(n_simulations + 1):(2 * n_simulations), ]

# --- 6. Verify Correlations (Optional, but highly recommended) ---
cat("\n--- Verification of Correlations ---\n")
cat("\nCorrelation matrix of original data:\n")
print(cor(original_inputs_df, method = "spearman")) # Use Spearman for rank correlations

cat("\nCorrelation matrix of X1_data_correlated:\n")
print(cor(X1_data_correlated, method = "spearman"))

cat("\nCorrelation matrix of X2_data_correlated:\n")
print(cor(X2_data_correlated, method = "spearman"))

# You can also check if the marginal distributions are preserved
# plot(ecdfs[[col_name]], main=paste("ECDF of", col_name))
# lines(ecdf(X1_data_correlated[[col_name]]), col="red") # Should be very similar

# --- 7. Perform Sobol' Analysis with Correlated Data ---
# Now use X1_data_correlated and X2_data_correlated
sobol_results_nles5_model <- soboljansen(
  model = calculate_NLES5_for_sobol,
  X1 = X1_data_correlated,
  X2 = X2_data_correlated,
  n = n_simulations,
  conf = 0.95
)

print(sobol_results_nles5_model)
