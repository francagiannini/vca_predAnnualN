
# Load packages
library(dtwclust)
library(sensitivity)
library(xgboost)
library(dplyr)
library(ggplot2)


# --- Adaptation for Sobol Sensitivity Analysis ---

# Now, when you define the NLES5_surrogate_model function for Sobol,
# it needs to perform the same one-hot encoding on the input samples (X)
# that the `sobol` function generates.

# You will need to manage the mapping of original variable names to their one-hot encoded versions
# in your `factor_ranges` for Sobol.

# Revised NLES5_surrogate_model function:
NLES5_surrogate_model_xgb <- function(X_sobol_raw) {
  # X_sobol_raw matrix where columns correspond to original input variables
  # (e.g., NS, na, Mau, CU, etc.) with values within their defined ranges.

  # Convert to a data frame
  X_sobol_df <- as.data.frame(X_sobol_raw)
  colnames(X_sobol_df) <- names(ohe_data %>% select(-PUdvaskF)) # Use original input names

  # Convert to factors
  for (col in categorical_vars) {
    if (col %in% names(X_sobol_df)) {
      X_sobol_df[[col]] <- as.factor(X_sobol_df[[col]])
      # Ensure factor levels are consistent with training data (if any issues)
      # e.g., levels(X_sobol_df[[col]]) <- levels(sasoutput_fil[[col]])
    }
  }

  # Apply the SAME dummy encoding transformation
  X_sobol_ohe <- data.frame(predict(dummy_model, newdata = X_sobol_df))

  # Make predictions using your trained XGBoost model
  predictions <- predict(xgb_surrogate_model, xgb.DMatrix(as.matrix(X_sobol_ohe)))
  return(predictions)
}


# --- Define Input Ranges/Distributions for Sobol (CRITICAL UPDATE) ---
# Here we name factor because it emulates a factorial design
# 'Factor_ranges' for Sobol should refer to the original input variables
# and their ranges. The `NLES5_surrogate_model_xgb` should handle the encoding internally.

# Re-define factor_ranges based on original numerical values or factor levels
# IMPORTANT: Use the distinct values for categorical variables to define min/max.
categorical_vars <- c("Mau", "Mfu", "Vau", "Vfu", "SoilG")

factor_ranges_ohe <- list()
for (col_name in names(sasoutput_fil %>% select(-PUdvaskF))) {
  if (col_name %in% categorical_vars) {
    # For categorical variables, min/max should be the range of their numeric levels
    # Ensure they are coerced to numeric if they are factors
    unique_levels <- as.numeric(levels(sasoutput_fil[[col_name]]))
    factor_ranges_ohe[[col_name]] <- c(min(unique_levels), max(unique_levels))
  } else {
    # For continuous variables
    factor_ranges_ohe[[col_name]] <- c(min(sasoutput_fil[[col_name]]), max(sasoutput_fil[[col_name]]))
  }
}

min_vals_ohe <- unlist(lapply(factor_ranges_ohe, function(x) x[1]))
max_vals_ohe <- unlist(lapply(factor_ranges_ohe, function(x) x[2]))

k_ohe <- length(min_vals_ohe)
N_sobol <- 1000 # Number of initial points for the Sobol sequence

# Run the Sobol analysis with the OHE-aware surrogate model
system.time({
  sobol_results_ohe <- sobol(
    model = NLES5_surrogate_model_xgb,
    X1 = matrix(runif(N_sobol * k_ohe), nrow = N_sobol, ncol = k_ohe, dimnames = list(NULL, names(factor_ranges_ohe))),
    X2 = matrix(runif(N_sobol * k_ohe), nrow = N_sobol, ncol = k_ohe, dimnames = list(NULL, names(factor_ranges_ohe))),
    n = N_sobol,
    min = min_vals_ohe,
    max = max_vals_ohe,
    conf = 0.95
  )
})

print(sobol_results_ohe)

# (Continue with plotting Sobol indices as before)





