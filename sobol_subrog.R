# Load packages
library(tidyverse)
library(readxl)
library(caret)
library(xgboost)
library(sensobol)

# --- ASSUMPTION: All necessary setup from previous steps has been run ---
# This means the following objects MUST be available in your R environment:
# 1. `sasoutput_fil` (your data.frame with original variables, factors converted)
# 2. `dummy_model` (the `dummyVars` object for one-hot encoding)
# 3. `categorical_vars` (character vector of your original categorical variable names)
# 4. `factor_levels_map` (list mapping categorical var names to their levels)
# 5. `xgb_surrogate_model` (your trained XGBoost model)
# 6. `original_input_names` (character vector of original input variable names from `sasoutput_fil`)
# 7. `model_expected_ohe_feature_names` (character vector of OHE feature names from your trained XGBoost model)

# If any of these are not defined, you MUST run the preceding code sections to create them.

# --- 1. Define Input Ranges/Distributions for `sensobol` ---
# `sensobol` expects a data.frame with min/max for each parameter or a list.
# Let's create a data.frame for the ranges.
input_bounds <- data.frame(
  parameters = original_input_names,
  min = NA, # Placeholder
  max = NA  # Placeholder
)

# Populate min/max values
for (i in seq_along(original_input_names)) {
  col_name <- original_input_names[i]
  if (col_name %in% categorical_vars) {
    unique_levels <- as.numeric(levels(sasoutput_fil[[col_name]]))
    input_bounds$min[i] <- min(unique_levels, na.rm = TRUE)
    input_bounds$max[i] <- max(unique_levels, na.rm = TRUE)
  } else {
    input_bounds$min[i] <- min(sasoutput_fil[[col_name]], na.rm = TRUE)
    input_bounds$max[i] <- max(sasoutput_fil[[col_name]], na.rm = TRUE)
  }
}

# The number of model evaluations (`N`) for `sensobol`
# `sensobol` internally handles the total number of runs based on `N` and `calc_second_order`.
# A higher N gives more accurate results. Start with 1000 or 2000 for decent results.
N_sensobol <- 1000

# --- 2. Generate the Sobol design matrix using `sensobol` ---
# `sobol_matrices()` generates the input matrices for the model evaluations.
# It uses quasi-random sequences like the Sobol sequence for efficient space-filling.
# We'll use the "sobol" type for the matrices.
# `calc_second_order = TRUE` calculates second-order interactions (which increases runs to (p+2)*N + p*(p-1)/2 * N if you want it).
# For just S1 and ST, calc_second_order=FALSE is sufficient and costs (p+2)*N.
sobol_design <- sobol_matrices(
  N = N_sensobol,
  params = input_bounds$parameters, # Correct argument name
  type = "QRN", # Generates values between 0 and 1
  order = "first" # Generates matrices for first and total order indices
)

# `sobol_design` is now a data.frame ready for model evaluation.
# It has N * (p+2) rows (where p is number of parameters) if calc_second_order=FALSE.

# --- 3. Evaluate the XGBoost surrogate model for each row in the design matrix ---
cat("\nEvaluating XGBoost surrogate model for sensobol design...\n")
system.time({
  # Your anonymous model function, adapted slightly for sensobol's output.
  # sensobol_matrices creates a data.frame, so X_sobol_raw will be a data.frame.
  sensobol_model_output <- function(X_sobol_df_input) {
    # X_sobol_df_input is a data.frame from `sobol_design`

    # Ensure column names are correct (sensobol handles this if `parameters` in `sobol_matrices` are set)
    # They should match `original_input_names`

    # Convert sampled columns back to factors, ensuring consistent levels
    for (col in categorical_vars) {
      if (col %in% names(X_sobol_df_input)) {
        X_sobol_df_input[[col]] <- factor(X_sobol_df_input[[col]], levels = factor_levels_map[[col]])
      }
    }

    # Apply the SAME one-hot encoding transformation using the stored dummy_model
    X_sobol_ohe_matrix <- predict(dummy_model, newdata = X_sobol_df_input)

    # Handle cases where predict.dummyVars might return a single row as a vector
    if (!is.matrix(X_sobol_ohe_matrix) && is.vector(X_sobol_ohe_matrix)) {
      X_sobol_ohe_matrix <- t(as.matrix(X_sobol_ohe_matrix))
    }

    # Ensure the OHE matrix columns match the model's expected features
    temp_matrix <- matrix(0, nrow = nrow(X_sobol_ohe_matrix), ncol = length(model_expected_ohe_feature_names),
                          dimnames = list(NULL, model_expected_ohe_feature_names))
    common_cols <- intersect(colnames(X_sobol_ohe_matrix), model_expected_ohe_feature_names)
    temp_matrix[, common_cols] <- X_sobol_ohe_matrix[, common_cols]
    X_sobol_ohe_matrix <- temp_matrix

    # Make predictions using your trained XGBoost model
    predictions <- predict(xgb_surrogate_model, xgb.DMatrix(as.matrix(X_sobol_ohe_matrix)))
    return(predictions)
  }

  # Apply the model function to the generated design matrix
  sensobol_outputs <- sensobol_model_output(sobol_design)
})

# --- 4. Calculate Sobol indices using `sensobol` ---
# `sobol_indices()` takes the design matrix and the model outputs.
sobol_results_sensobol <- sobol_indices(
  Y = sensobol_outputs,
  N = N_sensobol,
  parameters = input_bounds$parameters,
  type = "A", # Matches the type used in sobol_matrices()
  calc_second_order = FALSE # Must match calc_second_order in sobol_matrices()
)

print(sobol_results_sensobol)

# --- 5. Visualization of Sobol Indices from `sensobol` ---
# `sensobol` stores results in a data.frame ready for ggplot.
# It automatically provides S1 (first-order) and ST (total-order) indices.
# It's usually in a data.frame format, often accessed directly.

# Plotting `sensobol` results is typically straightforward using its internal plot function
# or by directly accessing the data.frame `sobol_results_sensobol`.

# Example of manual plotting (similar to your previous approach):
plot_sobol_data_sensobol <- sobol_results_sensobol %>%
  filter(
    !is.na(S1) | !is.na(ST) # Ensure there are values
  ) %>%
  select(parameter, S1, S1_conf, ST, ST_conf) %>%
  # S1_conf and ST_conf are usually the margin of error, not min/max CI.
  # So, S1_low = S1 - S1_conf, S1_high = S1 + S1_conf
  mutate(
    S1_low = S1 - S1_conf,
    S1_high = S1 + S1_conf,
    ST_low = ST - ST_conf,
    ST_high = ST + ST_conf
  ) %>%
  tidyr::gather(key = "Index_Type", value = "Value", S1, ST) %>%
  mutate(
    Conf_Low = ifelse(Index_Type == "S1", S1_low, ST_low),
    Conf_High = ifelse(Index_Type == "S1", S1_high, ST_high),
    Variable = fct_reorder(parameter, Value, .desc = TRUE) # Order by value
  )

ggplot(plot_sobol_data_sensobol, aes(x = Variable, y = Value, fill = Index_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = Conf_Low, ymax = Conf_High),
                position = position_dodge(width = 0.8), width = 0.2) +
  labs(title = "Sobol Sensitivity Indices (First-Order and Total-Order) with XGBoost Surrogate (sensobol)",
       x = "Input Variable",
       y = "Sensitivity Index Value",
       fill = "Index Type") +
  scale_fill_manual(values = c("S1" = "steelblue", "ST" = "coral")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
