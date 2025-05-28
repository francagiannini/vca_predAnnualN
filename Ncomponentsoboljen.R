# Nitrogen Model Sensitivity Analysis in R (Sobol' Indices)

# This script defines the nitrogen model and performs global sensitivity analysis
# using Sobol' indices via the 'sensitivity' package.

# Install and load the 'sensitivity' package if you haven't already:
# install.packages("sensitivity")
library(sensitivity)
library(tidyverse)
library(patchwork)

# 1. Define the Nitrogen Model Function
# The model calculates N based on various nitrogen parameters and input variables.
# N = βt NT + βCS MNCS + βCA MNCA + βudb MNudb + βm1 (M1+M2)/2 + βf0 F0 +
#     βf1 (F1+F2)/2 + βg0 G0 + βm1 (G1+G2)/2
#
# This function will be passed directly to the 'sobol' function as the 'model' argument.
# It expects a matrix 'X' where each column is an input variable and each row is a sample.
calculate_N_for_sobol <- function(X) {
  # Default Beta Parameter Values (treated as fixed for this sensitivity analysis)
  beta_t = 0.456793
  beta_CS = 0.04957
  beta_g0 = 0.014099
  beta_CA = 0.157044
  beta_udb = 0.016314
  beta_f0 = 0.038245
  beta_m1_M = 0.026499
  beta_m1_G = 0.0265
  beta_f1 = 0.025499

  # Extract input variables from the matrix X by column name
  NT = X[, "NT"]
  MNCS = X[, "MNCS"]
  G0 = X[, "G0"]
  MNCA = X[, "MNCA"]
  MNudb = X[, "MNudb"]
  F0 = X[, "F0"]
  M1 = X[, "M1"]
  M2 = X[, "M2"]
  G1 = X[, "G1"]
  G2 = X[, "G2"]
  F1 = X[, "F1"]
  F2 = X[, "F2"]

  # Calculate the N value for each row (sample)
  N_value <- (
    beta_t * NT +
      beta_CS * MNCS +
      beta_CA * MNCA +
      beta_udb * MNudb +
      beta_m1_M * ((M1 + M2) / 2) +
      beta_f0 * F0 +
      beta_f1 * ((F1 + F2) / 2) +
      beta_g0 * G0 +
      beta_m1_G * ((G1 + G2) / 2)
  )
  return(N_value)
}

# 2. Define Input Variable Names and their Ranges for Sampling
# These ranges are crucial for the global sensitivity analysis.
# Adjust them based on the actual variability observed in your dataset.
input_names <- c("NT", "MNCS", "G0", "MNCA", "MNudb", "F0", "M1", "M2", "G1", "G2", "F1", "F2")

# Define the min and max values for each input variable's uniform distribution
# These are the same ranges used in the previous successful Sobol' runs.
input_ranges <- list(
  NT = c(50, 150),
  MNCS = c(20, 80),
  G0 = c(5, 20),
  MNCA = c(10, 40),
  MNudb = c(5, 20),
  F0 = c(10, 30),
  M1 = c(1, 10),
  M2 = c(1, 10),
  G1 = c(2, 12),
  G2 = c(2, 12),
  F1 = c(3, 15),
  F2 = c(3, 15)
)

# Extract min and max values into vectors for the sobol() call
min_vals <- unlist(lapply(input_ranges, function(x) x[1]))
max_vals <- unlist(lapply(input_ranges, function(x) x[2]))

# 3. Perform Global Sensitivity Analysis using Sobol' Indices

cat("--- Global Sensitivity Analysis (Sobol' Indices) ---\n")

# Number of Monte Carlo simulations per variable (N in Sobol' terminology)
# A higher N means more accurate results but longer computation time.
# N = 2000 is a good starting point. Total runs will be N * (2k + 2)
n_simulations <- 2000

# Generate X1 and X2 matrices with explicit column naming
X1_data <- as.data.frame(lapply(input_names, function(var_name) {
  runif(n_simulations, min = input_ranges[[var_name]][1], max = input_ranges[[var_name]][2])
}))
colnames(X1_data) <- input_names # Explicitly set column names

X2_data <- as.data.frame(lapply(input_names, function(var_name) {
  runif(n_simulations, min = input_ranges[[var_name]][1], max = input_ranges[[var_name]][2])
}))
colnames(X2_data) <- input_names # Explicitly set column names

# Perform the Sobol' sensitivity analysis directly using the 'sobol' function.
# The 'model' argument is now your 'calculate_N_for_sobol' function.
# X1 and X2 are generated within the 'sobol' call, using the specified 'min' and 'max' ranges.
# Setting 'order = 2' will calculate first-order, total-order, and second-order indices.
sobol_results_n_model <- soboljansen(
  model = calculate_N_for_sobol,
  X1 = X1_data, # Use the explicitly named X1_data
  X2 = X2_data, # Use the explicitly named X2_data
  n = n_simulations,
  conf = 0.95#,
  #order = 2 # Request calculation of second-order indices
)

# Print the results
cat("Sobol' First-Order Sensitivity Indices (S1):\n")
print(sobol_results_n_model$S)
cat("\n")

cat("Sobol' Total-Order Sensitivity Indices (ST):\n")
print(sobol_results_n_model$T)
cat("\n")

# To access second-order indices (S2), if 'order = 2' was set in the initial sobol() call:
if (!is.null(sobol_results_n_model$S2)) {
  cat("Sobol' Second-Order Sensitivity Indices (S2):\n")
  print(sobol_results_n_model$S2)
  cat("\n")
} else {
  cat("Second-Order Sensitivity Indices (S2) were not calculated or are not available.\n")
  cat("Ensure 'order = 2' is set in the initial sobol() call for sample generation.\n\n")
}

# --- 4. Visualization of Sobol Indices ---
# Access the first-order indices
first_order_indices <- data.frame(
  Variable = rownames(sobol_results_n_model$S),
  S1 = sobol_results_n_model$S[, "original"],
  S1_conf_low = sobol_results_n_model$S[, "min. c.i."],
  S1_conf_high = sobol_results_n_model$S[, "max. c.i."]
)

# Access the total-order indices
total_order_indices <- data.frame(
  Variable = rownames(sobol_results_n_model$T),
  ST = sobol_results_n_model$T[, "original"],
  ST_conf_low = sobol_results_n_model$T[, "min. c.i."],
  ST_conf_high = sobol_results_n_model$T[, "max. c.i."]
)

# Combine for easier visualization
sobol_df <- dplyr::left_join(first_order_indices, total_order_indices, by = "Variable")

# Prepare data for plotting
# Ensure 'tidyverse' is loaded for 'ggplot2' and 'tidyr::gather', 'dplyr::mutate', 'forcats::fct_reorder'
# If not already loaded from your initial script:
# library(tidyverse) # or library(ggplot2); library(tidyr); library(dplyr); library(forcats)

plot_sobol_data <- sobol_df |>
  tidyr::gather(key = "Index_Type", value = "Value", S1, ST) |>
  dplyr::mutate(
    Conf_Low = ifelse(Index_Type == "S1", S1_conf_low, ST_conf_low),
    Conf_High = ifelse(Index_Type == "S1", S1_conf_high, ST_conf_high),
    Variable = forcats::fct_reorder(Variable, Value, .desc = TRUE) # Order by value
  )

# Create the plot
ggplot(plot_sobol_data, aes(x = Variable, y = Value, fill = Index_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = Conf_Low, ymax = Conf_High),
                position = position_dodge(width = 0.8), width = 0.2) +
  labs(title = "Sobol Sensitivity Indices (First-Order and Total-Order) for N Model",
       x = "Input Variable",
       y = "Sensitivity Index Value",
       fill = "Index Type") +
  scale_fill_manual(values = c("S1" = "steelblue", "ST" = "coral")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### from sampled dist ----


# Generate X1 and X2 data frames by sampling directly from 'sasoutput_fil'
# We use replace = TRUE to allow for sampling with replacement, as n_simulations
# might be larger than the number of rows in sasoutput_fil.
# na.omit() is used to remove any NA values from the columns before sampling,
# as 'sample' cannot handle NAs directly.


asignationN <- c(
  # N
  "MNCS"="NS",
  "MNCA"="na",
  "M1"="nlevelMin1",
  "M2"="nlevelMin2",
  "F0"="NlevelFix0",
  "F1"="NlevelFix1",
  "F2"="NlevelFix2",
  "G0"="NlevelGod0",
  "G1"="NlevelGod1",
  "G2"="NlevelGod2",
  "MNudb"="Nudb",
  "NT"="TN"
)

sasoutput_fil_N <- rename(sasoutput_fil, all_of(asignationN))

X1_data <- as.data.frame(lapply(input_names, function(var_name) {
  if (var_name %in% colnames(sasoutput_fil_N)) {
    # Sample from the observed values of the variable, omitting NAs
    sample(na.omit(sasoutput_fil_N[[var_name]]), size = n_simulations, replace = TRUE)
  } else {
    stop(paste0("Column '", var_name, "' not found in 'sasoutput_fil_N'. Please check input_names."))
  }
}))
colnames(X1_data) <- input_names # Explicitly set column names

X2_data <- as.data.frame(lapply(input_names, function(var_name) {
  if (var_name %in% colnames(sasoutput_fil_N)) {
    # Sample from the observed values of the variable, omitting NAs
    sample(na.omit(sasoutput_fil_N[[var_name]]), size = n_simulations, replace = TRUE)
  } else {
    stop(paste0("Column '", var_name, "' not found in 'sasoutput_fil_N'. Please check input_names."))
  }
}))
colnames(X2_data) <- input_names # Explicitly set column names

# Perform the Sobol' sensitivity analysis directly using the 'sobol' function.
# The 'model' argument is now your 'calculate_N_for_sobol' function.
# X1 and X2 are generated within the 'sobol' call, using the specified 'min' and 'max' ranges.
# Setting 'order = 2' will calculate first-order, total-order, and second-order indices.
sobol_results_n_model <- soboljansen(
  model = calculate_N_for_sobol,
  X1 = X1_data, # Use the explicitly named X1_data
  X2 = X2_data, # Use the explicitly named X2_data
  n = n_simulations,
  conf = 0.95#,
  #order = 2 # Request calculation of second-order indices
)

# Print the results
cat("Sobol' First-Order Sensitivity Indices (S1):\n")
print(sobol_results_n_model$S)
cat("\n")

cat("Sobol' Total-Order Sensitivity Indices (ST):\n")
print(sobol_results_n_model$T)
cat("\n")

# To access second-order indices (S2), if 'order = 2' was set in the initial sobol() call:
if (!is.null(sobol_results_n_model$S2)) {
  cat("Sobol' Second-Order Sensitivity Indices (S2):\n")
  print(sobol_results_n_model$S2)
  cat("\n")
} else {
  cat("Second-Order Sensitivity Indices (S2) were not calculated or are not available.\n")
  cat("Ensure 'order = 2' is set in the initial sobol() call for sample generation.\n\n")
}

# --- 4. Visualization of Sobol Indices ---
# Access the first-order indices
first_order_indices <- data.frame(
  Variable = rownames(sobol_results_n_model$S),
  S1 = sobol_results_n_model$S[, "original"],
  S1_conf_low = sobol_results_n_model$S[, "min. c.i."],
  S1_conf_high = sobol_results_n_model$S[, "max. c.i."]
)

# Access the total-order indices
total_order_indices <- data.frame(
  Variable = rownames(sobol_results_n_model$T),
  ST = sobol_results_n_model$T[, "original"],
  ST_conf_low = sobol_results_n_model$T[, "min. c.i."],
  ST_conf_high = sobol_results_n_model$T[, "max. c.i."]
)

# Combine for easier visualization
sobol_df <- dplyr::left_join(first_order_indices, total_order_indices, by = "Variable")

# Prepare data for plotting
# Ensure 'tidyverse' is loaded for 'ggplot2' and 'tidyr::gather', 'dplyr::mutate', 'forcats::fct_reorder'
# If not already loaded from your initial script:
# library(tidyverse) # or library(ggplot2); library(tidyr); library(dplyr); library(forcats)

plot_sobol_data <- sobol_df |>
  tidyr::gather(key = "Index_Type", value = "Value", S1, ST) |>
  dplyr::mutate(
    Conf_Low = ifelse(Index_Type == "S1", S1_conf_low, ST_conf_low),
    Conf_High = ifelse(Index_Type == "S1", S1_conf_high, ST_conf_high),
    Variable = forcats::fct_reorder(Variable, Value, .desc = TRUE) # Order by value
  )


sobolplotN <-
  plot_sobol_data |> mutate(Input=recode(Variable,!!!asignationN)) |>
  mutate(Input=recode(Input,!!!translate)) |>
  ggplot(aes(x = Input, y = Value, fill = Index_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = Conf_Low, ymax = Conf_High),
                position = position_dodge(width = 0.8), width = 0.2) +
  labs(title = "Sensitivity estimates",
       x = "Input",
       y = "Sensitivity index",
       fill = "Sobol' index") +
  scale_fill_manual(values = c("S1" = "steelblue", "ST" = "coral"),
                    labels = c("first-order", "total-order")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        plot.title = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 9))



uncplot <-
  ggplot(data.frame(y = sobol_results_n_model$y), aes(x = y)) +
  geom_histogram(binwidth = diff(range(sobol_results_n_model$y))/30, # Auto-adjust binwidth
                 fill = "darkseagreen", color = "white", alpha = 0.8) +
  labs(title = "Uncertainty estimates",
       x = "N ouput",
       y = "Absolute frequency") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 9))



uncplot + sobolplotN +plot_layout(ncol = 2,
                                  widths = c(0.5,1.5)
                                  #heights = c(1.2,1)
                                  )

Npoint3 <-
  uncplot + sobolplotN +plot_layout(ncol = 1,
                                  #widths = c(0.5,1.5)
                                  heights = c(1,1.5)
)

ggsave("Npoint3.png",Npoint3,
       width = 183, height =183*.8,
       units = "mm")

## Sampled input distribution ----


# --- 1. Prepare Sampled Data (from X1_data) ---

sampled_data_long <- X1_data |>
  pivot_longer(
    cols = everything(), # Selects all columns (your parameters)
    names_to = "Parameter", # New column for parameter names
    values_to = "Value"     # New column for sampled values
  ) |>
  mutate(Source = "Sampled (sensitivity analysis)") # Label this data set

# --- 2. Prepare Original Data (from sasoutput_fil_N) ---

original_data_long <- sasoutput_fil_N |>
  select(all_of(input_names)) |> # Select only the 12 relevant parameter columns
  pivot_longer(
    cols = everything(),
    names_to = "Parameter",
    values_to = "Value"
  ) |>
  na.omit() |> # Remove any NA values, similar to how sampling handled them
  mutate(Source = "Original (dataset)") # Label this data set

# --- 3. Combine both data frames into a single data frame ---
# This allows ggplot2 to plot both sets of densities in the same plot.
combined_plot_data <- bind_rows(sampled_data_long, original_data_long)


# --- 4. Create the faceted density plot with both curves ---
dist_samplvs.orig <-
  combined_plot_data |> mutate(Input=recode(Parameter,!!!asignationN)) |>
  mutate(Input=recode(Input,!!!translate)) |>
  ggplot(aes(x = Value, color = Source,
             fill = Source)) +
  geom_density(alpha = 0.2, linewidth = 0.8) +
  facet_wrap( ~ Input, scales = "free", ncol = 4,
              labeller = labeller(Parameter=facet_labels_expr)) +
  # Add informative titles and labels.
  labs(
    #title = "Sampled vs. Observed input distributions",
    x = "Input value",
    y = "Density",
    color = "Source",fill = "Source") +
  scale_color_manual(
    values = c(
      "Original (dataset)" = "darkgreen",
      "Sampled (sensitivity analysis)" = "steelblue")) +
  scale_fill_manual(
    values = c(
      "Original (dataset)" = "darkgreen",
      "Sampled (sensitivity analysis)" = "steelblue") ) +
  theme_light() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5),
    strip.text = element_text(size = 8,),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9))

dist_samplvs.orig

ggsave("dist_samplvs.orig.png",
       dist_samplvs.orig,
       units = "mm",
       width = 183,
       height = 150)

## Correlation -----
# 1. Calculate the correlation matrix

correlation_matrix <- cor(X1_data)

# 2. Prepare the correlation matrix for ggplot2

correlation_df <- as.data.frame(correlation_matrix)
correlation_df$Var1 <- rownames(correlation_df) # Add a column for the first variable in the pair

correlation_long <- correlation_df |>
  pivot_longer(
    cols = -Var1, # Select all columns except 'Var1'
    names_to = "Var2", # New column for the second variable in the pair
    values_to = "Correlation" # New column for the correlation value
  ) |>
  # Ensure the variables are factors with consistent levels for proper ordering on axes
  mutate(
    Var1 = factor(Var1, levels = colnames(correlation_matrix)),
    Var2 = factor(Var2, levels = colnames(correlation_matrix))
  )

# 3. Create the correlation heatmap plot
plot_correlation_heatmap <-
  ggplot(correlation_long, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile(color = "white") +
   geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
   scale_fill_gradient2(low = "darkred", mid = "white", high = "darkblue",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlation\nCoefficient") +
  # Add labels and a title to the plot.
  labs(
    title = "Correlation Matrix of Sampled Input Parameters",
    x = "", # No x-axis label as parameter names serve this purpose
    y = ""  # No y-axis label as parameter names serve this purpose
  ) +
  # Apply a minimal theme for a clean look.
  theme_minimal() +
  # Further customize theme elements for aesthetics and readability.
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9), # Rotate x-axis labels
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5), # Center title
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    legend.position = "right",          # Position of the color legend
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  coord_fixed() # Ensures that tiles are square, making the heatmap visually proportionate.

print(plot_correlation_heatmap)

# Display the correlation plot
library(GGally)

plot_correlation_scatterplot_matrix <- ggpairs(
  data = X1_data,
  columns = 1:ncol(X1_data),
  title = "Scatterplot Matrix of Sampled Input Parameters",

  # --- Customize the Upper Triangle (Correlations) ---
  upper = list(
    continuous = wrap("cor",
                      method = "spearman",
                      size = 3.5,
                      color = "black")
  ),

  # --- Customize the Lower Triangle (Scatterplots) ---
  lower = list(
    continuous = wrap("points",
                      alpha = 0.3,
                      size = 0.5,
                      color = "steelblue")
  ),

  # --- Customize the Diagonal (Distributions) ---
  # diag = list(
  #   continuous = wrap("density",
  #                     alpha = 0.6,
  #                     fill = "lightblue",
  #                     color = "darkblue")
  # ),

  progress = FALSE
) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 9, face = "bold")
  )

# Display the generated plot
print(plot_correlation_scatterplot_matrix)

