# Nitrogen Model Sensitivity Analysis in R (Sobol' Indices)

# This script defines the nitrogen model and performs global sensitivity analysis
# using Sobol' indices via the 'sensitivity' package.

# Install and load the 'sensitivity' package if you haven't already:
# install.packages("sensitivity")
library(sensitivity)
library(tidyverse)
library(patchwork)

# 1. Define the Nitrogen Model Function ----
# The model calculates N based on various nitrogen parameters and input variables.
# N = βt NT + βCS MNCS + βCA MNCA + βudb MNudb + βm1 (M1+M2)/2 + βf0 F0 +
#     βf1 (F1+F2)/2 + βg0 G0 + βm1 (G1+G2)/2
#
# This function will be passed directly to the 'sobol' function as the 'model' argument.
# It expects a matrix 'X' where each column is an input variable and each row is a sample.
calculate_N_for_sobol <- function(X) {
  # Default Beta Parameter Values (treated as fixed for this sensitivity analysis)
  beta_t = 0.456793   # Total N in topsoil (0-25 cm) (ton N)
  beta_CS = 0.04957    # Mineral N spring in harvest year
  beta_CA = 0.157044   # Mineral N autumn in harvest year
  beta_udb = 0.038245  # Mineral N from grazing in harvest year
  beta_m1_M = 0.026499 # Mineral N added to previous/pre-previous crop (M1+M2)
  beta_f0 = 0.016314   # Nitrogen fixation in harvest year
  beta_f1 = 0.026499    # Nitrogen fixation in previous/pre-previous crop (F1+F2)
  beta_g0 = 0.014099   # Organic N spring in harvest year
  beta_m1_G = 0.0265   # Organic N added to previous/pre-previous crop (G1+G2)
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

# 2. Sampling ----
# These ranges are crucial for the global sensitivity analysis.
# Adjust them based on the actual variability observed in your dataset.
input_N <- c("NT", "MNCS", "G0", "MNCA", "MNudb", "F0", "M1", "M2", "G1", "G2", "F1", "F2")

# 3. Perform Global Sensitivity Analysis using Sobol' Indices

# Number of Monte Carlo simulations per variable (N in Sobol' terminology)
# A higher N means more accurate results but longer computation time.
# N = 2000 is a good starting point. Total runs will be N * (2k + 2)
n_simulations <- 2000

sasoutput_fil_N <- sasoutput_fil_NLES5 |> dplyr::select(input_N)

# 1. Prepare the data for row-wise sampling
# Select only the necessary columns as defined in input_N

data_for_sampling <- sasoutput_fil_N[, input_N, drop = FALSE]

# Remove rows with NA values in any of the selected input columns
data_for_sampling_complete <- na.omit(data_for_sampling)

# 2. Generate X1_data by sampling rows
# Sample row indices with replacement
sampled_indices_X1 <- sample(1:nrow(data_for_sampling_complete), size = n_simulations, replace = TRUE)
X1_data <- data_for_sampling_complete[sampled_indices_X1, , drop = FALSE]
# Column names are inherited and match input_N

# 3. Generate X2_data by sampling rows (independently from X1)
# Sample row indices with replacement, independently from the first sample
sampled_indices_X2 <- sample(1:nrow(data_for_sampling_complete), size = n_simulations, replace = TRUE)
X2_data <- data_for_sampling_complete[sampled_indices_X2, , drop = FALSE]
# Column names are inherited and match input_N

# Ensure they are data frames (should be by default from subsetting)
X1_data <- as.data.frame(X1_data)
X2_data <- as.data.frame(X2_data)
dim(X1_data)==dim(X2_data)

# 3. Perform Sobol' Sensitivity Analysis ----
sobol_results_n_model <- soboljansen(
  model = calculate_N_for_sobol,
  X1 = X1_data, # Use the explicitly named X1_data
  X2 = X2_data, # Use the explicitly named X2_data
  n = n_simulations,
  conf = 0.95#,
  #order = 2 # Request calculation of second-order indices
)

ggplot(sobol_results_n_model)

# Print the results
print(sobol_results_n_model$S)

print(sobol_results_n_model$T)

# --- 4. Visualization of Sobol Indices ----
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
  plot_sobol_data |>
  ggplot(aes(x = Variable, y = Value, fill = Index_Type)) +
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


Npoint3 <-
  uncplot + sobolplotN +plot_layout(ncol = 1,
                                  #widths = c(0.5,1.5)
                                  heights = c(1,1.5)
)

ggsave("Npoint3.png",Npoint3,
       width = 183, height =183*.8,
       units = "mm")

# 5. Visualization of Obs vs. Sampled input distribution -----

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
  select(all_of(input_N)) |> # Select only the 12 relevant parameter columns
  pivot_longer(
    cols = everything(),
    names_to = "Parameter",
    values_to = "Value"
  ) |>
  na.omit() |> # Remove any NA values, similar to how sampling handled them
  mutate(Source = "Original (dataset)") # Label this data set

# --- Combine both data frames into a single data frame ---
# This allows ggplot2 to plot both sets of densities in the same plot.
combined_plot_data <- bind_rows(sampled_data_long, original_data_long)


# --- Create the faceted density plot with both curves ---
dist_samplvs.orig <-
  combined_plot_data |>
  ggplot(aes(x = Value, color = Source,
             fill = Source)) +
  geom_density(alpha = 0.2, linewidth = 0.8) +
  facet_wrap( ~ Parameter, scales = "free", ncol = 4,
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

## Correlation sampled -----
# 1. Calculate the correlation matrix

corr_sampl_matrix <- cor(X1_data)

# 2. Prepare the correlation matrix for ggplot2

corr_sampl_df <- as.data.frame(corr_sampl_matrix)
corr_sampl_df$Var1 <- rownames(corr_sampl_df) # Add a column for the first variable in the pair

corr_sampl_long <- corr_sampl_df |>
  pivot_longer(
    cols = -Var1, # Select all columns except 'Var1'
    names_to = "Var2", # New column for the second variable in the pair
    values_to = "Correlation" # New column for the correlation value
  ) |>
  # Ensure the variables are factors with consistent levels for proper ordering on axes
  mutate(
    Var1 = factor(Var1, levels = colnames(corr_sampl_matrix)),
    Var2 = factor(Var2, levels = colnames(corr_sampl_matrix))
  )

# 3. Create the correlation heatmap plot
plot_corr_sampl_heatmap <-
  ggplot(corr_sampl_long, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile(color = "white") +
   geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
   scale_fill_gradient2(low = "darkred", mid = "white", high = "darkblue",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlation\nCoefficient") +
  # Add labels and a title to the plot.
  labs(
    title = "Sampled",
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

print(plot_corr_sampl_heatmap)

# Correlation observed
# 1. Calculate the correlation matrix

corr_obs_matrix <- cor(sasoutput_fil_N)

# 2. Prepare the correlation matrix for ggplot2

corr_obs_df <- as.data.frame(corr_obs_matrix)
corr_obs_df$Var1 <- rownames(corr_obs_df) # Add a column for the first variable in the pair

corr_obs_long <- corr_obs_df |>
  pivot_longer(
    cols = -Var1, # Select all columns except 'Var1'
    names_to = "Var2", # New column for the second variable in the pair
    values_to = "Correlation" # New column for the correlation value
  ) |>
  # Ensure the variables are factors with consistent levels for proper ordering on axes
  mutate(
    Var1 = factor(Var1, levels = colnames(corr_obs_matrix)),
    Var2 = factor(Var2, levels = colnames(corr_obs_matrix))
  )

# 3. Create the correlation heatmap plot
plot_corr_obs_heatmap <-
  ggplot(corr_obs_long, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = "darkblue",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlation\nCoefficient") +
  # Add labels and a title to the plot.
  labs(
    title = "Observed",
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


plot_corr_obs_heatmap+
  plot_corr_sampl_heatmap+
  plot_layout(ncol = 2, widths = c(.5, .5),
              guides = "collect")

