# Nitrogen Leaching Model (NLES5) Sensitivity Analysis in R (Sobol' Indices)

# This script performs GSA for the NLES5 model
# using Sobol' indices via the 'sensitivity' package.

# Load necessary packages
library(sensitivity)
library(tidyverse)
library(patchwork)
library(readxl)
library(purrr)
library(GGally)

# Labeller function for facet plots, using the 'translate' vector
facet_labels_expr <- as_labeller(translate)

# --- 1. Define the NLES5 Model Wrapper Function for Sobol' ----
# This function takes a matrix 'X' (where each row is a sample and each column is an input variable)
# and returns a vector of the model's output (L in this case).
# 'nles5' function are already defined in  R environment from script fgkNLES5. If you can not find it run that script

calculate_NLES5_for_sobol <- function(X) {
  # Extract input variables from the matrix X by column name
  # Ensure all variables expected by the NLES5 component functions are present in X.
  Y = X[, "Y"]
  NT = X[, "NT"]
  MNCS = X[, "MNCS"]
  MNCA = X[, "MNCA"]
  MNudb = X[, "MNudb"]
  M1 = X[, "M1"]
  M2 = X[, "M2"]
  F0 = X[, "F0"]
  F1 = X[, "F1"]
  F2 = X[, "F2"]
  G0 = X[, "G0"]
  G1 = X[, "G1"]
  G2 = X[, "G2"]
  WC = X[, "WC"]

  M = X[, "M"]
  W = X[, "W"]
  MP = X[, "MP"]
  WP = X[, "WP"]

  CU = X[, "CU"]
  jbnr = X[, "jbnr"]

  AAa = X[, "AAa"]
  AAb = X[, "AAb"]
  APb = X[, "APb"]
  #p2 = X[, "p2"]
  #p3 = X[, "p3"]

  #EEA = X[, "EEA"]
  #Fdato = X[, "Fdato"]
  #EMA = X[, "EMA"]
  #ETS = X[, "ETS"]
  #EPJ = X[, "EPJ"]

  # Initialize a vector to store the primary model output (L) for each sample
  L_vec <- numeric(nrow(X))

  # Loop through each row (sample) in the input matrix X
  for (i in 1:nrow(X)) {
    # Calculate intermediate components using  provided functions
    Ntheta_val <- N_func(
      NT = NT[i], MNCS = MNCS[i], MNCA = MNCA[i], MNudb = MNudb[i],
      M1 = M1[i], M2 = M2[i], F0 = F0[i], F1 = F1[i], F2 = F2[i],
      G0 = G0[i], G1 = G1[i], G2 = G2[i], WC = WC[i]
    )

    C_val <- C_func(M = M[i], W = W[i], MP = MP[i], WP = WP[i])

    S_val <- S_func(CU = CU[i])

    P_val <- P_func(jbnr = jbnr[i], AAa = AAa[i], AAb = AAb[i], APb = APb[i])

    #Psas_val <- Psas_func(jbnr = jbnr[i], AAa = AAa[i], AAb = AAb[i], APb = APb[i], p2 = p2[i], p3 = p3[i])

    # Run the main nles5 function with the calculated components and other inputs
    nles5_results <- nles5(
      Y = Y[i],
      Ntheta = Ntheta_val,
      C = C_val,
      P = P_val,
      #Psas = Psas_val,
      S = S_val#,
      #EEA = EEA[i],
      #Fdato = Fdato[i],
      #EMA = EMA[i],
      #ETS = ETS[i],
      #EPJ = EPJ[i]
    )

    L_vec[i] <- ifelse(is.na(nles5_results$L), 0, ifelse(nles5_results$L > 0, nles5_results$L, 0)) #this is also like his in the SAS script this is something to be conccerned about
  }
  return(L_vec)
}

# --- 2. Prepare Input Data for sensitivity : Load and Map 'sasoutput' for Sampling ----

# Create a data frame 'sasoutput_fil_NLES5' that maps columns from
# 'sasoutput' to the expected input names of the NLES5 model.

sasoutput_fil_NLES5 <- data.frame(
  Y = sasoutput$Indb_aar,
  NT = sasoutput$TN,
  MNCS = sasoutput$NS,
  MNCA = sasoutput$na,
  MNudb = sasoutput$Nudb,
  M1 = sasoutput$nlevelMin1,
  M2 = sasoutput$nlevelMin2,
  F0 = sasoutput$NlevelFix0,
  F1 = sasoutput$NlevelFix1,
  F2 = sasoutput$NlevelFix2,
  G0 = sasoutput$NlevelGod0,
  G1 = sasoutput$NlevelGod1,
  G2 = sasoutput$NlevelGod2,
  WC = sasoutput$Vafgr_Kappa, # Should be 1 or 2 as per N_func logic

  M = sasoutput$Mau,
  W = sasoutput$Vau,
  MP = sasoutput$Mfu,
  WP = sasoutput$Vfu,

  CU = sasoutput$CU,
  jbnr = sasoutput$jbnr,

  AAa = sasoutput$d1,
  AAb = sasoutput$d2,
  APb = sasoutput$d3#,
  #p2 = sasoutput$p2,
  #p3 = sasoutput$p3,

  # These are fixed for this analysis, based on  prior defaults
  #EEA = 0,
  #Fdato = 1,
  #EMA = 0,
  #ETS = 0,
  #EPJ = 1/11
)

# Define the complete list of input parameter names for NLES5
input_names_nles5 <- colnames(sasoutput_fil_NLES5)

# Define n of the simulation runs note that is higher than the number of rows in sasoutput_fil_NLES5 so some observations will be there more than 1 time
n_simulations <- 1000

# 1. Prepare the data for row-wise sampling
# Select only the necessary columns as defined in input_N

data_for_sampling_complete <- na.omit(sasoutput_fil_NLES5[, input_names_nles5, drop = FALSE])

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

# --- 3. Perform Sobol' Sensitivity Analysis using 'soboljansen' ----

sobol_results_nles5_model <- soboljansen(
  model = calculate_NLES5_for_sobol, #  wrapper function for NLES5
  X1 = X1_data[complete.cases(X1_data),],                      # First set of input samples
  X2 = X2_data,                      # Second set of input samples
  n = n_simulations,                 # Number of simulations (as defined above)
  conf = 0.95                        # Confidence level for indices
  #,order = 2, nboot = 100
)

ggplot(sobol_results_nles5_model)

# Print the calculated Sobol' indices
print(sobol_results_nles5_model$S)
print(sobol_results_nles5_model$T)
# Note: Second-order indices (S2) are not typically computed by default with 'soboljansen'
# when X1 and X2 are pre-generated. If you need S2, the 'sobol' function might be more
# appropriate if  inputs were purely continuous and defined by ranges.

# 4. Visualization of Sobol' Indices and Uncertainty for NLES5 ----


# Extract first-order indices
first_order_indices_nles5 <- data.frame(
    Variable = rownames(sobol_results_nles5_model$S),
    S1 = sobol_results_nles5_model$S[, "original"],
    S1_conf_low = sobol_results_nles5_model$S[, "min. c.i."],
    S1_conf_high = sobol_results_nles5_model$S[, "max. c.i."]
  )

# Extract total-order indices
total_order_indices_nles5 <- data.frame(
  Variable = rownames(sobol_results_nles5_model$T),
  ST = sobol_results_nles5_model$T[, "original"],
  ST_conf_low = sobol_results_nles5_model$T[, "min. c.i."],
  ST_conf_high = sobol_results_nles5_model$T[, "max. c.i."]
)

# Combine both sets of indices into a single data frame
sobol_df_nles5 <- dplyr::left_join(first_order_indices_nles5, total_order_indices_nles5, by = "Variable")


# Invert the vector
translate_inverse <- c(
  "Y" = "Indb_aar",
  "MNCS" = "NS",
  "MNCA" = "na",
  "M1" = "nlevelMin1",
  "M2" = "nlevelMin2",
  "F0" = "NlevelFix0",
  "F1" = "NlevelFix1",
  "F2" = "NlevelFix2",
  "G0" = "NlevelGod0",
  "G1" = "NlevelGod1",
  "G2" = "NlevelGod2",
  "MNudb" = "Nudb",
  "NT" = "TN",
  "WC" = "Vafgr_Kappa",
  "M" = "Mau",
  "MP" = "Mfu",
  "W" = "Vau",
  "WP" = "Vfu",
  "AAa" = "d1",
  "AAb" = "d2",
  "APb" = "d3",
  "jbnr" = "SoilG",   # Third instance for "Soil group"
  "Clay" = "CU"
)



# Reshape data for ggplot2 and add descriptive labels for plotting
plot_sobol_data_nles5 <- sobol_df_nles5 |>
  tidyr::gather(key = "Index_Type", value = "Value", S1, ST) |>
  dplyr::mutate(
    Conf_Low = ifelse(Index_Type == "S1", S1_conf_low, ST_conf_low),
    Conf_High = ifelse(Index_Type == "S1", S1_conf_high, ST_conf_high),
    # Use the 'translate' vector to get readable labels for the input parameters
    Input = recode(Variable, !!!translate_inverse)) |>
    mutate(Component = recode(Input, !!!asignation)) |>
  arrange(desc(Component))

# --- Create the Sobol' Indices Plot ---
sobolplotNLES5 <-
  ggplot(plot_sobol_data_nles5, aes(x = reorder(Variable, Value), y = Value,
                                    group = Index_Type, fill = Component)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = Value, ymax = Conf_High),
                position = position_dodge(width = 0.8), width = 0.2) +
  labs(title = "Sensitivity Estimates for NLES5 Model",
       x = "Input Parameter",
       y = "Sensitivity Index",
       fill = "Sobol' Index") +
  scale_fill_manual(values = fixed_colors)+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

sobolplotNLES5


plot_sobol_data_nles5_ordered <- plot_sobol_data_nles5 %>%
  left_join(
    plot_sobol_data_nles5 %>%
      filter(Index_Type == "ST") %>%
      select(Variable, ST_Value_for_sort = Value) %>%
      distinct(), # Ensures one ST value per Variable for sorting
    by = "Variable"
  ) %>%
  # If an Variable somehow doesn't have an ST value, assign -Inf to sort it at the bottom
  mutate(ST_Value_for_sort = if_else(is.na(ST_Value_for_sort), -Inf, ST_Value_for_sort)) %>%
  # Arrange by Component, then by ST value (descending), then by Variable name as a final tie-breaker
  arrange(Component, desc(ST_Value_for_sort), Variable) %>%
  # Create the 'Variable' factor with levels in the new desired order
  mutate(Variable = factor(Variable, levels = unique(.$Variable)))

sobolplotNLES5_vertical_ordered <-
  ggplot(plot_sobol_data_nles5_ordered, # Use the new ordered data
         aes(
           x = Variable, # 'Input' is now an ordered factor
           y = Value,#ifelse(Value>0,Value,0),
           fill = Component,
           alpha = Index_Type
         )) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.9),
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(ymin = Value, #ifelse(Value>0,Value,0),
        ymax = ifelse(Value>0,Conf_High,Conf_Low), group = Index_Type), # Corrected ymin
    position = position_dodge(width = 0.9),
    width = 0.25,
    color = "gray60",
    linewidth = 0.5
  ) +
  scale_fill_manual(values = fixed_colors) + # Your custom colors
  scale_alpha_manual(
    values = c("S1" = 0.25, "ST" = 1.0),
    name = "Sobol' Index",
    labels = c("S1 (First-order)", "ST (Total-order)")
  ) +
  coord_flip() + # This makes the plot vertical!
  labs(
    title = "Sensitivity estimates",
    x = "Input parameter", # This label will now appear on the Y-axis
    y = "Index Value", # This label will now appear on the X-axis
    fill = "Component"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = rel(0.9)), # Was axis.text.x, now for Input labels
    axis.text.x = element_text(size = rel(0.9)), # Was axis.text.y, now for Value labels
    # No angle needed for axis.text.y (Input labels) if they fit well
    # If Input labels are long, you might adjust their properties or plot margins
    axis.title.y = element_text(size = rel(1.0), margin = margin(r = 10)), # Was axis.title.x
    axis.title.x = element_text(size = rel(1.0), margin = margin(t = 10)), # Was axis.title.y
    plot.title = element_text(
      size = rel(1.2),
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "right",
    legend.title = element_text(size = rel(1.0), face = "bold"),
    legend.text = element_text(size = rel(0.9)),
    panel.grid.major.y = element_blank(), # Removes horizontal grid lines (along parameters)
    panel.grid.major.x = element_line(color = "gray90"), # Keeps vertical grid lines (along values)
    panel.grid.minor = element_blank()
  ) +
  guides(alpha = guide_legend(override.aes = list(fill = 'gray50')))

 print(sobolplotNLES5_vertical_ordered)

 ggsave("sobolplotNLES5_vertical_ordered.png",
        sobolplotNLES5_vertical_ordered,
        width = 183, height = 240, units = "mm", dpi = 300)
 #

#
#  sobolplotNLES5_vertical_ordere_point <-
#    ggplot(plot_sobol_data_nles5_ordered, # Use the new ordered data
#           aes(
#             x = Variable, # 'Input' is now an ordered factor
#             y = Value,
#             fill = Component,
#             color = Component,
#             alpha = Index_Type,
#             shape = Index_Type # Add shape aesthetic for points
#           )) +
#    geom_point(
#      stat = "identity",
#      position = position_dodge(width = 0.9),
#      linewidth = 0.5,
#      size = 4 # Adjust point size for better visibility
#    ) +
#    geom_errorbar(
#      aes(ymin = Conf_Low, ymax = Conf_High, group = Index_Type), # Corrected ymin
#      position = position_dodge(width = 0.9),
#      width = 0.25,
#      color = "gray50",
#      linewidth = 0.5
#    ) +
#    scale_fill_manual(values = fixed_colors) + # Your custom colors
#    scale_color_manual(values = fixed_colors) +
#    scale_alpha_manual(
#      values = c("S1" = 0.4, "ST" = 1.0),
#      name = "Sobol' Index",
#      labels = c("S1 (First-order)", "ST (Total-order)")
#    ) +
#    coord_flip() + # This makes the plot vertical!
#    labs(
#      title = "Sensitivity estimates NLES5 Model",
#      x = "Input parameter", # This label will now appear on the Y-axis
#      y = "Sensitivity Index Value", # This label will now appear on the X-axis
#    ) +
#    theme_minimal(base_size = 10) +
#    theme(
#      axis.text.y = element_text(size = rel(0.9)), # Was axis.text.x, now for Input labels
#      axis.text.x = element_text(size = rel(0.9)), # Was axis.text.y, now for Value labels
#      # No angle needed for axis.text.y (Input labels) if they fit well
#      # If Input labels are long, you might adjust their properties or plot margins
#      axis.title.y = element_text(size = rel(1.0), margin = margin(r = 10)), # Was axis.title.x
#      axis.title.x = element_text(size = rel(1.0), margin = margin(t = 10)), # Was axis.title.y
#      plot.title = element_text(
#        size = rel(1.2),
#        hjust = 0.5,
#        face = "bold"
#      ),
#      legend.position = "right",
#      legend.title = element_text(size = rel(1.0), face = "bold"),
#      legend.text = element_text(size = rel(0.9)),
#      panel.grid.major.y = element_blank(), # Removes horizontal grid lines (along parameters)
#      panel.grid.major.x = element_line(color = "gray90"), # Keeps vertical grid lines (along values)
#      panel.grid.minor = element_blank()
#    ) +
#    guides(alpha = guide_legend(override.aes = list(fill = 'gray50')))
#
#  print(sobolplotNLES5_vertical_ordere_point)
#
#
# # --- Create the Uncertainty Plot (Histogram of Model Output) ---
# # This shows the distribution of the NLES5 model's output (L_nuar).
# uncplotNLES5 <-
#   ggplot(data.frame(y = sobol_results_nles5_model$y), aes(x = y)) +
#   geom_histogram(binwidth = diff(range(sobol_results_nles5_model$y, na.rm = TRUE))/30, # Dynamic binwidth
#                  fill = "darkseagreen", color = "white", alpha = 0.8) +
#   labs(title = "Uncertainty Estimates for NLES5 Output",
#        x = "NLES5 Output(L)",
#        y = "Absolute Frequency") +
#   theme_minimal() +
#   theme(plot.title = element_text(size = 12, hjust = 0.5),
#         axis.title = element_text(size = 10),
#         axis.text = element_text(size = 9))
#
# # --- Combine Plots using 'patchwork' ---
# combined_nles5_plots <- uncplotNLES5 + sobolplotNLES5 +
#   plot_layout(ncol = 1, widths = c(0.5, 1.5)) # Adjust layout for better presentation
#
# print(combined_nles5_plots)
#
# # Save the combined plot to a file
# ggsave("NLES5_Sensitivity_Uncertainty.png", combined_nles5_plots,
#        width = 250, height = 150, units = "mm", dpi = 300)
#
#   ## 5. Comparison of Sampled vs. Original Input Distributions ----
#
#   # This section helps you visualize if the sampling process for sensitivity analysis
#   # (from  'sasoutput' data) accurately represents the original distributions.
#
#   # --- 1. Prepare Sampled Data (from X1_data) in Long Format ---
sampled_data_long_nles5 <- X1_data |>
  pivot_longer(
    cols = everything(), # Selects all input parameter columns
    names_to = "Parameter", # Column for parameter names
    values_to = "Value"     # Column for sampled values
  ) |>
  mutate(Source = "Sampled") # Label this dataset

# --- 2. Prepare Original Data (from sasoutput_fil_NLES5) in Long Format ---
original_data_long_nles5 <- sasoutput_fil_NLES5 |>
  select(all_of(input_names_nles5)) |> # Select only the relevant input columns
  pivot_longer(
    cols = everything(),
    names_to = "Parameter",
    values_to = "Value"
  ) |>
  na.omit() |> # Remove any NA values from original data for fair comparison
  mutate(Source = "Observed") # Label this dataset

# --- 3. Combine both data frames into a single data frame ---
combined_plot_data_nles5 <- bind_rows(sampled_data_long_nles5, original_data_long_nles5)

# --- 4. Create the Faceted Density Plot ---
# This plot shows the density distribution for each input parameter,
# comparing the sampled data (used for GSA) against  original dataset.
dist_samplvs_orig_nles5_cont <-
  combined_plot_data_nles5 |>
  filter(!Parameter %in% c("M", "MP", "W", "WP", "WC", "jbnr","Y")) |>
  mutate(Component = recode(recode(Parameter, !!!translate_inverse), !!!asignation)) |>
  arrange(Value) |>
  ggplot(aes(x = Value, linetype = Source, fill = Component, colour = Component)) +
  geom_density(alpha = 0.2, linewidth = 0.8) + # 'colour' and 'linetype' apply to the line, 'fill' to the area
  facet_wrap(~ Parameter , scales = "free", ncol = 4) +
  labs(
    x = "Input Value",
    y = "Density",
    title = "Continous inputs" # Added a main title example
  ) +
  # Use the 'name' argument in scales for legend titles
  scale_fill_manual(
    name = "Component Group", # Legend title for fill
    values = fixed_colors
  ) +
  scale_colour_manual(
    name = "Component Group", # Identical name to merge with fill legend
    values = fixed_colors
  ) +
  scale_linetype_manual(
    name = "Source",    # Legend title for linetype
    values = c("Observed" = "solid", "Sampled" = "dashed"), # Example linetypes
    labels = c("Observed", "Sampled") # Optional: Shorter labels for the legend
  ) +
  theme_light() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), # Adjusted title size
    strip.text = element_text(size = 9, face = "bold"), # Made strip text a bit larger and bold
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.position = "right",
    legend.box = "vertical", # Arranges multiple legends vertically
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

print(dist_samplvs_orig_nles5_cont)

data_for_plot_cat <- combined_plot_data_nles5 |>
  filter(Parameter %in% c("M", "MP", "W", "WP", "WC", "jbnr")) |>
  mutate(
    Value = as.factor(Value), # Treat 'Value' as a factor for categorical plotting
    Component_temp = recode(recode(Parameter, !!!translate_inverse), !!!asignation) # Temp component for ordering
  ) |>
  # Calculate N for each Parameter-Source-Value combination
  count(Parameter, Component_temp, Source, Value, name = "N") |>
  # Calculate total N for each Parameter-Source (for proportions)
  group_by(Parameter, Component_temp, Source) |>
  mutate(Proportion = N / sum(N)) |>
  ungroup() |>
  # Create the interaction term for x-axis labels. Using \n for a line break if possible.
  mutate(
    x_axis_label = interaction(Source, Value,  sep = "-", drop = TRUE),
    Component = Component_temp # Final component column
  ) |>
  # Order data for consistent factor levels if Parameter is made a factor for faceting
  arrange(Component, Parameter, Source, Value) |>
  # Ensure Parameter is a factor with levels in the arranged order for facet_wrap
  mutate(Parameter = factor(Parameter, levels = unique(Parameter)))

dist_samplvs_orig_nles5_cat <-
  data_for_plot_cat |>
  ggplot(aes(
    x = x_axis_label,         # Use the pre-created interaction term for x-axis
    y = Proportion,           # Use the calculated Proportion
    fill = Component,
    colour = Component,       # For bar outline
    #linetype = Source         # For bar outline style
    ,alpha=Source
  )) +
  geom_col(linewidth = 0.5) + # Using geom_col. Adjusted alpha for better visibility.
  # linewidth for the border.
  facet_wrap(~ Parameter, scales = "free_x", ncol = 2) + # Use "free_x" as category names/counts differ
  scale_y_continuous(
    name = "Relative Frequency ",
    labels = scales::percent_format(accuracy = 1) # Format y-axis as percentage
  ) +
  labs(
    x = " ", # Updated x-axis label
    title = "Categorical inputs" # Example title
    # Legend titles are best set in their respective scales below
  ) +
  scale_fill_manual(name = "Component Group", values = fixed_colors) +
  scale_colour_manual(name = "Component Group", values = fixed_colors) + # For bar outlines

  theme_light() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(), # Often good to remove vertical grid lines for categorical plots
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 9, face = "bold"), # Facet labels
    axis.title = element_text(size = 11),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7), # ROTATED X-AXIS LABELS & smaller
    axis.text.y = element_text(size = 9),
    legend.position = "right", # No legend for this plot
  )

print(dist_samplvs_orig_nles5_cat)


combined_plots_layout <- dist_samplvs_orig_nles5_cat +
  dist_samplvs_orig_nles5_cont +
  plot_layout(
    guides = 'collect', # Collects common legends
    heights = c(1, 3)    # Example: give equal height; adjust as needed
    # You can also use ncol=1 for a single column layout if preferred
  )


# Save the distribution comparison plot
ggsave("NLES5_dist_samplvs_orig.png",
       dist_samplvs_orig_nles5_cont,
       units = "mm",
       width = 183, # Adjusted width for more facets
       height = 250, # Adjusted height
       dpi = 300)


## 6. Correlation Analysis of Sampled Input Parameters

  # Understanding the correlations between  input parameters is crucial,
  # as strong correlations can affect the interpretation of Sobol' indices.

  # 1. Calculate the Pearson correlation matrix for the sampled input data (X1_data)
corr_sampled_matrix_nles5 <- cor(X1_data)

# 2. Prepare the correlation matrix for a heatmap visualization with ggplot2
corr_sampled_df_nles5 <- as.data.frame(corr_sampled_matrix_nles5)
corr_sampled_df_nles5$Var1 <- rownames(corr_sampled_df_nles5) # Add a column for the first variable

corr_sapled_long_nles5 <- corr_sampled_df_nles5 |>
  pivot_longer(
    cols = -Var1, # Select all columns except 'Var1'
    names_to = "Var2", # New column for the second variable
    values_to = "Correlation" # New column for the correlation value
  ) |>
  # Ensure variables are factors with consistent levels for proper ordering on axes
  mutate(
    Var1 = factor(Var1, levels = colnames(corr_sampled_matrix_nles5)),
    Var2 = factor(Var2, levels = colnames(corr_sampled_matrix_nles5))
  )

# 3. Create the correlation heatmap plot
plot_corr_sampled_heatmap_nles5 <-
  ggplot(corr_sapled_long_nles5, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(color = "white") + # Draws the colored tiles with white borders
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 2) + # Add correlation values as text
  scale_fill_gradient2(low = "darkred", mid = "white", high = "darkblue",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlation\nCoefficient") +
  labs(title = "Sampled",
       x = "", y = "") + # No axis labels as parameter names serve this purpose
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6), # Rotate x-axis labels
    axis.text.y = element_text(size = 6),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), # Center title
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  ) +
  coord_fixed() # Ensures square tiles

print(plot_corr_sampled_heatmap_nles5)


# 1. Calculate the Pearson correlation matrix for the sampled input data (X1_data)
corr_matrix_nles5 <- cor(sasoutput_fil_NLES5)

# 2. Prepare the correlation matrix for a heatmap visualization with ggplot2
corr_obs_df_nles5 <- as.data.frame(corr_matrix_nles5)
corr_obs_df_nles5$Var1 <- rownames(corr_obs_df_nles5) # Add a column for the first variable

corr_obs_long_nles5 <- corr_obs_df_nles5 |>
  pivot_longer(
    cols = -Var1, # Select all columns except 'Var1'
    names_to = "Var2", # New column for the second variable
    values_to = "Correlation" # New column for the correlation value
  ) |>
  # Ensure variables are factors with consistent levels for proper ordering on axes
  mutate(
    Var1 = factor(Var1, levels = colnames(corr_matrix_nles5)),
    Var2 = factor(Var2, levels = colnames(corr_matrix_nles5))
  )

# 3. Create the correlation heatmap plot
plot_corr_obs_heatmap_nles5 <-
  ggplot(corr_obs_long_nles5, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(color = "white") + # Draws the colored tiles with white borders
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 2) + # Add correlation values as text
  scale_fill_gradient2(low = "darkred", mid = "white", high = "darkblue",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlation\nCoefficient") +
  labs(title = "Observed",
       x = "", y = "") + # No axis labels as parameter names serve this purpose
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6), # Rotate x-axis labels
    axis.text.y = element_text(size = 6),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), # Center title
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  ) +
  coord_fixed() # Ensures square tiles

corr_heat_obsvssampl <-
plot_corr_sampled_heatmap_nles5+
  plot_corr_obs_heatmap_nles5+
  plot_layout(ncol = 2, widths = c(.5, .5),
              guides = "collect")

ggsave("corr_heat_obsvssamp.png",
       corr_heat_obsvssampl,
       units = "mm",
       width = 250, # Adjusted width for more facets
       height = 100, # Adjusted height
       dpi = 300)

