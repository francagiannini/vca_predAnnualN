# Load necessary packages
library(tidyverse)
library(readxl)
library(randomForest)
library(xgboost)
library(caret)

# 1. Data Loading ----
predus <- readRDS("pred_NLES5NUAR.RDS")

summary(predus$L)

# Create the density plot for L
density_plot_L <- ggplot(predus, aes(x = L)) +
  geom_density(fill = "#377ea2", color = "#377ea2", alpha = 0.5) + # Using a custom light orange color
  labs(
    #title = "Distribution of Predicted annual Leaching",
    x = "NLES5 predicted annual leaching (L)",
    y = "Density"
  ) +
  theme_minimal() + # Consistent with your previous plots
  theme(
    #plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
    #axis.title = element_text(face = "bold"), # Bold axis titles
    panel.grid.major = element_line(colour = "grey90")#, # Lighter grid lines
    #panel.grid.minor = element_blank(colour = "grey90") # Remove minor grid lines
  )+
  scale_x_continuous(breaks = seq(0, 300, 50), limits = c(0,263.65)) # Adjust x-axis breaks

print(density_plot_L)

# You can save the plot if needed
ggsave("density_plot_L.png",
       plot = density_plot_L,
       units = "mm",
       width = 150,
       height = 100,
       dpi = 300)

table(predus$SoilG)
# 1.1. Renaming and grouping ----
# Is d or p the percolation groups why some are negative?
# hist(predus$d2)
# hist(predus$p2)
# plot(predus$d2, predus$p2)

# 1.3 Select NLES5 inputs (i.e. variables) for analysis ----

predus_fil <- predus |>
  select(
    L,    # Predicted values (response for RF/XGBoost)
    NS, na, nlevelMin1, nlevelMin2, NlevelFix0, NlevelFix1, NlevelFix2,
    NlevelGod0, NlevelGod1, NlevelGod2, Nudb, TN, Vafgr_Kappa, # Nitrogen related terms
    Mau, Mfu, Vau, Vfu, # Crop related terms
    CU,             # Clay content (Soil)
    Indb_aar,     # Year trend
    d1,d2, d3,
    SoilG # Percolation/Drainage related terms
  ) |> mutate(SoilG= recode_factor(SoilG, "S"=0, "C"=1)) |>
  mutate(Mau=as.factor(Mau),
         Mfu=as.factor(Mfu),
         Vau=as.factor(Vau),
         Vfu=as.factor(Vfu),
         SoilG= as.factor(SoilG),
         )

# --- 2. Random Forest Surrogate Model (Optimized for High Fidelity i.e, over fitting) ----

set.seed(151190) # For reproducibility

rf_surrogate_model <- randomForest(
  L ~ .,
  data = predus_fil,
  ntree = 5000,      # Increase number of trees further (e.g., to 5000)
  mtry = ncol(predus_fil) - 1, # Set mtry to all predictors (p)
  nodesize = 1,      # Allow single-observation terminal nodes
  maxnodes = NULL,   # Let nodesize control tree depth (no explicit max)
  sampsize = nrow(predus_fil), # Use entire dataset for each tree (pure bagging mode)
  importance = TRUE  # Still useful for initial importance checks
)

print(rf_surrogate_model)

library(ranger)

rf_ranger_model <- ranger(
  formula = L ~ .,
  data = predus_fil,
  num.trees = 5000,
  mtry = ncol(predus_fil) - 1,
  min.node.size = 1,
  sample.fraction = 1, # Use entire dataset for each tree (pure bagging)
  importance = "impurity"
)

print(rf_ranger_model)

# Compare the two models
plot(rf_surrogate_model$predicted,rf_ranger_model$predictions)


cat("R-squared of High-Fidelity Random Forest on training data:",
    rf_surrogate_model$rsq[length(rf_surrogate_model$rsq)], "\n")

# Evaluate performance on the training data (should be very high)
rf_predictions_train <- predict(rf_surrogate_model, newdata = predus_fil)
rf_rmse_train <- sqrt(mean((predus_fil$L - rf_predictions_train)^2))
cat("RMSE of High-Fidelity Random Forest on training data:", rf_rmse_train, "\n")

# You might want to save this model if it's the one you'll use for Sobol
saveRDS(rf_surrogate_model, "rf_surrogate_model_max_fidelity.rds")


# --- 3. XGBoost Surrogate Model (Optimized for High Fidelity) ----

## --- Create a dummy variables for categories ---
# Create a dummyVars object

dummy_model <- dummyVars(~ ., sep="_",data = predus_fil)

# Apply the dummy encoding
ohe_data <- data.frame(predict(dummy_model, newdata = predus_fil))

set.seed(151190) # For reproducibility

# Prepare data for XGBoost (matrix format)
X_train <- as.matrix(ohe_data |> select(-L))
y_train <- ohe_data$L

# Define parameters for high fidelity
xgb_params_high_fidelity <- list(
  objective = "reg:squarederror", # Regression task
  eval_metric = "rmse",         # Root Mean Squared Error
  nrounds = 5000,               # Very high number of boosting rounds
  eta = 0.005,                   # Very small learning rate
  max_depth = 7,               # Deep trees to capture interactions
  subsample = 0.7,              # Not Use all data for each tree
  colsample_bytree = 0.8,       # Not Use all features for each tree
  min_child_weight = 0.8,         # Allow very small child weights
  gamma = 0,                    # No minimum loss reduction for splits
  lambda = 0,                   # No L2 regularization
  alpha = 0                     # No L1 regularization
)

xgb_surrogate_model <- xgb.train(
  params = xgb_params_high_fidelity,
  data = xgb.DMatrix(X_train, label = y_train),
  nrounds = xgb_params_high_fidelity$nrounds,
  verbose = 0 # Set to 1 to see training progress
)

# Evaluate performance on the training data (should be very high)
xgb_predictions_train <- predict(xgb_surrogate_model, xgb.DMatrix(X_train))
xgb_rmse_train <- sqrt(mean((y_train - xgb_predictions_train)^2))

cat("RMSE of High-Fidelity XGBoost on training data:", xgb_rmse_train, "\n")
cat("XGBoost R-squared (Training):", cor(xgb_predictions_train, ohe_data$L), "\n")

# You might want to save this model
saveRDS(xgb_surrogate_model, "xgb_surrogate_model_high_fidelity.rds")
xgb.save(xgb_surrogate_model, "xgb_surrogate_model_high_fidelity.model")
#xgb_surrogate_model <- readRDS("xgb_surrogate_model_high_fidelity.rds")

#xgb_surrogate_model <- xgb.load("xgb_surrogate_model_high_fidelity.model")

# --- Comparison ---
cat("\n--- Training Performance Comparison ---\n")
cat("RF R-squared (Training):", rf_surrogate_model$rsq[length(rf_surrogate_model$rsq)], "\n")
cat("XGBoost R-squared (Training):", cor(xgb_predictions_train, ohe_data$L), "\n")
cat("RF RMSE (Training):", rf_rmse_train, "\n")
cat("XGBoost RMSE (Training):", xgb_rmse_train, "\n")

plot(xgb_predictions_train, ohe_data$L)
abline(0,1, col="orange3")
plot(rf_surrogate_model$predicted, ohe_data$L)
plot(rf_ranger_model$predictions, ohe_data$L)
abline(0,1)
plot(rf_ranger_model$predictions, xgb_predictions_train)

# Create a data frame for observed values
df_observed <- data.frame(
  L = predus_fil$L,
  Type = "Predicted by NLES5"
)

# Create a data frame for predicted values
df_predicted <- data.frame(
  L = xgb_predictions_train,
  Type = "Predicted by XGBoost"
)

# Combine the two data frames
combined_df <- rbind(df_observed, df_predicted)

# Ensure 'Type' is a factor for proper plotting and legend ordering
combined_df$Type <- factor(combined_df$Type, levels = c("Predicted by XGBoost", "Predicted by NLES5"))

# 2. Create the combined density plot
combined_density_plot <- ggplot(combined_df, aes(x = L, fill = Type, color = Type)) +
  geom_density(alpha = 0.2, linewidth = 0.3) + # Use linewidth for the outline
  scale_fill_manual(values = c("Predicted by NLES5" = "darkgreen", "Predicted by XGBoost" = "steelblue")) + # Custom fill colors
  scale_color_manual(values = c("Predicted by NLES5" = "darkgreen", "Predicted by XGBoost" = "darkblue")) + # Custom outline colors
  labs(
    x = "Annual N leaching (L)",
    y = "Density",
    fill = "Data Type",
    color = "Data Type"
  ) +
  theme_minimal() +
  theme(
    #plot.title = element_text(hjust = 0.5, face = "bold"),
    #axis.title = element_text(face = "bold"),
    legend.position = "bottom", # Place legend at the bottom
    legend.title = element_blank(), # Remove legend title for brevity
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank()
  )

print(combined_density_plot)

ggsave("density_plot_L.png",
       plot = combined_density_plot,
       units = "mm",
       width = 150,
       height = 100,
       dpi = 300)

# 4. Results and visualization ----
# Variable Importance
## --- 4.1 Extract Feature Importance ----
# The xgb.importance function requires the feature names.
importance_matrix <- xgb.importance(
  feature_names = xgb_surrogate_model$feature_names, # Pass the column names of your training data matrix
  model = xgb_surrogate_model
)

# Print the importance matrix (data.table)
print(importance_matrix)

xgb.plot.importance(importance_matrix = importance_matrix, top_n = 55) # Plot top N features


prefixes <- c("Mau", "Mfu", "Vau", "Vfu","SoilG", "d1","d2","d3",
              "Indb_aar", "NS", "na", "nlevelMin1", "nlevelMin2",
              "NlevelFix0", "NlevelFix1", "NlevelFix2",
              "NlevelGod0", "NlevelGod1", "NlevelGod2",
              "Nudb", "TN", "CU", "Vafgr_Kappa")


importance_matrix_processed <- importance_matrix |>
  mutate(BaseFeature = str_extract(Feature, paste0(
    "^(",
    paste(prefixes, collapse = "|"),
    "|[A-Za-z]+)(?=\\d|_|$)"
  )))


## 4.2. Renaming and grouping ----
# Define the mapping for the crop prefixes
translate <- c(
  #𝐿 = 𝜏(𝑌 – 1991) + {(𝜇 + 𝜃𝑖𝑁 + 𝐶)𝜅}(𝑃 𝑆)p

  # Y
  "Indb_aar"="Y",

  # N
  "NS"="MNCS",
  "na"="MNCA",
  "nlevelMin1"="M1",
  "nlevelMin2"="M2",
  "NlevelFix0"="F0",
  "NlevelFix1"="F1",
  "NlevelFix2"="F2",
  "NlevelGod0"="G0",
  "NlevelGod1"="G1",
  "NlevelGod2"="G2",
  "Nudb"="MNudb",
  "TN"="NT",
  "Vafgr_Kappa"="WC",

  # C
  "Mau" = "M",
  "Mfu" = "MP",
  "Vau" = "W",
  "Vfu" = "WP",

  #P
  "d1"="AAa",
  "d2"="AAb",
  "d3"="APb",
  "SoilGS"="jbnr",
  "SoilGC"="jbnr",
  "SoilG"="jbnr",

  #S
  "CU"= "Clay")

asignation <- c(
  # Y
  "Indb_aar"="Y",

  # N
  "NS"="N",
  "na"="N",
  "nlevelMin1"="N",
  "nlevelMin2"="N",
  "NlevelFix0"="N",
  "NlevelFix1"="N",
  "NlevelFix2"="N",
  "NlevelGod0"="N",
  "NlevelGod1"="N",
  "NlevelGod2"="N",
  "Nudb"="N",
  "TN"="N",
  "Vafgr_Kappa"="N",

  # C
  "Mau" = "C",
  "Mfu" = "C",
  "Vau" = "C",
  "Vfu" = "C",

  #P
  "d1"="P",
  "d2"="P",
  "d3"="P",
  "p2"="P",
  "p3"="P",
  "SoilGS"="P",
  "SoilGC"="P",
  "SoilG"="P",

  #S
  "CU"="S"
)

## 4.3 Overal results Input contribution -----

importance_rf <- as.data.frame(rf_surrogate_model$importance) |>
  bind_cols('BaseFeature'=row.names(as.data.frame(rf_surrogate_model$importance)),
            'ranger'=(rf_ranger_model$variable.importance)/sum(rf_ranger_model$variable.importance)*100) |>
  mutate(Input = recode(BaseFeature, !!!translate)) |>
  mutate(IncMSE_rel = `%IncMSE` / sum(`%IncMSE`) * 100,
         IncNodePurity_rel=IncNodePurity/sum(IncNodePurity) * 100)


# Group by the identified crop prefixes and sum their contributions
contributions <- importance_matrix_processed |>
  mutate(Component = recode(BaseFeature, !!!asignation),
         Input = recode(BaseFeature, !!!translate))|>
  group_by(Input, Component) |>
  summarise(
    Total_Gain = sum(Gain)*100,
    Total_Cover = sum(Cover)*100,
    Total_Frequency = sum(Frequency)*100
  ) |>
  left_join(importance_rf, by = 'Input') |>
  mutate(Input = factor(Input, levels = ordered_input_levels),
         Component = factor(Component, levels = c("Y","N","C","S","P")))

# Define color palette components.
fixed_colors <- c(
  "S" = "#e45a3d",
  "N" = "#377eb8",
  "C" = "#3AA600",
  "Y" = "#FFAE00",
  "P" = "#985aa1"
)

# Visualization: Point Plot
pp_plot <-
  contributions |> pivot_longer(
    cols = c(
      #xgb
      Total_Gain,
      #Total_Cover,
      Total_Frequency,
      #rf
      ##IncMSE_rel,
      #IncNodePurity_rel,
      #ranger,
      # IncNodePurity_rel
    ),
    names_to = "metric",
    values_to = "contribution"
  ) |>
  arrange(Component)|>
  mutate(metric = factor(metric, levels = c("Total_Gain", "Total_Frequency"))) |>
  ggplot(aes(
    x = contribution ,
    y = fct_reorder(Input, desc(Component)),
    color = Component,
    fill= metric
  )) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_manual(values = fixed_colors) +
  labs(x = "Relative contribution", y = "Input", color = "Component") +
  geom_segment(aes(
    x = 0,
    xend = contribution,
    y = Input,
    yend = Input
  )) +
  facet_grid(. ~ metric) +
  scale_y_discrete(
    labels = function(x)
      gsub("", "", x)
  )

pp_plot

seggain <-
  contributions |>
  ggplot(aes(
    x = Total_Gain ,
    y = fct_reorder(Input, desc(Component)),
    color = Component
  )) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_manual(values = fixed_colors) +
  labs(x = "Relative contribution (%)", y = "Input", color = "Component") +
  geom_segment(aes(
    x = 0,
    xend = Total_Gain,
    y = Input,
    yend = Input
  )) +
  scale_y_discrete(
    labels = function(x)
      gsub("", "", x)
  )

segfreq <-
  contributions |>
  ggplot(aes(
    x = Total_Frequency ,
    y = fct_reorder(Input, desc(Component)),
    color = Component
  )) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_manual(values = fixed_colors) +
  labs(x = "Feauture frequency (%)", y = "Input", color = "Component") +
  geom_segment(aes(
    x = 0,
    xend = Total_Frequency,
    y = Input,
    yend = Input
  )) +
  scale_y_discrete(
    labels = function(x)
      gsub("", "", x)
  )


# Visualization: Stacked Bar Plot
bar_plot <-
  contributions |> pivot_longer(
    cols = c(
      #xgb
      Total_Gain,
      #Total_Cover,
      Total_Frequency,
      #rf
      IncMSE_rel,
      IncNodePurity_rel,
      #ranger,
      # IncNodePurity_rel
    ),
    names_to = "metric",
    values_to = "contribution"
  ) |>
  #mutate(metric = factor(metric, levels = c("Total_Gain", "Total_Frequency"))) |>
  ggplot(aes(y = contribution,
             x = 1,
             group = fct_reorder(Input, -desc(Component)))) +
  geom_bar(
    stat = "identity",
    position = "stack",
    aes(fill = Component),
    color = "#FFF",
    linewidth = 0.1
  ) +
  geom_text(
    aes(label = ifelse(contribution > 1.5, paste(gsub("", "", Input), ";", round(contribution, 1), " "), paste(""))
    ),
    position = position_stack(vjust = 0.5),
    hjust = 0.5,
    angle = 0,
    color = "white",
    size = 3
  ) +
  theme_light() +
  scale_fill_manual(values = fixed_colors) +
  labs(y = "Relative contribution (%)", x = "", fill = "Component") +
  scale_x_continuous(breaks = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  facet_grid(.~metric)

bar_plot

# second option for report


ptotgain <- contributions |>
  ggplot(aes(y = Total_Gain,
             x = 1,
             group = fct_reorder(Input, -desc(Component)))) +
  geom_bar(
    stat = "identity",
    position = "stack",
    aes(fill = Component),
    color = "#FFF",
    linewidth = 0.1
  ) +
  geom_text(
    aes(label = ifelse(Total_Gain > 2, paste(gsub("", "", Input), ";", round(Total_Gain, 1), ""), paste(" "))
    ),
    position = position_stack(vjust = 0.5),
    hjust = 0.5,
    angle = 0,
    color = "white",
    size = 3
  ) +
  theme_minimal() +
  scale_fill_manual(values = fixed_colors) +
  labs(y = " ",
       x = "", fill = "Component") +
  scale_x_continuous(breaks = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) # Total_Gain

pfrequency <- contributions |>
  ggplot(aes(y = Total_Frequency,
             x = 1,
             group = fct_reorder(Input, -desc(Component)))) +
  geom_bar(
    stat = "identity",
    position = "stack",
    aes(fill = Component),
    color = "#FFF",
    linewidth = 0.1,
    alpha = 3
  ) +
  geom_text(
    aes(label = ifelse(Total_Frequency > 1.5, paste(gsub("", "", Input), ";", round(Total_Frequency, 1), " "), paste(""))
    ),
    position = position_stack(vjust = 0.5),
    hjust = 0.5,
    angle = 0,
    color = "white",
    size = 3
  ) +
  theme_minimal() +
  scale_fill_manual(values = fixed_colors) +
  labs(y = "Feauture frequency (%)", x = "", fill = "Component") +
  scale_x_continuous(breaks = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) #total frequency

paddedgain <- contributions |>
  select(Component, Total_Gain) |>
  group_by(Component) |>
  summarise(contribution = sum(Total_Gain), .groups = "drop") |>
  mutate(
    metric = "Total_Gain_Component", # Label for the summary bar
    Input = Component,# Use Component as label
    percentage = contribution / sum(contribution) * 100
  ) |>
  ggplot(aes(y = percentage,
             x = 1,
             group = fct_reorder(Input, -desc(Component)))) +
  geom_bar(
    stat = "identity",
    position = "stack",
    aes(fill = Component),
    color = "#FFF",
    linewidth = 0.1
  ) +
  geom_text(
    aes(label = ifelse(percentage > 1.5, paste(gsub("", "", Input), ";", round(percentage, 1), " "), paste(""))
    ),
    position = position_stack(vjust = 0.5),
    hjust = 0.5,
    angle = 0,
    color = "white",
    size = 3
  ) +
  theme_minimal() +
  scale_fill_manual(values = fixed_colors) +
  labs(y = "Relative contribution (%)", x = "", fill = "Component") +
  scale_x_continuous(breaks = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) #total contr comp

## Figure 6 -----
# Combine with custom widths
gain_plot <-
  paddedgain +
  theme(
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  ptotgain +
  theme(axis.ticks.y = element_blank(),
        #axis.text.y = element_blank(),
        legend.position = "bottom") +
  seggain +
  theme(axis.ticks.y = element_blank(), legend.position = "none") +
  plot_layout(
    ncol = 4,
    widths = c(0.8,2.1,2)
  )+plot_annotation(tag_levels = 'a')

gain_plot

ggsave("gain_plot.png",
       gain_plot,
       units = "mm",
       width = 183,
       height = 210)

freq_plot <-
  segfreq +
  theme(axis.ticks.y = element_blank(),
        legend.position = "none") +
  pfrequency + theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )+
  plot_layout(
    ncol = 2,
    widths = c(1,1.5),
    axis_titles = 'collect_y'
  )

ggsave("freq_plot.png",
       freq_plot,
       units = "mm",
       width = 183,
       height = 200)


# Ensure necessary packages are loaded

# --- Assuming 'contributions' data frame and 'fixed_colors' are defined ---
# (from your previous code snippet)

# 1. Prepare the data for the dumbbell plot
# We need the data in a 'long' format for the points, and a 'wide' format for the segments
### Figure 7 -----
dumbbell_data_long <- contributions |>
  pivot_longer(
    cols = c(Total_Gain, Total_Frequency), # Select the two metrics you want to compare
    names_to = "metric",
    values_to = "contribution"
  ) |>
  arrange(Component) |>
  mutate(metric = factor(metric, levels = c("Total_Gain", "Total_Frequency"))) # Ensure order for legend

# To draw the connecting segments, we need the start (Gain) and end (Frequency) points
# on the same row for each Input variable.
dumbbell_data_wide <- dumbbell_data_long |>
  pivot_wider(
    names_from = metric,
    values_from = contribution
  ) |>
  select(Input, Component, Total_Gain, Total_Frequency) # Select relevant columns

# 2. Create the Dumbbell Plot
dumbbell_plot <-
  ggplot(dumbbell_data_long, aes(y = reorder(Input, desc(Input)), x = contribution)) +
  # First, add the connecting segments
  geom_segment(data = dumbbell_data_wide,
               aes(x = Total_Gain, xend = Total_Frequency, y = reorder(Input, desc(Input)), yend = reorder(Input, desc(Input))),
               color = "grey60", linewidth = 0.8) + # Grey segments to connect the points
  # Then, add the points for each metric (Gain and Frequency)
  geom_point(aes(color = Component, shape = metric), size = 4) + # Color by Component, Shape by Metric
  scale_color_manual(values = fixed_colors) + # Use your predefined component colors
  scale_shape_manual(values = c("Total_Gain" = 16, "Total_Frequency" = 17), # Use common point shapes (solid circle, solid triangle)
                     labels = c("Total_Gain" = "Gain", "Total_Frequency" = "Frequency")) + # Labels for shape legend
  labs(
    #title = "Relative contribution of input variables (Gain vs. Frequency)",
    x = "Relative contribution (%)",
    y = "Input",
    color = "Component", # Legend title for component colors
    shape = "Metric" # Legend title for metric shapes
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5), # Center and bold plot title
    #axis.title = element_text(face = "bold"), # Bold axis titles
    legend.position = "right", # Place legend at the bottom
    legend.box = "right", # Arrange legend items horizontally
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines for cleaner look
    panel.grid.minor = element_blank(), # Remove minor grid lines
    axis.text.y = element_text(angle = 0, hjust = 1) # Ensure Y-axis labels are readable
  )


# Print the plot
print(dumbbell_plot)

# You can save the plot if needed
ggsave("dumbbell_plot_contribution.png",
       plot = dumbbell_plot,
       width = 183,
       units= "mm",
       height = 150,
       dpi = 300)


# 📈 1. Gain (Most Commonly Used)
# Definition: The average improvement in accuracy brought by a feature to the branches it is used in.
# Interpretation: A higher gain means the feature is more important in reducing error.
# Use case: Best for understanding which features contribute most to model performance.
# Think of it as: "How much does this feature help the model make better decisions?"
#
# 🌿 2. Cover
# Definition: The average number of samples affected by splits on this feature.
# Interpretation: A higher cover means the feature is used in splits that affect more data points.
# Use case: Useful for understanding how broadly a feature is used across the dataset.
# Think of it as: "How many data points does this feature influence?"
#
# 🔁 3. Frequency (a.k.a. Weight)
# Definition: The number of times a feature is used in all trees.
# Interpretation: A higher frequency means the feature is used more often in the model.
# Use case: Good for identifying features that are consistently used, even if their individual impact is small.
# Think of it as: "How often does the model rely on this feature?"

# 🧠 Now, comparing with classical sensitivity indices (like in global sensitivity analysis):
#   Metric	XGBoost (Gain, Cover, Frequency)	Classical Sensitivity Indices (e.g., Sobol, FAST)
# Goal	Model-specific feature importance	Model-agnostic input-output sensitivity
# Scope	Local to model structure (tree splits)	Global across input space
# Gain	Measures improvement in loss	Similar to first-order Sobol index (main effect)
# Cover	Measures data coverage of splits	No direct analog, but loosely related to input density
# Frequency	Measures how often a feature is used	No direct analog, but reflects model reliance
# Interpretability	Easy to compute and visualize	More rigorous but computationally expensive

## 4.4. Sensitivity marginal behavior -----

library(pdp)

# Train a sample model (if needed)
#xgb_model <- xgboost(data = as.matrix(X), label = y, nrounds = 100, objective = "reg:squarederror", verbose = 0)

# Create a partial dependence plot for a feature
# pdp_result_xgb <- pdp::partial(
#   object = xgb_surrogate_model,
#   pred.var = 'Indb_aar',
#   train = as.data.frame(ohe_data |> dplyr::select(xgb_surrogate_model$feature_names)),
#   grid.resolution = 100
# )
#
# pdp_result_rf <- pdp::partial(
#   object = rf_surrogate_model,
#   pred.var = 'Indb_aar',
#   train = as.data.frame(predus_fil |> dplyr::select(row.names(rf_surrogate_model$importance))),
#   grid.resolution = 100
# )
#
# # Plot with ggplot2
# ggplot(pdp_result_xgb, aes(x = Indb_aar, y = yhat)) +
#   geom_point(color = "steelblue", size = 1) +
#   geom_line(color = "steelblue", size = .8) +
#   geom_smooth(method = "loess", se = TRUE, color = "steelblue", size = 0.5) +
#   labs(
#     title = "Partial Dependence Plot",
#     x = "Indb_aar",
#     y = "Predicted Response"
#   ) +
#   theme_minimal()
#
# ggplot(pdp_result_rf, aes(x = Indb_aar, y = yhat)) +
#   geom_point(color = "steelblue", size = 1) +
#   geom_line(color = "steelblue", size = .8) +
#   geom_smooth(method = "loess", se = TRUE, color = "steelblue", size = 0.5) +
#   labs(
#     title = "Partial Dependence Plot",
#     x = "Indb_aar",
#     y = "Predicted Response"
#   ) +
#  theme_minimal()

### Figure 8 ------
categorical_vars <- c("Mau", "Mfu", "Vau", "Vfu", "SoilG_1", "SoilG_2","SoilG" ,"Vafgr_Kappa")

# Function to compute PDP for both models
compute_pdp <- function(feature) {

  pdp_result_xgb <- pdp::partial(
    object = xgb_surrogate_model,
    pred.var = feature,
    train = as.data.frame(ohe_data |> dplyr::select(xgb_surrogate_model$feature_names)),
    grid.resolution = 20
  ) |> mutate(variable = feature, model = "XGBoost") |>
    rename(covariable = feature) |>
    select(variable, model, yhat,covariable)

  pdp_result_rf <- pdp::partial(
    object = rf_surrogate_model,
    pred.var = feature,
    train = as.data.frame(predus_fil |> dplyr::select(row.names(rf_surrogate_model$importance))),
    grid.resolution = 20
  )|> mutate(variable = feature, model = "rf")|>
    rename(covariable = feature) |>
    select(variable, model, yhat,covariable)

  rbind(pdp_result_xgb, pdp_result_rf)
}

# Compute PDPs for all features
pdp_cont <- lapply(setdiff(row.names(rf_surrogate_model$importance),
                           categorical_vars),
                   compute_pdp)
#saveRDS(pdp_cont, "pdp_cont.RDS")

#pdp_cont <- readRDS(pdp_cont.RDS)

# Plot all PDPs in a faceted plot
pdpplot <- do.call(rbind,pdp_cont) |>
  filter(model=="XGBoost")|>
  mutate(Component=recode(variable, !!!asignation),
         Input=recode(variable,!!!translate)) |>
  mutate(Component = factor(Component, levels = c("Y","N","C","S","P")),
         Input = factor(Input, levels = ordered_input_levels)) |>
  filter(Input!= "G2")|>
  ggplot(aes(x = covariable, y = yhat, color=Component#, linetype = model
             )) +
  geom_line(size = 0.5) +
  geom_point(size = .5) +
  geom_smooth(method = "loess", se = TRUE, size = 0.8) +
  scale_color_manual(values = fixed_colors)+
  facet_wrap(~ reorder(Input, Component), ncol = 4, scales = "free_x")+
  #facet_grid(variable~model , scales = "free") +
  labs(
    x = "Feature Value",
    y = "Predicted Response"
  ) +
  theme_light() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10,),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9))

pdpplot

ggsave("pdp.png",
       plot = pdpplot,
       width = 183,
       units= "mm",
       height = 220,
       dpi = 300)

## --- 4.5. Partial Dependence Plots for Categorical Variables ----


feature <- "d2"  # or any other feature name

# Create the sequence of values for the feature
feature_vals <- seq(
  min(ohe_data[[feature]], na.rm = TRUE),
  max(ohe_data[[feature]], na.rm = TRUE),
  length.out = 50
)

# Create the prediction grid with correct column names
pred_grid <- expand.grid(
  SoilG_1 = factor(c(0, 1)),
  dummy = feature_vals
)

# Rename 'dummy' to the actual feature name
names(pred_grid)[names(pred_grid) == "dummy"] <- feature

# Now call pdp::partial
pdp_cat <- pdp::partial(
  object = xgb_surrogate_model,
  pred.var = c(feature, "SoilG_1"),
  train = as.data.frame(ohe_data |> dplyr::select(xgb_surrogate_model$feature_names)),
  pred.grid = pred_grid
)

pdp_cat |>
  ggplot(aes(x = d2, y = yhat, color = as.factor(SoilG_1))) +
  geom_point(size = 1) +
  geom_smooth(method = "loess", se = TRUE, size = 0.5) +
  geom_line() +
  labs(
    title = paste("Partial Dependence Plot for", feature),
    x = feature,
    y = "Predicted Response"
  ) +
  theme_minimal()


pdp3 <- lapply(c("d1", "d2", "d3", "Indaar"),
               pdp_point3)

# # --- 5. Shapley Values for Feature Importance only xgb---
# # Load specific packages for SHAP values
# library(shapr)      # For calculating SHAP values
# library(shapviz)    # For visualizing SHAP values
#
# # Calculates SHAP values for the model using the preprocessed training data.
# shap_values_xgb <-
#   shapviz(xgb_surrogate_model,
#           X_pred =data.matrix(ohe_data |>
#                                 select(xgb_surrogate_model$feature_names)))
#
#
# # Visualize SHAP values
# sv_importance(shap_values_xgb, show_numbers = TRUE)
# sv_importance(shap_values_xgb, kind = "beeswarm")
#
# sv_dependence(shap_values_xgb, v = setdiff(row.names(rf_surrogate_model$importance),
#                                            categorical_vars), color_var = "SoilG")
#
#
# sv_dependence(shap_values_xgb, v = "d1", color_var = "SoilG_1")
#
#
# sv_interaction(shap_values_xgb, max_display = 20, show_numbers = TRUE)
#
#

cv_results_robust <- sasoutput_fil_NLES5 |>
  summarise(across(everything(),
                   # Using abs() on the mean for a more stable calculation
                   ~ (sd(.x, na.rm = TRUE) / abs(mean(.x, na.rm = TRUE))) * 100
  ))

print(cv_results_robust)
