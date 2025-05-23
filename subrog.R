# Load necessary packages
library(tidyverse)
library(readxl)
library(randomForest)
library(xgboost)
library(caret)

# 1. Data Loading ----
sasoutput <- read_excel("Scenarier20190909B4_found0325.xls")
table(sasoutput$SoilG)
# 1.1. Renaming and grouping ----
# Is d or p the percolation groups why some are negative?
hist(sasoutput$d2)
hist(sasoutput$p2)
plot(sasoutput$d2, sasoutput$p2)

# 1.3 Select NLES5 inputs (i.e. variables) for analysis ----

sasoutput_fil <- sasoutput |>
  select(
    PUdvaskF,    # Predicted values (response for RF/XGBoost)
    NS, `NA`, nlevelMin1, nlevelMin2, NlevelFix0, NlevelFix1, NlevelFix2,
    NlevelGod0, NlevelGod1, NlevelGod2, Nudb, TN, # Nitrogen related terms
    Mau, Mfu, Vau, Vfu, # Crop related terms
    CU,             # Clay content (Soil)
    Indb_aar,     # Year trend
    d1,d2, d3,
    SoilG # Percolation/Drainage related terms
  ) |> mutate(SoilG= recode_factor(SoilG, "S"=0, "C"=1)) |>
  rename('na' = `NA`) |>
  mutate(Mau=as.factor(Mau),
         Mfu=as.factor(Mfu),
         Vau=as.factor(Vau),
         Vfu=as.factor(Vfu),
         SoilG= as.factor(SoilG))

# --- 2. Random Forest Surrogate Model (Optimized for High Fidelity i.e, over fitting) ----

set.seed(151190) # For reproducibility

rf_surrogate_model <- randomForest(
  PUdvaskF ~ .,
  data = sasoutput_fil,
  ntree = 5000,      # Increase number of trees further (e.g., to 5000)
  mtry = ncol(sasoutput_fil) - 1, # Set mtry to all predictors (p)
  nodesize = 1,      # Allow single-observation terminal nodes
  maxnodes = NULL,   # Let nodesize control tree depth (no explicit max)
  sampsize = nrow(sasoutput_fil), # Use entire dataset for each tree (pure bagging mode)
  importance = TRUE  # Still useful for initial importance checks
)

print(rf_surrogate_model)

cat("R-squared of High-Fidelity Random Forest on training data:",
    rf_surrogate_model$rsq[length(rf_surrogate_model$rsq)], "\n")

# Evaluate performance on the training data (should be very high)
rf_predictions_train <- predict(rf_surrogate_model, newdata = sasoutput_fil)
rf_rmse_train <- sqrt(mean((sasoutput_fil$PUdvaskF - rf_predictions_train)^2))
cat("RMSE of High-Fidelity Random Forest on training data:", rf_rmse_train, "\n")

# You might want to save this model if it's the one you'll use for Sobol
saveRDS(rf_surrogate_model, "rf_surrogate_model_max_fidelity.rds")


# --- 3. XGBoost Surrogate Model (Optimized for High Fidelity) ----

## --- Create a dummy variables for categories ---
# Create a dummyVars object

dummy_model <- dummyVars(~ ., sep = '_', data = sasoutput_fil)

# Apply the dummy encoding
ohe_data <- data.frame(predict(dummy_model, newdata = sasoutput_fil))

set.seed(151190) # For reproducibility

# Prepare data for XGBoost (matrix format)
X_train <- as.matrix(ohe_data |> select(-PUdvaskF))
y_train <- ohe_data$PUdvaskF

# Define parameters for high fidelity
xgb_params_high_fidelity <- list(
  objective = "reg:squarederror", # Regression task
  eval_metric = "rmse",         # Root Mean Squared Error
  nrounds = 5000,               # Very high number of boosting rounds
  eta = 0.005,                   # Very small learning rate
  max_depth = 5,               # Deep trees to capture interactions
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
cat("XGBoost R-squared (Training):", cor(xgb_predictions_train, ohe_data$PUdvaskF), "\n")

# You might want to save this model
#saveRDS(xgb_surrogate_model, "xgb_surrogate_model_high_fidelity.rds")
#xgb.save(xgb_surrogate_model, "xgb_surrogate_model_high_fidelity.model")
#xgb_surrogate_model <- readRDS("xgb_surrogate_model_high_fidelity.rds")

# --- Comparison ---
cat("\n--- Training Performance Comparison ---\n")
cat("RF R-squared (Training):", rf_surrogate_model$rsq[length(rf_surrogate_model$rsq)], "\n")
cat("XGBoost R-squared (Training):", cor(xgb_predictions_train, ohe_data$PUdvaskF), "\n")
cat("RF RMSE (Training):", rf_rmse_train, "\n")
cat("XGBoost RMSE (Training):", xgb_rmse_train, "\n")

plot(xgb_predictions_train, ohe_data$PUdvaskF)
abline(0,1, col="orange3")
plot(rf_surrogate_model$predicted, ohe_data$PUdvaskF)
abline(0,1)


# 4. Results and visualization ----
# Variable Importance
# --- 4.1 Extract Feature Importance ----
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
                   "Nudb", "TN", "CU")


importance_matrix_processed <- importance_matrix |>
  mutate(BaseFeature = str_extract(Feature, paste0(
    "^(",
    paste(prefixes, collapse = "|"),
    "|[A-Za-z]+)(?=\\d|_|$)"
  )))


# 4.2. Renaming and grouping ----
# Define the mapping for the crop prefixes
translate <- c(
  #𝐿 = 𝜏(𝑌 – 1991) + {(𝜇 + 𝜃𝑖𝑁 + 𝐶)𝜅}(𝑃 𝑆)p

  # Y
  "Indb_aar"="Year",

  # N
  "NS"="N mineral Spring",
  "na"="N mineral Autum",
  "nlevelMin1"="N mineral prev. year",
  "nlevelMin2"="N mineral 2nd prev. year",
  "NlevelFix0"="N Fixation",
  "NlevelFix1"="N Fixation prev. year",
  "NlevelFix2"="N Fixation 2nd prev. year",
  "NlevelGod0"="N Organic",
  "NlevelGod1"="N Organic prev. year",
  "NlevelGod2"="N Fixation 2nd prev. year",
  "Nudb"="N grazing",
  "TN"="Total Soil N",

  # C
  "Mau" = "Main crop",
  "Mfu" = "Main crop prev. year",
  "Vau" = "Winter crop",
  "Vfu" = "Winter crop prev. year",

  #P
  "d1"="Percolation April to August",
  "d2"="Percolation September to March",
  "d3"="Percolation September to March prev. year",
  "SoilGS"="Soil texture",
  "SoilGC"="Soil texture",

  #S
  "CU"="Clay")

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

  # C
  "Mau" = "C",
  "Mfu" = "C",
  "Vau" = "C",
  "Vfu" = "C",

  #P
  "d1"="P",
  "d2"="P",
  "d3"="P",
  "SoilGS"="P",
  "SoilGC"="P",

  #S
  "CU"="S"
)

# 4.3 Overal results Input contribution -----

# Group by the identified crop prefixes and sum their contributions
contributions <- importance_matrix_processed |>
  mutate(Component = recode(BaseFeature, !!!asignation),
           Input = recode(BaseFeature, !!!translate))|>
    group_by(Input, Component) |>
    summarise(
      Total_Gain = sum(Gain),
      Total_Cover = sum(Cover),
      Total_Frequency = sum(Frequency)
    ) |>
    arrange(desc(Total_Frequency))



# Define color palette components.
fixed_colors <- c(
  "S" = "#e45a3d",
  "N" = "#377eb8",
  "C" = "#3AA600",
  "Y" = "#FFAE00",
  "P" = "#985aa1"
)

# Visualization: Point Plot
pp_plot <- contributions |> pivot_longer(
  cols = c(Total_Gain, Total_Cover, Total_Frequency),
  names_to = "metric",
  values_to = "contribution"
) |>
  ggplot(aes(
    x = contribution ,
    y = fct_rev(Input),
    color = Component
  )) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_manual(values = fixed_colors) +
  labs(x = "Relative contribution", y = "Input", color = "Component") +
  geom_segment(aes(
    x = 0,
    xend = contribution,
    y = fct_rev(Input),
    yend = fct_rev(Input)
  )) +
  facet_wrap( ~ metric) +
  scale_y_discrete(
    labels = function(x)
      gsub("", "", x)
  )

pp_plot

# Visualization: Stacked Bar Plot
bar_plot <-
  contributions |> pivot_longer(
    cols = c(Total_Gain, Total_Cover, Total_Frequency),
    names_to = "metric",
    values_to = "contribution"
  ) |> ggplot(aes(y = contribution,
                            x = 1, group = Input)) +
  geom_bar(
    stat = "identity",
    position = "stack",
    aes(fill = Component),
    color = "#FFF",
    linewidth = 0.1
  ) +
  geom_text(
    aes(label = ifelse(
      metric >= 1, paste(gsub(" ", "", Input), ";", round(contribution, 2), " "), ""
    )),
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
  facet_wrap(~metric)

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

# 4.4. Sensitivity marginal behavior -----

library(pdp)

# Assume you have a trained XGBoost model and training data
# X: feature matrix, y: target vector

# Train a sample model (if needed)
xgb_model <- xgboost(data = as.matrix(X), label = y, nrounds = 100, objective = "reg:squarederror", verbose = 0)

# Create a partial dependence plot for a feature
pdp_result <- partial(
  object = xgb_surrogate_model,
  pred.var = "Clay",
  train = as.data.frame(ohe_data),
  grid.resolution = 50
)

 # Plot with ggplot2
ggplot(pdp_result, aes(x = feature1, y = yhat)) +
  geom_line(color = "steelblue", size = 1.2) +
  labs(
    title = "Partial Dependence Plot",
    x = "Feature 1",
    y = "Predicted Response"
  ) +
  theme_minimal()

