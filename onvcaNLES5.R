library(tidyverse)
library(ggthemes)
library(ggpubr)
library(patchwork)
library(readxl)
library(lme4) # For mixed-effects models
library(VCA)  # For variance component analysis
library(randomForest)

# 1. Data Loading ----
# from output sas fitting nles5 model
sasoutput <- read_excel("Scenarier20190909B4_found0325.xls")

# Inspect the some of the terms to understand its distribution.
table(sasoutput$SoilG)

colnames(sasoutput)

#This analysis aimed to understand the variance contribution of input variables
#in the NLES5 prediction of annual leaching values.
#We prepossessed by categorizing continuous numerical variables into
#discrete groups using a probabilistic approach based on quantiles generating tree levels in each variable.
#A linear mixed-effects model was then employed to assess the relative influence
# of the four components of NLES5: #
# nitrogen-related inputs,
#crop types,
#soil characteristics (clay content),
#temporal trends (year),
#and percolation/drainage on the predicted leaching.
# Variance component analysis (VCA) was conducted to quantify
#the proportion of variance attributed to each factor, with the results visualized through point and stacked bar plots.
#These visualizations provided a  overview of the relative importance of each factor.
# The analysis, while informative, has a important limitation because
#it doesn't account for correlations or interactions between the input variables.
# We might over- or underestimate the importance of certain variables
# that we know are strongly correlated or interact with others.
# As a second alternative we have build a "surrogated model" some people call this a "meta model"
# i.e., a model of a model where the response variable is the annual leaching predicted by NLES5
# and the explanatory variables the inputs. For this we performed a random forest because
# while it is not a sensitivity analysis in the strict sense,
# the algorithm can identify how correlations in inputs can affect the response
# better than the classical VCA.
# The differences in both methods can be talking about
# the complex level of interaction in the inputs,
# and the presence of non-linear patterns,
# the Random Forest surrogate model might provide a more accurate representation of variable importance.
# But for understanding the nature of predictor-outcome relationships
# and have data that fits linear assumptions, the VCA offers better interpretability
# Again, a complete sensitivity analysis would require a more accurate approach.


# 2. Uncertainity of predicted and observed Nleaching ----
scatter <- ggplot(sasoutput, aes(x = PUdvaskF, y = Udvask)) +
  geom_point(alpha=0.3, size = 2) + # Add a color aesthetic for the points
  geom_abline(slope = 1, intercept = 0, color = "coral") +
  geom_rug() +
  labs(x = "Predicted Annual Leaching (PUdvaskF)",
       y = "Observed Annual Leaching (Udvask)") + # Set the legend title
  theme_pubr() +
  theme(legend.position = c(0.15, 0.9))

ggplot(sasoutput, aes(x = PUdvaskF, fill = "Marginal Density")) +
  geom_density(alpha = 0.4) +
  theme_classic2() +
  theme(legend.position = "none") +
  labs(x = "Predicted Annual Leaching (PUdvaskF)",
       y = "Density") +
  scale_fill_manual(values = c("Marginal Density" = "coral")) # Set a fill color

dens_x <- ggplot(sasoutput, aes(x = PUdvaskF, fill = "Marginal Density")) +
  geom_density(alpha = 0.4) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("Marginal Density" = "coral")) # Set a fill color

dens_y <- ggplot(sasoutput, aes(x = Udvask, fill = "Marginal Density")) +
  geom_density(alpha = 0.4) +
  theme_void() +
  theme(legend.position = "none") +
  coord_flip() +
  scale_fill_manual(values = c("Marginal Density" = "lightcoral")) # Set a fill color

dens_x  + plot_spacer() + scatter + dens_y +
  plot_layout(widths = c(3, 1), heights = c(1, 3))

# Calculate summary statistics of the predicted leaching
cat("\nSummary Statistics of Predicted Leaching:\n")
print(summary(sasoutput$PUdvaskF))
cat("\nQuantiles:\n")
print(quantile(sasoutput$PUdvaskF, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))

# 3. Sensitivity ----
# 3.1 Select relevant columns(i.e. variables) for analysis----
sasoutput_fil <- sasoutput |>
  select(
    #Udvask,     # Observed values
    PUdvaskF,   # Predicted values
    NS, `NA`, nlevelMin1, nlevelMin2, NlevelFix0, NlevelFix1, NlevelFix2,
    NlevelGod0, NlevelGod1, NlevelGod2, Nudb, TN, # Nitrogen related terms
    Mau, Mfu, Vau, Vfu, # Crop related terms
    CU,          # Clay content (Soil)
    Indb_aar,    # Year trend
    Drain, SoilG # Percolation/Drainage related terms
  ) |> rename('na' = `NA`)

# 3.2. Categorization Function ----
categorize_numeric_probabilistic <- function(df) {
  #' Categorizes numeric columns using quantiles, handling low-variance columns.
  #'
  #' @param df A data frame.
  #' @return A data frame with categorized numeric columns.
  numeric_cols <- df |> select(where(is.numeric)) |> names()

  for (col in numeric_cols) {
    quantiles <- quantile(df[[col]], probs = c(0, 0.33, 0.66, 1), na.rm = TRUE)

    # Handle low-variance columns by creating equidistant categories.
    if (length(unique(quantiles)) == 1) {
      min_val <- min(df[[col]], na.rm = TRUE)
      max_val <- max(df[[col]], na.rm = TRUE)
      breaks <- seq(min_val, max_val, length.out = 4)
      df <- df |>
        mutate(
          !!paste0(col, "_cat") := cut(
            .data[[col]],
            breaks = breaks,
            labels = paste0(col, "_Eq", 1:3),
            include.lowest = TRUE
          )
        )
    } else {
      # Categorize using quantiles.
      df <- df |>
        mutate(
          !!paste0(col, "_cat") := case_when(
            .data[[col]] >= quantiles[1] &
              .data[[col]] < quantiles[2] ~ paste0(col, "_Q1"),
            .data[[col]] >= quantiles[2] &
              .data[[col]] < quantiles[3] ~ paste0(col, "_Q2"),
            .data[[col]] >= quantiles[3] &
              .data[[col]] <= quantiles[4] ~ paste0(col, "_Q3"),
            TRUE ~ as.character(.data[[col]]) # Handle NAs
          )
        )
    }

    # Handle single-value categories by creating 0/non-0 categories.
    if (length(unique(df[[paste0(col, "_cat")]][!is.na(df[[paste0(col, "_cat")]])])) == 1) {
      df <- df |>
        mutate(!!paste0(col, "_cat") := case_when(.data[[col]] == 0 ~ "0", TRUE ~ "Non-0"))
    }
  }
  return(df)
}

# Apply categorization function to the selected data.
categorized_sasoutput <- categorize_numeric_probabilistic(sasoutput_fil)

# Select categorized columns for the model.
sasoutput_fil_cat <- categorized_sasoutput |>
  select(
    PUdvaskF,
    NS_cat,
    na_cat,
    nlevelMin1_cat,
    nlevelMin2_cat,
    NlevelFix0_cat,
    NlevelFix1_cat,
    NlevelFix2_cat,
    NlevelGod0_cat,
    NlevelGod1_cat,
    NlevelGod2_cat,
    Nudb_cat,
    TN_cat,
    # Categorized Nitrogen terms
    Mau,
    Mfu,
    Vau,
    Vfu,
    # Crop terms
    CU_cat,
    # Categorized Clay
    Indb_aar_cat,
    # Categorized Year trend
    Drain_cat,
    SoilG # Categorized Percolation/Drainage
  )

# 3.3. VCA Mixed-Effects Model ----
# Fit a linear mixed-effects model to assess the variance contribution of different factors.
res_VCA <- lme4::lmer(
  PUdvaskF ~ 1 +
    (1 | TN_cat) + (1 | NS_cat) + (1 | na_cat) +
    (1 | nlevelMin1_cat) + (1 | nlevelMin2_cat) +
    (1 | NlevelFix0_cat) + (1 | NlevelFix1_cat) + (1 | NlevelFix2_cat) +
    (1 | NlevelGod0_cat) + (1 | NlevelGod1_cat) + (1 | NlevelGod2_cat) +
    (1 | Nudb_cat) +
    (1 | Mau) + (1 | Mfu) + (1 | Vau) + (1 | Vfu) +
    (1 | CU_cat) +
    (1 | Indb_aar_cat) +
    (1 | Drain_cat) + (1 | SoilG),
  na.action = na.omit,
  REML = TRUE,
  data = sasoutput_fil_cat
)

# Summarize the model results.
summary(res_VCA)

# Variance Component Analysis
# Extract variance components from the model.
vca <- as.data.frame(lme4::VarCorr(res_VCA))

# Group input variables into expert-defined categories.
expert_groups <- c(
  "NS_cat" = "Ninp",
  "na_cat" = "Ninp",
  "nlevelMin1_cat" = "Ninp",
  "nlevelMin2_cat" = "Ninp",
  "NlevelFix0_cat" = "Ninp",
  "NlevelFix1_cat" = "Ninp",
  "NlevelFix2_cat" = "Ninp",
  "NlevelGod0_cat" = "Ninp",
  "NlevelGod1_cat" = "Ninp",
  "NlevelGod2_cat" = "Ninp",
  "Nudb_cat" = "Ninp",
  "TN_cat" = "Ninp",
  "Mau" = "Crop",
  "Mfu" = "Crop",
  "Vau" = "Crop",
  "Vfu" = "Crop",
  "CU_cat" = "Soil",
  "Indb_aar_cat" = "Trend",
  "Drain_cat" = "Percolation",
  "SoilG" = "Percolation",
  "Residual" = "Residual"
)

# Calculate variance proportions and organize the data for plotting.
vca_rel <- vca |>
  group_by(grp) |>
  summarise(VCA = vcov / sum(vca$vcov) * 100) |>
  arrange(VCA, grp) |>
  mutate(component = expert_groups[grp]) |>
  group_by(component) |>
  mutate(ordered_grp = fct_reorder(grp, -VCA)) |>
  ungroup()

# 3.4. Random Forest Model ----
set.seed(151190) # For reproducibility

rf_model <- randomForest(
  PUdvaskF ~ .,
  data = sasoutput_fil,
  importance = TRUE,
  ntree = 500 # Adjust as needed
)

# 4. Results and visualization ----
# Variable Importance
importance_df <- as.data.frame(importance(rf_model))
importance_df$Variable <- rownames(importance_df)
importance_df <- importance_df |>
  group_by(Variable) |>
  summarise(Total_IncMSE = sum(`%IncMSE`)) |>
  mutate(RF = Total_IncMSE / sum(Total_IncMSE) * 100)

merged_data <- vca_rel  |>
  mutate(Variable = str_replace(grp, "_cat", "")) %>%
  left_join(importance_df, by = "Variable")


merged_long <-
  merged_data |> pivot_longer(
    cols = c(VCA, RF),
    names_to = "method",
    values_to = "rel_contribution"
  )

# Define color palette components.
fixed_colors <- c(
  "Soil" = "#e45a3d",
  "Ninp" = "#377eb8",
  "Crop" = "#3AA600",
  "Trend" = "#FFAE00",
  "Percolation" = "#985aa1"
)

# Visualization: Point Plot
pp_plot <- merged_long |>
  ggplot( aes(
  x = rel_contribution,
  y = fct_rev(ordered_grp),
  color = component
)) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_manual(values = fixed_colors) +
  labs(x = "Relative contribution (%)", y = "Input", color = "Component") +
  geom_segment(aes(
    x = 0,
    xend = rel_contribution,
    y = fct_rev(ordered_grp),
    yend = fct_rev(ordered_grp)
  ))+
  facet_wrap(~method)+
  scale_y_discrete(labels = function(x) gsub("_cat", "", x))

# Visualization: Stacked Bar Plot
bar_plot <-
merged_long |> ggplot(aes(y = rel_contribution,
                          x = 1, group = ordered_grp)) +
  geom_bar(
    stat = "identity",
    position = "stack",
    aes(fill = component),
    color = "#FFF",
    linewidth = 0.1
  ) +
  geom_text(
    aes(label = ifelse(
      rel_contribution >= 1, paste(gsub("_cat", "", grp), ";", round(rel_contribution, 2), "%"), ""
    )),
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
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  facet_wrap(~method)


# Combine Plots

# Your Point Plot
pp_plot_no_legend <- merged_long |>
  ggplot(aes(
    x = rel_contribution,
    y = fct_rev(ordered_grp),
    color = component
  )) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_manual(values = fixed_colors) +
  labs(x = "Relative contribution (%)", y = "Input") +
  geom_segment(aes(
    x = 0,
    xend = rel_contribution,
    y = fct_rev(ordered_grp),
    yend = fct_rev(ordered_grp)
  )) +
  facet_wrap(~method) +
  scale_y_discrete(labels = function(x) gsub("_cat", "", x)) +
  theme(legend.position = "none")

# Stacked Bar Plot
bar_plot_no_legend <-
  merged_long |> ggplot(aes(y = rel_contribution,
                            x = 1, group = ordered_grp)) +
  geom_bar(
    stat = "identity",
    position = "stack",
    aes(fill = component),
    color = "#FFF",
    linewidth = 0.1
  ) +
  geom_text(
    aes(label = ifelse(
      rel_contribution >= 1, paste(gsub("_cat", "", grp), ";", round(rel_contribution, 2), "%"), ""
    )),
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
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "none") + # Remove legend here as well
  facet_wrap(~method)

# Extract the legend from one of the plots (e.g., pp_plot)
common_legend <- get_legend(
  bar_plot_no_legend + theme(legend.position = "bottom",
                             legend.direction = "horizontal")
)

combined_plot <-
  (pp_plot_no_legend + bar_plot_no_legend) / common_legend +
  plot_layout(heights = c(1, .1), widths = c(1, 1)) # Adjust heights as needed

print(combined_plot)


combined_plot <-
  pp_plot_no_legend / bar_plot_no_legend / common_legend +
  plot_layout(heights = c(1, 1, 0.2)) # Adjust heights as needed

print(combined_plot)
