library(tidyverse)
# Output estimation from Nles5 Scenarier20190829B4.rtf (N-normer) Model S37. Kørt 20190920
#Fastholdte parametre i følgende modeller trend=-0.11080


parameter <- c(
  "Teta2", "bf0", "bm1", "bf1", "bg0", "deltas1", "deltas2", "nys2",
  "g2", "g3", "g4", "g5", "g6", "g8", "g9", "g10", "g11", "g12", "g13",
  "g14", "gv2", "gv3", "gv4", "gv5", "gv6", "gv7", "gv8", "l2", "l3",
  "l4", "l12", "lv2", "lv3", "lv4", "lv5", "lv6", "lv7", "lv8", "lv9",
  "lv10", "my", "b1", "b2", "b5", "b10", "eta3", "delta1", "delta2", "ny2"
)

estimate <- c(
  1.2051, 0.0163, 0.0265, 0.0137, 0.0141, 0.00119, 0.00111, 0.000856,
  -6.7438, -7.2794, -13.4935, -17.4783, -11.1924, -0.6400, 3.5344,
  -7.3187, -1.2484, 19.5236, -6.2293, -2.8664, -2.0553, -0.4557,
  -15.9592, -3.7918, -14.5962, -1.0486, -21.0597, 2.8466, 0.6643,
  1.1600, 4.0497, 9.7042, 10.6009, 9.3544, 13.2411, 5.4827, -1.5721,
  7.4126, 7.3958, 10.9751, 25.9582, 0.0291, 0.1570, 0.0382, 0.4568,
  0.00185, 0.000798, 0.000745, 0.000638
)

std_error <- c(
  0.1107, 0.00553, 0.00612, 0.00544, 0.00880, 0.000437, 0.000306, 0.000163,
  2.7253, 3.0893, 4.1828, 5.3884, 4.8207, 3.1956, 2.9728, 2.2949, 8.0491,
  9.7452, 9.1540, 3.2005, 1.1849, 1.4986, 2.6737, 1.7005, 2.5691, NA, 6.2083,
  1.0308, 2.0000, 9.8377, 1.5589, 2.8636, 3.4470, 2.9015, 5.1011, 3.0945,
  2.9633, NA, 7.9762, 9.3177, 4.3418, 0.00704, 0.0343, 0.0111, 0.2022,
  0.00456, 0.000233, 0.000180, 0.000144
)

conf_lower <- c(
  0.9881, 0.00547, 0.0145, 0.00308, -0.00316, 0.000337, 0.000508, 0.000537,
  -12.0885, -13.3379, -21.6966, -28.0457, -20.6464, -6.9071, -2.2956,
  -11.8193, -17.0338, 0.4119, -24.1815, -9.1431, -4.3791, -3.3947,
  -21.2027, -7.1267, -19.6345, NA, -33.2351, 0.8250, -3.2580, -18.1330,
  0.9926, 4.0884, 3.8409, 3.6642, 3.2372, -0.5860, -7.3835, NA, -8.2467,
  -7.2983, 17.4433, 0.0153, 0.0899, 0.0166, 0.0602, -0.00709, 0.000341,
  0.000392, 0.000356
)

conf_upper <- c(
  1.4222, 0.0272, 0.0385, 0.0244, 0.0314, 0.00205, 0.00171, 0.00117,
  -1.3991, -1.2208, -5.2904, -6.9110, -1.7384, 5.6271, 9.3644, -2.8182,
  14.5369, 38.6353, 11.7229, 3.4103, 0.2684, 2.4832, -10.7156, -0.4569,
  -9.5579, NA, -8.8844, 4.8682, 4.5865, 20.4531, 7.1068, 15.3201,
  17.3610, 15.0447, 23.2450, 11.5514, 4.2393, NA, 23.0382, 29.2485,
  34.4730, 0.0429, 0.2242, 0.0599, 0.8534, 0.0108, 0.00126, 0.00110,
  0.000920
)

# Create the data frame
org_param <- data.frame(
  Parameter = parameter,
  Estimate = estimate,
  Std.Error = std_error,
  Conf.Lower = conf_lower,
  Conf.Upper = conf_upper
) |> mutate(Estimation = "KK19")


asignation <- c(
  "Teta2"="Nitrogen",

  "bf0"="Nitrogen", "bm1"="Nitrogen", "bf1"="Nitrogen", "bg0"="Nitrogen",
  "b1"="Nitrogen", "b2"="Nitrogen", "b5"="Nitrogen", "b10"="Nitrogen",

  "deltas1"="Percolation-Soil", "deltas2"="Percolation-Soil", "nys2"="Percolation-Soil",
  "delta1"="Percolation-Soil", "delta2"="Percolation-Soil", "ny2"="Percolation-Soil",

  "g2"= "Crop", "g3"= "Crop", "g4"= "Crop", "g5"= "Crop", "g6"= "Crop",
  "g8"= "Crop", "g9"= "Crop", "g10"= "Crop", "g11"= "Crop", "g12"= "Crop", "g13"= "Crop",
  "g14"= "Crop", "gv2"= "Crop", "gv3"= "Crop", "gv4"= "Crop",
  "gv5"= "Crop", "gv6"= "Crop", "gv7"= "Crop", "gv8"= "Crop",

  "l2"= "Crop", "l3"= "Crop", "l4"= "Crop", "l12"= "Crop",
  "lv2"= "Crop", "lv3"= "Crop", "lv4"= "Crop", "lv5"= "Crop", "lv6"= "Crop",
  "lv7"= "Crop", "lv8"= "Crop", "lv9"= "Crop", "lv10"= "Crop",

  "my"= "Crop",

  "eta3"= "Percolation-Soil"
)

# # Define color palette components.
# fixed_colors <- c(
#   "Soil" = "#e45a3d",
#   "Nitrogen" = "#377eb8",
#   "Crop" = "#3AA600",
#   "Year" = "#FFAE00",
#   "Percolation" = "#985aa1"
# )

# --- Second estimation 25  ----

# Create the data frame
# new_param_df <- read.csv(
#  # "C:/Users/au710823/OneDrive - Aarhus universitet/NyMarkmodel/NLES5_SAS/reproducedfgk25.csv",
#   header = TRUE, stringsAsFactors = FALSE) |>
#   rename(
#       "Std.Error"=StdErr ,
#       "Conf.Lower"=LowerCL,
#       "Conf.Upper"=UpperCL
#   ) |> mutate(Estimation = "fgk25")

new_param_df <- readxl::read_excel("C:/Users/au710823/OneDrive - Aarhus universitet/NyMarkmodel/NLES5_SAS/fgk25.XLS",
  sheet = "Estimater") |>
  rename(
    "Std.Error"=StdErr ,
    "Conf.Lower"=LowerCL,
    "Conf.Upper"=UpperCL
  ) |> mutate(Estimation = "fgk25")


# --- Combine Data Frames ---

combined_param <- bind_rows(org_param, new_param_df) |>
  mutate(Component = recode(Parameter, !!!asignation))

# --- Create Forest Plot ---

# Filter out parameters with missing standard errors or confidence intervals
# and remove fixed parameters that are 0 (e.g., g1, gv1, l1)
combined_param_filtered <- combined_param %>%
  filter(!is.na(Std.Error)) %>%
  filter(Estimate != 0)

# Create the plot using ggplot2
# The plot displays the estimate with confidence intervals for each parameter,
# separated by model.
ggplot(combined_param_filtered, aes(x = Parameter, y = Estimate, color = Estimation)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Conf.Lower, ymax = Conf.Upper), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_flip() + # Flip the coordinates to make it a forest plot
  labs(
    title = "Comparison of Parameter Estimates (Original vs. New Model)",
    y = "Estimate",
    x = "Parameter"
  ) +
  theme_minimal() +
  facet_wrap(~ Component, scales = "free") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_blank()
  )

### Error

new_pred <- readxl::read_excel("C:/Users/au710823/OneDrive - Aarhus universitet/NyMarkmodel/NLES5_SAS/fgk25.XLS",
                                   sheet = "PredB4")

head(new_pred)


new_pred |> ggplot(aes(x=PUdvaskF, y=Udvask, colour = as.character(SoilG) )) +
  geom_point(alpha=.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "Comparison of Predicted vs. Observed Values",
    x = "Predicted Values",
    y = "Observed Values"
  ) +
  scale_color_grey()+
  theme_minimal()

new_pred |> ggplot(aes(x=RUdvaskF)) +
  geom_histogram(alpha=.7) +
  scale_color_grey()+
  theme_minimal()


#MSE RMSE
new_pred |> #group_by(Mau) |>
  summarise(
  MSE = mean((Udvask - PUdvaskF)^2),
  RMSE = sqrt(mean((Udvask - PUdvaskF)^2)),
  RMSE_rel = sqrt(mean((Udvask - PUdvaskF)^2))/mean(Udvask)*100)


new_pred |> #group_by(Mau) |>
  summarise(
    MSE = mean((sqrt(Udvask) - sqrt(PUdvaskF))^2),
    RMSE = sqrt(mean((sqrt(Udvask) - sqrt(PUdvaskF))^2)),
    RMSE_rel = sqrt(mean((sqrt(Udvask) - sqrt(PUdvaskF))^2))/mean(sqrt(Udvask))*100)
