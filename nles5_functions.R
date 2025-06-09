# L function NLES5 -----
nles5 <- function(Y,#Temporal tendency
                  Ntheta,# Nitrogen
                  C, # Crop
                  P,# Perc
                  Psas,# Perc sas code
                  S, # Soil
                  # EEA, Fdato, EMA, ETS, EPJ, # NUAR coef
                  tao = -0.1108, # Parameter for time effect
                  mu = 23.51,# Parameter for base level se =4.3418
                  k = 1.5,# Parameter for scaling C and N
                  rho = 1.085 # Scaling factor accounting for bias transformation this is a problem
                  ) {
                  # Calculate the nitrogen value using the N function with provided parameters
                  #𝐿 = 𝜏(𝑌 – 1991) + {(𝜇 + 𝜃𝑖𝑁 + 𝐶)^𝜅}(𝑃 𝑆)p

                  L <-  tao * (Y - 1991) + ((mu + Ntheta + C)^k) * (P * S) * rho

                  L_percSAS <-  tao * (Y - 1991) + ((mu + Ntheta + C)^k) * (Psas * S) * rho

                  #Lwr <- (tao * (Y - 1991) + ((mu + Ntheta + C)^k)*(P * S)) * rho #this is how it is in SAS

                  # Calculate the nitrogen value using the N function with provided parameters
                  #𝐿NUAR = (𝜏(𝑌 − 1991) + {(𝜇 + 𝜃𝑖𝑁 + 𝐶)𝜅}(𝑃 𝑆)p)(1 - EEA Fdato - EMA - ETS)(1- EPJ)
                  #L_nuar <- L * (1 - EEA*Fdato - EMA - ETS) * (1 - EPJ)

                  return(list(
                    L = L,
                    #Lwr,
                    #L_nuar=L_nuar,
                    Ntheta = Ntheta,
                    C = C,
                    P = P,
                    Psas = Psas,
                    S = S,
                    L_percSAS = L_percSAS
                  ))

                  }

# N component ----
# 1.1 Define the Nitrogen Model Function
# The model calculates N based on various nitrogen parameters and input variables.
# N = βt NT + βCS MNCS + βCA MNCA + βudb MNudb + βm1 (M1+M2)/2 + βf0 F0 +
#     βf1 (F1+F2)/2 + βg0 G0 + βm1 (G1+G2)/2

N_func <- function(
    #Inputs
  NT, # Total N in topsoil (0-25 cm) (ton N)
  MNCS, # Mineral N spring in harvest year
  MNCA, # Mineral N autumn in harvest year
  MNudb, # Mineral N from grazing in harvest year
  M1, M2, # N mineral prev years
  F0, F1, F2, #N from N fixation
  G0, G1, G2, #N from Organic fertilization
  WC,  # autum vegetation

  #define param
  beta_t = 0.456793,    # Total N in topsoil (0-25 cm) (ton N)
  beta_CS = 0.04957,    # Mineral N spring in harvest year
  beta_CA = 0.157044,   # Mineral N autumn in harvest year
  beta_udb = 0.038245,  # Mineral N from grazing in harvest year
  beta_m1_M = 0.026499, # Mineral N added to previous/pre-previous crop (M1+M2)
  beta_f0 = 0.016314,   # Nitrogen fixation in harvest year
  beta_f1 = 0.026499,    # Nitrogen fixation in previous/pre-previous crop (F1+F2)
  beta_g0 = 0.014099,   # Organic N spring in harvest year
  beta_m1_G = 0.0265,   # Organic N added to previous/pre-previous crop (G1+G2)
  theta_2 = 1.205144 # Factor for N in winter crops (WC=2) or not (WC=1), default is 1.205144

) {


  # Calculate the N value based on the provided parameters and input variables
  N <- (
    beta_t * NT +
      beta_CS * MNCS +
      beta_CA * MNCA +
      beta_udb * MNudb +
      beta_m1_M * (M1 + M2) / 2 +
      beta_f0 * F0 +
      beta_f1 * (F1 + F2) / 2 +
      beta_g0 * G0 +
      beta_m1_G * (G1 + G2) / 2
  )
  Ntheta <- ifelse(WC==1,N, N * theta_2)
  return(Ntheta)
}

# C component ----


C_func <- function(M, W, MP, WP) {
  # Parameter tables
  M_p <- c(
    M1 = 0,              # Vintersæd (Winter cereals)
    M2 = -6.744,         # Vårsæd (Spring cereals)
    M3 = -7.279,         # Bælgsæd-korn blanding (Legume-cereal mix)
    M4 = -13.493,        # Græs og kløvergræs (Grass and clover grass)
    M5 = -17.478,        # Frøgræs (Seed grass)
    M6 = -11.192,        # Brak (Fallow)
    M8 = -0.640,         # Sukkerroer og foderroer (Sugar beets and fodder beets)
    M9 = 3.534,          # Majshelsæd og kartofler (Maize silage and potatoes)
    M10 = -7.319,         # Vinterraps (Winter rapeseed)
    M11 = -1.248,        # Vintersæd efter græs, kløvergræs, frøgræs og brak (Winter cereals after grass, clover, seed grass, or fallow)
    M12 = 19.524,        # Majshelsæd efter græs, kløvergræs, frøgræs og brak (Maize silage after grass, clover, seed grass, or fallow)
    M13 = -6.229,        # Vårsæd efter græs, kløvergræs, frøgræs og brak (Spring cereals after grass, clover, seed grass, or fallow)
    M14 = -2.866         # Bælgsæd og vårraps (Legumes and spring rapeseed)
  )

  W_p <- c(
    W1 = 0,              # Vintersæd (Winter cereals)
    W2 = -2.055,         # Bar jord (Bare soil)
    W3 = -0.456,         # Bar jord efter majshelsæd og kartofler (Bare soil after maize silage and potatoes)
    W4 = -15.959,        # Efterafgrøder, undersået græs og brak (Catch crops, undersown grass, and fallow)
    W5 = -3.792,         # Ukrudt og spildkornsplanter (Weeds and volunteer cereals)
    W6 = -14.596,        # Græs, kløvergræs, vinterraps, roer, frøgræs og udlægsafgrøder (Grass, clover, winter rapeseed, beets, seed grass, and undersown crops)
    W8 = -21.060,         # Vintersæd efter græs og frøgræs (Winter cereals after grass and seed grass)
    W9 = -1.049         # Græs og kløvergræs pløjet sent efterår eller vinter (Grass and clover plowed late autumn or winter)
  )

  MP_p <- c(
    MP1 = 0,             # Vintersæd (Winter cereals)
    MP2 = 2.847,         # Andre afgrøder end vintersæd, græs, kløvergræs, frøgræs og brak (Other crops)
    MP3 = 0.664,         # Græs, kløvergræs, frøgræs og brak (Grass, clover, seed grass, and fallow)
    MP4 = 1.160,          # Vår- og vinterafgrøder efter græs, kløvergræs, frøgræs og brak (Spring/winter crops after grass, clover, seed grass, or fallow)
    MP12 = 2.847
  )

  WP_p <- c(
    WP1 = 0,             # Vintersæd (Winter cereals)
    WP2 = 9.704,         # Bar jord og spildkorn (Bare soil and volunteer cereals)
    WP3 = 10.601,        # Græs og kløvergræs (Grass and clover)
    WP4 = 9.354,         # Efterafgrøder (Catch crops)
    WP5 = 13.241,        # Frøgræs og brak (Seed grass and fallow)
    WP6 = 5.483,         # Sukkerroer, foderroer og hamp (Sugar beets, fodder beets, and hemp)
    WP7 = -1.572,        # Bar jord efter majshelsæd og kartofler (Bare soil after maize silage and potatoes)
    WP8 = 7.413,         # Vinterraps (Winter rapeseed)
    WP9 = 7.396,         # Bar jord eller vintersæd efter græs, kløvergræs, frøgræs eller brak ompløjet *forår* (Bare soil or winter cereals after spring-plowed grass, clover, seed grass, or fallow)
    WP10 = 10.975        # Bar jord eller vintersæd efter græs, kløvergræs, frøgræs eller brak ompløjet *efterår* (Bare soil or winter cereals after autumn-plowed grass, clover, seed grass, or fallow)
  )

  #Look up values
  C <- M_p[[paste("M", M, sep = "")]] + #paste("M",M[[1]], sep = "")
    W_p[[paste("W", W, sep = "")]] +

    MP_p[[paste("MP", MP, sep = "")]] +
    WP_p[[paste("WP", WP, sep = "")]]

  return(C)
}

# S component ----
S_func <-  function(P_ler=0.001849, CU){

  S <- exp(P_ler*CU)

  return(S)
}

# P component documentation ----
P_func <- function(jbnr,
                   AAa, #d1 April–August current year
                   AAb, #d2 Sept–March current year
                   APb, #d3   Sept–March preceding year

                   # Parameters for sandy soil (JB1 + JB3)
                   delta1s = 0.001194,     # April–August percolation on sandy soil
                   delta2s = 0.001107,     # Sept–March percolation on sandy soil
                   noo2s = 0.000856,  # Sept–March percolation in preceding year on sandy soil

                   # Parameters for loamy soil (all other JB types)
                   delta1c = 0.000798,     # April–August percolation on loamy soil
                   delta2c = 0.000745,     # Sept–March percolation on loamy soil
                   noo2c = 0.000638  # Sept–March percolation in preceding year on loamy soil
) {

  # Compute P based on soil type
  if (jbnr <= 3) {
    # Sandy soil formula
    P <- (1 - exp(-delta1s * AAa - delta2s * AAb)) * exp(-noo2s * APb)
  } else {
    # Loamy soil formula
    P <- (1 - exp(-delta1c * AAa - delta2c * AAb)) * exp(-noo2c * APb)
  }

  return(P)
}

# P component SAS script ----
Psas_func <- function(jbnr,
                      AAa, #d1 April–August current year
                      AAb, #d2 Sept–March current year
                      APb, #d3   Sept–March preceding year
                      p2,
                      p3,

                      # Parameters for sandy soil (JB1 + JB3)
                      delta1s = 0.001194,     # April–August percolation on sandy soil
                      delta2s = 0.001107,     # Sept–March percolation on sandy soil
                      noo2s = 0.000856,  # Sept–March percolation in preceding year on sandy soil

                      # Parameters for loamy soil (all other JB types)
                      delta1c = 0.000798,     # April–August percolation on loamy soil
                      delta2c = 0.000745,     # Sept–March percolation on loamy soil
                      noo2c = 0.000638  # Sept–March percolation in preceding year on loamy soil
) {

  if (jbnr <= 3) {
    # Sandy soil formula for Psas
    Psas <- (1 - exp(-delta1s * AAa - delta2s * AAb - delta2s * APb)) * exp(-noo2s * p2 - noo2s * p3)
  } else {
    # Loamy soil formula for Psas
    Psas <-  (1 - exp(-delta1c * AAa - delta2c * AAb - delta2c * APb)) * exp(-noo2c * p2 - noo2c * p3)
  }

  return(Psas)
}

# Some translation vectors to match documentation -----

# Invert the vector
translate <- c(
  #𝐿 = 𝜏(𝑌 – 1991) + {(𝜇 + 𝜃𝑖𝑁 + 𝐶)𝜅}(𝑃 𝑆)p

  # Y
  "Indb_aar"="Y",

  # N
  "NS"="MN_CS",
  "na"="MN_CA",
  "nlevelMin1"="M1",
  "nlevelMin2"="M2",
  "NlevelFix0"="F0",
  "NlevelFix1"="F1",
  "NlevelFix2"="F2",
  "NlevelGod0"="G0",
  "NlevelGod1"="G1",
  "NlevelGod2"="G2",
  "Nudb"="MN_udb",
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
  "SoilGS"="Soil group",
  "SoilGC"="Soil group",
  "SoilG"="Soil group",

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
  "jbnr" = "SoilG",
  "Clay" = "CU"
)




# Exploring the impact of each input ----
sasoutput <- read_excel("Scenarier20190909B4_found0325.xls")|>
  rename('na' = `NA`)

sasoutput_fil_NLES5 <- sasoutput |>
  select(
    PUdvaskF,    # Predicted values (response for RF/XGBoost)
    NS, na, nlevelMin1, nlevelMin2, NlevelFix0, NlevelFix1, NlevelFix2,
    NlevelGod0, NlevelGod1, NlevelGod2, Nudb, TN, Vafgr_Kappa, # Nitrogen related terms
    Mau, Mfu, Vau, Vfu, # Crop related terms
    CU,             # Clay content (Soil)
    Indb_aar,     # Year trend
    d1,d2, d3, p2,p3,
    SoilG # Percolation/Drainage related terms
  ) |> mutate(SoilG= recode_factor(SoilG, "S"=0, "C"=1)) |>
  mutate(Mau=as.factor(Mau),
         Mfu=as.factor(Mfu),
         Vau=as.factor(Vau),
         Vfu=as.factor(Vfu),
         SoilG= as.factor(SoilG))


head(sasoutput_fil_NLES5)

# --- Helper function to get the mode for categorical variables ---
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


# --- Calculate mean/mode for all input columns to be used as fixed values ---
fixed_values <- list()
for (col_name in names(sasoutput_fil_NLES5)) {
  # Identify categorical columns that use direct indices or distinct values
  if (col_name %in% c("WC", "M", "W", "MP", "WP", "jbnr")) {
    fixed_values[[col_name]] <- get_mode(sasoutput_fil_NLES5[[col_name]])
  } else {
    fixed_values[[col_name]] <- mean(sasoutput_fil_NLES5[[col_name]])
  }
}

fixed_values

# --- Perform Marginal Effect Analysis ---
# Initialize a list to store results for each varied variable
marginal_effects_results <- list()

# List of all input variables (columns) to iterate through
input_variables_to_test <- names(sasoutput_fil_NLES5)

# Loop through each input variable to study its marginal effect
for (var_to_vary in input_variables_to_test) {

  # Create a temporary data frame, initialized with fixed (mean/mode) values for all columns
  # This ensures that for each row in temp_df, all variables except 'var_to_vary' are constant.
  temp_df <- as.data.frame(lapply(fixed_values, function(x) rep(x, nrow(sasoutput_fil_NLES5))))
  names(temp_df) <- names(fixed_values)

  # Replace the column for the current 'var_to_vary' with its actual observed values
  # from the original dataset.
  temp_df[[var_to_vary]] <- sasoutput_fil_NLES5[[var_to_vary]]

  # Calculate derived inputs for the nles5 function based on the modified temp_df
  # Initialize vectors to store these derived values
  Ntheta_vals <- numeric(nrow(temp_df))
  C_vals <- numeric(nrow(temp_df))
  P_vals <- numeric(nrow(temp_df))
  S_vals <- numeric(nrow(temp_df))
  Psas_vals <- numeric(nrow(temp_df))

  # Iterate through each row of the temporary dataframe to calculate the derived components.
  # This is necessary because helper functions might not be fully vectorized for all cases
  # (e.g., paste() for C_func, conditional logic in P_func).
  for (i in 1:nrow(temp_df)) {
    row_data <- temp_df[i, ]

    Ntheta_vals[i] <- N_func(
      NT = row_data$NT, MNCS = row_data$MNCS, MNCA = row_data$MNCA, MNudb = row_data$MNudb,
      M1 = row_data$M1, M2 = row_data$M2, F0 = row_data$F0, F1 = row_data$F1, F2 = row_data$F2,
      G0 = row_data$G0, G1 = row_data$G1, G2 = row_data$G2, WC = row_data$WC
    )

    C_vals[i] <- C_func(
      M = row_data$M, W = row_data$W, MP = row_data$MP, WP = row_data$WP
    )

    P_vals[i] <- P_func(
      jbnr = row_data$jbnr, AAa = row_data$AAa, AAb = row_data$AAb, APb = row_data$APb
    )

    S_vals[i] <- S_func(
      CU = row_data$CU
    )

    Psas_vals[i] <- Psas_func(
      jbnr = row_data$jbnr, AAa = row_data$AAa, AAb = row_data$AAb, APb = row_data$APb,
      p2 = row_data$p2, p3 = row_data$p3
    )
  }

  # Now, run the main nles5 function with the derived inputs for the current 'var_to_vary'
  nles5_output <- nles5(
    Y = temp_df$Y,
    Ntheta = Ntheta_vals,
    C = C_vals,
    P = P_vals,
    Psas = Psas_vals,
    S = S_vals
  )

  # Store the results in a data frame
  results_df <- data.frame(
    varied_variable = var_to_vary,                 # Name of the variable being varied
    original_value = sasoutput_fil_NLES5[[var_to_vary]], # The actual values of the varied variable
    L = nles5_output$L,
    L_percSAS = nles5_output$L_percSAS,
    Ntheta = nles5_output$Ntheta,
    C = nles5_output$C,
    P = nles5_output$P,
    S = nles5_output$S,
    Psas = nles5_output$Psas
  )
  # Add this results data frame to our list
  marginal_effects_results[[var_to_vary]] <- results_df
}
ordered_input_levels <- c(
  # C components (from your previous 'C' group Mau, Mfu, Vau, Vfu)
  "M", "MP", "W", "WP",
  # N components (from your previous 'N' group NS, na, nlevel*, Nlevel*, Nudb, TN, Vafgr_Kappa)
  "MNCS", "MNCA", "M1", "M2", "F0", "F1", "F2", "G0", "G1", "G2", "MNudb", "NT", "WC",
  # P components (from your previous 'P' group d1, d2, d3, p2, p3, SoilGS, SoilGC, SoilG)
  # Note: p2, p3, SoilGS, SoilGC were in the original P group, but not in translate_inverse's values.
  # If those are expected to be present, you'll need to update translate_inverse and this order.
  # Assuming translate_inverse's values are exhaustive for the new Input names:
  "AAa", "AAb", "APb","p2","p3",
  "jbnr", # Corresponds to SoilG
  # S components (from your previous 'S' group CU)
  "CU",
  # Y components (from your previous 'Y' group Indb_aar)
  "Y"
)

ordered_component_levels <- c("Y", "C", "N", "P", "S")
# Combine all results into a single comprehensive data frame
final_results_df <- do.call(rbind, marginal_effects_results) |>
  mutate(Input = recode(varied_variable, !!!translate_inverse)) |>
  mutate(Component = recode(Input, !!!asignation)) |>  # Add a column to indicate the model used
  mutate(
    # Convert 'Input' and 'Component' to factors with the desired order
    varied_variable = factor(varied_variable, levels = ordered_input_levels),
    Component = factor(Component, levels = ordered_component_levels)
  )

# --- Prepare Data for Plotting ---
# Pivot the results to a long format suitable for ggplot2's facet_wrap.
# This makes it easy to plot multiple output components against the varied input.
plot_data <- final_results_df %>%
  pivot_longer(
    cols = c(L, L_percSAS),#, Ntheta, C, P, S, Psas), # Columns representing the outputs/components
    names_to = "model",               # New column for the name of the output/component
    values_to = "predicted_value"                # New column for the predicted value
  )

fixed_color_variable <- c(
  # Y
  "Indb_aar"="#FFAE00",

  # N
  "NS" = "#66B3EE",
  "na" = "#377eb8",
  "nlevelMin1" = "#2A5B8A",
  "nlevelMin2" = "#1B3C5B",
  "NlevelFix0" = "#66B3EE",
  "NlevelFix1" = "#377eb8",
  "NlevelFix2" = "#2A5B8A",
  "NlevelGod0" = "#1B3C5B",
  "NlevelGod1" = "#66B3EE",
  "NlevelGod2" = "#377eb8",
  "Nudb" = "#2A5B8A",
  "TN" = "#1B3C5B",
  "Vafgr_Kappa" = "#66B3EE",

  # C
  "Mau" = "#6BBF33",
  "Mfu" = "#3AA600",
  "Vau" = "#287300",
  "Vfu" = "#3AA600" ,

  #P
  "d1" = "#BE9FD2",
  "d2" = "#985AA1",
  "d3" = "#73417A",
  "p2" = "#4D2C52",
  "p3" = "#BE9FD2",
  "SoilGS" = "#985AA1",
  "SoilGC" = "#73417A",
  "SoilG" = "#4D2C52",

  #S
  "CU"="#e45a3d"
)

fixed_color_variable_inverse <- c(
  "Y" = "#FFAE00",
  "MNCS" = "#66B3EE",
  "MNCA" = "#377eb8",
  "M1" = "#2A5B8A",
  "M2" = "#1B3C5B",
  "F0" = "#66B3EE",
  "F1" = "#377eb8",
  "F2" = "#2A5B8A",
  "G0" = "#1B3C5B",
  "G1" = "#66B3EE",
  "G2" = "#377eb8",
  "MNudb" = "#2A5B8A",
  "NT" = "#1B3C5B",
  "WC" = "#66B3EE",
  "M" = "#6BBF33",
  "MP" = "#3AA600",
  "W" = "#287300",
  "WP" = "#3AA600",
  "AAa" = "#BE9FD2",
  "AAb" = "#985AA1",
  "APb" = "#73417A",
  "p2" = "#4D2C52",
  "p3" = "#BE9FD2",
  "jbnr" = "#4D2C52",
  "CU" = "#e45a3d"
)
library(ggrepel)

# --- Plotting: Faceted Scatter Plot ---
# Create the ggplot object.
# Each point represents a predicted output value for a specific original input value,
# with all other inputs fixed at their mean/mode.
p <- plot_data |> filter(Component %in% c(N, P)) |>
  ggplot( aes(x = original_value, y = predicted_value,
                                  colour = varied_variable,
                           group = varied_variable)) +
  geom_point() + # Scatter points with some transparency
  # Add a LOESS smooth line to visualize trends. 'se = FALSE' removes standard error shading.
  #geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
  geom_line()+
  labs(
    x = "Original Input Variable Value (others Fixed at Mean/Mode)",
    y = "Predicted Output Value"
  ) +
  # Create a grid of plots, one for each combination of varied input and output component.
  # 'scales = "free_x"' allows x-axis limits to vary per facet, crucial for different input ranges.
  # 'ncol' can be adjusted to control the number of columns in the plot grid.
  facet_grid(model~ Component, scales = "free_x") +
  scale_color_manual(values = fixed_color_variable_inverse)+
  theme_minimal() + # Use a minimal theme for a clean look
  # geom_text_repel(
  #                 aes(label = varied_variable), # What text to display (the combined label)
  #                 size = 2.5, # Adjust font size as needed
  #                 box.padding = 0.5, # Padding around the label box
  #                 point.padding = 0.5, # Padding around the point
  #                 min.segment.length = 0.1, # Minimum length of the segment connecting label to point
  #                 force = 1, # Force of repulsion between overlapping labels
  #                 max.overlaps = Inf, # Allow all labels to be drawn (use cautiously for very dense plots)
  #                 show.legend = FALSE # Do not show a legend for these labels themselves
  # )+
  theme(
    strip.text = element_text(size = 8), # Adjust font size for facet labels
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6), # Rotate and resize x-axis labels for readability
    axis.text.y = element_text(size = 7), # Resize y-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold the main title
  )+
  guides(color = guide_legend(ncol = 1))
# Print the plot
print(p)
