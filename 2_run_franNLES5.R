library(tidyverse)
library(ggthemes)
library(ggpubr)
library(patchwork)
library(readxl)
library(purrr)
# This script runs the NLES5 model using parameters from sas data frame.

sasoutput <- read_excel("Scenarier20190909B4_found0325.xls") |>
  rename("na"="NA")

# Sample Data Frame for Parameters
# Ensure the column names exactly match the arguments expected by your functions.

param_table <- data.frame(
  Y = sasoutput$Indb_aar,    # Year of main crop harvest
  NT = sasoutput$TN,         # Total N in topsoil (0-25 cm) (ton N)
  MNCS = sasoutput$NS,       # Mineral N spring in harvest year
  MNCA = sasoutput$na,       # Mineral N autumn in harvest year (now correctly named 'MNCA')
  MNudb = sasoutput$Nudb,    # Mineral N from grazing in harvest year
  M1 = sasoutput$nlevelMin1, # N mineral from previous crops (1st year)
  M2 = sasoutput$nlevelMin2, # N mineral from previous crops (2nd year)
  F0 = sasoutput$NlevelFix0, # N from N fixation current year
  F1 = sasoutput$NlevelFix1, # N from N fixation 1st previous year
  F2 = sasoutput$NlevelFix2, # N from N fixation 2nd previous year
  G0 = sasoutput$NlevelGod0, # N from Organic fertilization current year
  G1 = sasoutput$NlevelGod1, # N from Organic fertilization 1st previous year
  G2 = sasoutput$NlevelGod2, # N from Organic fertilization 2nd previous year
  WC = sasoutput$Vafgr_Kappa,     # Factor for N in winter crops (WC=2) or not (WC=1)

  M = sasoutput$Mau,         # Main Crop type for C_func
  W = sasoutput$Vau,         # Winter crop type for C_func
  MP = sasoutput$Mfu,        # Previous crop type for C_func
  WP = sasoutput$Vfu,        # Previous winter crop type for C_func

  CU = sasoutput$CU,         # Clay content topsoil for S_func (in g/kg)
  jbnr = sasoutput$jbnr,     # Soil type for P_func (1-3 for sandy, >3 for loamy)

  AAa = sasoutput$d1,        # April–August current year precipitation
  AAb = sasoutput$d2,        # Sept–March current year precipitation
  APb = sasoutput$d3#,     # Sept–March preceding year precipitation
  #p2 = sasoutput$p2,       # Perc other
  #p3 = sasoutput$p3       # Perc other


  #EEA = c(0.45, 0.45, 0), # relative effect of cover crops on leaching
  #Fdato = c(1, 1, 1), # relative effect for the sowing date of the cover crop/harvest date of the main crop
  #EMA = c(0.05, 0.05, 0.05), # relative effect of intercrops
  #ETS = c(0.2, 0.2, 0.2), # relative effect of early winter sowing
  #EPJ = c(0.04, 0.04, 0.04) # relative effect of precision farming
)

print(head(param_table))

# Function to run NLES5 for a single row of parameters
# This function takes all the columns of param_table as individual arguments
run_nles5_pmap <- function(Y,
                           NT, MNCS, MNCA, MNudb, M1, M2, F0, F1, F2, G0, G1, G2, WC,
                           M, W, MP, WP,
                           CU,
                           jbnr, AAa, AAb, APb#, #p2, p3,
                           #EEA = 0, Fdato=1, EMA=0, ETS=0, EPJ=1/11
                           ) {
  #browser()  # Pause execution for debugging

  # Calculate Ntheta using N_func
  Ntheta <- N_func(
    NT = NT,
    MNCS = MNCS,
    MNCA = MNCA,
    MNudb = MNudb,
    M1 = M1,
    M2 = M2,
    F0 = F0,
    F1 = F1,
    F2 = F2,
    G0 = G0,
    G1 = G1,
    G2 = G2,
    WC = WC
  )

  # Calculate C using C_func. M, W, MP, WP are now numeric from param_table,
  # and C_func constructs the lookup keys (e.g., "M" + M).
  C <- C_func(
    M = M,
    W = W,
    MP = MP,
    WP = WP
  )

  # Calculate S using S_func
  S <- S_func(
    CU = CU
  )

  # Calculate P using P_func
  P <- P_func(
    jbnr = jbnr,
    AAa = AAa,
    AAb = AAb,
    APb = APb
  )

  # Calculate P using P_func
  # Psas <- Psas_func(
  #   jbnr = jbnr,
  #   AAa = AAa,
  #   AAb = AAb,
  #   APb = APb,
  #   p2 = p2,
  #   p3 = p3
  # )

  # Run the main nles5 function with all necessary parameters
  L_results <- nles5(
    Y = Y,
    Ntheta = Ntheta,
    C = C,
    P = P,
    #Psas = Psas,
    S = S,
    #EEA = EEA,
    #Fdato = Fdato,
    #EMA = EMA,
    #ETS = ETS,
    #EPJ = EPJ
  )

  # Return results as a data frame row

#output nles5(): return(list(L=L,
  # Lwr=Lwr,
  # #L_nuar=L_nuar,
  # Ntheta=Ntheta,
  # C=C,
  # P=P,
  # Psas=Psas,
  # S=S,
  # L_percSAS=L_percSAS
  # ))

  return(data.frame(
    L = L_results$L,
    #Lwr = L_results$Lwr,
    #L_nuar = L_results$L_nuar,
    Ntheta = L_results$Ntheta,
    C = L_results$C,
    P = L_results$P,
    #Psas = L_results$Psas,
    S = L_results$S#,
    #L_percSAS = L_results$L_percSAS
  ))

  }



# --- Run the model for each row ---
results_list <- pmap(param_table, run_nles5_pmap)

# Combine the list of data frames (each representing a row of results)
final_results <- bind_rows(results_list) |>
  bind_cols(sasoutput)

saveRDS(final_results, file = "pred_NLES5NUAR.RDS")

# final_results |> ggplot(aes(x = PUdvaskF, y = L)) +
#   geom_point(col= "darkblue", alpha=0.2)  +
#   geom_point(aes(y=L_percSAS), col="darkgreen", alpha=0.2)+
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   scale_y_continuous(
#     sec.axis = sec_axis(
#       transform =  ~ .,
#       name = "L_percSAS" # Name for the secondary Y-axis
#     ))+
#   theme(legend.position = "bottom")
#
# final_results |> pivot_longer(cols = c(L,L_percSAS, PUdvaskF),
#                               names_to = "variable",
#                               values_to = "Annual leaching") |>
#   ggplot(aes(x = `Annual leaching`, fill = variable, col = variable)) +
#   geom_density(alpha = 0.1) +
#   theme_minimal() +
#   scale_color_manual(values = c("darkblue","darkgreen", "darkred"))+
#   scale_fill_manual(values = c("darkblue","darkgreen", "darkred"))+
#   theme(legend.position = "bottom")
#
#
# final_results |> ggplot(aes(x = Nitrogen*ifelse(Vafgr_Kappa==1,1, 1.205144), y = Ntheta)) +
#   geom_point(col= "darkblue", alpha=0.2)  +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   theme(legend.position = "bottom")
#
#
#
# final_results |> ggplot(aes(x = Drain, y = P ))+
#   geom_point(col= "purple", alpha=0.2)  +
#   geom_point(aes(y=Psas), col="pink", alpha=0.2)+
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
#   scale_y_continuous(
#     sec.axis = sec_axis(
#       transform =  ~ .,
#       name = "P_percSAS" # Name for the secondary Y-axis
#     ))+
#   theme_minimal() +
#   theme(legend.position = "bottom")

