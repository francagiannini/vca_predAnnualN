# L function NLES5 -----
nles5 <- function(Y, #Temporal tendency
                  Ntheta, # Nitrogen
                  C, # Crop
                  P, # Perc
                  # Psas, # Perc sas code
                  S, # Soil
                  # EEA=0, Fdato=1, EMA=0.05, ETS=0.05, EPJ=1/11, # NUAR coef
                  tao = -0.1108,# Parameter for time effect
                  mu = 23.51,# Parameter for base level se =4.3418
                  k = 1.5,# Parameter for scaling C and N
                  rho = 1.085 # Scaling factor accounting for bias transformation this is a problem
) {
  # Calculate the nitrogen leaching using the N function with provided parameters
  #𝐿 = 𝜏(𝑌 – 1991) + {(𝜇 + 𝜃𝑖𝑁 + 𝐶)^𝜅}(𝑃 𝑆)p

  L <-  tao * (Y - 1991) + ((mu + Ntheta + C)^k)*((P * S) * rho)

  # Calculate L using the percolation function as found in SAS code
  # L_percSAS <-  tao * (Y - 1991) + ((mu + Ntheta + C)^k)*(Psas * S) * rho
  # Calculate L by rho affecting entire expression as in SAS code
  Lwr <- (tao * (Y - 1991) + ((mu + Ntheta + C)^k)*(P * S)) * rho

  #𝐿NUAR = (𝜏(𝑌 − 1991) + {(𝜇 + 𝜃𝑖𝑁 + 𝐶)𝜅}(𝑃 𝑆)p)(1 - EEA Fdato - EMA - ETS)(1- EPJ)
  # L_nuar <- L * (1 - EEA*Fdato - EMA - ETS) * (1 - EPJ)

  return(list(L=L,
              Lwr=Lwr,
              #L_nuar=L_nuar,
              Ntheta=Ntheta,
              C=C,
              P=P,
              #Psas=Psas,
              S=S#,

              #L_percSAS=L_percSAS
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
  beta_CS = 0.049570,    # Mineral N spring in harvest year
  beta_CA = 0.157044,   # Mineral N autumn in harvest year
  beta_udb = 0.038245,  # Mineral N from grazing in harvest year
  beta_m1_M = 0.026499, # Mineral N added to previous/pre-previous crop (M1+M2)
  beta_f0 = 0.016314,   # Nitrogen fixation in harvest year
  beta_f1 = 0.026499,    # Nitrogen fixation in previous/pre-previous crop (F1+F2)
  beta_g0 = 0.014099,   # Organic N spring in harvest year
  beta_m1_G = 0.026499,   # Organic N added to previous/pre-previous crop (G1+G2)
  theta_2 = 1.205144 # Factor for N in winter crops (WC=2) or not (WC=1), default is 1.205144

) {


  # Calculate the N value based on the provided parameters and input variables
  N <- (
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
S_func <-  function(P_ler=-0.001849, CU){

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

