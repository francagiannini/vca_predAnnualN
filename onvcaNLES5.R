library(tidyverse)
library(ggthemes)
library(ggpubr)
library(readxl)

sasoutput <- read_excel("Scenarier20190909B4_found0325.xls")


sasoutput_fil <- sasoutput |>
  select(Udvask, # Observed
         PUdvaskF, # Predicted
         # Nterm
         NS,#1
         `NA`,#2
         nlevelMin1,
         nlevelMin2,
         NlevelFix0,
         NlevelFix1,
         nlevelFix2,
         NlevelGod0,
         NlevelGod1,
         nlevelGod2,
         Nudb,
         TN
         #Crop
         Mau, #Main crop
         Mfu, #Main crop prevous year
         Vau, #Winter crop
         Vfu, #Winter crop previous year

         #Soil
         Ler1,#Clay

         #Trend (year)
         Indb_aar

         #Percolation, it is drain but for sandt and loamy soils diferentiated
         Drain

         )
