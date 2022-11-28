############### data-cleaning ##################
# this file takes raw files and wrangles them so that they are in a format approporiate for the Pacific Salmon Explorer
################################################

### Loading packages ###

library(dplyr)
library(tibble)
library(scales)
library(ggplot2)
library(readxl)
library(reshape2)
library(stringr)

#### Juvenile survey data ####

##read in files where we have juvenile survey data ##
# Directory to draw data from #
dir.data <- "~/Salmon Watersheds Dropbox/Eric Hertz/X Drive/1_PROJECTS/1_Active/Steelhead/3_Data_Analysis/steelhead-status/data/raw-population-indicator-data"
dir.results <- "~/Salmon Watersheds Dropbox/Eric Hertz/X Drive/1_PROJECTS/1_Active/Steelhead/3_Data_Analysis/steelhead-status/output/pse-data"
setwd(dir.data)

col.juv <- read.csv("columbia_2021_OBMEP_juv.csv", header=TRUE, na.string="")
vimibc.juv <- read.csv("vimi_2021_blackcreekfence_juvenilesteelhead_derek_022022.csv", header=TRUE)
vimikeogh.juv <- read.csv("vimi_2021_Keogh_instream_database_2021_11_16_trevordavies.csv", header=TRUE, na.string="")
vimichek.juv <- read.csv("VIMI_2021_Korman&Schick_Cheakamus_Adult_and_Juvenile_SH_juveniledata.csv", header=TRUE, na.string="")
fraserseton.juv <- read.csv("Fraser_2020_Buchanan_setonriver_BCHydro_juv.csv", header=TRUE, na.string="")
skeenakispiox.juv <- read.csv("Skeena_2002_Williamson_Kispiox_juv.csv", header=TRUE, na.string="")
skeenaother.juv <- read.csv("Skeena_1993_1994_Bustard_Kitwanga_Morice_Sustut_Zymoetz_juv.csv", header=TRUE, na.string="")
vimiother.juv <- read.csv("VIMI_2020_McCulloch_VI_Rivers_juv.csv", header=TRUE, na.string="")
colother.juv <- read.csv("columbia_2021_OBMEP_juv.csv", header=TRUE, na.string="")
fraserother.juv <- read.csv("Fraser_2015_Decker_Thompson_juv.csv", header=TRUE, na.string="")


# create blank dataframe
Names <- c("Conservation Unit Name",	"Species",	"Year",	"absolute_abundance",	"Latitude",	"Longitude",	
           "Enumeration Location Name/Description",	"Enumeration Method","density_abundance")

df.juv = data.frame(matrix(vector(), 0, 9,
                       dimnames=list(c(), Names)),
                stringsAsFactors=F)

### Black Creek
vimibc.juv <- vimibc.juv %>%
  rename (Year = year,
          absolute_abundance = total_juvenile_steelhead
          ) %>%
  select(., Year, absolute_abundance)

df.juv <- 
  bind_rows(df.juv,vimibc.juv)

df.juv$Species <- "Steelhead"
df.juv$Latitude <- 49.850802
df.juv$Longitude <- -125.100656
df.juv$Enumeration.Method <- "Fence"
df.juv$Conservation.Unit.Name <- "East Vancouver Island Winter"
df.juv$Enumeration.Location.Name.Description <- "Black Creek"
#df.juv$density_abundance <- "NA"


# Keough
vimikeogh.juv <- vimikeogh.juv %>%
  rename (Year = year,
          absolute_abundance = N_wild,
          Species = species
  ) %>%
  select(., Year, absolute_abundance, Species,life_stage)

vimikeogh.juv <-
 subset (vimikeogh.juv, Species == "WST")  %>%
  subset (., life_stage == "juvenile (smolt)") %>%
  select(., Year, absolute_abundance, Species)

vimikeogh.juv$Conservation.Unit.Name <- strrep("East Vancouver Island Winter",1)
vimikeogh.juv$Species <- strrep("Steelhead",1)
vimikeogh.juv$Latitude <- 50.67593152
vimikeogh.juv$Longitude <- -127.3500607
vimikeogh.juv$Enumeration.Location.Name.Description <- "Keogh River"
vimikeogh.juv$Enumeration.Method <- "Fence"
#vimikeogh.juv$density_abundance <- "NA"

df.juv <- 
  bind_rows(df.juv,vimikeogh.juv)

## chekamus
vimichek.juv <- vimichek.juv %>%
  rename (Enumeration.Location.Name.Description = River,
          absolute_abundance = Abundance
  ) %>%
  select(., Enumeration.Location.Name.Description, Year, absolute_abundance, Age)%>%
  subset (., Age == "1+") %>%
  select(., Enumeration.Location.Name.Description, Year, absolute_abundance)

vimichek.juv$Conservation.Unit.Name <- strrep("South Coast Winter",1)
vimichek.juv$Species <- strrep("Steelhead",1)
vimichek.juv$Latitude <- 49.79119771
vimichek.juv$Longitude <- -123.1672473
vimichek.juv$Enumeration.Method <- "Electrofishing"
#vimichek.juv$density_abundance <- "NA"
vimichek.juv$Longitude[vimichek.juv$Enumeration.Location.Name.Description=="Brohm"] <- -123.120242
vimichek.juv$Latitude[vimichek.juv$Enumeration.Location.Name.Description=="Brohm"] <- 49.806974


vimichek.juv$absolute_abundance <- as.numeric(vimichek.juv$absolute_abundance) 

df.juv <- 
  bind_rows(df.juv,vimichek.juv)


## seton
fraserseton.juv <- fraserseton.juv %>%
  rename (Enumeration.Location.Name.Description = stream_name,
          Year = year,
          absolute_abundance = juvenile_abundance_estimate
  ) %>%
  select(., Enumeration.Location.Name.Description, Year, absolute_abundance)

fraserseton.juv$Conservation.Unit.Name <- strrep("Mid Fraser Summer",1)
fraserseton.juv$Species <- strrep("Steelhead",1)
fraserseton.juv$Latitude <- 50.67081909
fraserseton.juv$Longitude <- -121.9533358
fraserseton.juv$Enumeration.Method <- "Electrofishing"
#fraserseton.juv$density_abundance <- "NA"



df.juv <- 
  bind_rows(df.juv,fraserseton.juv)


## kispiox
skeenakispiox.juv <- skeenakispiox.juv %>%
  rename (Enumeration.Location.Name.Description = stream_name,
          Year = year,
          density_abundance = density_fry_per_100m2
  ) %>%
  select(., Enumeration.Location.Name.Description, Year, density_abundance)

skeenakispiox.juv$Conservation.Unit.Name <- strrep("Kispiox",1)
skeenakispiox.juv$Species <- strrep("Steelhead",1)
skeenakispiox.juv$Latitude <- 55.55680749
skeenakispiox.juv$Longitude <- -127.9224601
skeenakispiox.juv$Enumeration.Method <- "Electrofishing"
#skeenakispiox.juv$absolute_abundance <- "NA"


df.juv <- 
  bind_rows(df.juv,skeenakispiox.juv)


## other skeena
skeenaother.juv <- skeenaother.juv %>%
  rename (Enumeration.Location.Name.Description = stream_name,
          Year = year,
          density_abundance = density_parr_per_100m2,
          Conservation.Unit.Name = cu_name,
          Enumeration.Method = method,
          Latitude = lat,
          Longitude = long
  )  %>%
  select(., Enumeration.Location.Name.Description, Year, density_abundance, Conservation.Unit.Name,Enumeration.Method,Longitude,Latitude)

skeenaother.juv$Species <- strrep("Steelhead",1)


df.juv <- 
  bind_rows(df.juv,skeenaother.juv)


## other vimi
vimiother.juv <- vimiother.juv %>%
  rename (Enumeration.Location.Name.Description = stream_name,
          Year = year,
          density_abundance = density_avg_fry_per_100m2_adj,
          Conservation.Unit.Name = cu_name,
          Enumeration.Method = method,
          Latitude = lat,
          Longitude = long
  )  %>%
  select(., Enumeration.Location.Name.Description, Year, density_abundance, Conservation.Unit.Name,Enumeration.Method,Longitude,Latitude)

vimiother.juv$Species <- strrep("Steelhead",1)


df.juv <- 
  bind_rows(df.juv,vimiother.juv)


## other columbia
colother.juv <- colother.juv %>%
  rename (Enumeration.Location.Name.Description = survey_stream_name,
          Year = year,
          absolute_abundance = age1_pop_est,
          Conservation.Unit.Name = cu_name,
          Enumeration.Method = method,
          Latitude = lat,
          Longitude = long
  )  %>%
  select(., Enumeration.Location.Name.Description, Year, absolute_abundance, Conservation.Unit.Name,Enumeration.Method,Longitude,Latitude)

colother.juv$Species <- strrep("Steelhead",1)


df.juv <- 
  bind_rows(df.juv,colother.juv)


## other fraser
fraserother.juv <- fraserother.juv %>%
  rename (Enumeration.Location.Name.Description = stream_name,
          Year = year,
          absolute_abundance = age1_standing_stock,
          Conservation.Unit.Name = cu_name,
          Enumeration.Method = method,
          Latitude = lat,
          Longitude = long
  )  %>%
  select(., Enumeration.Location.Name.Description, Year, absolute_abundance, Conservation.Unit.Name,Enumeration.Method,Longitude,Latitude)

fraserother.juv$Species <- strrep("Steelhead",1)
fraserother.juv$absolute_abundance <- as.numeric(fraserother.juv$absolute_abundance) 


df.juv <- 
  bind_rows(df.juv,fraserother.juv)



### write csv
setwd(dir.results)
write.csv(df.juv, "Steelhead_juvenile_20221121.csv", row.names=FALSE)


