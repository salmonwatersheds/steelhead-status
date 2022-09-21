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
dir.data <- "~/Salmon Watersheds Dropbox/Eric Hertz/X Drive/1_PROJECTS/1_Active/Steelhead/3_Data_Analysis/steelhead-status/data/raw-data"
dir.results <- "~/Salmon Watersheds Dropbox/Eric Hertz/X Drive/1_PROJECTS/1_Active/Steelhead/3_Data_Analysis/steelhead-status/output/pse-data"
setwd(dir.data)

setwd(dir.data)

col.juv <- read.csv("columbia_2021_OBMEP_ELECTROJUVENILEESTIMATE_202206.csv", header=TRUE, na.string="")
vimibc.juv <- read.csv("vimi_2021_blackcreekfence_juvenilesteelhead_derek_022022.csv", header=TRUE)
vimikeogh.juv <- read.csv("vimi_2021_Keogh_instream_database_2021_11_16_trevordavies.csv", header=TRUE, na.string="")
vimichek.juv <- read.csv("VIMI_2021_Korman&Schick_Cheakamus_Adult_and_Juvenile_.csv", header=TRUE, na.string="")

# create blank dataframe
Names <- c("Conservation Unit Name",	"Species",	"Year",	"Abundance",	"Latitude",	"Longitude",	"Enumeration Location Name/Description",	"Enumeration Method")

df.juv = data.frame(matrix(vector(), 0, 8,
                       dimnames=list(c(), Names)),
                stringsAsFactors=F)

### Black Creek
vimibc.juv <- vimibc.juv %>%
  rename (Year = year,
          Abundance = total_juvenile_steelhead
          ) %>%
  select(., Year, Abundance)

df.juv <- 
  bind_rows(df.juv,vimibc.juv)

df.juv$Species <- "Steelhead"
df.juv$Latitude <- 49.850802
df.juv$Longitude <- -125.100656
df.juv$Enumeration.Method <- "Fence"
df.juv$Conservation.Unit.Name <- "East Vancouver Island Winter"
df.juv$Enumeration.Location.Name.Description <- "Black Creek"



# Keough
vimikeogh.juv <- vimikeogh.juv %>%
  rename (Year = year,
          Abundance = N_wild,
          Species = species
  ) %>%
  select(., Year, Abundance, Species,life_stage)

vimikeogh.juv <-
 subset (vimikeogh.juv, Species == "WST")  %>%
  subset (., life_stage == "juvenile (smolt)") %>%
  select(., Year, Abundance, Species)

vimikeogh.juv$Conservation.Unit.Name <- strrep("East Vancouver Island Winter",1)
vimikeogh.juv$Species <- strrep("Steelhead",1)
vimikeogh.juv$Latitude <- 50.67593152
vimikeogh.juv$Longitude <- -127.3500607
vimikeogh.juv$Enumeration.Location.Name.Description <- "Keogh River"
vimikeogh.juv$Enumeration.Method <- "Fence"

df.juv <- 
  bind_rows(df.juv,vimikeogh.juv)

## chekamus
vimichek.juv <- vimichek.juv %>%
  rename (Enumeration.Location.Name.Description = River
  ) %>%
  select(., Enumeration.Location.Name.Description, Year, Abundance, Age)%>%
  subset (., Age == "1+") %>%
  select(., Enumeration.Location.Name.Description, Year, Abundance)

vimichek.juv$Conservation.Unit.Name <- strrep("South Coast Winter",1)
vimichek.juv$Species <- strrep("Steelhead",1)
vimichek.juv$Latitude <- 49.79119771
vimichek.juv$Longitude <- -123.1672473
vimichek.juv$Enumeration.Method <- "Electrofishing"
vimichek.juv$Longitude[vimichek.juv$Enumeration.Location.Name.Description=="Brohm"] <- -123.120242
vimichek.juv$Latitude[vimichek.juv$Enumeration.Location.Name.Description=="Brohm"] <- 49.806974


vimichek.juv$Abundance <- as.numeric(vimichek.juv$Abundance) 

df.juv <- 
  bind_rows(df.juv,vimichek.juv)

### write csv
setwd(dir.results)
write.csv(df.juv, "Steelhead_juvenile_20220907.csv", row.names=FALSE)


library(openxlsx)
library("rio")
xls <- dir(pattern = "xlsx")
created <- mapply(convert, xls, gsub("xlsx", "csv", xls))

## Familiarization with dataset

file.info("~/YourDirectoryHere/file_name.csv")$size

#an initial look at the data frame
str(df)


## Check for structural errors

#Variable labels

df <- df %>% rename(employees = How.many.employees.does.your.company.or.organization.have.)

colnames(df)[2]

## Check for data irregularities 

## Document data versions and changes made 