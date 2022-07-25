## Script to join two working files - one with steelhead CU/stream names and one with coordinates for the streams 

#load packages
library(stringr)
library(dplyr)
library(readr)
library(purrr)

#for all csvs in this folder, join by the column "stream name"
data_join <- list.files(path = "data/cu-data/steelhead-cu-2022", # Identify all CSV files
                        pattern = "*.csv", full.names = TRUE) %>% lapply(read_csv) %>%                              # Store all files in list
  reduce(full_join, by = "stream_name")

# Full-join data sets into one data set 
data_join  
# Print data to RStudio console
print(data_join)
#produce an output file that is the cu-stream lookup file with coordinates 
write.csv(data_join, "data/cu-data/2_steelhead_cu_population_lookup_2022.csv", row.names = FALSE)
