library(tidyverse)
library(tidyr)
library(dplyr)

sh_cu_stream <- read.csv("data/cu-data/2_steelhead_cu_population_lookup_2022.csv")

#remove unnecassary columns
df = subset(sh_cu_stream, select = -c(watershed_code, parkinson_cu, tautz_cu, phylogenetic_group, parkinson2005cu, COL_H, coord_certainty, comments_alternate, coord_source))

# add species (spp) column with all values SH 
df['sppqualifiedid'] = 'SH'

df['Indicator'] = "NA"


write.csv(df, "output/pse-data/sh_streams_decoder_Aug.2022.csv", row.names = FALSE)
