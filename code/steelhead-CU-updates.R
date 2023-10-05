###############################################################################
#
# Pacific Salmon Explorer data manipulation and analyses for steelhead CUs
# Author: Eric Hertz
# Date: May 16, 2023
#
###############################################################################

# This R code manipulates, analyzes and outputs data for the Pacific Salmon 
# Foundation's Pacific Salmon Explorer (www.salmonexplorer.ca), specifically 
# the following datasets for each steelhead Conservation Unit (CU):
#    > Number 1, 3 and 4 which include data on observed and estimated spawner abundance, 
#      catch and total exploitation rates. 
#    > Number 5 which includes data on estimated recruitment and spawner abundance by 
#      brood year.
#    > Number 101-103, 202 and 279 which include data on current abundance, biological 
#      benchmarks and status, and trends in abundance over time.


###############################################################################
# Datasets 1, 3 and 4
###############################################################################

#------------------------------------------------------------------------------
# 1. Load required packages
#------------------------------------------------------------------------------

source("code/required_packages.R")

#------------------------------------------------------------------------------
# 2. Load data and create matrices to store output
#------------------------------------------------------------------------------

date <- "May162023"# set date to include in output file name

raw.counts <- read.delim(file="Data/NuSEDS_escapement_data_Steelhead_Sep_2022.txt", header=TRUE)# load NuSEDS data

escape <- read.delim(file="Data/Steelhead_TRTC_2022-Sep.txt", header=TRUE)# load reconstructed CU escapement estimates

cu_master.all.spp <- read.csv("Data/Steelhead_CU_decoder_Sep-2022.csv",  header=TRUE)# load master decoder file

end.yr <- 2022 # set the end year for all data 
n.cu <- nrow(cu_master.all.spp) # EA: set the number of CUs to use for dataframe dimensions

#species <- c("CM","CK","CO","PKE","PKO","SEL","SER")
species <- c("Pink","Coho","Chinook","Chum","Sockeye-River","Sockeye-Lake","SH")# vector of species considered

# lglcounts matrix for storing data related to reconstructed spawner abundance
# Note: this is now referred to as "estimated abudnance" in the PSE; should change code to reflect this.
lglcounts <- matrix(nrow=n.cu, ncol=length(seq(1950, end.yr))+3)
lglcounts[,1] <- as.character(unique(cu_master.all.spp$cuid))
colnames(lglcounts) <- c("cu",seq(1950, end.yr),"location","species")

# NuSEDS count matrix for storing data related to observed spawner abundance
nusedscounts <- matrix(nrow=n.cu, ncol=length(seq(1950, end.yr))+3)
nusedscounts[,1] <- as.character(unique(cu_master.all.spp$cuid))
colnames(nusedscounts) <- c("cu",seq(1950, end.yr),"location","species")

#NuSEDS stream count matrix for storing data related to observed spawner abundance 
#from individual streams
stream.counts <- matrix(NA,1,end.yr-1946)
colnames(stream.counts) <- c("CUID","species","location",seq(1950,end.yr))

# totalrun matrix for storing data related to total run size
totalrun <- matrix(nrow=n.cu, ncol=length(seq(1950, end.yr))+3)
totalrun[,1] <- as.character(unique(cu_master.all.spp$cuid))
colnames(totalrun) <- c("cu",seq(1950, end.yr),"location","species")

# cdncatch matrix for storing data related to Canadian harvest
cdncatch <- matrix(nrow=n.cu, ncol=length(seq(1950, end.yr))+3)
cdncatch[,1]<-as.character(unique(cu_master.all.spp$cuid))
colnames(cdncatch)<- c("cu",seq(1950, end.yr),"location","species")

# uscatch matrix for storing data related to US harvest
uscatch <- matrix(nrow=n.cu, ncol=length(seq(1950, end.yr))+3)
uscatch[,1] <- as.character(unique(cu_master.all.spp$cuid))
colnames(uscatch) <- c("cu",seq(1950, end.yr),"location","species")

# exploit matrix for storing data related to exploitation rates
exploit <- matrix(nrow=n.cu,ncol=length(seq(1950, end.yr))+3)
exploit[,1] <- as.character(unique(cu_master.all.spp$cuid))
colnames(exploit) <- c("cu",seq(1950, end.yr),"location","species")

#------------------------------------------------------------------------------
# 3. Loop through each species and CU to populate the matrices
#------------------------------------------------------------------------------
for (a in species) {# loop through each species
  cu_master <- subset(cu_master.all.spp,spp==a)
  cu_id <- as.vector(cu_master$cuid)	
  
  for (h in cu_id){# loop through each CU
    
    spp <- as.character(cu_master$spp[cu_master$cuid==h])# ID species under consideration
    
    # subset master data files for the CU under consideration
    LGL <- subset(escape, StatArea.CU == as.character(cu_master$trtc_cu[cu_master$cuid==h]))
    NuSEDS <- subset(raw.counts,CU_findex == as.character(cu_master$ncc_cu[cu_master$cuid==h]))
    
    # extract required LGL escapement data and manipulate it
    LGL.counts <- as.numeric(gsub(",","", LGL$TE)) 
    LGL.cdncatch <- as.numeric(gsub(",","", LGL$CDN.Harvest))  
    LGL.uscatch <- as.numeric(gsub(",","", LGL$Total.Harvest)) - LGL.cdncatch
    LGL.totalrun <- as.numeric(gsub(",","", LGL$Total.Run)) 
    LGL.exploit <- as.numeric(gsub(",","", LGL$Total.ER)) 
    year <- LGL$Year
    LGL.year.count <- cbind(year,LGL.counts, LGL.cdncatch, LGL.uscatch, LGL.totalrun, LGL.exploit)
    
    # extract NuSEDS escapement data and manipulate it
    all.counts.streams<-
      cbind(NuSEDS[,6],NuSEDS[,76:(end.yr-1874)]) # NOTE dimensions set based on exact format of file
    nuseds.streams <- cbind(rep(as.character(cu_master$cuid[cu_master$cuid==h]) ,dim(all.counts.streams)[1]),rep(a,dim(all.counts.streams)[1]), all.counts.streams)
    colnames(nuseds.streams) <- c("CUID","species","location",seq(1950,end.yr))
    stream.counts <- rbind(stream.counts, nuseds.streams)
    
    all.counts<-NuSEDS[,76:(end.yr-1874)]# NOTE dimensions are based on exact format of ncc_streams_escapement file
    NuSEDS.counts<-as.numeric(colSums (all.counts, na.rm = TRUE, dims = 1))
    NuSEDS.counts[NuSEDS.counts==0] = NA
    year<-c(1950:end.yr)
    NuSEDS.year.count<-cbind(NuSEDS.counts,year)
    
    # merge LGL and NuSEDS escapement estimates into single file
    LGL.NuSEDS.escape<-merge(NuSEDS.year.count,LGL.year.count,all.x=TRUE)
    
    # remove years with questionable reconstructed spawner estimates (based on recommendations from Karl English)
    # if(spp=="Chinook"){
    # 	LGL.full<-LGL.NuSEDS.escape
    # 	LGL.NuSEDS.escape[1:35,3]<-NA
    # }
    
    #if(spp=="Sockeye-Lake"){
    #		LGL.full<-LGL.NuSEDS.escape
    #		LGL.NuSEDS.escape[1:10,3:6]<-NA
    
    #}
    
    # if(spp=="Pink"){	
    # 	if(h=="Nass-Portland-Observatory"){ 
    # 		LGL.NuSEDS.escape[c(TRUE,FALSE),2] <- NA				
    # 		}
    # 	}
    
    # populate lglcounts matrix
    lglcounts[lglcounts[,1]==h,2:74] <- round(LGL.NuSEDS.escape[,3],digits=0)		
    lglcounts[lglcounts[,1]==h,75] <- as.character(cu_master$cuid[cu_master$cuid==h])
    lglcounts[lglcounts[,1]==h,76] <- as.character(cu_master$species[cu_master$cuid==h])
    
    # populate nusedscounts matrix
    nusedscounts[nusedscounts[,1]==h,2:74] <- round(LGL.NuSEDS.escape[,2], digits = 0)		
    nusedscounts[nusedscounts[,1]==h,75] <- as.character(cu_master$cuid[cu_master$cuid==h])
    nusedscounts[nusedscounts[,1]==h,76] <- as.character(cu_master$species[cu_master$cuid==h])
    
    # populate totalrun matrix
    totalrun[totalrun[,1]==h,2:74] <- round(LGL.NuSEDS.escape[,6], digits = 0)		
    totalrun[totalrun[,1]==h,75] <- as.character(cu_master$cuid[cu_master$cuid==h])
    totalrun[totalrun[,1]==h,76] <- as.character(cu_master$species[cu_master$cuid==h])
    
    # populate cdncatch matrix
    cdncatch[cdncatch[,1]==h,2:74] <- round(LGL.NuSEDS.escape[,4], digits=0)		
    cdncatch[cdncatch[,1]==h,75] <- as.character(cu_master$cuid[cu_master$cuid==h])
    cdncatch[cdncatch[,1]==h,76] <- as.character(cu_master$species[cu_master$cuid==h])
    
    # populate uscatch matrix
    uscatch[uscatch[,1]==h,2:74] <- round(LGL.NuSEDS.escape[,5], digits = 0)		
    uscatch[uscatch[,1]==h,75] <- as.character(cu_master$cuid[cu_master$cuid==h])
    uscatch[uscatch[,1]==h,76] <- as.character(cu_master$species[cu_master$cuid==h])
    
    # populate exploit matrix
    exploit[lglcounts[,1]==h,2:74] <- LGL.NuSEDS.escape[,7]		
    exploit[lglcounts[,1]==h,75] <- as.character(cu_master$cuid[cu_master$cuid==h])
    exploit[lglcounts[,1]==h,76] <- as.character(cu_master$species[cu_master$cuid==h])
  }
}

#------------------------------------------------------------------------------
# 4. Output datasets for the Pacific Salmon Explorer
#------------------------------------------------------------------------------

# output dataset 1 part 1 (CU level spawner abundance)
ds.1.1.a <- as.data.frame(nusedscounts)
ds.1.1.a[is.na(ds.1.1.a)] <- NA
ds.1.1 <- melt(ds.1.1.a, id.vars = c("cu","location","species"),variable.name = "year",value.name = "NuSEDS counts by CU")

ds.1.2.a <- as.data.frame(lglcounts)
ds.1.2.a[is.na(ds.1.2.a)] <- NA
ds.1.2 <- melt(ds.1.2.a, id.vars = c("cu","location","species"),variable.name = "year",value.name = "LGL counts")

ds.1.25.a <- as.data.frame(totalrun)
ds.1.25.a[is.na(ds.1.25.a)] <- NA
ds.1.25 <- melt(ds.1.25.a, id.vars = c("cu","location","species"),variable.name = "year",value.name = "Total run")

ds.1 <- cbind(ds.1.2[2:5],ds.1.1[,5],ds.1.25[,5])
colnames(ds.1) <- c("CUID","Species","Year","LGL counts","NuSEDS counts by CU","Total run")

# Re-format dataframe per Katy's requests #
#alt.spp <- c("Chum", "Chinook", "Coho", "Pink (Even)", "Pink (Odd)", "Lake Sockeye", "River Sockeye")
#spp.table <- as.data.frame(cbind(species, alt.spp)) # Species name decoder table

#ds.1$Species <- as.character(ds.1$Species) # Change from factor to character 

#for (k in 1:length(ds.1$CUID)) { spp = ds.1$Species[k]; ds.1$Species[k] = as.character(spp.table[spp.table$species==spp,]$alt.spp) }

write.csv(ds.1,file=paste("output/pse-data/dataset_1part1.",date,".csv",sep=""))

# output dataset 1 part 2 (stream level spawner abundance)
stream.id <- read.delim(file="Data/sh_streams_decoder_Apr212023.txt",header=TRUE)
# For the purposes of data compilation, filter out streams that are not in NuSEDS #
# Streams where "Indicator" column is NA are not in NuSEDs 
#stream.id <- stream.id[(which(!is.na(stream.id$Indicator))),]
stream.id <- cbind(stream.id[,1], stream.id[,5:6])
colnames(stream.id) <- c("streamid","streamname","CUID")

ds.1.31.a <- as.data.frame(stream.counts)
colnames(ds.1.31.a) <- c("CUID","Species","streamname",seq(1950,end.yr))
ds.1.31.a[is.na(ds.1.31.a)] <- NA
ds.1.31.a <- ds.1.31.a[-1,] 

ds.1.31 <- melt(ds.1.31.a, id.vars = c("CUID","Species","streamname"),variable.name = "year",value.name = "NuSEDS counts by stream")
ds.1 <- merge(ds.1.31,stream.id,by=c("streamname","CUID"))
ds.1 <- ds.1 [complete.cases(ds.1 ),]

# Re-format dataframe per Katy's requests #
ds.1 <- ds.1[,c(6,3,4,5,1,2)]

#ds.1$Species <- as.character(ds.1$Species) # Change from factor to character 

#for (k in 1:length(ds.1$CUID)) { spp = ds.1$Species[k]; ds.1$Species[k] = as.character(spp.table[spp.table$species==spp,]$alt.spp) }

write.csv(ds.1,file=paste("output/pse-data/dataset_1part2.",date,".csv",sep=""))

# output dataset 3 
ds.3.24.a <- as.data.frame(cdncatch)
ds.3.24.a[is.na(ds.3.24.a)] <- NA
ds.3.24 <- melt(ds.3.24.a, id.vars = c("cu","location","species"),variable.name = "year",value.name = "CDN catch")

ds.3.23.a <- as.data.frame(uscatch)
ds.3.23.a[is.na(ds.3.23.a)] <- NA
ds.3.23 <- melt(ds.3.23.a, id.vars = c("cu","location","species"),variable.name = "year",value.name = "US catch")

ds.3 <- merge(ds.3.24, ds.3.23,by=c("cu","location","species","year"))
ds.3  <- ds.3[,-2]
colnames(ds.3) <- c("CUID","Species","Year","CDN catch","US catch")
write.csv(ds.3,file=paste("output/pse-data/dataset_3.",date,".csv",sep=""))

# output dataset 4
ds.4.a <- as.data.frame(exploit)
ds.4.a[is.na(ds.4.a)] <- NA
ds.4 <- melt(ds.4.a, id.vars = c("cu","location","species"),variable.name = "year",value.name = "Total exploitation rate")
ds.4 <- ds.4[,-2]
colnames(ds.4) <- c("CUID","Species","Year","Total exploitation rate")
write.csv(ds.4,file=paste("output/pse-data/dataset_4.",date,".csv",sep=""))

###############################################################################
# Dataset 5
###############################################################################

#------------------------------------------------------------------------------
# 1. Load required packages
#------------------------------------------------------------------------------

source("code/required_packages.R")

#------------------------------------------------------------------------------
# 2. Load data and create matrices to store output
#------------------------------------------------------------------------------
date <- "Nov212022"# set date to include in output file name

brood_tables<- read.delim(file="Data/Steelhead_age_2022_Sep.txt", header=TRUE)# load brood table data

#cu_master.all.spp <- read.delim("Data/Nass_conservationunits.txt")# load master decoder file

#cu_master.all.spp <- subset(cu_master.all.spp, use == "TRUE")# drop CUs that are not considered any further

end.yr <- 2017 # set the end year for all data 
n.cu <- nrow(cu_master.all.spp) # set number of CUs for dataframe dimensions

#species <- c("CM",  "CK",  "CO",  "PKE", "PKO", "SEL", "SER")
species <- c("Pink","Coho","Chinook","Chum","Sockeye-River","Sockeye-Lake","SH")# vector of species considered

# spawn matrix to store data related to spawner abundance
spawn <- matrix(nrow=n.cu,ncol=length(seq(1950, end.yr))+3)
spawn[,1]<-as.character(unique(cu_master.all.spp$cuid))
colnames(spawn)<- c("cu",seq(1950, end.yr),"location","species")

# recruit matrix to store data related to recruitment
recruit <- matrix(nrow=n.cu,ncol=length(seq(1950, end.yr))+3)
recruit[,1]<-as.character(unique(cu_master.all.spp$cuid))
colnames(recruit)<- c("cu",seq(1950, end.yr),"location","species")

#------------------------------------------------------------------------------
# 3. Loop through each species and CU to populate the matrices
#------------------------------------------------------------------------------

for (a in species) {# loop through each species
  cu_master <- subset(cu_master.all.spp,spp==a)
  cu_id <- as.vector(cu_master$cuid)  
  
  for (h in cu_id){# loop through each CU
    
    spp <- as.character(cu_master$spp[cu_master$cuid==h])# ID species under consideration
    
    # subset master data files for the CU under consideration
    LGL_brood<-subset(brood_tables, StatArea.CU == as.character(cu_master$trtc_cu[cu_master$cuid==h]))
    years<-seq(1950,end.yr,1)
    dim(years) <- c(68,1)
    colnames(years)<-"BroodYear"
    LGL_brood <-merge(years, LGL_brood,all.x=TRUE) 
    
    # drop incomplete brood years
    if (spp=="Chinook"){ 
      for (ii in 1:68){
        if(is.na(LGL_brood[ii,7])){
          LGL_brood[ii,12]<-NA}
        if(is.na(LGL_brood[ii,8])){
          LGL_brood[ii,12]<-NA}
        if(is.na(LGL_brood[ii,9])){
          LGL_brood[ii,12]<-NA}
        if(is.na(LGL_brood[ii,10])){
          LGL_brood[ii,12]<-NA}
      }
    }
    
    if (spp=="Chum"){
      for (ii in 1:68){
        if(is.na(LGL_brood[ii,7])){
          LGL_brood[ii,12]<-NA}
        if(is.na(LGL_brood[ii,8])){
          LGL_brood[ii,12]<-NA}
        if(is.na(LGL_brood[ii,9])){
          LGL_brood[ii,12]<-NA}
      }  	
    }
    
    if (spp=="Sockeye-Lake"){
      if(h=="Meziadin"){
        for (ii in 1:68){
          if(is.na(LGL_brood[ii,8])){
            LGL_brood[ii,12]<-NA}
          if(is.na(LGL_brood[ii,9])){
            LGL_brood[ii,12]<-NA}
          if(is.na(LGL_brood[ii,10])){
            LGL_brood[ii,12]<-NA}
        }
      }else{
        for (ii in 1:68){
          if(is.na(LGL_brood[ii,8])){
            LGL_brood[ii,12]<-NA}
          if(is.na(LGL_brood[ii,9])){
            LGL_brood[ii,12]<-NA}
          if(is.na(LGL_brood[ii,10])){
            LGL_brood[ii,12]<-NA}
        }
      }				
    }
    
    if (spp=="Coho"){
      for (ii in 1:68){
        if(is.na(LGL_brood[ii,7])){
          LGL_brood[ii,12]<-NA}
        if(is.na(LGL_brood[ii,8])){
          LGL_brood[ii,12]<-NA}
      }		
    }
    
    if (spp=="Sockeye-River"){
      for (ii in 1:65){
        if(is.na(LGL_brood[ii,7])){
          LGL_brood[ii,12]<-NA}
        if(is.na(LGL_brood[ii,8])){
          LGL_brood[ii,12]<-NA}
        if(is.na(LGL_brood[ii,9])){
          LGL_brood[ii,12]<-NA}
        if(is.na(LGL_brood[ii,10])){
          LGL_brood[ii,12]<-NA}
        if(is.na(LGL_brood[ii,11])){
          LGL_brood[ii,12]<-NA}
      }
    }
    
    # populate spawn matrix with escapement
    spawn[spawn[,1]==h,2:69] <- round(LGL_brood[,5], digits = 0)		
    spawn[spawn[,1]==h,70] <- as.character(cu_master$cuid[cu_master$cuid==h])
    spawn[spawn[,1]==h,71] <- as.character(cu_master$species[cu_master$cuid==h])
    
    # populate recruit matrix with recruits
    recruit[recruit[,1]==h,2:69] <- round(LGL_brood[,12], digits = 0)		
    recruit[recruit[,1]==h,70] <- as.character(cu_master$cuid[cu_master$cuid==h])
    recruit[recruit[,1]==h,71] <- as.character(cu_master$species[cu_master$cuid==h])
  }
}

#------------------------------------------------------------------------------
# 4. Output datasets for the Pacific Salmon Explorer
#------------------------------------------------------------------------------

# output dataset 5; note that log recuits per spawner, kalman filter estimates and residuals from stock-recruitment relationship are not presented at this time
ds.5.26.a <- as.data.frame(spawn)
ds.5.26.a [is.na(ds.5.26.a )] <- NA
ds.5.26 <- melt(ds.5.26.a , id.vars = c("cu","location","species"),variable.name = "year", value.name = "Spawners")

ds.5.27.a <- as.data.frame(recruit)
ds.5.27.a [is.na(ds.5.27.a )] <- NA
ds.5.27 <- melt(ds.5.27.a , id.vars = c("cu","location","species"),variable.name = "year", value.name = "Recruits")

ds.5 <- merge(ds.5.26, ds.5.27,by=c("cu","location","species","year"))
ds.5$KF_alpha <- NA 
ds.5$lnRS <- NA
ds.5$Rickery_resida <- NA
ds.5 <- ds.5[,-1]

colnames(ds.5) <- c("CUID","Species","Year","Spawners","Recruits","KF_alpha","lnRS","Ricker_resid")

# calculate recruits per spawner and remove brood years with recruits-per-spawner greater than 30 
ds.5$R_S<-as.numeric(ds.5$Recruits)/as.numeric(ds.5$Spawners)
ds.5$R_S [is.na(ds.5$R_S )] <- 0
subset(ds.5,R_S >30)# to track which observation are being removed
# ds.5 <- subset(ds.5,R_S <30)
# ds.5 <- within(ds.5,rm(R_S))

write.csv(ds.5,file=paste("output/pse-data/dataset_5.",date,".csv",sep=""))

###############################################################################
# Datasets 101-103, 202 and 279
###############################################################################

#------------------------------------------------------------------------------
# 1. Load required packages and functions
#------------------------------------------------------------------------------

source("code/required_packages.R")
source("code/cu_status_functions.R")# functions for quantifying CU status

#------------------------------------------------------------------------------
# 2. Load data and create matrices to store output
#------------------------------------------------------------------------------
date <- "May162023"

raw.counts <- read.delim(file="Data/NuSEDS_escapement_data_Steelhead_Sep_2022.txt", header=TRUE)# load NuSEDS data
# EA: any problem with using age table instead/as well?
escape <- read.delim(file="Data/Steelhead_TRTC_2022-Sep.txt", header=TRUE)# load reconstructed CU escapement estimates

cu_master.all.spp <- read.csv("Data/Steelhead_CU_decoder_Sep-2022.csv", header=TRUE)# load master decoder file

cu_master.all.spp <- subset(cu_master.all.spp, use %in% c("TRUE","ABUND"))# drop CUs that are not considered any 

end.yr <- 2022 # set the end year for all data (used for estimating current abundance)
bench.yr <- 2022 # if you want to estimate benchmarks based on a subset of the time series, set alternative end year here. Otherwise, set to same as end.yr

species <- c("Pink","Coho","Chinook","Chum","Sockeye-River","Sockeye-Lake","SH") # vector of species considered

n.cu <- length((unique(cu_master.all.spp$cuid))) # Number of CUs that were included in analysis
# summary matrix for storing data related to summary statistics
summary.out <- matrix(NA,n.cu,7)
summary.out[,1] <-  as.character(unique(cu_master.all.spp$cuid))
colnames(summary.out) <- c("cu","species","no. survey","no. indicator","max spw", "min spw", "gen length")

# running average matrix for storing data related to running average abundance estimates
run.avg.escape <- matrix(nrow=n.cu,ncol=length(seq(1950, end.yr))+3)
run.avg.escape[,1]<-as.character(unique(cu_master.all.spp$cuid))
colnames(run.avg.escape)<- c("cu",seq(1950, end.yr),"location","species")

# trend matrix to store trend parameters
trends<- matrix(nrow=n.cu, ncol=6)
trends[,1]<-as.character(unique(cu_master.all.spp$cuid))
colnames(trends) <- c("cu","Intercept","Slope","Percent_Change","location","species")

# benchmark matrix to store benchmarks
stat.met.bch <- matrix(NA,n.cu,22)
stat.met.bch[,1] <-  as.character(unique(cu_master.all.spp$cuid))
colnames(stat.met.bch) <- c("cu","curr_spw","curr_vs_lt","trend","trend_lower"
                            ,"trend_upper", "exp","sgen","sgen_lower","sgen_upper"
                            ,"smsy","smsy_lower", "smsy_upper", "25%_spw","75%_spw"
                            ,"15%_hab_smax","55%_hab_smax","uopt","uopt_lower"
                            ,"uopt_upper","location","species")

# status matrix to store data on status
curr_status <- matrix(NA,n.cu,28) # 28 = number of columns (shouldn't change)
curr_status[,1] <-  as.character(unique(cu_master.all.spp$cuid))
colnames(curr_status) <- c("cu","sr_red","sr_yellow","sr_green","sr_red_prob"
                           ,"sr_yellow_prob", "sr_green_prob","hist_red","hist_yellow"
                           ,"hist_green","hab_red","hab_yellow","hab_green","ratio_red"
                           ,"ratio_yellow","ratio_green","trend_red","trend_yellow"
                           ,"trend_green","trend_red_prob", "trend_yellow_prob"
                           ,"trend_green_prob","exp_red","exp_green","exp_red_prob"
                           ,"exp_green_prob","location","species")

# abundance matrix to store data on current abundance
curr_abun <- matrix(NA,n.cu,7)
curr_abun[,1] <-  as.character(unique(cu_master.all.spp$cuid))
colnames(curr_abun) <- c("cu","abund","start_year","end_year","rule","location","CUID")

#------------------------------------------------------------------------------
# 3. Loop through each species and CU to populate the matrices
#------------------------------------------------------------------------------

# Geometric mean function #
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=TRUE) / length(x[!is.na(x)]))
}

for (a in species) { # loop through each species
  cu_master <- subset(cu_master.all.spp,spp==a)
  
  cu_id <- as.vector(cu_master$cuid)  
  
  for (h in cu_id){ # loop through each CU
    
    spp <- as.character(cu_master$spp[cu_master$cuid==h]) #ID species under consideration
    
    #subset master data files for the CU under consideration
    LGL <- subset(escape, StatArea.CU == as.character(cu_master$trtc_cu[cu_master$cuid==h]))
    NuSEDS <- subset(raw.counts,CU_findex == as.character(cu_master$ncc_cu[cu_master$cuid==h]))
    
    #extract required LGL escapement data and manipulate it
    LGL.counts<-as.numeric(gsub(",","", LGL$TE))/1000
    year<-LGL$Year
    LGL.year.count<-cbind(LGL.counts,year)
    
    # extract NuSEDS escapement data and manipulate it
    all.counts<-NuSEDS[,76:(end.yr-1874)]# dimensions based on NuSEDs file (will need to be re-coded if data are changed)
    NuSEDS.counts<-as.numeric(colSums (all.counts, na.rm = TRUE, dims = 1))/1000
    NuSEDS.counts[NuSEDS.counts==0] = NA
    year<-c(1950:end.yr)
    NuSEDS.year.count<-cbind(NuSEDS.counts,year)
    
    # merge LGL and NuSEDS escapement estimates into single file
    LGL.NuSEDS.escape<-merge(NuSEDS.year.count,LGL.year.count,all.x=TRUE)
    
    # Set reference objects for row corresponding to last year of LGL.NuSEDS.escape data frame #
    end.yr.row <- which(LGL.NuSEDS.escape$year==end.yr)
    bench.yr.row <- which(LGL.NuSEDS.escape$year==bench.yr) # Set this if you want to subset years used to generate benchmarks.

    #-----------------------------------------#
    # calculate most recent abundance estimate#
    #-----------------------------------------#
    
    # if CU has no data on spawner abundance then current abundance = "Data deficient" and it is "Tier 4"
    if (length(which(!is.na(LGL.NuSEDS.escape[,2])))==0){
      if (is.na(sum(LGL.NuSEDS.escape[,3],na.rm=T))==FALSE){
        curr_abun[which(curr_abun[,1] == h),2] <- "Data deficient"
        curr_abun[which(curr_abun[,1] == h),3] <- NA
        curr_abun[which(curr_abun[,1] == h),4] <- NA
        curr_abun[which(curr_abun[,1] == h),5] <- 4
      }
    }
    # if CU has data on spawner abundance over the last generation for the CU as a whole (i.e., reconstucted abundance) then "Tier 1"
    if (length(which(!is.na(LGL.NuSEDS.escape[,3])))>0){  
      gen <- cu_master$gen_length[cu_master$cuid==h]
      curab <- LGL.NuSEDS.escape[((end.yr.row+1)-gen):end.yr.row,3]
      if(is.nan(gm_mean(curab,na.rm=T))==F){
        curr_abun[which(curr_abun[,1] == h),2] <- round(gm_mean(curab,na.rm=T)*1000)
        curr_abun[which(curr_abun[,1] == h),3] <- LGL.NuSEDS.escape[((end.yr.row+1)-gen),1]
        curr_abun[which(curr_abun[,1] == h),4] <- LGL.NuSEDS.escape[end.yr.row,1]
        curr_abun[which(curr_abun[,1] == h),5] <- 1
      }else{ 
        # if CU has data on CU level spawner abundance but not over the last generation then "Tier 2"
        NonNAindex <- which(!is.na(LGL.NuSEDS.escape[,3]))
        LastYear <- max(NonNAindex)
        FirstYear <- LastYear-(cu_master$gen_length[cu_master$cuid==h]-1)	
        if (a %in% c("PKO","PKE")) {FirstYear <- LastYear}
        window <- LGL.NuSEDS.escape[FirstYear:LastYear,]
        curr_abun[which(curr_abun[,1] == h),2] <- round(gm_mean(LGL.NuSEDS.escape[FirstYear:LastYear,3],na.rm=T)*1000)
        curr_abun[which(curr_abun[,1] == h),3] <- window[min(which(!is.na(window[,3]))),1]
        curr_abun[which(curr_abun[,1] == h),4] <- LGL.NuSEDS.escape[LastYear,1]
        curr_abun[which(curr_abun[,1] == h),5] <- 2	
      }
    }else{ 
      # if CU has data on spawner abundance from spawner surveys but not reconstructions then "Tier 3"
      if (length(which(!is.na(LGL.NuSEDS.escape[,2])))>0){
        NonNAindex <- which(!is.na(LGL.NuSEDS.escape[,2]))
        LastYear <- max(NonNAindex)
        FirstYear <- LastYear-(cu_master$gen_length[cu_master$cuid==h]-1)	
        window <- LGL.NuSEDS.escape[FirstYear:LastYear,]
        curr_abun[which(curr_abun[,1] == h),2] <- round(gm_mean(LGL.NuSEDS.escape[FirstYear:LastYear,2],na.rm=T)*1000)
        curr_abun[which(curr_abun[,1] == h),3] <- window[min(which(!is.na(window[,2]))),1]
        curr_abun[which(curr_abun[,1] == h),4] <- LGL.NuSEDS.escape[LastYear,1]
        curr_abun[which(curr_abun[,1] == h),5] <- 3			
      }
    }
    curr_abun[curr_abun[,1]==h,6] <- as.character(cu_master$cuid[cu_master$cuid==h])
    curr_abun[curr_abun[,1]==h,7] <- as.character(cu_master$species[cu_master$cuid==h])
    
    #------------------------------------#
    # calculate 25th and 50 percentiles#
    #------------------------------------#  
    
    stat.met.bch [which(stat.met.bch[,1] == h),14] <- quantile(LGL.NuSEDS.escape[1:bench.yr.row,3],probs=0.25,na.rm=T)*1000 
    stat.met.bch [which(stat.met.bch[,1] == h),15] <- quantile(LGL.NuSEDS.escape[1:bench.yr.row,3],probs=0.5,na.rm=T)*1000 
    
    #------------------------------------------------------#
    # calculate smoothed running avg. escapement and output#
    #------------------------------------------------------#
    #generate smoothed time series of escapement
    gen <-cu_master$gen_length[cu_master$cuid==h]	
    # running average function with right alignment (i.e., average for 2015 is 2012, 2013, 
    # 2014 and 2015 if generation length is 4 years)
    smooth.esc.1 <- rollapply(log(LGL.NuSEDS.escape$LGL.counts*1000),width=gen,by=1,FUN=gm_mean,na.rm=TRUE,align="right",fill=NA) 
    smooth.esc <- cbind(LGL.NuSEDS.escape,smooth.esc.1)
    smooth.esc[smooth.esc=="NaN"] = NA
    
    if(spp %in% c("PKE", "PKO")){
      smooth.esc.p <-log(merge(NuSEDS.year.count,LGL.year.count,all.x=TRUE)*1000)[,3]
    }
    
    #insert smoothed running average into matrix
    if (sum(smooth.esc[,4],na.rm=TRUE)>0) {
      if(spp %in% c("PKE", "PKO")){ smooth.p <- smooth.esc.p
      smooth.p[is.na(smooth.p)] <- -99999 }
      smooth.other <- smooth.esc
      smooth.other[is.na(smooth.other)] <- -99999
      ifelse(spp %in% c("PKE", "PKO"),run.avg.escape[run.avg.escape[,1]==h,2:(end.yr.row+1)] <- round(smooth.p, digits = 2), run.avg.escape[run.avg.escape[,1]==h,2:(end.yr.row+1)] <- round(smooth.other[,4], digits = 2))
    }
    
    run.avg.escape[run.avg.escape[,1]==h,(end.yr.row+2)] <- as.character(cu_master$cuid[cu_master$cuid==h])
    run.avg.escape[run.avg.escape[,1]==h,(end.yr.row+3)] <- as.character(cu_master$species[cu_master$cuid==h])
    
    start.yr <- subset(LGL.NuSEDS.escape, select = c("year","LGL.counts"))
    start.yr <- na.omit(start.yr)
    start.yr <- head(start.yr,1)
    start.yr <- start.yr$year
    #--------------------------------------------------------------------------#
    #Fit Bayesian linear regression to smothed abundance to estimate time trend#
    #--------------------------------------------------------------------------#
    
    if (sum(smooth.esc[,4],na.rm=TRUE)>0) { 
      model = paste("	
                    model {
                    for (i in 1:N){
                    y[i] ~ dnorm(y.hat[i], tau)
                    y.hat[i] <- a + b * x[i]
                    }
                    a ~ dnorm(0, .0001)
                    b ~ dnorm(0, .0001)
                    tau <- pow(sigma, -2)
                    sigma ~ dunif(0, 100)
                    }")
      
      cat(model, file = "linear.escape.bug")
      
      my.mod.full = jags(data = list("x"= (LGL.NuSEDS.escape$year-start.yr),
                                     "y"= smooth.esc[,4],"N"= length(smooth.esc[,4])), 
                         inits = , parameters.to.save= c("a","b"), model.file="linear.escape.bug", 
                         n.chains = 1, n.burnin = 5000, n.thin = 5, n.iter = 100000, DIC = TRUE)
      
      ifelse(spp %in% c("PKO", "PKE"), start<-length(smooth.esc[,4])-((3*gen)-1),start<-length(smooth.esc[,4])-9)
      start.full<-min(which(!is.na(smooth.esc[,4])))
      end.full<-length(smooth.esc[,4])-max(which(!is.na(smooth.esc[,4])))
      out.full<-as.mcmc(my.mod.full)
      raw.change.full<-((exp(out.full[[1]][,2]*(length(start.full:(length(smooth.esc[,4])-end.full)))))-1)*100	  	
      trends[trends[,1]==h,2]<-round(median(out.full[[1]][,1]), digits = 2)
      trends[trends[,1]==h,3]<-round(median(out.full[[1]][,2]), digits = 4)
      trends[trends[,1]==h,4]<-round(median(raw.change.full), digits = 0)		
    }
    
    trends[trends[,1]==h,5] <- as.character(cu_master$cuid[cu_master$cuid==h])
    trends[trends[,1]==h,6] <- as.character(cu_master$species[cu_master$cuid==h])
    
    #---------------#
    #Run status code#
    #---------------#
    
    # set current abundance and exploitation estimates
    AvgEsc <- gm_mean(LGL.NuSEDS.escape[(length(LGL.NuSEDS.escape[,3]) -(cu_master$gen_length[cu_master$cuid==h]-1)):length(LGL.NuSEDS.escape[,3]),3], na.rm=T)*1000
    if(a %in% c("PKE","PKO")) {AvgEsc <- tail(LGL[!is.na(LGL$TE),]$TE, 1)} 
    AvgER <- mean(tail(LGL$Total.ER, cu_master$gen_length[cu_master$cuid==h]), na.rm=T) 
    if(a %in% c("PKE","PKO")){ AvgER <- tail(LGL[!is.na(LGL$Total.ER),]$Total.ER, 1) } 
    
    i<-cu_master$sr_code[cu_master$cuid==h]
    
    # import alpha and beta posterior estimates from heirarchical Bayesian stock-recruitment analysis for each species and then derive stock recruitment based benchmarks against which to evaluate status  
    if (spp=="Chinook"){
      if (is.na(i)==F){
        d2<-read.table(file="Data/bayesian_posterior_samples/CK.post.out",header=T)
        nsims<- dim(d2)[1]
        i<-cu_master$sr_code[cu_master$cuid==h]
        ii<- which(names(d2)==paste("a.",i,".",sep=""));jj<- which(names(d2)==paste("b.",i,".",sep=""))
        bench<- round(getProd(d2[,ii],d2[,jj],AvgEsc,AvgER),digits=2)
        #if(i==1){bench[7,] <- NA}# make NA for Ecstall becasue no current abundance estimate
        #if(i==7){bench[7,] <- NA}# make NA for Upper Bulkley becasue no current abundance estimate
      }else{
        bench <- matrix(NA,10,3)
      }
    }
    
    if (spp=="Chum"){
      if (is.na(i)==F){
        d2<-read.table(file="Data/bayesian_posterior_samples/CM.post.out",header=T)
        nsims<- dim(d2)[1]
        i<-cu_master$sr_code[cu_master$cuid==h]
        ii<- which(names(d2)==paste("a.",i,".",sep=""));jj<- which(names(d2)==paste("b.",i,".",sep=""))
        bench<- round(getProd(d2[,ii],d2[,jj],AvgEsc,AvgER),digits=2)
      }else{
        bench <- matrix(NA,10,3)
      }
    }
    
    if (spp=="Sockeye-Lake"){
      if (is.na(i)==F){
        d2<-read.table(file="Data/bayesian_posterior_samples/SX.post.out",header=T)
        nsims<- dim(d2)[1]
        i<-cu_master$sr_code[cu_master$cuid==h]
        ii<- which(names(d2)==paste("a.",i,".",sep=""));jj<- which(names(d2)==paste("b.",i,".",sep=""))
        bench<- round(getProd(d2[,ii],d2[,jj],AvgEsc,AvgER),digits=2)
        #if(i==9){bench[7,] <- NA}# make NA for Johnson becasue no current abundance estimate		
      }else{
        bench <- matrix(NA,10,3)
      }
    }
    
    if (spp=="Pink") {
      if (is.na(i)==F){
        d2<-read.table(file="Data/bayesian_posterior_samples/PK.post.out",header=T)
        nsims<- dim(d2)[1]
        i<-cu_master$sr_code[cu_master$cuid==h]
        ii<- which(names(d2)==paste("a.",i,".",sep=""));jj<- which(names(d2)==paste("b.",i,".",sep=""))
        bench<- round(getProd(d2[,ii],d2[,jj],AvgEsc,AvgER),digits=2)
      }else{
        bench <- matrix(NA,10,3)
      }
    }
    
    if (spp=="Coho"){
      if (is.na(i)==F){
        d2<-read.table(file="Data/bayesian_posterior_samples/CO.post.out",header=T)
        nsims<- dim(d2)[1]
        i<-cu_master$sr_code[cu_master$cuid==h]
        ii<- which(names(d2)==paste("a.",i,".",sep=""));jj<- which(names(d2)==paste("b.",i,".",sep=""))
        bench<- round(getProd(d2[,ii],d2[,jj],AvgEsc,AvgER),digits=2)
      }else{
        bench <- matrix(NA,10,3)
      }
    }
    
    if (spp=="SH"){ # EA: this should change if there are SER CUs to include. 
      bench <- matrix(NA,10,3)
      
      
      
    }
    
    #-----------------#
    #set status colors#
    #-----------------#
    
    G<-"#009900"
    Y<-"#FFFF00"
    Rr<-"#CC0000"
    W<-"#FFFFFF"
    
    #---------------------------------------------------------------------------------------------#	
    #estimate status by metric (only percentile and SR benchmarks are estimated in current verion)#
    #---------------------------------------------------------------------------------------------#	
    
    COSEWIC<- matrix(nrow=2, ncol=3)# trend over last three generations status  assesment
    exploit<- matrix(nrow=2, ncol=2)# exploitation		
    ltm<- matrix(nrow=1, ncol=3)# long-term mean	
    habitat<- matrix(nrow=1, ncol=3)# habitat
    
    hist.escape<- matrix(nrow=1, ncol=3)# historic escapement  	
    ifelse(AvgEsc < quantile(LGL.NuSEDS.escape[,3],probs=0.25,na.rm=T)*1000, hist.escape[1,1]<-Rr,hist.escape[1,1]<-W)
    if(AvgEsc > quantile(LGL.NuSEDS.escape[,3],probs=0.25,na.rm=T)*1000){if(AvgEsc < quantile(LGL.NuSEDS.escape[,3],probs=0.5,na.rm=T)*1000){hist.escape[1,2]<-Y}else{hist.escape[1,2]<-W}}else{hist.escape[1,2]<-W}
    ifelse( AvgEsc > quantile(LGL.NuSEDS.escape[,3],probs=0.5,na.rm=T)*1000, hist.escape[1,3]<-G,hist.escape[1,3]<-W)
    
    SR<- matrix(nrow=2, ncol=3)# stock-recruitment
    ifelse(bench[7,1] >0, SR[1,1]<-Rr,SR[1,1]<-W)
    ifelse(bench[7,2]>0, SR[1,2]<-Y,SR[1,2]<-W)
    ifelse(bench[7,3]>0, SR[1,3]<-G,SR[1,3]<-W)
    
    ifelse(bench[7,1] >0, SR[2,1]<-bench[7,1],SR[2,1]<-"0")
    ifelse(bench[7,2] >0, SR[2,2]<-bench[7,2] ,SR[2,2]<-"0")
    ifelse(bench[7,3] >0, SR[2,3]<-bench[7,3],SR[2,3]<-"0")
    
    #---------------------#
    #benchmarks and status#
    #---------------------#
    
    curr_status[which(curr_status[,1] == h),2:26] <-c(SR[1,],SR[2,],hist.escape,habitat,ltm,COSEWIC[1,], COSEWIC[2,], exploit[1,],exploit[2,])
    stat.met.bch[,1] <-  as.character(unique(cu_master.all.spp$cuid))
    stat.met.bch [which(stat.met.bch[,1] == h),2] <- ifelse(is.na(AvgEsc),"NA", round(AvgEsc,digits=0))
    stat.met.bch [which(stat.met.bch[,1] == h),7] <- ifelse(is.na(AvgER),"NA", round(AvgER,digits=2))
    if(is.na(bench[1,1])==TRUE){stat.met.bch [which(stat.met.bch[,1] == h),8:13] <- "NA" }else{ stat.met.bch [which(stat.met.bch[,1] == h),8:13] <- as.numeric(c(bench[2,1:3],bench[3,1:3])) }
    
    curr_status[which(curr_status[,1] == h),27] <- as.character(cu_master$cuid[cu_master$cuid==h])
    curr_status[which(curr_status[,1] == h),28] <- as.character(cu_master$species[cu_master$cuid==h])
    
    stat.met.bch [which(stat.met.bch[,1] == h),21] <- as.character(cu_master$cuid[cu_master$cuid==h])
    stat.met.bch [which(stat.met.bch[,1] == h),22] <- as.character(cu_master$species[cu_master$cuid==h])
  }
}


#------------------------------------------------------------------------------
# 4. Output datasets for the Pacific Salmon Explorer
#------------------------------------------------------------------------------

# output dataset 101
ds.101 <- cbind(curr_status[,27],curr_status[,28],curr_status[,-1])
ds.101 <- ds.101[,-28]
ds.101 <- ds.101[,-28]
colnames(ds.101)[1:2] <- c("CUID", "Species")
write.csv(ds.101,file=paste("output/pse-data/dataset_101.",date,".csv",sep=""))

# output dataset 102
ds.102 <- cbind(stat.met.bch[,21], stat.met.bch[,22], stat.met.bch[,-1])
ds.102 <- ds.102[,-22]
ds.102 <- ds.102[,-22]
colnames(ds.102)[1:2] <- c("CUID", "Species")
write.csv(ds.102,file=paste("output/pse-data/dataset_102.",date,".csv",sep=""))

# output dataset 103
ds.103.a <- as.data.frame(run.avg.escape)
ds.103.a [ds.103.a == "-99999"] <- NA
ds.103 <- melt(ds.103.a, id.vars = c("cu","location","species"),variable.name = "year", value.name = "avgescapelog")
ds.103  <- ds.103[,-1]
colnames(ds.103) <- c("CUID","Species","Year","avgescapelog")
write.csv(ds.103,file=paste("output/pse-data/dataset_103.",date,".csv",sep=""))

# output dataset 202
ds.202.a <- as.data.frame(trends)
ds.202 <- cbind(trends[,5], trends[,6], trends[,-1])
ds.202 <- ds.202[,-6]
ds.202 <- ds.202[,-6]
colnames(ds.202)[1:2] <- c("CUID", "Species")
write.csv(ds.202,file=paste("output/pse-data/dataset_202.",date,".csv",sep=""))

# output dataset 279
ds.279 <- as.data.frame(curr_abun)
ds.279 <- cbind(ds.279[,6], ds.279[,7], ds.279[,-1])
ds.279 <- ds.279[,-7]
ds.279 <- ds.279[,-7]
colnames(ds.279) <- c("CUID", "Species","Abundance","Start year","End year","Rule")
write.csv(ds.279,file=paste("output/pse-data/dataset_279.",date,".csv",sep=""))