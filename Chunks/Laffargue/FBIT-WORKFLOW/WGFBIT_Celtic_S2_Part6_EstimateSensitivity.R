#####################################
######### WGFBIT - October 2019 ##########
##################################### >)))°>
######## SCRIPT dataset compilation - Celtic/BoB Eco-region (WGFBIT)
##
### R version 3.6.2

#################################
##
### PART R - SENSITIVITY Estimate
##
#################################


## SELECT SURVEY DATASET
SURVEY <- "EVHOE"  # EVHOE, DEMERSALES, ...

SEL.FishVAR <- "SurfSARmean"

## SELECT FISHING THRESHOLD
SEL.FISHthresh <- T # T: Selection of data only under fishing level threshold (or F)
FISHthresh <- 0.5 # Fishing SAR threshold
##

LIST.ECOREG <- "Bay of Biscay and the Iberian Coast"


### LOAD FILES
## ALREADY CREATED DATA FILES from PREVIOUS SCRIPTS
# Names of Files created from previous script parts
if(T){
SEL.FILEnames <- c(
"WGFBIT_Celtic_S2_BIOdataset.RData",
"WGFBIT_Celtic_S2_BIOLOGICAL_BTA.RData",
"WGFBIT_Celtic_S2_FISHINGdataset.RData",
"WGFBIT_Celtic_S2_HABITATdataset.RData")

SEL.DIRname <- DIR.OUPUT[[paste(1)]]

for(FILE in SEL.FILEnames){
load(file.path(SEL.DIRname,FILE))
}
}


### LOAD MODELLING RESULTS
# Specific from biological dta input (survey's specific)

SEL.FILEnames2 <- c(
"WGFBIT_Celtic_S2_RASTERstack.RData",
"RESULT_GLMMmodels.RData",
"Coefficients_Bdata.RData")

if(SEL.FISHthresh){SEL.DIRname2 <- file.path(SEL.DIRname,SURVEY,paste("FishPresThresh",FISHthresh,sep="_"))}else{SEL.DIRname2 <- file.path(SEL.DIRname,SURVEY,paste("FishPresThresh","None",sep="_"))}

for(FILE in SEL.FILEnames2){
load(file.path(SEL.DIRname2,FILE))}


### LOAD STANDARD GRIDS
GRIDS <- list()
TEMP.FILES <- list.files(DIR.SOURCE[[paste("5b")]])
TEMP.FILES <- TEMP.FILES[grep("grid_sensitivity",TEMP.FILES)]
for(FILE in TEMP.FILES){
GRIDS[[paste(FILE)]] <- loadRData(file=file.path(DIR.SOURCE[[paste("5b")]],FILE))
rm(FILE)}

rm(TEMP.FILES)



### PREDICT IMPACT
#### DERIVED FROM FBIT SCRIPT

#### !!! WARNING ####
#### !!! WARNING ####
#### !!! WARNING use UNIQUE DEPLETION RATE FOR TEST ONLY !!!!! ####
###################################################################
RASTERstack$Fd=0.1*RASTERstack[[SEL.FishVAR]]
#### !!! WARNING ####
#### !!! WARNING ####
#### !!! WARNING ####
###################################################################



###################################################################
###################################################################
## Compute Relative Benthic status from RBS function
# Computed by cell !!! TIME CONSUMING !!!
RASTERstack$RBS <- NA
for(i in which(is.na(getValues(RASTERstack$intercept))==F)){
RASTERstack$RBS[i] <- RBS(RASTERstack$Fd[i],RASTERstack$slope[i],RASTERstack$intercept[i])}

RASTERstackproj <- stack(projectRaster(RASTERstack,crs=CRS(HABITATdataset$EUNISproj)))

##
###################################################################
###################################################################



###################################################################
###################################################################
## TEST PLOT
plot(RASTERstack[[c(SEL.FishVAR,"medLong","RBS")]])

TEMP.EXTENT <- merge(extent(HABITATdataset$RegionLim$BoundariesDefinition$Celtic),extent(HABITATdataset$RegionLim$BoundariesDefinition$BoB))

TEMP.RASTER <- crop(RASTERstackproj,TEMP.EXTENT)
for(NAME in c(SEL.FishVAR,paste(SEL.FishVAR,"Thresh",sep="_"),"medLong","RBS")){
windows()
plot(TEMP.RASTER[[NAME]])}

if(F){TEMP.RASTER[["RBS"]][is.na(TEMP.RASTER[["RBS"]])] <- -1
TEMP.RASTER[["RBS"]][TEMP.RASTER[["RBS"]]>0] <- NA}

## Result correlation tests
if(F){
TEMP.TESTsel <- c(SEL.FishVAR,"RBS","medLong")
TEMP.TEST <- values(TEMP.RASTER[[TEMP.TESTsel]])
TEMP.TEST <- TEMP.TEST[is.na(rowSums(TEMP.TEST))==F,]
cor(TEMP.TEST)}
##
###################################################################
###################################################################



##########################
## Save results
if(file.exists(SEL.DIRname2)){
setwd(SEL.DIRname2)
save(RASTERstack,file="WGFBIT_Celtic_S2_RASTERstack_RBS.RData")
save(RASTERstackproj,file="WGFBIT_Celtic_S2_RASTERstackproj_RBS.RData")}
##########################
##########################
##########################




###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################





### TO COMPLETE FROM HERE
### TO COMPLETE FROM HERE
### TO COMPLETE FROM HERE
### TO COMPLETE FROM HERE
### FBIT SOURCE SCRIPT = Habitatstatefishing.R


#### TO COMPLETE BY GEAR !!!!!
#### TO COMPLETE BY GEAR !!!!!
#### TO COMPLETE BY GEAR !!!!!
#### TO COMPLETE BY GEAR !!!!!
#### TO COMPLETE BY GEAR !!!!!

  ## Extract fisheries dataset
  for(YEAR in names(FISHINGdataset)){}
  YEAR <- 2009
  FisheriesMet <- FISHINGdataset[[paste(YEAR)]][["DBF"]]

  FisheriesMet$csquares <- FisheriesMet$c_square

##################################################################
# calculate state for each year

# loop each year and calculate state 
  state_year <- c()
  ccname <- c()

  FisheriesMet <- FisheriesMet[is.na(match(FisheriesMet$csquares,Region@data$csquares))==F,] # Select fisheries dataset in the Region of interest

  FisheriesMet<-cbind(FisheriesMet, Region@data[match(FisheriesMet$csquares,Region@data$csquares), c("intercept","slope")])
  
  
  # depletion rates is based on  Hiddink et al. PNAS 2017 / Rijnsdorp et al. ICES 2020
  
    Parameter.DEPLETION <- 
  data.frame(rbind(c("DRB_MOL",0.200),c("OT_CRU",0.100),c("OT_MIX_CRU",0.100),c("OT_MIX_CRU_DMF",0.100),c("OT_DMF",0.026),c("OT_MIX",0.074),c("OT_MIX_DMF_BEN",0.074),c("OT_MIX_DMF_PEL",0.074),c("OT_SPF",0.009),c("SDN_DMF",0.009),c("SSC_DMF",0.016),c("TBB_CRU",0.060),c("TBB_DMF",0.140),c("TBB_MOL",0.060)),stringsAsFactors=F)
    colnames(Parameter.DEPLETION) <- c("GEAR","Parameter")
    Parameter.DEPLETION$Parameter <- as.numeric(Parameter.DEPLETION$Parameter)

 
  for (i in 1: length(Period)){
    loopdata <- FisheriesMet
    Depl_DRB_MOL          <- 0.200 * loopdata[,paste("DRB_MOL_surface_sar",Period[i],sep="_")]
    Depl_OT_CRU           <- 0.100 * loopdata[,paste("OT_CRU_surface_sar",Period[i],sep="_")]
    Depl_OT_MIX_CRU       <- 0.100 * loopdata[,paste("OT_MIX_CRU_surface_sar",Period[i],sep="_")]
    Depl_OT_MIX_CRU_DMF   <- 0.100 * loopdata[,paste("OT_MIX_CRU_DMF_surface_sar",Period[i],sep="_")]
    Depl_OT_DMF           <- 0.026 * loopdata[,paste("OT_DMF_surface_sar",Period[i],sep="_")]
    Depl_OT_MIX           <- 0.074 * loopdata[,paste("OT_MIX_surface_sar",Period[i],sep="_")]
    Depl_OT_MIX_DMF_BEN   <- 0.074 * loopdata[,paste("OT_MIX_DMF_BEN_surface_sar",Period[i],sep="_")]
    Depl_OT_MIX_DMF_PEL   <- 0.074 * loopdata[,paste("OT_MIX_DMF_PEL_surface_sar",Period[i],sep="_")]
    Depl_OT_SPF           <- 0.009 * loopdata[,paste("OT_SPF_surface_sar",Period[i],sep="_")]
    Depl_SDN_DMF          <- 0.009 * loopdata[,paste("SDN_DMF_surface_sar",Period[i],sep="_")]
    Depl_SSC_DMF          <- 0.016 * loopdata[,paste("SSC_DMF_surface_sar",Period[i],sep="_")]
    Depl_TBB_CRU          <- 0.060 * loopdata[,paste("TBB_CRU_surface_sar",Period[i],sep="_")]
    Depl_TBB_DMF          <- 0.140 * loopdata[,paste("TBB_DMF_surface_sar",Period[i],sep="_")]
    Depl_TBB_MOL          <- 0.060 * loopdata[,paste("TBB_MOL_surface_sar",Period[i],sep="_")]
    



###

#####













###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################