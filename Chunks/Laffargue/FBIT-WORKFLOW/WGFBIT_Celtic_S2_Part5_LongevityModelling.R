#####################################
######### WGFBIT - October 2019 ##########
##################################### >)))°>
######## SCRIPT dataset compilation - Celtic/BoB Eco-region (WGFBIT)
##
### R version 3.6.2

#################################
##
### PART R - LONGEVITY modelling
##
#################################


## SELECT SURVEY DATASET
SURVEY <- "EVHOE"  # EVHOE, DEMERSALES, ...

SEL.VARfishing <- 'SurfSARyearly'

####
## SELECT FISHING THRESHOLD
SEL.FISHthresh <- T # T: Selection of data only under fishing level threshold (or F)
FISHthresh <- 0.01 # Fishing SAR threshold

if(SEL.FISHthresh==F){FISHthresh <- "None"}
##
####

### LOAD FILES
## ALREADY CREATED DATA FILES from PREVIOUS SCRIPTS
# Names of Files created from previous script parts
if(T){
SEL.FILEnames <- c(
"WGFBIT_Celtic_S2_BIOdataset.RData",
"WGFBIT_Celtic_S2_BIOLOGICAL_BTA.RData",
"WGFBIT_Celtic_S2_FISHINGdataset.RData",
"WGFBIT_Celtic_S2_HABITATdataset.RData",
"WGFBIT_Celtic_S2_RASTERstack.RData")

SEL.DIRname <- DIR.OUPUT[[paste(1)]]

for(FILE in SEL.FILEnames){
load(file.path(SEL.DIRname,FILE))
}
}

##### CREATE OUTPUT FOLDER
if(SEL.FISHthresh==F){if(file.exists(file=file.path(DIR.OUPUT[[paste(1)]],SURVEY))==F){dir.create(file.path(DIR.OUPUT[[paste(1)]],SURVEY))}}else{
if(file.exists(file=file.path(DIR.OUPUT[[paste(1)]],SURVEY,paste("FishPresThresh",FISHthresh,sep="_")))==F){dir.create(file.path(DIR.OUPUT[[paste(1)]],SURVEY,paste("FishPresThresh",FISHthresh,sep="_")))}}
#####





### DATA SOURCE FROM STANDARD FORMAT DATASET
TEMP.DATAbio <- BIOdataset$FORMATED[[SURVEY]]$DATA_BTA
TEMP.DATAstation <- BIOdataset$FORMATED[[SURVEY]]$STATION

### Cumulated biomass computation
BioVAR <- "Biomass"
TEMP.DATAbio[,BioVAR] <- as.numeric(TEMP.DATAbio[,BioVAR])
TEMP.BiomassTrait <- TEMP.DATAbio[,grep(SEL.TRAIT,colnames(TEMP.DATAbio))]*TEMP.DATAbio[,BioVAR]

# Traits Biomasse / Station
TEMP.BiomassTrait2 <- data.frame("Station"=TEMP.DATAbio$Station,TEMP.BiomassTrait)[is.na(rowSums(TEMP.BiomassTrait))==F,] # Supprime NAs
TEMP.BiomassTrait2 <- aggregate(.~Station, data=TEMP.BiomassTrait2,FUN="sum")
rownames(TEMP.BiomassTrait2) <-  TEMP.BiomassTrait2[,1]
TEMP.BiomassTrait2 <-  TEMP.BiomassTrait2[,-1]

# traits Biomass /Station cumulated sum
TEMP.BiomassTrait3 <- t(apply(TEMP.BiomassTrait2, 1, cumsum))

### Scale cum biomass 0 to 1 ## MODIF PLAFF TO CHECK !!!!
#TEMP.BiomassTrait3 <- TEMP.BiomassTrait3/TEMP.BiomassTrait3[,ncol(TEMP.BiomassTrait3)]
### MODIF PLAFF TO CHECK !!!!

## 



# get longevity categories separate for each station 
Station <-rep(rownames(TEMP.BiomassTrait3),ncol(TEMP.BiomassTrait3)-1)

Cumb <- NULL
for(i in c(1:(ncol(TEMP.BiomassTrait3)-1))){Cumb <- c(Cumb,TEMP.BiomassTrait3[,i])}
Longevity <-c(rep(1,nrow(TEMP.BiomassTrait3)),rep(3,nrow(TEMP.BiomassTrait3)),rep(10,nrow(TEMP.BiomassTrait3)))  
fulldat   <- data.frame(Station, Cumb,Longevity) 


#### REVOIR CI-DESSOUS !!!!
fulldat$ll <- log1p(fulldat$Longevity)

#fulldat$ll <-fulldat$Longevity ### MODIF PLAFF TO CHECK !!!!

# add a small number to values very close to 0 and 1 
for (i in 1:(nrow(fulldat))){
  if (fulldat$Cumb[i] < 1e-3){ fulldat$Cumb[i] <- 1e-3}
  if (fulldat$Cumb[i] > 0.999){fulldat$Cumb[i] <- 0.999}
}   


####                                                                                                                                                                                                    
####
####

fulldat2 <- merge(fulldat, TEMP.DATAstation@data, by="Station")
#head(fulldat2)





############################
##### LONGEVITY MODELLING
######

RESULT.GLMMmodels <- list() # Object for results strorage


#########################
# FILTERING LOW FISHING PRESSURE LOCATIONS ONLY
if(SEL.FISHthresh){fulldat2 <- fulldat2[which(fulldat2[,SEL.VARfishing]<=FISHthresh),]
RESULT.GLMMmodels[["FISHthreshold"]] <- FISHthresh

#####
#####
##### CORRIGER/COMPLETER
#### REVOIR CALCUL [,SEL.VARfishing] dans RASTER #####
RASTERstack$SurfSARmean_Thresh <- RASTERstack$SurfSARmean
values(RASTERstack$SurfSARmean_Thresh)[values(RASTERstack$SurfSARmean_Thresh)>FISHthresh] <- NA
#####
#####
#####
#####
#####

}else{RESULT.GLMMmodels[["FISHthreshold"]] <- "None"

#####
#####
##### CORRIGER/COMPLETER
#####
#####
RASTERstack$SurfSARmean_Thresh <- RASTERstack$SurfSARmean}
#####
#####
#####
#####

#########################
##


## MODELS DEFINITION


## Without trawling
SEL.VAR <- c("Raster_Depth","Raster_Chl","Raster_Temp","Raster_ALLenergy","Raster_Substrate") # 
LIST.VAR <- list()
SEL.MODELS <- list()
for(i in c(1:length(SEL.VAR))){TEMP <- combn(SEL.VAR,i)
for(j in c(1:dim(TEMP)[2])){
LIST.VAR[[length(LIST.VAR)+1]] <- c("Station","Cumb","ll",as.vector(TEMP[,j]))
SEL.MODELS <- c(SEL.MODELS,as.formula(paste("Cumb","~",paste("ll",paste(as.vector(TEMP[,j]),collapse="+"),"(1 | Station)",sep="+"))))}
}

RESULT.GLMMmodels[["VAR"]] <- SEL.VAR
RESULT.GLMMmodels[["MODELS"]] <- SEL.MODELS


## With trawling
if(F){LIST.VAR <- list()
SEL.MODELS <- list()
for(i in c(1:length(SEL.VAR))){TEMP <- combn(SEL.VAR,i)
for(j in c(1:dim(TEMP)[2])){
LIST.VAR[[length(LIST.VAR)+1]] <- c("Station","Cumb","ll","SurfaceSARmean",as.vector(TEMP[,j]))
SEL.MODELS <- c(SEL.MODELS,as.formula(paste("Cumb","~",paste("ll","SurfaceSARmean",paste(as.vector(TEMP[,j]),collapse="+"),"(1 | Station)",sep="+"))))}
}
} # T/F

# TO COMPLETE !!!!
# TO COMPLETE !!!!
# TO COMPLETE !!!!
# TO COMPLETE !!!!

##



## Glmm Models computation
RESULT.GLMMmodels[["FIT"]] <- list()

RESULTS.AIC <- NULL
RESULTS.SINGULAR <- NULL
for(i in 1:length(SEL.MODELS)){
print(paste(i,"/",length(SEL.MODELS)))
TEMP.data <- na.omit(fulldat2[,LIST.VAR[[i]]])

if(nrow(TEMP.data)>0)
{

TEMP.MODEL <- try(lme4::glmer(SEL.MODELS[[i]], data=TEMP.data, family=binomial),silent =T)
if(class(TEMP.MODEL)=="glmerMod"){RESULT.GLMMmodels[["FIT"]][[i]] <- TEMP.MODEL
RESULTS.AIC <- c(RESULTS.AIC,AIC(RESULT.GLMMmodels[["FIT"]][[i]]))
RESULTS.SINGULAR <- c(RESULTS.SINGULAR,isSingular(RESULT.GLMMmodels[["FIT"]][[i]]))}else{
print("NO FIT")
RESULT.GLMMmodels[["FIT"]][[i]] <- NA
RESULTS.AIC <- c(RESULTS.AIC,NA)
RESULTS.SINGULAR <- c(RESULTS.SINGULAR,NA)}
}


rm(TEMP.data, i)}

RESULT.GLMMmodels[["AIC"]] <- RESULTS.AIC
RESULT.GLMMmodels[["SINGULAR"]] <- RESULTS.SINGULAR


## SELECT BEST MODEL
# Singular GLM fit not considered
if(length(which(RESULTS.SINGULAR==F))>0){
TEMP.SEL <-  which(RESULTS.SINGULAR==F)[which.min(RESULTS.AIC[which(RESULTS.SINGULAR==F)])]
}else{
TEMP.SEL <- which.min(RESULTS.AIC)} # Model selection criteria # or indicate i e.g. 31 #

MODELminAIC.Model <- SEL.MODELS[[TEMP.SEL]]
MODELminAIC.Coef <- fixef(RESULT.GLMMmodels[["FIT"]][[TEMP.SEL]])

RESULT.GLMMmodels[["BESTmodel"]] <- MODELminAIC.Model
RESULT.GLMMmodels[["BESTmodelCoef"]] <- MODELminAIC.Coef 

modcoeff  <-  MODELminAIC.Coef

## Compute medlong from selected model coefficients
# prepare grid specific information to predict longevity at a certain location
RES <- qlogis(0.5)- MODELminAIC.Coef[["(Intercept)"]]
for(NAME in na.omit(names(MODELminAIC.Coef)[match(SEL.VAR,names(MODELminAIC.Coef))])){
RES <- RES-MODELminAIC.Coef[NAME]*RASTERstack[[sub(".*_","",NAME)]]
rm(NAME)}
medLong <- exp(RES/MODELminAIC.Coef["ll"])
rm(RES)

# Attribute medLong to Raster stack
RASTERstack$medLong <- medLong
var_SPDF <- RASTERstack
var_SPDF$medLong <- medLong
#var_SPDF$medLong <- RES/MODELminAIC.Coef["ll"] ### MODIF PLAFF TO CHECK !!!!


# to predict impact we will not use median longevity but the longevity distribution
  # hence estimate the slope and intercept for each gridcell
RASTERstack$slope <- modcoeff["ll"]  # slope of binomial model

intercept <- modcoeff[["(Intercept)"]]
for(NAME in na.omit(names(modcoeff)[match(SEL.VAR,names(modcoeff))])){
intercept <- intercept+MODELminAIC.Coef[NAME]*RASTERstack[[sub(".*_","",NAME)]]
rm(NAME)}

RASTERstack$intercept <- intercept



#### SUMMARY & TEST PLOT

#### TO COMPLETE ####
#### TO COMPLETE ####
#### TO COMPLETE ####
#### TO COMPLETE ####

# Selected models
TEMP.MODELsel <- which(is.na(RESULT.GLMMmodels$SINGULAR)==F&RESULT.GLMMmodels$SINGULAR==F)
TEMP.MODELSvalid <- cbind(TEMP.MODELsel,as.character(RESULT.GLMMmodels$MODELS[TEMP.MODELsel]),
                          round(as.numeric(RESULT.GLMMmodels$AIC[TEMP.MODELsel]),digits=3))

TEMP.MODELSvalid <- TEMP.MODELSvalid[order(as.numeric(TEMP.MODELSvalid[,3])),]

#### TO COMPLETE ####
#### TO COMPLETE ####
#### TO COMPLETE ####
#### TO COMPLETE ####
#### TO COMPLETE ####

####



windows()
plot(var_SPDF)
plot(var_SPDF[["medLong"]])

## Save results
setwd(file.path(DIR.OUPUT[[paste(1)]],SURVEY,paste("FishPresThresh",FISHthresh,sep="_")))
save(RASTERstack,file="WGFBIT_Celtic_S2_RASTERstack.RData")
save(RESULT.GLMMmodels,file="RESULT_GLMMmodels.RData")
save(modcoeff,file="Coefficients_Bdata.RData")  

setwd(file.path(DIR.OUPUT[[paste(1)]],SURVEY))
save(RASTERstack,file="WGFBIT_Celtic_S2_RASTERstack.RData")
save(RESULT.GLMMmodels,file="RESULT_GLMMmodels.RData")
save(modcoeff,file="Coefficients_Bdata.RData")  









#### SCRIPT END ####
##############################################################################################################################################################################################################
##############################################################################################################################################################################################################
##############################################################################################################################################################################################################
##############################################################################################################################################################################################################
##############################################################################################################################################################################################################
##############################################################################################################################################################################################################








### TO COMPLETE FROM HERE
### TO COMPLETE FROM HERE
### TO COMPLETE FROM HERE
### TO COMPLETE FROM HERE

## Carto longevity classes
if(F){
summary(var_SPDF[["medLong"]])
cuts=c(0, 5, 10,12,15,18,20,25,30) #set breaks
pal <- colorRampPalette(c("white","black"))
blueorange <- c("#feb24c","#fed976","#ffffcc","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494")
pal <- colorRampPalette(c("green","red"))
plot(var_SPDF[["medLong"]],breaks=cuts, col = blueorange)







##############
########## FROM SOURCE SCRIPT 
######### NS_CS_longevity.R
if(F){#### script to calculate longevity composition of benthic communities in North Sea and Celtic Sea
#### following relationships by Rijnsdorp et al (Ecol App 2018)

  #### calculate median longevity and fraction of biomass in <3 3-10 and more than 10 years
  b0 <- -5.682973    ### intercept
  b1 <-  3.502197    ### ln(longevity)
  b2 <- -0.082575    ## ln(trawling)
  b3 <-  0.021062    ### mud
  b4 <-  0.018634    ### gravel
  b5 <-  0.042260    ## ln(shear stress)
  b6 <- -0.118155    ## ln(trawling):ln(shear stress)
  b7 <- -0.019554    ## ln(longevity):Gravel

  Mud <- Region@data$Mud
  Gravel <- Region@data$Gravel
  Stress <- log(Region@data$Shearstress+0.01)
  Trawling <- log(rep(0.01,nrow(Region@data)))

  medlong<-exp((logit(0.5)-b0-b3*Mud-b4*Gravel-b5*Stress-b6*Trawling*Stress-b2*Trawling)/(b1+b7*Gravel))
  medlong[medlong<1]  <- 1
  Region@data$medlong<-medlong

  ll<-log(1)
  Region@data$Lone<-expit(b0+b7*Gravel*ll+b3*Mud+b4*Gravel+b5*Stress+b6*Trawling*Stress+b2*Trawling+b1*ll)
  ll<-log(3)
  Region@data$Lthree<-expit(b0+b7*Gravel*ll+b3*Mud+b4*Gravel+b5*Stress+b6*Trawling*Stress+b2*Trawling+b1*ll)
  ll<-log(10)
  Region@data$Lten<-expit(b0+b7*Gravel*ll+b3*Mud+b4*Gravel+b5*Stress+b6*Trawling*Stress+b2*Trawling+b1*ll)

  rm(list = c('b0','b1','b2','b3','b4','b5','b6','b7','Mud','Gravel','Stress','Trawling','ll','medlong'))
} # T/F  


############# TRASH
############# TRASH
############# TRASH
############# TRASH



######################################################
######################################################
######################################################
###### PREPARATION DATA JOCHEN(ILVO) 
# Read BTS data ####
if(F){
dirpath <- "C:/Users/jdepestele/OneDrive - ILVO/gitr/follow_others/ICES/DATRAS_BTS"
input <- paste(dirpath,"input", sep="/")

setwd(dirpath)

bent <- readRDS(paste0(input,"/HH_and_HL_combined_BTSdata_IrishSea7a_CelticSea7fg.RDS"))

bent <- bent[which(bent$phylum %in% c("Arthropoda","Mollusca","Echinodermata","Annelida",
                                      "Cnidaria","Porifera","Platyhelminthes","Sipuncula",
                                      "Priapulida","Nemertea","Acanthocephala")),]
bent <- bent[which(!bent$class %in% c("Cephalopoda")),]
bent <- bent[!is.na(bent$Valid_Aphia),]
# SpatialPointsDataframe
# coords <- bent[,c("HaulLong","HaulLat")]
# crs = CRS("+proj=longlat +datum=WGS84")
# sbent <- SpatialPointsDataFrame(coords      = coords,
#                                 data        = bent,
#                                 proj4string = crs)
# class(sbent)

# # Make a map using sf
# coordinates(bent) <- ~ShootLong+ShootLat
# proj4string(bent) <- CRS("+proj=longlat +datum=WGS84")
# 
# bent <- st_as_sf(bent,coords = 1:2)
# bent <- bent %>% 
#   st_transform(crs = 4326)
# bent_map <- get_stamenmap(
#   bbox = unname(st_bbox(bent)),
#   zoom = 5, maptype = 'toner-lite', source = 'stamen'
# ) %>% ggmap()
# 
# bent_map + geom_sf(data = bent, 
#                       aes(col = factor(Year)),alpha = 0.25, show.legend = 'point', inherit.aes = F)
# 

# Read the trait matrix from Pascal ####
load(file="//clo.be/Home/Home_d1/jdepestele/Documents/000_data/biota/traits_benth/WGFBIT_Celtic_S2_BIOLOGICAL_BTA.RData")
# Traits Matrix (a lot of traits)
head(BIOLOGICAL_BTA$MEGAFAUNA$MATRIX.FUNCnew)
# Longevity data only
BIOLOGICAL_BTA$MEGAFAUNA$MATRIX.FUNCnew[,c("LO1","LO2","LO3","LO4")]
# You will also find all taxonomic informations (extracted from worms)
BIOLOGICAL_BTA$MEGAFAUNA$MATRIX.FUNCtaxa

Ldat <- setDT(BIOLOGICAL_BTA$MEGAFAUNA$MATRIX.FUNCnew[,c("LO1","LO2","LO3","LO4")])
Ldat$valid_name <- row.names(BIOLOGICAL_BTA$MEGAFAUNA$MATRIX.FUNCnew[,c("LO1","LO2","LO3","LO4")])
setkey(Ldat, valid_name)

spdat <- setDT(BIOLOGICAL_BTA$MEGAFAUNA$MATRIX.FUNCtaxa)
setkey(spdat, valid_name)

Ldat <- spdat[Ldat]

# # Benthic traits from the ICES tutorial:
# LdatNS <- read.csv("//clo.be/Home/Home_d1/jdepestele/Documents/000_data/biota/traits_benth/Benthic_Data_tutorial.csv",header=T,sep=";")
# LdatNS$genus <- sapply(strsplit(LdatNS$Nomen," "), "[", 1)
# unique(LdatNS$genus)


# Merge longevity traits with BTS taxa ####
setDT(bent)
setkey(bent, Valid_Aphia)
setkey(Ldat, AphiaID)

nrow(bent[!Valid_Aphia %in% Ldat$AphiaID])/nrow(bent) # 11% lost data
# unique(bent[!bent$Valid_Aphia %in% Ldat$AphiaID]$phylum)
# unique(bent[!bent$Valid_Aphia %in% Ldat$AphiaID]$family)
unique(bent[!bent$Valid_Aphia %in% Ldat$AphiaID]$genus)

Ldat  <- Ldat[,c("AphiaID","LO1","LO2","LO3","LO4")]
Lbent <- Ldat[bent]

names(Lbent)


summary(Lbent$HLNoAtLngt)
summary(Lbent$TotalNo)
summary(Lbent$CatCatchWgt)

Lbent[,c("HLNoAtLngt","TotalNo","CatCatchWgt")] <- lapply(Lbent[,c("HLNoAtLngt","TotalNo","CatCatchWgt")], 
                                                        as.numeric)

Lbent_wgt <- Lbent[,c("StNo","HaulNo","Year","Area_27","HaulLat","HaulLong","AphiaID","genus","scientificname","sweptarea",
                  "CatCatchWgt","SubWgt",
                  "LngtClass","HLNoAtLngt","TotalNo","LO1","LO2","LO3","LO4")]
Lbent_wgt <- Lbent_wgt %>%
  group_by(StNo, HaulNo, Year, Area_27,
           HaulLat,HaulLong,scientificname, sweptarea,LO1,LO2,LO3,LO4) %>%
  dplyr::summarize(CatCatchWgt_kg=mean(CatCatchWgt,na.rm=T)/1000,
                   HLNoAtLngt=sum(HLNoAtLngt,na.rm=T),
                   TotalNo=mean(TotalNo,na.rm=T))


# STANDARDISE PER SWEPT AREA
Lbent_wgt$kg_per_sqm <- Lbent_wgt$CatCatchWgt_kg / Lbent_wgt$sweptarea

# remove NAs
Lbent_wgt = Lbent_wgt[!is.na(Lbent_wgt$kg_per_sqm),]

nrow(Lbent_wgt)

summary(Lbent_wgt$HaulLat)
summary(Lbent_wgt$HaulLong)
Lbent_wgt <- Lbent_wgt[Lbent_wgt$HaulLong>-8 & Lbent_wgt$HaulLong< -2,]

DT <- setDT(Lbent_wgt)
DT <- DT[, .(L1=mean(LO1,na.rm=T),L1_3=mean(LO2,na.rm=T),L3_10=mean(LO3,na.rm=T),L10=mean(LO4,na.rm=T),kg_per_sqm=mean(kg_per_sqm,na.rm=T)),
   by=.(StNo,HaulNo,Area_27,HaulLat,HaulLong)]  #Year,
Lbent_wgt <- DT

# Make a map using sf
sfLbent_wgt <- as.data.frame(Lbent_wgt)
coordinates(sfLbent_wgt) <- ~HaulLong+HaulLat
proj4string(sfLbent_wgt) <- CRS("+proj=longlat +datum=WGS84")

sfLbent_wgt <- st_as_sf(sfLbent_wgt,coords = 1:2)
sfLbent_wgt <- sfLbent_wgt %>%
  st_transform(crs = 4326)

bent_map <- get_stamenmap(
  bbox = c(-9, 49.5, -2, 55.5),# unname(st_bbox(sfLbent_wgt)),
  zoom = 7, maptype = 'terrain', source = 'stamen'
) %>% ggmap()

sfLbent_wgt$L1cat <- cut(sfLbent_wgt$L1, breaks = c(-1e-10,0.25,0.5,0.75,1))
sfLbent_wgt$L1_3cat <- cut(sfLbent_wgt$L1_3, breaks = c(-1e-10,0.25,0.5,0.75,1))
sfLbent_wgt$L3_10cat <- cut(sfLbent_wgt$L3_10, breaks = c(-1e-10,0.25,0.5,0.75,1))
sfLbent_wgt$L10cat <- cut(sfLbent_wgt$L10, breaks = c(-1e-10,0.25,0.5,0.75,1))
# sfLbent_wgt$L10cat <- cut(sfLbent_wgt$L10, breaks = unique(quantile(sfLbent_wgt$L10,probs=seq.int(0,1, by=1/5),na.rm=T)))
bent_map + geom_sf(data = sfLbent_wgt[!is.na(sfLbent_wgt$L1),],
                      aes(col = as.factor(L1cat), size = L1),
                   alpha = 0.25, show.legend = 'point', inherit.aes = F)
bent_map + geom_sf(data = sfLbent_wgt[!is.na(sfLbent_wgt$L1_3),],
                   aes(col = as.factor(L1_3cat), size = L1_3),
                   alpha = 0.25, show.legend = 'point', inherit.aes = F)
bent_map + geom_sf(data = sfLbent_wgt[!is.na(sfLbent_wgt$L3_10),],
                   aes(col = as.factor(L3_10cat), size = L3_10),
                   alpha = 0.25, show.legend = 'point', inherit.aes = F)
bent_map + geom_sf(data = sfLbent_wgt,
                   aes(col = as.factor(L10cat), size = L10),
                   alpha = 0.25, show.legend = 'point', inherit.aes = F)
bent_map + geom_sf(data = sfLbent_wgt,
                   aes(col = kg_per_sqm, size = L10),
                   alpha = 0.25, show.legend = 'point', inherit.aes = F)



# open benthic data
Lbent_wgt <- as.data.frame(Lbent_wgt)
Lbent_wgt[,c("L1","L1_3","L3_10","L10")] <- Lbent_wgt$kg_per_sqm * Lbent_wgt[,c("L1","L1_3","L3_10","L10")] # multiply biomass with fuzzy coded trait data
# Lbent_wgt$sample_ID <- paste(Lbent_wgt$Year, Lbent_wgt$StNo, Lbent_wgt$HaulNo, sep="__")
Lbent_wgt$sample_ID <- paste(Lbent_wgt$StNo, Lbent_wgt$HaulNo,Lbent_wgt$HaulLat,Lbent_wgt$HaulLong, sep="__")

# summarize benthic data per sample ID and calculate the fraction
namesCol <- c("kg_per_sqm","L1","L1_3","L3_10","L10")
statdat <- aggregate(Lbent_wgt[, namesCol], by= list(Lbent_wgt$sample_ID), FUN=function(x){sum(x, na.rm=T)})
names(statdat) <- c("ID","Biomass","L1","L1_3","L3_10","L10")


# now link to environmental conditions

LAEA <- "+proj=longlat +datum=WGS84"

Hauls <-  unique(Lbent_wgt[,c("sample_ID", "HaulLong", "HaulLat")])

Coords <- as.data.frame(rgdal::project(cbind(Hauls$HaulLong, Hauls$HaulLat), LAEA))
names(Coords) <- c("Lon_proj", "Lat_proj")
LonData <- cbind.data.frame(statdat, Coords)

LonData_SPDF <- LonData
coordinates(LonData_SPDF) <- ~Lon_proj+Lat_proj
spplot(LonData_SPDF[,"L10"]) #, at=c(seq(0,10,0.5),seq(11,51,2))


#We load Env variables

fondo_stack <- stack("//clo.be/Home/Home_d1/jdepestele/Documents/4_ICES/2021_FBIT/From_José/TraitAnalysis/AllEnvLayers.tif")
names(fondo_stack) <- c("Depth", "Chl", "Temp", "Currents", "Waves", "AllEnergy", "Substrate")
plot(fondo_stack)


e <- as(extent(-8, -2, 49, 55.5), 'SpatialPolygons')
crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"
r <- projectRaster(fondo_stack, crs=LAEA)
NWWenv <- crop(r, e)
plot(NWWenv)


# now link to environmental conditions
#we extract values from raster
EnvVarHauls <- as.data.frame(raster::extract(NWWenv, Coords))

DatabyHaul_Selected <- cbind.data.frame(LonData, EnvVarHauls)
DatabyHaul_Selected <- na.omit(DatabyHaul_Selected)
head(DatabyHaul_Selected)


###
# TRAWLING INTENSITY IS NOT INCLUDE WHILE IT SHOULD!
pathdir <- "C:/Users/jdepestele/OneDrive - ILVO/gitr/fbit_benthis_nationaal_vv"
input <- paste(pathdir,"/0_input/",sep="")
sar <- readRDS(paste0(input,"SARsurf_EUmet_20092017.RDS"))
head(sar)
library(sfdSAR)
sar$lat <- csquare_lat(sar$squares)
sar$lon <- csquare_lon(sar$squares)
LonLatProj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'"
TemplateRast <- raster(xmn=-26,xmx=12,ymn=32,ymx=70, res=c(0.05,0.05), crs=LonLatProj)
sar$depletion <- sar$depl_rate * sar$SARsurf
SpData <- sar[,c("year","lat","lon","depletion")]
coordinates(SpData) <- c("lon", "lat")

Year_vect <- unique(SpData$year)
Myrast <- list()

for (i in 1:length(Year_vect)){
  DataByYear <- SpData[SpData$year==Year_vect[i],]
  Myrast[[i]] <- rasterize(DataByYear, TemplateRast, field="depletion", method="ngb")}

AllYears <- do.call(stack, Myrast)
names(AllYears) <- as.character(Year_vect)
plot(AllYears)
plot(fondo_stack) # different axes from AllYears !!

AllYears_proj <- projectRaster(AllYears, fondo_stack$Depth, method="ngb")
AllYears_proj <- mask(AllYears_proj, fondo_stack$Depth)
plot(AllYears_proj)

#Mean study period
meanF <- stackApply(AllYears_proj, indices =  rep(1,nlayers(AllYears_proj)), fun = "mean", na.rm = T)
maxF <- stackApply(AllYears_proj, indices =  rep(1,nlayers(AllYears_proj)), fun = "max", na.rm = T)

meanF_proj <- projectRaster(meanF, crs=LAEA)
plot(meanF_proj)

DatabyHaul_Selected$TrawlingEffort <- raster::extract(meanF_proj, cbind(DatabyHaul_Selected$Lon_proj,DatabyHaul_Selected$Lat_proj))
summary(DatabyHaul_Selected$TrawlingEffort)
DatabyHaul_Selected[is.na(DatabyHaul_Selected$TrawlingEffort),]$TrawlingEffort <- 0

###

e <- as(extent(-8, -2, 49, 55.5), 'SpatialPolygons')
crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"
r2 <- projectRaster(meanF_proj, crs=LAEA)
NWWsar <- crop(r2, e)
names(NWWsar) <- "Depletion"
plot(NWWsar)

nlayers(NWWenv)
plot(NWWenv)

NWWenvd <- stack(NWWenv, NWWsar)
nlayers(NWWenvd)
plot(NWWenvd)
}
######################
######################
######################
######################
######################




# With trawling:
if(F){
mod1   <-  glmer(Cumb ~ ll + Depth + Chl +Temp+AllEnergy +Substrate + TrawlingEffort + (1 | Uniq), data=fulldat, family=binomial)
mod2   <-  glmer(Cumb ~ ll + Depth + Chl +Temp+AllEnergy +Substrate  + (1 | Uniq), data=fulldat, family=binomial)
mod3   <-  glmer(Cumb ~ ll + Depth + Chl +Temp+AllEnergy + TrawlingEffort + (1 | Uniq), data=fulldat, family=binomial)
mod4   <-  glmer(Cumb ~ ll + Depth + Chl +Temp +Substrate + TrawlingEffort + (1 | Uniq), data=fulldat, family=binomial)
mod5   <-  glmer(Cumb ~ ll + Depth + Chl +AllEnergy +Substrate + TrawlingEffort + (1 | Uniq), data=fulldat, family=binomial)
mod6   <-  glmer(Cumb ~ ll + Depth + Temp+AllEnergy +Substrate + TrawlingEffort + (1 | Uniq), data=fulldat, family=binomial)
mod7   <-  glmer(Cumb ~ ll + Chl +Temp+AllEnergy +Substrate + TrawlingEffort + (1 | Uniq), data=fulldat, family=binomial)
mod8   <-  glmer(Cumb ~ ll + Depth + Chl +Temp+AllEnergy  + (1 | Uniq), data=fulldat, family=binomial)
mod9   <-  glmer(Cumb ~ ll + Depth + Chl +Temp + TrawlingEffort + (1 | Uniq), data=fulldat, family=binomial)
mod10   <-  glmer(Cumb ~ ll + Depth + Chl + AllEnergy + TrawlingEffort + (1 | Uniq), data=fulldat, family=binomial)
mod11   <-  glmer(Cumb ~ ll + Depth + Temp  + AllEnergy + TrawlingEffort + (1 | Uniq), data=fulldat, family=binomial)
mod12   <-  glmer(Cumb ~ ll + Chl + Temp  + AllEnergy + TrawlingEffort + (1 | Uniq), data=fulldat, family=binomial)

AICdf <- AIC(mod1,mod2,mod3, mod4,mod5,mod6, mod7, mod8, mod9, mod10, mod11, mod12)
AICdf[order(AICdf[,2]),]
AICdf[which(AICdf[,2]==min(AICdf[,2])),]

# models give a singular fit --> the random effect is very small (but you can argue that it, in principle, has to be included)

modcoeff  <-  fixef(mod11)

# models give a singular fit --> the random effect is very small (but you can argue that it, in principle, has to be included)
modcoeff_mod8  <-  fixef(mod8)
modcoeff_mod8 
# model 8
# (Intercept)            ll         Depth           Chl          Temp     AllEnergy 
# -1.958235e+01  1.057818e+01 -1.702439e-02 -1.657814e-02  7.409814e-02  5.530291e-04 
coef_int_mod8 <- modcoeff[1]
coef_ll_mod8  <- modcoeff[2]
coef_Depth_mod8 <- modcoeff[3]
coef_Chl_mod8 <- modcoeff[4]
coef_Temp_mod8 <- modcoeff[5]
coef_AllEnergy_mod8 <- modcoeff[6]

var_SPDF_mod8 <- NWWenvd
var_SPDF_mod8$medLong <- exp((qlogis(0.5)- coef_int_mod8 - coef_Depth_mod8*var_SPDF_mod8$Depth - coef_Temp_mod8*var_SPDF_mod8$Temp
                         - coef_AllEnergy_mod8*var_SPDF_mod8$AllEnergy - coef_Chl_mod8*var_SPDF_mod8$Chl)/ coef_ll_mod8)
windows()
plot(var_SPDF_mod8)
plot(var_SPDF_mod8[["medLong"]])

summary(var_SPDF_mod8[["medLong"]])
cuts=c(0, 4, 4.5, 5, 5.5, 6, 10) #set breaks
pal <- colorRampPalette(c("white","black"))
blueorange <- c("#feb24c","#fed976","#ffffcc","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494")
pal <- colorRampPalette(c("green","red"))
plot(var_SPDF_mod8[["medLong"]],breaks=cuts, col = blueorange)
# plot(var_SPDF[["medLong"]],breaks=cuts, col = pal(length(cuts)))


# model 11
modcoeff_mod11  <-  fixef(mod11)
modcoeff_mod11 
# model 11
# (Intercept)             ll          Depth           Temp      AllEnergy TrawlingEffort 
# -1.955549e+01   1.061971e+01  -2.229183e-02   3.985256e-02   5.545446e-04   7.317996e-02 

coef_int_mod11 <- modcoeff[1]
coef_ll_mod11  <- modcoeff[2]
coef_Depth_mod11 <- modcoeff[3]
coef_Temp_mod11 <- modcoeff[4]
coef_AllEnergy_mod11 <- modcoeff[5]
coef_Depl_mod11 <- modcoeff[6]

var_SPDF_mod11 <- NWWenvd
var_SPDF_mod11$medLong <- exp((qlogis(0.5)- coef_int_mod11 - coef_Depth_mod11*var_SPDF_mod11$Depth - 
                                 coef_Temp_mod11*var_SPDF_mod11$Temp - coef_AllEnergy_mod11*var_SPDF_mod11$AllEnergy - 
                                 coef_Depl_mod11*var_SPDF_mod11$Depletion)/ coef_ll_mod11)
windows()
plot(var_SPDF_mod11)
plot(var_SPDF_mod11[["medLong"]])

summary(var_SPDF_mod11[["medLong"]])
cuts=c(3, 4, 4.5, 5, 5.5, 6.5) #set breaks
pal <- colorRampPalette(c("white","black"))
blueorange <- c("#feb24c","#fed976","#ffffcc","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494")
pal <- colorRampPalette(c("green","red"))
plot(var_SPDF_mod11[["medLong"]],breaks=cuts, col = blueorange)
# plot(var_SPDF[["medLong"]],breaks=cuts, col = pal(length(cuts)))


        # NWWenvsar <- fulldat
        # var_SPDF <- NWWenvsar
        # 
        # coordinates(var_SPDF) <- ~Lon_proj + Lat_proj
        # 
        # # Get ready for mapping
        # coords <- NWWenvsar[,c("Lon_proj","Lat_proj")]
        # crs = CRS("+proj=longlat +datum=WGS84")
        # 
        # NWWenvsar$medLong <- exp((qlogis(0.5)- coef_int - coef_Depth*NWWenvsar$Depth - coef_Temp*NWWenvsar$Temp
        #                          - coef_SAR*NWWenvsar$TrawlingEffort - coef_AllE*NWWenvsar$AllEnergy)/ coef_ll)
        # summary(NWWenvsar$medLong)
        # NWWenvsar$medLongcat <- cut(NWWenvsar$medLong, 
        #                             breaks = unique(quantile(NWWenvsar$medLong,probs=seq.int(0,1, by=1/5),na.rm=T)))
        # 
        # coordinates(NWWenvsar) <- ~Lon_proj+Lat_proj
        # proj4string(NWWenvsar) <- CRS("+proj=longlat +datum=WGS84")
        # 
        # sfNWWenvsar <- st_as_sf(NWWenvsar,coords = 1:2)
        # sfNWWenvsar <- sfNWWenvsar %>%
        #   st_transform(crs = 4326)
        # sfNWWenvsar_map <- get_stamenmap(
        #   bbox = c(-9, 49.5, -2, 55.5),
        #   zoom = 5, maptype = 'terrain', source = 'stamen'
        # ) %>% ggmap()
        # 
        # windows()
        # sfNWWenvsar_map + geom_sf(data = sfNWWenvsar,
        #                       aes(col = medLongcat),alpha = 0.25, show.legend = 'point', inherit.aes = F)



sfLbent_wgt$L10cat <- cut(sfLbent_wgt$L10, breaks = unique(quantile(sfLbent_wgt$L10,probs=seq.int(0,1, by=1/5),na.rm=T)))
windows()
bent_map + geom_sf(data = sfLbent_wgt,
                   aes(col = as.factor(L10cat), size = L10),
                   alpha = 0.25, show.legend = 'point', inherit.aes = F)

sfLbent_wgt$L3_10cat <- cut(sfLbent_wgt$L3_10, breaks = unique(quantile(sfLbent_wgt$L3_10,probs=seq.int(0,1, by=1/5),na.rm=T)))
windows()
bent_map + geom_sf(data = sfLbent_wgt,
                   aes(col = as.factor(L3_10cat), size = L3_10),
                   alpha = 0.25, show.legend = 'point', inherit.aes = F)

sfLbent_wgt$L1_3cat <- cut(sfLbent_wgt$L1_3, breaks = unique(quantile(sfLbent_wgt$L1_3,probs=seq.int(0,1, by=1/5),na.rm=T)))
windows()
bent_map + geom_sf(data = sfLbent_wgt,
                   aes(col = as.factor(L1_3cat), size = L1_3),
                   alpha = 0.25, show.legend = 'point', inherit.aes = F)

sfLbent_wgt$L1cat <- cut(sfLbent_wgt$L1, breaks = unique(quantile(sfLbent_wgt$L1,probs=seq.int(0,1, by=1/5),na.rm=T)))
windows()
bent_map + geom_sf(data = sfLbent_wgt,
                   aes(col = as.factor(L1cat), size = L1),
                   alpha = 0.25, show.legend = 'point', inherit.aes = F)
                   
                   }
                 }