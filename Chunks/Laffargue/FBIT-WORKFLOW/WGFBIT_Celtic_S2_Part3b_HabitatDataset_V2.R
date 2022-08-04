#####################################
######### WGFBIT - October 2019 ##########
##################################### >)))°>
######## SCRIPT dataset compilation - Celtic/BoB Eco-region (WGFBIT)
##
### R version 3.5.1
###############
####


######## TO MODIFY
### Utilize Jose extraction script and envirmnemental dataset format (RASTER)
#######


#################################
### SUMMARY
#################################


#################
### Load files ##
#################

## Formatted biological dataset
load(file=file.path(DIR.OUPUT[[paste(1)]],paste(SEL.ScriptName,"BIOdataset.RData",sep=""))) 

##

### HABITAT DATASET
## Load previously created file
if(file.exists(file.path(DIR.OUPUT[[paste(1)]],paste(SEL.ScriptName,"HABITATdataset.RData",sep="")))){
load(file=file.path(DIR.OUPUT[[paste(1)]],paste(SEL.ScriptName,"HABITATdataset.RData",sep="")))}else{
HABITATdataset <- list()}


######################################################
### LOAD STANDARD C-SQUARE GRIDS
GRIDS <- list()
TEMP.FILES <- list.files(DIR.SOURCE[[paste("5b")]])
TEMP.FILES <- TEMP.FILES[grep("grid_sensitivity",TEMP.FILES)]
for(FILE in TEMP.FILES){
GRIDS[[paste(FILE)]] <- loadRData(file=file.path(DIR.SOURCE[[paste("5b")]],FILE))
rm(FILE)}

GRIDS[["WHOLE"]] <- rbind(GRIDS[[1]],GRIDS[[2]])
###
###


#################
#################
#################


#################################
### SUMMARY
#################################

### 1.2 - Data selection sources: Files names

# 1.2.1 - EUNIS/MSFD
SOURCEFILE.EUNIS <-  "EUSeaMap2016_WesternWaters.shp"
#SOURCEFILE.EUNIS <-  "EuSeaMap2019_BoB_CS.shp"

# 1.2.2 - Depth
SOURCEFILE.DEPTH <- "GebcoDATA_2014_RASTER_TOTAL.RData"

# 1.2.3 - Fishing
REFname.FISHING <-  "OSPAR_intensity_total"
SOURCEFILE.FISHING <- list.files(DIR.SOURCE[[paste(4)]])[grep(REFname.FISHING,list.files(DIR.SOURCE[[paste(4)]]))]


      ### PARAMETER

# Definition of standard geographic coordinates format      
SEL.STANDARDproj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
     
      ### 2.3.2 - Region selection
      
      ## 2.3.2.1 - Define polygon boundaries for subregion
      

gt <- (GridTopology(c(-16, 36), c(0.05, 0.05), c(360, 560))) # c(long, lat), c(cellsize long, lat), c(nb of grids long, lat)
 
SELECT.BOUNDARIES <- list()

     ### Boundaries to check /define
SELECT.BOUNDARIES[["Celtic"]] <- c("LongMin"=-13,"LongMax"=-1,"LatMin"=48,"LatMax"=52)
SELECT.BOUNDARIES[["BoB"]] <- c("LongMin"=-10,"LongMax"=-1,"LatMin"=43,"LatMax"=48)


## CREATION of EXTENDED BOUNDARIES FROM WHOLE BIOLOGICAL DATASET
TEMP.LON <- NULL
TEMP.LAT <- NULL
for(SURVEY in names(BIOdataset[["FORMATED"]])){
print(paste("Match stations vs habitats types for",SURVEY,"dataset"))                         
        # 2.3.4.1 - Format survey coordinates
TEMP.DATA <- as.data.frame(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]])
TEMP.DATA$x <- as.numeric(as.vector(TEMP.DATA[["Longitude"]]))
TEMP.DATA$y <- as.numeric(as.vector(TEMP.DATA[["Latitude"]]))

TEMP.LON <- c(min(c(TEMP.LON,TEMP.DATA$x)),max(c(TEMP.LON,TEMP.DATA$x)))
TEMP.LAT <- c(min(c(TEMP.LAT,TEMP.DATA$y)),max(c(TEMP.LAT,TEMP.DATA$y)))}

coordinates(TEMP.DATA) <- ~ x+y     
proj4string(TEMP.DATA) <- SEL.STANDARDproj


SELECT.BOUNDARIES[["WHOLE_Extended"]] <- c("LongMin"=floor(min(TEMP.LON)*10)/10,"LongMax"=ceiling(max(TEMP.LON)*10)/10,"LatMin"=floor(min(TEMP.LAT)*10)/10,"LatMax"=ceiling(max(TEMP.LAT)*10)/10)
SELECT.BOUNDARIES[["Celtic_Extended"]] <- c("LongMin"=-16,"LongMax"=-1,"LatMin"=48,"LatMax"=ceiling(max(TEMP.LAT)*10)/10)
SELECT.BOUNDARIES[["BoB_Extended"]] <- c("LongMin"=-16,"LongMax"=-1,"LatMin"=floor(min(TEMP.LAT)*10)/10,"LatMax"=48)

rm(TEMP.LON,TEMP.LAT)
     ###




#### TEMPORARY ENVIRONMENT FILES FROM JOSE SCRIPT#####

## Load raster file and attribution to each Surveys stations
SEL.RASTER <- T # "raster" : raster stack from Jose Script / "raw" : .shp from previous script (loading HABITATdataset)
if(SEL.RASTER){                                                                                                       
RASTERstack <- stack(file.path(DIR.OUPUT[[paste(1)]],paste("AllEnvLayers.tif",sep=""))) # Load raster stack from Jose Script
names(RASTERstack) <- c("Depth", "Chl","Temp", "Currents", "waves","ALLenergy","Substrate" )


for(SURVEY in names(BIOdataset[["FORMATED"]])){
TEMP.DATA <- as.data.frame(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]])
TEMP.DATA$x <- as.numeric(as.vector(TEMP.DATA[["Longitude"]]))
TEMP.DATA$y <- as.numeric(as.vector(TEMP.DATA[["Latitude"]]))
coordinates(TEMP.DATA) <- ~ x+y     
proj4string(TEMP.DATA) <- SEL.STANDARDproj

for(VAR in  names(RASTERstack)){
print(paste("Raster", VAR, "to", SURVEY))
TEMP <- RASTERstack[[VAR]]
TEMP2 <- extract(TEMP, TEMP.DATA)
BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]][[paste("Raster",VAR,sep="_")]] <- TEMP2
rm(TEMP,TEMP2)}
rm(SURVEY)}
}

## Save New dataset
 # Save R object
if(F){
TEMP.DIR <- file.path(DIR.OUPUT[[paste(1)]])
if(file.exists(TEMP.DIR)==F){dir.create(TEMP.DIR)}
save(BIOdataset,file=file.path(TEMP.DIR,paste(SEL.ScriptName,"BIOdataset.RData",sep="")))
save(RASTERstack,file=file.path(TEMP.DIR,paste(SEL.ScriptName,"RASTERstack.RData",sep="")))
rm(TEMP.DIR)}
  #

##


####




  


######################################
   ### 2.3 - HABITAT DATASET "EUNIS" (MSFD broad habitat types)
######################################

#EUNIS.POL                   <-  rgdal::readOGR(dsn=file.path(DIR.SOURCE[[paste(6)]],SOURCEFILE.EUNIS)) #,proj4string=SEL.EUNISproj)           rgdal::readOGR

EUNIS.POL <- shapefile(file.path(DIR.SOURCE[[paste(6)]],SOURCEFILE.EUNIS)) 
EUNIS.POL <- gBuffer(EUNIS.POL, byid=TRUE, width=0)

LIST.EUNIScode <- unique(EUNIS.POL@data[,c("HAB_TYPE","Biozone","Substrate")])
LIST.EUNIScode$NAME <- paste(LIST.EUNIScode[["Biozone"]],LIST.EUNIScode[["Substrate"]])

SEL.EUNISproj <- proj4string(EUNIS.POL) # CRS("+proj=longlat +ellps=WGS84") # EUNIS CRS

HABITATdataset[["EUNISproj"]] <- SEL.EUNISproj

      ## 2.3.1 - FILTER / MODIFY EUNIS HABITATS TYPES
       ### WARNING TRANSFORMED MSFD SEDIMENT TYPES
if(T){
EUNIS.POL@data[,"MSFD_BH17"][EUNIS.POL@data[,"MSFD_BH17"]=="Upper bathyal sediment or Upper bathyal rock and biogenic reef"] <- "Upper bathyal sediment"}

###
##
           ## Harmonize habitat types names
     SEL.HABtypeSource <- c("HAB_TYPE","MSFD_BH17")
     SEL.HABtype <- c("EUNIS","MSFD")
     
     colnames(EUNIS.POL@data)[match(SEL.HABtypeSource,colnames(EUNIS.POL@data))] <- SEL.HABtype
           ##

      # ALL subregions
TEMP.REGIONtot <- NULL
for(BOUND in names(SELECT.BOUNDARIES)){
TEMP.REGIONtot[["LONG"]] <- c(TEMP.REGIONtot[["LONG"]],as.vector(SELECT.BOUNDARIES[[BOUND]][c("LongMin","LongMax")]))
TEMP.REGIONtot[["LAT"]] <- c(TEMP.REGIONtot[["LAT"]],as.vector(SELECT.BOUNDARIES[[BOUND]][c("LatMin","LatMax")]))
rm(BOUND)}
TEMP.REGIONtot[["LONG"]] <- range(TEMP.REGIONtot[["LONG"]])
TEMP.REGIONtot[["LAT"]] <- range(TEMP.REGIONtot[["LAT"]])
SELECT.BOUNDARIES[[paste(names(SELECT.BOUNDARIES),collapse=".")]] <- c("LongMin"=min(TEMP.REGIONtot[["LONG"]]),"LongMax"=max(TEMP.REGIONtot[["LONG"]]),"LatMin"=min(TEMP.REGIONtot[["LAT"]]),"LatMax"=max(TEMP.REGIONtot[["LAT"]]))
rm(TEMP.REGIONtot)

HABITATdataset[["RegionLim"]][["BoundariesDefinition"]] <- SELECT.BOUNDARIES

       
      #  2.3.2.3 -  Intersection of Eunis shapefiles with polygons boundaries
      
      # Subregion Polygon

HABITATdataset[["SubregionPolygon"]] <- list()      
HABITATdataset[["SubregionEUNIS"]] <- list()


### RUN FIRST INTERSECTION FOR THE MAXIMAL SPATIAL EXTENT !!!!
## Intersect next extent from this first one
LongMin=NULL;LatMin=NULL;LongMax=NULL;LatMax=NULL
for(SEL.BOUND in names(HABITATdataset[["RegionLim"]][["BoundariesDefinition"]]))
{LongMin=min(c(LongMin,SELECT.BOUNDARIES[[SEL.BOUND]]["LongMin"]))
       LatMin=min(c(LatMin,SELECT.BOUNDARIES[[SEL.BOUND]]["LatMin"]))
       LongMax=max(c(LongMax,SELECT.BOUNDARIES[[SEL.BOUND]]["LongMax"]))
       LatMax=max(c(LatMax,SELECT.BOUNDARIES[[SEL.BOUND]]["LatMax"]))}
# Create crop polygon
if(T){
Pol1=rbind(c(LongMin,LatMin),c(LongMin,LatMax),c(LongMax,LatMax),c(LongMax,LatMin),c(LongMin,LatMin))
Pols1=Polygons(list(Polygon(Pol1)),"Pols1")
TEMP.BoundariesMAX=SpatialPolygons(list(Pols1))
proj4string(TEMP.BoundariesMAX) <- proj4string(EUNIS.POL)
HABITATdataset[["SubregionPolygon"]][["WHOLE"]] <- TEMP.BoundariesMAX}
rm(LongMin, LongMax, LatMin, LatMax)


temp.intersect <- gIntersects(EUNIS.POL, HABITATdataset[["SubregionPolygon"]][["WHOLE"]],byid=T)
EUNIS.POLselectRegion <- gIntersection(EUNIS.POL, HABITATdataset[["SubregionPolygon"]][["WHOLE"]],byid=T)
TEMP.POLlist <- sub("Pols1 ","",names(EUNIS.POLselectRegion))
TEMP.DATAslot <- data.frame(EUNIS.POL@data[temp.intersect[1,],])

TEMP.PolyIDs <- sapply(slot(EUNIS.POLselectRegion, "polygons"), function(x) slot(x, "ID"))
rownames(TEMP.DATAslot) <- TEMP.PolyIDs
EUNIS.POLselectRegion <- SpatialPolygonsDataFrame(EUNIS.POLselectRegion,TEMP.DATAslot)

HABITATdataset[["SubregionEUNIS"]][["WHOLE"]] <- EUNIS.POLselectRegion

#gArea(sp3)


 
## TO COMPLETE

###

### Delimit each habitat datasets for each selected bound
##

for(SEL.BOUND in names(HABITATdataset[["RegionLim"]][["BoundariesDefinition"]]))
{print(paste("EUNIS map selection for subregion",SEL.BOUND))
       LongMin=SELECT.BOUNDARIES[[SEL.BOUND]]["LongMin"]
       LatMin=SELECT.BOUNDARIES[[SEL.BOUND]]["LatMin"]
       LongMax=SELECT.BOUNDARIES[[SEL.BOUND]]["LongMax"]
       LatMax=SELECT.BOUNDARIES[[SEL.BOUND]]["LatMax"]

## TO CHECK / COMPLETE
# spTransform(CRS(prj)) %>%

## Create crop polygon
if(T){
Pol1=rbind(c(LongMin,LatMin),c(LongMin,LatMax),c(LongMax,LatMax),c(LongMax,LatMin),c(LongMin,LatMin))
Pols1=Polygons(list(Polygon(Pol1)),"Pols1")
TEMP.Boundaries2=SpatialPolygons(list(Pols1))
proj4string(TEMP.Boundaries2) <- proj4string(HABITATdataset[["SubregionEUNIS"]][["WHOLE"]])

HABITATdataset[["SubregionPolygon"]][[SEL.BOUND]] <- TEMP.Boundaries2}

#spTransform(CRS(prj))

temp.intersect <- gIntersects(HABITATdataset[["SubregionEUNIS"]][["WHOLE"]], HABITATdataset[["SubregionPolygon"]][[SEL.BOUND]],byid=T)
EUNIS.POLselectRegion <- gIntersection(HABITATdataset[["SubregionEUNIS"]][["WHOLE"]], HABITATdataset[["SubregionPolygon"]][[SEL.BOUND]],byid=T)
TEMP.POLlist <- sub("Pols1 ","",names(EUNIS.POLselectRegion))
TEMP.DATAslot <- data.frame(HABITATdataset[["SubregionEUNIS"]][["WHOLE"]]@data[as.vector(temp.intersect),])

TEMP.PolyIDs <- sapply(slot(EUNIS.POLselectRegion, "polygons"), function(x) slot(x, "ID"))
rownames(TEMP.DATAslot) <- TEMP.PolyIDs
EUNIS.POLselectRegion <- SpatialPolygonsDataFrame(EUNIS.POLselectRegion,TEMP.DATAslot)

HABITATdataset[["SubregionEUNIS"]][[SEL.BOUND]] <- EUNIS.POLselectRegion

rm(LongMin,LatMin,LongMax,LatMax,SEL.BOUND,EUNIS.POLselectRegion,temp.intersect,TEMP.POLlist,TEMP.DATAslot,TEMP.PolyIDs,TEMP.Boundaries2,Pols1,Pol1)}       #
       

 
 
     ### 2.3.3 -  Surface by subregion and habitat type
     # preliminary projection of the dataset ...
     
HABITATdataset[["SubregionHABITATsurf"]] <- list()

     for(SEL.BOUND in names(HABITATdataset[["RegionLim"]][["BoundariesDefinition"]])){
     print(paste("Surface calculation by habitat type for sub-region",SEL.BOUND))
HAB.SURFACE <- list()
EUNIS.POLselectRegion <- HABITATdataset[["SubregionEUNIS"]][[SEL.BOUND]]

     TEMP.SURF <- spTransform(EUNIS.POLselectRegion, CRS("+proj=utm +zone=30 +datum=WGS84"))
     EUNIS.POLselectRegion@data["PolySurface"] <-  rgeos::gArea(TEMP.SURF,byid=TRUE)
     
          ##
     
     TEMP.FORMULA <- as.formula(paste("PolySurface~",paste(SEL.HABtype,collapse="+")))
     HAB.SURFACE[["ALL"]] <- aggregate(TEMP.FORMULA,EUNIS.POLselectRegion@data,FUN="sum")
     rm(TEMP.FORMULA)
     
          for(HAB in SEL.HABtype){
     TEMP.FORMULA <- as.formula(paste("PolySurface~",HAB))
     HAB.SURFACE[[HAB]] <- aggregate(TEMP.FORMULA,HAB.SURFACE[["ALL"]],FUN="sum")
     HAB.SURFACE[[HAB]] <- HAB.SURFACE[[HAB]][order(HAB.SURFACE[[HAB]]$PolySurface,decreasing=T),]
     rm(TEMP.FORMULA,HAB)}

     HABITATdataset[["SubregionHABITATsurf"]][[SEL.BOUND]] <- HAB.SURFACE

rm(SEL.BOUND,HAB.SURFACE,EUNIS.POLselectRegion,TEMP.SURF)}
     
     #


      



   # PDF report    
SEL.PDF2=F
if(SEL.PDF2==T)
{setwd(paste(DIR.OUPUT[[paste(1)]]))
pdf(file=file.path(DIR.OUPUT[[paste(1)]],paste(paste(SEL.ScriptName,"HABITATdataset",sep="_"), ".pdf",sep="")))
   #


for(SEL.BOUND in names(HABITATdataset[["RegionLim"]][["BoundariesDefinition"]])){

par(mfrow=c(1,1),mar=c(10,6,1,1))
for(HAB in SEL.HABtype){barplot(HABITATdataset[["SubregionHABITATsurf"]][[SEL.BOUND]][[HAB]]$PolySurface,names.arg=HABITATdataset[["SubregionHABITATsurf"]][[SEL.BOUND]][[HAB]][[HAB]],las=2,ylab="Total surface",cex.names=0.8)
title(main=paste(HAB,"habitats for Subregion",SEL.BOUND),cex.main=0.6)}                                                                           


# Test plot
if(T){

par(mfrow=c(1,1),mar=c(1,1,1,1))
plot(HABITATdataset[["SubregionEUNIS"]][[SEL.BOUND]])
plot(rworldmap::getMap(resolution = "high"), add=T,col="darkgrey",lty=0) # Coastline
title(main=paste(SEL.BOUND, "boundaries for analysis"),cex.main=0.6)
plot(HABITATdataset[["SubregionPolygon"]][[SEL.BOUND]],add=T,lwd=2)
}
}

# END PDF REPORT
if(T){for(i in 1:length(dev.list())){dev.off()
rm(i)}}
} # End pdf report
#
       


#######################################################
#######################################################
#######################################################
#######################################################




####### NEEDED PREVIOUS SCRIPT PART TO BE RUN BEFORE (at least initiation script)
####### NEEDED following object to be loaded
##




        ### 2.3.4 - CHECK/Describe sampling stations OVER EUNIS or MSFD habitats
            
for(SURVEY in names(BIOdataset[["FORMATED"]])){
print(paste("Match stations vs habitats types for",SURVEY,"dataset"))                         
        # 2.3.4.1 - Format survey coordinates
TEMP.DATA <- as.data.frame(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]])
TEMP.DATA$x <- as.numeric(as.vector(TEMP.DATA[["Longitude"]]))
TEMP.DATA$y <- as.numeric(as.vector(TEMP.DATA[["Latitude"]]))
coordinates(TEMP.DATA) <- ~ x+y     
proj4string(TEMP.DATA) <- SEL.EUNISproj
#



        # 2.3.4.2 -  Attribute EUNIS/MSFD habitats to stations
        # Select EUNIS/MSFD polygons
        TEMP1 <- sp::over(TEMP.DATA,HABITATdataset[["SubregionEUNIS"]][[grep("\\.",names(HABITATdataset[["RegionLim"]][["BoundariesDefinition"]]))]])
        TEMP.DATA$EUNIS <- TEMP1$EUNIS
        TEMP.DATA$HAB_SURF <- TEMP1$Shape_Area
        TEMP.DATA$MSFD <- TEMP1$MSFD

        # 2.3.4.3 -  Attribute subregions to stations
        for(SEL.BOUND in names(HABITATdataset[["RegionLim"]][["BoundariesDefinition"]])){
        TEMP.DATA[[SEL.BOUND]] <- sp::over(TEMP.DATA,HABITATdataset[["SubregionPolygon"]][[SEL.BOUND]])
        rm(SEL.BOUND)}
        #
        BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]] <- TEMP.DATA
        #
        rm(TEMP1)
        
rm(SURVEY,TEMP.DATA)} # End Loop by survey
         
 # Save R object
if(F){TEMP.DIR <- file.path(DIR.OUPUT[[paste(1)]])
if(file.exists(TEMP.DIR)==F){dir.create(TEMP.DIR)}
save(BIOdataset,file=file.path(TEMP.DIR,paste(SEL.ScriptName,"BIOdataset.RData",sep="")))
rm(TEMP.DIR)}
  #





if(F)
{# SELECT METADATA within the selected polygons limit
TEMP.intersect <- as.vector(gIntersects(SELECT.Boundaries,TEMP.DATA,byid=T))
TEMP.DATA <- TEMP.DATA[TEMP.intersect,]
rm(TEMP.intersect)
#
}


######################################
 ### 2.4 - Depth dataset
######################################


    ## 2.4.1 - Attribute depths to stations

DATA.DEPTH <- loadRData(file.path(DIR.SOURCE[[paste("5a")]],SOURCEFILE.DEPTH))

for(SURVEY in names(BIOdataset[["FORMATED"]])){
BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]$Depth <- extract(DATA.DEPTH, BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]], method='simple', buffer=NULL)
rm(SURVEY)}


 # Save R object
if(F){
TEMP.DIR <- file.path(DIR.OUPUT[[paste(1)]])
if(file.exists(TEMP.DIR)==F){dir.create(TEMP.DIR)}
save(BIOdataset,file=file.path(TEMP.DIR,paste(SEL.ScriptName,"BIOdataset.RData",sep="")))
rm(TEMP.DIR)}
  #


  ###
  
  
  
  
  
######################################  
 ### 2.5 - Fishing pressure dataset
######################################

 SEL.VARfishing <- 'SurfSAR'
 FISHINGdataset <- list()
 
 SEL.COORDINATESnames <- c("mid_lon","mid_lat")
 
 for(FILE in SOURCEFILE.FISHING){
 YEAR <- as.numeric(gsub("\\D", "", FILE))
 if(is.na(YEAR)){YEAR <- FILE} 
 FISHINGdataset[[paste(YEAR)]] <- list()
 FISHINGdataset[[paste(YEAR)]][["DBF"]] <- foreign::read.dbf(file.path(DIR.SOURCE[[paste(4)]],FILE))
 TEMP0 <- data.frame(FISHINGdataset[[paste(YEAR)]][["DBF"]][,c(SEL.COORDINATESnames,SEL.VARfishing)])
 colnames(TEMP0)[match(SEL.COORDINATESnames,colnames(TEMP0))] <- c("MidLon","MidLat")
 coordinates(TEMP0) <- ~MidLon+MidLat
 ext <- extent(TEMP0)
 TEMP.raster <- raster(ncols=length(unique(TEMP0$MidLon)), nrows=length(unique(TEMP0$MidLat)),ext) 
 FISHINGdataset[[paste(YEAR)]][["RASTER"]] <- rasterize(TEMP0, TEMP.raster, SEL.VARfishing, fun=max)
 crs(FISHINGdataset[[paste(YEAR)]][["RASTER"]]) <- SEL.STANDARDproj # Apply standard projection defined above : To CHECK / MAKE it MORE GENERIC !!!!
 rm(TEMP0,TEMP.raster,FILE,YEAR)}
 
  ## 2.5.2 - Attribute FISHING variable to stations
  
  for(SURVEY in names(BIOdataset[["FORMATED"]])){
TEMP.FISHING <- NULL
  for(NAME in names(FISHINGdataset)){
TEMP.FISHING <- cbind(TEMP.FISHING,extract(FISHINGdataset[[NAME]][["RASTER"]], BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]], method='simple', buffer=NULL))
rm(NAME)}

### CORRECT FOR FALSE 0 (NO DATA ~ 0)
for(j in 1:ncol(TEMP.FISHING)){
TEMP.ZERO <- which(is.na(TEMP.FISHING[,j]))
if(length(TEMP.ZERO)>0){
for(i in TEMP.ZERO){
if(sum(TEMP.FISHING[i,],na.rm=T)>0){TEMP.FISHING[,j][i] <- 0}
rm(i)}
}
rm(j)}
###

colnames(TEMP.FISHING) <- paste(SEL.VARfishing,names(FISHINGdataset),sep="")
TEMP.FISHING <- cbind(TEMP.FISHING,rowMeans(TEMP.FISHING,na.rm=T))

colnames(TEMP.FISHING)[ncol(TEMP.FISHING)] <- paste(SEL.VARfishing,"mean",sep="")



###
#### ADAPT FISHING PRESSURE TO BIO DATASET YEAR
## TO COMPLETE / ADAPT  /MODIFY 
## TO COMPLETE / ADAPT  /MODIFY 
## TO COMPLETE / ADAPT  /MODIFY 
## Simple RULE = YEAR of SAMPLING
## TO COMPLETE / ADAPT  /MODIFY 
## TO COMPLETE / ADAPT  /MODIFY 
## TO COMPLETE / ADAPT  /MODIFY 

TEMP.FishYEARs <- sub(SEL.VARfishing,"",colnames(TEMP.FISHING))

TEMP.FISHING <- cbind(TEMP.FISHING,NA)
colnames(TEMP.FISHING)[ncol(TEMP.FISHING)] <- paste(SEL.VARfishing,"yearly",sep="")

for(YEAR in unique(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]$Year)){
TEMP.SEL <- which(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]$Year==YEAR)
TEMP.COL <- which(colnames(TEMP.FISHING)==paste(SEL.VARfishing,YEAR,sep=""))
if(length(TEMP.COL)==1){
TEMP.FISHING[,paste(SEL.VARfishing,"yearly",sep="")][TEMP.SEL] <- TEMP.FISHING[,TEMP.COL][TEMP.SEL]}
rm(TEMP.COL,TEMP.SEL)}

####
###


BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]] <- cbind(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]],TEMP.FISHING)

rm(SURVEY)}



  ## 2.5.3 - Add fishing pressure data to raster stack
  ##### COMPUTE FISHING SAR RASTER / Add to environmental Raster

for(i in 1:length(FISHINGdataset)){FISHINGdataset[[i]]$RASTERproj <- projectRaster(FISHINGdataset[[i]]$RASTER,RASTERstack, method="bilinear")}


  ##### TO COMPLETE /CORRECT : MODIFY ATTRIBUTION OF Fishing SAR depending on Bio sampling YEAR
  ###
  ### Compute mean fishing SAR
  
TEMP <- FISHINGdataset[[1]]$RASTERproj
TEMP2 <- FISHINGdataset[[1]]$RASTERproj
for(i in 2:length(FISHINGdataset)){
TEMP2 <- stack(TEMP2,resample(FISHINGdataset[[i]]$RASTERproj,FISHINGdataset[[1]]$RASTERproj))
}

TEMP3 <- mean(TEMP2)
names(TEMP3) <-  paste(SEL.VARfishing,"mean",sep="") 
RASTERstack <-  stack(TEMP3,resample(RASTERstack,TEMP3))
####
#####


############################################
#### Adapt RASTERstack to C-Square grid 
#### TO MODIFY BEFORE THIS STEP ??? !!!!
####

RASTERstackproj <- projectRaster(RASTERstack,crs=crs(GRIDS[["WHOLE"]]))
TEMP <- rasterToPoints(RASTERstackproj, fun=NULL, spatial=T)
TEMP.OVER <- over(GRIDS[["WHOLE"]],TEMP,returnList=T)
TEMP.DATA <- NULL
for(i in 1:length(TEMP.OVER)){
TEMP.DATA <- rbind(TEMP.DATA,lapply(lapply(TEMP.OVER[[i]],FUN="as.numeric"),FUN="median"))}

TEMP@data[as.numeric(sub("g","",names(TEMP.OVER))),]
values(RASTERstackproj[[1]])

###
############################################
############################################



  
  
  ##
  
  
  
  ## 2.5.3 - Save R object
TEMP.DIR <- file.path(DIR.OUPUT[[paste(1)]])
if(file.exists(TEMP.DIR)==F){dir.create(TEMP.DIR)}
save(FISHINGdataset,file=file.path(TEMP.DIR,paste(SEL.ScriptName,"FISHINGdataset.RData",sep="")))
save(BIOdataset,file=file.path(TEMP.DIR,paste(SEL.ScriptName,"BIOdataset.RData",sep="")))
save(RASTERstack,file=file.path(TEMP.DIR,paste(SEL.ScriptName,"RASTERstack.RData",sep="")))
save(RASTERstackproj,file=file.path(TEMP.DIR,paste(SEL.ScriptName,"RASTERstackproj.RData",sep="")))
save(HABITATdataset,file=file.path(DIR.OUPUT[[paste(1)]],paste(SEL.ScriptName,"HABITATdataset.RData",sep="")))
  #

rm(TEMP.DIR)
  #




  ### 2.6 - C-Square dataset

## TO COMPLETE
## TO COMPLETE
## TO COMPLETE
## TO COMPLETE
## TO COMPLETE
 

  ###
 
 

###################################################
###################################################
###################################################
###################################################
###################################################
###################################################

