#####################################
######### WGFBIT - October 2019 ##########
##################################### >)))°>
######## SCRIPT dataset compilation - Celtic/BoB Eco-region (WGFBIT)
##
### R version 3.6.2

#################################
##
### PART 4 - PRELIMINARY tests
##
#################################


#################################
### SUMMARY
#################################

 

###################################################
#### 3 - PRELIMINARY tests
###################################################

# Habitat types to be tested
SEL.HABtype <- c("EUNIS","MSFD","Raster_Substrate") 
SEL.COMPONENT <- "MEGAFAUNA"

# Select BTA dataset*

## Select longevity traits tye to be considered for indicator computation
SEL.TRAIT <- c("LO") # Select list of traits to considere (!!!! TO CHECK : LO or LG !!!!!)

## Source/outputs Directories
# Defined in initiation script : WGFBIT_Celtic_S2_Part0_Initiation.R 


##

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


   # SELECTION YEARS
   #SEL.YEARS <- c(2008:max(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]$Year))
   #BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]] <- BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]][is.na(match(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]$Year,SEL.YEARS))==F,]
   #BIOdataset[["FORMATED"]][[SURVEY]][["DATA"]] <-   BIOdataset[["FORMATED"]][[SURVEY]][["DATA"]][is.na(match(BIOdataset[["FORMATED"]][[SURVEY]][["DATA"]]$Station,unique(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]$Station)))==F,]
   #



##


## datasets coverage, habitat, species , traits range ...

## TO COMPLETE
## TO COMPLETE
## TO COMPLETE
# aggregate by components / gear !!!!
#for(GEAR in unique(BIOLOGICAL_METADATA[["ALL"]]$SamplingGearCAT)) #

## TO COMPLETE
## TO COMPLETE
## TO COMPLETE


 ######################################################
 ### 3.0 - TOTAL Number of stations / Fishing pressure threhold and biomass estimates by longevity categories
 ######################################################
  

#TESTS
if(F){TEMP <- unique((BIOdataset$FORMATED$EVHOE$DATA_BTA[,c("valid_name","LO1")]))
length(which(is.na(TEMP$LO1)))/nrow(TEMP)

TEMP.MAP1 <- BIOdataset$FORMATED$EVHOE$DATA_BTA[BIOdataset$FORMATED$EVHOE$DATA_BTA$Year<2017,]

plot(unique(TEMP.MAP1[,c("Longitude", "Latitude")]))
plot(rworldmap::getMap(resolution = "high"), add=T,col="darkgrey",lty=0)

TEMP.MAP1[,c("LO1w","LO2w","LO3w","LO4w")] <- TEMP.MAP1$Biomass*TEMP.MAP1[,c("LO1","LO2","LO3","LO4")]
TEMP.MAP2 <- TEMP.MAP1[is.na(TEMP.MAP1$LO1w)==F,]
TEMP.MAP2 <- aggregate(cbind(LO1w,LO2w,LO3w,LO4w)~Station+Latitude+Longitude,data=TEMP.MAP2,FUN="sum")

TEMP.MAP2[,c("LO1w","LO2w","LO3w","LO4w")] <- log1p(TEMP.MAP2[,c("LO1w","LO2w","LO3w","LO4w")])
TEMP.MAP2$Short <- rowSums(TEMP.MAP2[,c("LO1w","LO2w","LO3w")])
TEMP.MAP2$High <- rowSums(TEMP.MAP2[,c("LO4w","LO3w")])

TEMP.MAP2$Longitude <- as.numeric(TEMP.MAP2$Longitude)
TEMP.MAP2$Latitude <- as.numeric(TEMP.MAP2$Latitude)


##### PLOT BY QUANTILES
TEMP.VAR <- "High" # "High", "LO4w"


            #### WHOLE POINTS
TEMP.MAP3 <- TEMP.MAP2
TEMP.CUT <- cut(TEMP.MAP3[,TEMP.VAR],breaks=quantile(TEMP.MAP3[,TEMP.VAR],probs = seq(0, 1, 0.20)),include.lowest =T)
TEMP.SIZE <- match(TEMP.CUT,levels(TEMP.CUT))/1.5
TEMP.MAP3[,c("QuantLow","QuantUp")] <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", TEMP.CUT) ),
      upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", TEMP.CUT) ))

TEMP.MAP3[,c("x", "y")] <- TEMP.MAP3[,c("Longitude", "Latitude")]
coordinates(TEMP.MAP3) <- c("x", "y")
crs(TEMP.MAP3) <- CRS("+proj=longlat +ellps=WGS84") #####!!! WARNING DEPENDANT ON SOURCE COORDINATES !!!!

par(mfrow=c(1,1),mar=c(2.5,2.5,1,1)) 
plot(TEMP.MAP3@data[,c("Longitude", "Latitude")],type="n")
plot(rworldmap::getMap(resolution = "high"), add=T,col="darkgrey",lty=0)
sp::plot(TEMP.MAP3,pch=20,cex=TEMP.SIZE,col=adjustcolor("darkgreen", alpha.f=0.3),add=T)
legend(x="bottomleft",legend=levels(TEMP.CUT),pt.bg=adjustcolor("darkgreen", alpha.f=0.3),pch=20,pt.cex=sort(unique(TEMP.SIZE)))

            #### WITH low SAR filter
FISHthresh <- 0.5 # SAR threshold
TEMP.MAP4 <- TEMP.MAP3[na.omit(match(BIOdataset$FORMATED$EVHOE$STATION$Station[which(BIOdataset$FORMATED$EVHOE$STATION$SurfSARyearly<=FISHthresh)],TEMP.MAP3$Station)),]

par(mfrow=c(1,1),mar=c(2.5,2.5,1,1))
plot(TEMP.MAP4@data[,c("Longitude", "Latitude")],type="n")
plot(rworldmap::getMap(resolution = "high"), add=T,col="darkgrey",lty=0)
plot(TEMP.MAP4,pch=20,cex=TEMP.SIZE,col=adjustcolor("darkgreen", alpha.f=0.3),add=T)
legend(x="bottomleft",legend=levels(TEMP.CUT),pt.bg=adjustcolor("darkgreen", alpha.f=0.3),pch=20,pt.cex=sort(unique(TEMP.SIZE)))

}



 ######################################################
 ### 3.1 - Number of stations for each habitat TYPES
 ######################################################
 
       ## 3.1.1 - covered Habitats
if(SEL.PDF==T)
{setwd(paste(DIR.OUPUT[[paste(1)]]))
pdf(file=file.path(DIR.OUPUT[[paste(1)]],paste(paste(SEL.ScriptName,"HABITATcoverage",sep="_"), ".pdf",sep="")))}
 
CoveredHabitat <- list()

for(SURVEY in names(BIOdataset[["FORMATED"]])){
CoveredHabitat[[SURVEY]] <- list()

if(length(na.omit(as.vector(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]@data[,"Year"])))>0){
TEMP.MetadataSummaryTABLE <- table(data.frame(as.vector(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]@data[,"Station"]),as.vector(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]@data[,"Year"])))}else{
TEMP.MetadataSummaryTABLE <- table(data.frame(as.vector(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]@data[,"Station"])))}

par(mfrow=c(2,2),mar=c(4,4,2,2))

plot(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]],col="darkblue")
plot(rworldmap::getMap(resolution = "high"), add=T,col="darkgrey",lty=0)
title(main=paste(SURVEY," (",paste(range(unique(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]@data[,"Year"])),collapse="-"),")",sep=""),cex.main=0.9)

if(length(na.omit(ncol(TEMP.MetadataSummaryTABLE)))>0){TEMP.BARPLOT <- colSums(TEMP.MetadataSummaryTABLE)}else{TEMP.BARPLOT <- sum(TEMP.MetadataSummaryTABLE)}

barplot(TEMP.BARPLOT,main=paste("Number of stations covered by",SURVEY),las=2,ylab="Number of stations",cex.main=0.6)
## 


par(mfrow=c(2,2),mar=c(8,4,2,0))
for(TYPE in SEL.HABtype)
{
     
TEMP.FORMULA <- as.formula(paste("HAB_SURF","~",TYPE,sep=""))

TEMP.BIO1 <- BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]@data

if(length(na.omit(TEMP.BIO1[,TYPE]))>0&length(na.omit(TEMP.BIO1[,"HAB_SURF"]))){

TEMP.BIOsurf <- aggregate(TEMP.FORMULA, data=TEMP.BIO1,FUN="sum")
TEMP.BIOsurfORDER <- as.vector(TEMP.BIOsurf[,TYPE][order(TEMP.BIOsurf$HAB_SURF,decreasing=T)])


TEMP.BIO1$Freq <- 1
TEMP.BIO2 <- table(as.data.frame(TEMP.BIO1)[,c("Freq",TYPE)])
TEMP.BIOmissing <- TEMP.BIO2[,TEMP.BIO2==0]
TEMP.BIO2 <- TEMP.BIO2[,TEMP.BIO2>0]
TEMP.BIO2 <- TEMP.BIO2[match(TEMP.BIOsurfORDER,names(TEMP.BIO2))]

if(T){#par(mfrow=c(1,1),mar=c(10,5,2,2))
barplot(TEMP.BIO2,names.arg=colnames(TEMP.BIO2),las=2,xlab=paste(TYPE,"habitat"),ylab="Number of stations",cex.names=0.7)
title(main=paste("Covered habitats from",SURVEY),cex.main=0.8)}

if(F){for(SEL.BOUND in names(HABITATdataset[["RegionLim"]][["BoundariesDefinition"]])){
CoveredHabitat[[SURVEY]][[SEL.BOUND]] <- list()
CoveredHabitat[[SURVEY]][[SEL.BOUND]][[TYPE]] <- names(TEMP.BIO2)
CoveredHabitat[[SURVEY]][[SEL.BOUND]][[TYPE]] <- cbind(CoveredHabitat[[SURVEY]][[SEL.BOUND]][[TYPE]],rainbow(n=length(CoveredHabitat[[SURVEY]][[SEL.BOUND]][[TYPE]])))
rm(SEL.BOUND)}
}} #else{print(paste("WARNING no",TYPE,"DATA for",SURVEY,"/",SEL.BOUND))}
rm(TYPE)}
rm(SURVEY)}




##

      ## 3.1.2 - Map / sampling coverage
 
for(SEL.BOUND in names(HABITATdataset[["RegionLim"]][["BoundariesDefinition"]])){

         #
         par(mfrow=c(1,1))
         plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
         text(5.5, 8, paste("Map of samples coverage", sep=""), cex=1.2)
         text(5.5, 7, paste("Region=",SEL.BOUND, sep=" "), cex=1.1)         
        #
 

for(SURVEY in names(BIOdataset[["FORMATED"]])){
for(TYPE in SEL.HABtype)
{

print(paste("Map of sampling dataset coverage for",SEL.BOUND,SURVEY,TYPE))

# TEMP.BIO <- BIOLOGICAL_METADATA[["ALL"]][BIOLOGICAL_METADATA[["ALL"]]$SamplingGearCAT==GEAR,]
if(F){TEMP.EUNIS.POL <- EUNIS.POL[is.na(match(EUNIS.POL@data[,c(TYPE)],CoveredHabitat[[SURVEY]][[SEL.BOUND]][[TYPE]][,1]))==F,]
TEMP.EUNIS.POL@data["COLOR"] <- CoveredHabitat[[SURVEY]][[SEL.BOUND]][[TYPE]][,2][match(TEMP.EUNIS.POL@data[,c(TYPE)],CoveredHabitat[[SURVEY]][[SEL.BOUND]][[TYPE]][,1])]
# plot(TEMP.BIO[,c("longitude","latitude")],col="white")

plot(HABITATdataset[["SubregionPolygon"]][[SEL.BOUND]])
plot(TEMP.EUNIS.POL,col=unlist(TEMP.EUNIS.POL@data["COLOR"]),add=T,lty=0)
points(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]][,c("Longitude","Latitude")],pch=19,cex=0.5)
plot(rworldmap::getMap(resolution = "high"), add=T,col="darkgrey",lty=0)
title(main=paste())}

if(F){plot(data.frame(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]][,c("Longitude","Latitude")])[,c("Longitude","Latitude")],type="n")
plot(TEMP.EUNIS.POL,col=unlist(TEMP.EUNIS.POL@data["COLOR"]),add=T,lty=0)
points(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]][,c("Longitude","Latitude")],pch=19,cex=0.5)
plot(rworldmap::getMap(resolution = "high"), add=T,col="darkgrey",lty=0)
title(main=paste())}


# By subregion
TEMP.EUNIS.POL <- HABITATdataset[["SubregionEUNIS"]][[SEL.BOUND]][is.na(match(HABITATdataset[["SubregionEUNIS"]][[SEL.BOUND]]@data[,c(TYPE)],CoveredHabitat[[SURVEY]][[SEL.BOUND]][[TYPE]][,1]))==F,]
TEMP.EUNIS.POL@data["COLOR"] <- CoveredHabitat[[SURVEY]][[SEL.BOUND]][[TYPE]][,2][match(TEMP.EUNIS.POL@data[,c(TYPE)],CoveredHabitat[[SURVEY]][[SEL.BOUND]][[TYPE]][,1])]

# plot for selected subregion
par(mfrow=c(1,1),mar=c(1,1,1,1))
plot(HABITATdataset[["SubregionEUNIS"]][[SEL.BOUND]])
plot(rworldmap::getMap(resolution = "high"), add=T,col="darkgrey",lty=0) # Coastline
plot(HABITATdataset[["SubregionPolygon"]][[SEL.BOUND]],add=T,lwd=2)
#plot(TEMP.EUNIS.POL,col=unlist(TEMP.EUNIS.POL@data["COLOR"]),add=T,lty=0)
points(BIOdataset[["FORMATED"]][[SURVEY]]$STATION,col="black")
title(main=paste("Covered", TYPE,"habitat from",SURVEY, "sampling","for",SEL.BOUND,"subregion"))

#


rm(TYPE)}
rm(SURVEY)} #
rm(SEL.BOUND)} # 


# END PDF REPORT
if(T){for(i in 1:length(dev.list())){dev.off()
rm(i)}}



##
if(F){
for(name in names(BIOLOGICAL_METADATA))
{
setwd(DIR.SOURCE[[paste(6)]])
write.table(BIOLOGICAL_METADATA[[name]],file=paste("CORR",sub(".csv","",name),".csv",sep=""),sep=";",row.names = F)
rm(name)}}
##



 
 ### 3.2 - BTA matrix coverage for species lists by main habitats

for(SURVEY in names(BIOdataset[["FORMATED"]])){

 if(SEL.PDF==T)
{setwd(paste(DIR.OUPUT[[paste(1)]]))
pdf(file=file.path(DIR.OUPUT[[paste(1)]],paste(paste(SEL.ScriptName,"BTAcoverage",SURVEY,sep="_"), ".pdf",sep="")))}


for(REGION in names(HABITATdataset[["RegionLim"]][["BoundariesDefinition"]])){

for(TYPE in SEL.HABtype)
{




         #
         par(mfrow=c(1,1))
         plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
         text(5.5, 8, paste("Distribution of biomass in longevity classes", sep=""), cex=1.2)
         text(5.5, 7, paste("Region=",REGION, sep=" "), cex=1.1)         
         text(5.5, 6, paste("Habitat type=",TYPE, sep=" "), cex=1.1)
         text(5.5, 5, paste("Survey=",SURVEY, sep=" "), cex=1.1)
        #

TEMP.name <- FILE.NAME[["BTA"]][,1][FILE.NAME[["BTA"]][,1]==SEL.COMPONENT]

#BIOLOGICAL_BTA[[SEL.COMPONENT]]

TEMP.DATA <- BIOdataset[["FORMATED"]][[SURVEY]]$DATA

if(nrow(TEMP.DATA)>3){

if(F){TEMP.DATA <- TEMP.DATA[is.na(match(TEMP.DATA$TaxaSOURCE,BIOdataset[["FORMATED"]][[SURVEY]][["FilteredSpeciesList"]]$TaxaSOURCE))==F,]} # apply species filter (already done before)

TEMP.DATA$Biomass <- as.numeric(TEMP.DATA$Biomass)

TEMP.STATION <- BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]][which(BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]][[REGION]]==1),] # apply subregion filter

# station/espece Matrix
TEMP.DATAmatrix <- table(TEMP.DATA[,c("Station","valid_name")])


TEMP.BTA <- BIOLOGICAL_BTA[[SEL.COMPONENT]][["MATRIX.FUNCnew"]]
TEMP.BTA <- TEMP.BTA[,which(is.na(match(gsub("([^A-Za])+", "", x = colnames(TEMP.BTA)),SEL.TRAIT))==F)]
TEMP.BTA <- TEMP.BTA[rowSums(TEMP.BTA)>0,] # Suppress null rows
SEL.LISTtraits <- colnames(TEMP.BTA)

# ADD HABITAT informaton to dataset
TEMP.DATA[[TYPE]] <- BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]][[TYPE]][match(TEMP.DATA$Station, BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]$Station)]
#

# ADD BTA information at the orginal taxa level in the biological dataset
TEMP.DATA$BTA <- 0
TEMP.DATA$BTA[is.na(match(TEMP.DATA$valid_name,rownames(TEMP.BTA)))==F] <- 1
TEMP.DATA <- cbind(TEMP.DATA,TEMP.BTA[match(TEMP.DATA$valid_name,rownames(TEMP.BTA)),])
TEMP.DATA$Merge_to_traits[TEMP.DATA$BTA==1] <- TEMP.DATA$valid_name[TEMP.DATA$BTA==1]
#


                   # identify missing species
                   TEMP.FORMULA <- as.formula(paste("Biomass~valid_name+",TYPE,sep=""))
                   TEMP.MISSING2 <- TEMP.DATA[TEMP.DATA$BTA==0,]
                   TEMP.MISSING2 <- aggregate(TEMP.FORMULA,data=TEMP.MISSING2,FUN="sum")
                   
                    if(F){
                   setwd(setwd(DIR.SOURCE[[paste(6)]]))
                   save(TEMP.MISSING2,file=paste("TEMP_MISSING",GEAR,TYPE,".Rdata",sep=""))}
                   #

         ## 3.2.1 -  Species richness coverage by habitat type
         
TEMP.richness <- table(TEMP.DATA[,c(TYPE,"valid_name")])
TEMP.richness[TEMP.richness>0] <- 1
TEMP.richnessMISS <- TEMP.richness
TEMP.richnessMISS[,is.na(match(colnames(TEMP.richnessMISS), unique(TEMP.DATA$valid_name[TEMP.DATA$BTA==0])))] <- 0

TEMP.BARPLOT <- as.matrix(rbind(rowSums(TEMP.richness-TEMP.richnessMISS),rowSums(TEMP.richnessMISS)))
TEMP.BARPLOT <- TEMP.BARPLOT[,colSums(TEMP.BARPLOT)>0] # keep only habitats with data
rownames(TEMP.BARPLOT) <- c("Covered", "Missing")

TEMP.BARPLOTprop <- t(t(TEMP.BARPLOT)/rowSums(t(TEMP.BARPLOT)))
rownames(TEMP.BARPLOTprop) <- c("Covered", "Missing")


par(mfrow=c(1,1),mar=c(17,5,2,2))
#barplot(colSums(TEMP.BARPLOT),ylab="Taxa richness",main=paste(TYPE,"habitats"),las=2)
barplot(TEMP.BARPLOT,las=2,ylab="Taxa richness",
legend.text = T,args.legend = list(x = "topleft", bty = "n"))
title(main=paste("Species richness coverage from longevity matrix",paste("(",SEL.COMPONENT," ",SURVEY,")",sep="")),cex.main=0.8,sub=TYPE)


par(mfrow=c(1,1),mar=c(17,5,2,2))
barplot(TEMP.BARPLOTprop,las=2,ylab="Taxa richness",
legend.text = T,args.legend = list(x = "bottom", bty = "n"))
title(main=paste("Species richness coverage from longevity matrix",paste("(",SEL.COMPONENT," ",SURVEY,")",sep="")),cex.main=0.6)

                                                                                                                                                                      



         ## 3.2.2 - Biomass coverage by habitat type

                            # GLOBAL identification of missing species
                   TEMP.FORMULA <- as.formula(paste("Biomass~",TYPE,sep=""))
                   TEMP.MISSINGmass1 <- aggregate(TEMP.FORMULA,data=TEMP.DATA,FUN="sum")
                   TEMP.MISSINGmass2 <- aggregate(TEMP.FORMULA,data=TEMP.DATA[TEMP.DATA$BTA==1,],FUN="sum")
                   TEMP.MISSINGmass  <- merge(TEMP.MISSINGmass1,TEMP.MISSINGmass2,by=TYPE)
                   TEMP.MISSINGmass  <- cbind(TEMP.MISSINGmass, TEMP.MISSINGmass[,2]-TEMP.MISSINGmass[,3])
                   colnames(TEMP.MISSINGmass)[c(2:4)] <- c("Total","Covered","Missing")


TEMP.MATRIXmass <- t(as.matrix(TEMP.MISSINGmass[,c("Covered","Missing")]))
colnames(TEMP.MATRIXmass) <- as.vector(TEMP.MISSINGmass[[TYPE]])
TEMP.MATRIXmassprop <- t(t(TEMP.MATRIXmass)/rowSums(t(TEMP.MATRIXmass)))


par(mfrow=c(1,1),mar=c(17,5,2,2))
barplot(TEMP.MATRIXmassprop,las=2,ylab="Biomass proportion (%)",
legend.text = T,args.legend = list(x = "topright", bty = "n"))
title(main=paste("Biomass coverage from longevity matrix",paste("(",SEL.COMPONENT," ",SURVEY,")",sep="")),cex.main=0.8,sub=TYPE)


         # by sampling station
                   TEMP.FORMULA <- as.formula(paste("Biomass~",paste("Station",TYPE,sep="+"),sep=""))
                   TEMP.MISSINGmass1 <- aggregate(TEMP.FORMULA,data=TEMP.DATA,FUN="sum")
                   TEMP.MISSINGmass2 <- aggregate(TEMP.FORMULA,data=TEMP.DATA[TEMP.DATA$BTA==1,],FUN="sum")
                   TEMP.MISSINGmass  <- merge(TEMP.MISSINGmass1,TEMP.MISSINGmass2,by=c("Station",TYPE))
                   TEMP.MISSINGmass  <- cbind(TEMP.MISSINGmass, TEMP.MISSINGmass[,3]-TEMP.MISSINGmass[,4])
                   colnames(TEMP.MISSINGmass)[c(3:5)] <- c("Total","Covered","Missing")
                   TEMP.MISSINGmass$prop <- TEMP.MISSINGmass$Covered/TEMP.MISSINGmass$Total

                   TEMP.FORMULA2 <- as.formula(paste("prop~",TYPE,sep=""))
                   boxplot(TEMP.FORMULA2,data=TEMP.MISSINGmass,las=2,ylab="Covered biomass proportion")
                   title(main=paste("Biomass coverage by station from longevity matrix",paste("(",SEL.COMPONENT," ",SURVEY," dataset",")",sep="")),cex.main=0.8) 




 ### 3.3 - Biomass Distribution / Longevity classes & Habitat
 

   ## 3.3.1 - For all stations with biological observations
TEMP.BIOMASS  <- TEMP.DATA[TEMP.DATA$BTA==1,] # Select stations with longevity information
TEMP.BIOMASS <- data.frame(TEMP.BIOMASS[,TYPE],TEMP.BIOMASS[,SEL.LISTtraits]*TEMP.BIOMASS$Biomass)
colnames(TEMP.BIOMASS)[1] <- TYPE


            # biomass distribution FOR ALL the dataset
             LONGprop.HAB <- NULL
             for(HAB in unique(TEMP.BIOMASS[,TYPE]))
             {TEMP.BIOMASS2 <- TEMP.BIOMASS[TEMP.BIOMASS[,TYPE]==HAB,]
             LONGprop.HAB <- rbind(LONGprop.HAB,colSums(TEMP.BIOMASS2[,SEL.LISTtraits],na.rm = T))
             rm(HAB)}
             row.names(LONGprop.HAB) <- unique(TEMP.BIOMASS[,TYPE])
             LONGprop.HABrel <- LONGprop.HAB/rowSums(LONGprop.HAB)

         
             # Order / Habitat mean depth
             TEMP.formulaDepth <- as.formula(paste("Depth~",TYPE))
             TEMP.DepthOrder <- aggregate(TEMP.formulaDepth,data=BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]],FUN="mean")    
             TEMP.DepthOrder <- TEMP.DepthOrder[order(TEMP.DepthOrder$Depth,decreasing=T),]
             LONGprop.HABrel <- LONGprop.HABrel[match(TEMP.DepthOrder[[TYPE]],rownames(LONGprop.HABrel)),]
             rm(TEMP.formulaDepth)
             #


 
             # BARPLOT

                          TEMP.COL <- grey(0:(nrow(t(LONGprop.HAB))-1) / (nrow(t(LONGprop.HAB))-1))
                          TEMP.NAMES <- rownames(t(LONGprop.HAB))
                          TEMP.NAMES <- sub("X","",TEMP.NAMES)
                          TEMP.NAMES <- sub("w","",TEMP.NAMES)
                          TEMP.NAMES <- sub("\\.","-",TEMP.NAMES)
                          TEMP.NAMES[TEMP.NAMES=="-10years"] <- ">10years"
                          TEMP.NAMES[TEMP.NAMES=="-1year"] <- "<1year"
                          TEMP.NAMES <- sub("year"," year",TEMP.NAMES)

             par(mfrow=c(1,1),mar=c(10,5,4,2),mai=c(1,1,1,1))
             
             if(F){barplot(t(LONGprop.HAB),las=2,main=paste("Biomass distribution for each habitat and longevity class (",SEL.COMPONENT,")",sep=""),col=TEMP.COL)}
                         
             barplot(t(LONGprop.HABrel),las=2,main=paste("Biomass distribution & longevity class (",SEL.COMPONENT,"/",SURVEY, " dataset",")",sep=""),sub=paste("Habitat:",TYPE)
             ,col=TEMP.COL,cex.names=0.6,space=0.3,ylab="Biomass proportion",cex.main=0.8)
            legend("topright",rev(TEMP.NAMES),fill=rev(TEMP.COL))



 TEMP.BIOMASSsta  <- TEMP.DATA[TEMP.DATA$BTA==1,] # Select stations with longevity information
 TEMP.BIOMASSsta <- data.frame(TEMP.BIOMASSsta[,c("Station",TYPE)],TEMP.BIOMASSsta[,SEL.LISTtraits]*TEMP.BIOMASSsta$Biomass)
 TEMP.BIOMASSsta[,"Station"] <- as.character(TEMP.BIOMASSsta[,"Station"])
 TEMP.BIOMASSsta[,TYPE] <- as.character(TEMP.BIOMASSsta[,TYPE])
 TEMP.BIOMASSsta[,SEL.LISTtraits] <- lapply(TEMP.BIOMASSsta[,SEL.LISTtraits],"as.character")
 TEMP.BIOMASSsta[,SEL.LISTtraits] <- lapply(TEMP.BIOMASSsta[,SEL.LISTtraits],"as.numeric")
 TEMP.BIOMASSsta <- aggregate(as.formula(paste(paste("cbind(",paste(SEL.LISTtraits,collapse=","),")",sep=""),"~",paste("Station", TYPE,sep="+"))),data= TEMP.BIOMASSsta,FUN="sum")
 
 
 
 TEMP.BIOMASSsta <- reshape2::melt(TEMP.BIOMASSsta, measure.vars=SEL.LISTtraits,id.vars=c(TYPE))
 colnames(TEMP.BIOMASSsta)[colnames(TEMP.BIOMASSsta)==TYPE] <- "Habitat" 

 # Order habitat by depth
 TEMP.BIOMASSsta <- TEMP.BIOMASSsta[order(match(TEMP.BIOMASSsta$Habitat,TEMP.DepthOrder[[TYPE]])),] 
 #
 
 library(ggplot2)
 if(F){ggplot2::ggplot(data = TEMP.BIOMASSsta, aes(x=Habitat, y=value)) + geom_boxplot(aes(fill=variable))}


 # Boxplot all
 box2 <- ggplot(data = TEMP.BIOMASSsta, aes(x=Habitat, y=value)) + 
             geom_boxplot(aes(fill=variable))
 plot(box2 + facet_wrap( ~ Habitat, scales="free"))
 
 
 # Boxplot by habitat
 box1 <- ggplot(data = TEMP.BIOMASSsta, aes(x=Habitat, y=value)) + 
             geom_boxplot(aes(fill=variable),outlier.shape = NA)+scale_y_continuous(limits = quantile(TEMP.BIOMASSsta$value, c(0.1, 0.9)))
 plot(box1 + facet_wrap( ~ Habitat, scales="free"))

 rm(box2,box1,TEMP.BIOMASSsta,TEMP.COL,TEMP.NAMES,TEMP.DepthOrder)
        
        
       

        ## 3.3.2 - fuzzy pca to check for species/traits distribution for each habitat type ...
        


 ## TO ADD / COMPLETE
 ## TO ADD / COMPLETE
 ## TO ADD / COMPLETE
 ## TO ADD / COMPLETE
 # fuzzy pca to check for species/traits distribution for each habitat type ...
TEMP.BIOMASSfuzzy  <- TEMP.DATA[TEMP.DATA$BTA==1,] # Select stations with longevity information


## Distrib TOTAL

TEMP.df <- TEMP.BIOMASSfuzzy[,SEL.LISTtraits]
TEMP.colblocks <- table(gsub('[[:digit:]]+', '', colnames(TEMP.df)))
TEMP.roww <- TEMP.BIOMASSfuzzy[,"Biomass"]

TEMP.fuzzy <-  ade4::prep.fuzzy.var(TEMP.df, TEMP.colblocks, row.w = TEMP.roww)
TEMP.FCA <- ade4::dudi.fca(TEMP.fuzzy, scannf = FALSE, nf = 2)
TEMP.FPCA <- ade4::dudi.fpca(TEMP.fuzzy, scannf = FALSE, nf = 2)

## plot
ade4::scatter(TEMP.FCA, csub = 3, clab.moda = 1)
if(F){ade4::scatter(TEMP.FPCA, csub=3, clab.moda = 2)}
if(F){biplot(TEMP.FPCA)}
if(F){screeplot(TEMP.FPCA, npcs = 20, type = c("barplot"))}
#
rm(TEMP.BIOMASSfuzzy)

## DISTRIB par HABITAT

TEMP.DATA[[TYPE]][is.na(TEMP.DATA[[TYPE]])] <- "NONE"

for(HAB in unique(TEMP.DATA[[TYPE]])){

TEMP.BIOMASSfuzzy  <- TEMP.DATA[TEMP.DATA$BTA==1&TEMP.DATA[,TYPE]==HAB,] # Select stations with longevity information

TEMP.df <- TEMP.BIOMASSfuzzy[,SEL.LISTtraits]
TEMP.df <- TEMP.df[,colSums(TEMP.df)!=0]

TEMP.colblocks <- table(gsub('[[:digit:]]+', '', colnames(TEMP.df)))
TEMP.roww <- TEMP.BIOMASSfuzzy[,"Biomass"]

TEMP.fuzzy <-  ade4::prep.fuzzy.var(TEMP.df, TEMP.colblocks, row.w = TEMP.roww)
TEMP.FCA <- ade4::dudi.fca(TEMP.fuzzy, scannf = FALSE, nf = 2)
TEMP.FPCA <- ade4::dudi.fpca(TEMP.fuzzy, scannf = FALSE, nf = 2)

## plot
ade4::scatter(TEMP.FCA, csub = 3, clab.moda = 1)
title(sub=paste(TYPE,HAB))
if(F){ade4::scatter(TEMP.FPCA, csub=3, clab.moda = 2)}
if(F){biplot(TEMP.FPCA)}
if(F){screeplot(TEMP.FPCA, npcs = 20, type = c("barplot"))}
#
rm(TEMP.BIOMASSfuzzy,TEMP.df)

rm(HAB)}



 #
 ## TO ADD / COMPLETE
 ## TO ADD / COMPLETE
 ## TO ADD / COMPLETE
 ## TO ADD / COMPLETE
        


        
        ## 3.3.3 - Against fishing pressure / with NULL or low fishing pressure
 
 
 ### Test for coverage of stations with low or NULL fishing pressure
# SELthreshold.FISHnull <- ?????
###

SEL.VARfishing <- "SurfaceSARmean"        
        # ADD FISHING information to dataset
        TEMP.DATA[[paste(SEL.VARfishing,"mean",sep="")]] <- BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]@data[[paste(SEL.VARfishing,"mean",sep="")]][match(TEMP.DATA$Station, BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]]@data$Station)]

        for(HAB in unique(TEMP.DATA[[TYPE]])){


 ### TO COMPLETE FROM HERE
          ### TO COMPLETE FROM HERE
               ### TO COMPLETE FROM HERE
                    ### TO COMPLETE FROM HERE
                         ### TO COMPLETE FROM HERE
                              ### TO COMPLETE FROM HERE


        
        # HAB <- "Upper bathyal sediment"
        TEMP.BIOMASSsta2  <- TEMP.DATA[TEMP.DATA[[TYPE]]==HAB&TEMP.DATA$BTA==1,] # Select HAB stations with longevity information
        
        if(length(na.omit(TEMP.BIOMASSsta2[[TYPE]]))>0){
        
        TEMP.BIOMASSsta2 <- data.frame(TEMP.BIOMASSsta2[,c("Station",TYPE,paste(SEL.VARfishing,"mean",sep=""))],TEMP.BIOMASSsta2[,SEL.LISTtraits]*TEMP.BIOMASSsta2$Biomass)
        TEMP.BIOMASSsta2[,"Station"] <- as.character(TEMP.BIOMASSsta2[,"Station"])
        TEMP.BIOMASSsta2[,TYPE] <- as.character(TEMP.BIOMASSsta2[,TYPE])
        TEMP.BIOMASSsta2[,SEL.LISTtraits] <- lapply(TEMP.BIOMASSsta2[,SEL.LISTtraits],"as.character")
        TEMP.BIOMASSsta2[,SEL.LISTtraits] <- lapply(TEMP.BIOMASSsta2[,SEL.LISTtraits],"as.numeric")
        TEMP.BIOMASSsta2 <- aggregate(as.formula(paste(paste("cbind(",paste(SEL.LISTtraits,collapse=","),")",sep=""),"~",paste("Station", TYPE,paste(SEL.VARfishing,"mean",sep=""),sep="+"))),data= TEMP.BIOMASSsta2,FUN="sum")

        if(nrow(TEMP.BIOMASSsta2)>3){

        # Fishing effort classes
        TEMP.QUANT <- quantile(TEMP.BIOMASSsta2[,paste(SEL.VARfishing,"mean",sep="")],probs = seq(0, 1, 0.10))
        TEMP.QUANT[1] <- 0
        TEMP.QUANT <- unique(TEMP.QUANT)
        TEMP.BIOMASSsta2$quant <- cut(TEMP.BIOMASSsta2[,paste(SEL.VARfishing,"mean",sep="")], TEMP.QUANT)
     
        par(mfrow=c(2,2), mar=c(6,4,2,2))
        for(VAR in SEL.LISTtraits){
        boxplot(as.formula(paste(VAR,"~","quant")),data= TEMP.BIOMASSsta2,outline=F,las=2,main=paste(VAR),ylab="Biomass")}                     
        # glm fit
        #par(mfrow=c(2,2), mar=c(6,4,2,2))
        SEL.family <- "gaussian" # "gaussian", "poisson"
        SEL.LIMp <- 0.05
                # Select/Transfo data
        TEMP.glmData <- TEMP.BIOMASSsta2[,c(paste(SEL.VARfishing,"mean",sep=""),SEL.LISTtraits)]
        if(F){TEMP.glmData <- as.data.frame(apply(TEMP.glmData, 2, function(i) log1p(i)/max(log1p(i))))}
        if(T){TEMP.glmData <- as.data.frame(apply(TEMP.glmData, 2, function(i) log1p(i)))}


        # Glm & plot
        par(mfrow=c(1,1), mar=c(6,4,2,2))
        plot(range(TEMP.glmData[,paste(SEL.VARfishing,"mean",sep="")]),range(as.vector(unlist(TEMP.glmData[,SEL.LISTtraits]))),type="n",xlab=paste(SEL.VARfishing,"mean",sep=""),ylab="Biomass")
        title(main=paste(TYPE,HAB))
        
        TEMP.COLOR <- rainbow(length(SEL.LISTtraits))
        legend("topright", SEL.LISTtraits, pch = 20,col=TEMP.COLOR,bty="n")
        
        for(VAR in SEL.LISTtraits){
        
        TEMP.glmdata2 <- TEMP.glmData[,c(paste(SEL.VARfishing,"mean",sep=""),VAR)]
       
        
        TEMP.glm <- glm(as.formula(paste(VAR,"~",paste(SEL.VARfishing,"mean",sep=""))),data=TEMP.glmdata2,family = SEL.family)
        #summary(TEMP.glm)
        
        TEMP.xvar <-  seq(0,ceiling(max(TEMP.glmdata2[,paste(SEL.VARfishing,"mean",sep="")])),0.01)
        TEMP.yvar <- predict(TEMP.glm,list(SurfaceSARmean=TEMP.xvar),type="response")
        points(TEMP.glmdata2[,c(paste(SEL.VARfishing,"mean",sep=""),VAR)],pch=20,col=TEMP.COLOR[which(SEL.LISTtraits==VAR)])
        
        TEMP.p <- summary(TEMP.glm)$coefficients[,4][2] # "Pr(>|z|
        if(TEMP.p<SEL.LIMp){lines(TEMP.xvar,TEMP.yvar,col=TEMP.COLOR[which(SEL.LISTtraits==VAR)],lwd=2)}
        rm(VAR)}
        } # End test number of stations for HAB
        }else{       #
         par(mfrow=c(1,1))
         plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
         text(5.5, 8, paste("NO DATA for HAB",HAB, sep=""), cex=1.5)
        #
} # Test Nb Data
        rm(HAB,TEMP.BIOMASSsta2)}

#

        
        ##
        
        
}else{

         #
         par(mfrow=c(1,1))
         plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
         text(5.5, 8, paste("NO DATA", sep=""), cex=1.5)
        #
}  # Test empty dataset

rm(TYPE)} # End loop habitat type
rm(REGION)} # End loop region delimitation

# END PDF REPORT
if(T){for(i in 1:length(dev.list())){dev.off()
rm(i)}}

rm(SURVEY)}
      
 ### TO COMPLETE
### TO COMPLETE
### TO COMPLETE

     ##


#### TRASH #####
#### TRASH #####
#### TRASH #####
#### TRASH #####
#### TRASH #####


#### TRASH #####
#### TRASH #####
#### TRASH #####
#### TRASH #####
#### TRASH #####


#### 4 - Formatting and writing files for each FBIT steps
 

SEL.STDfield.BENTHIC_DATA <- c("STATION_ID","Merge_to_traits","Biomass",BTA.NAMEstd[,2]) #sample_ID	Nomen	Merge_to_traits	sample_m2	Biomass	L1	L1_3	L3_10	L10
SEL.STDfilename.BENTHIC_DATA <- c("Benthic_Data")
SEL.STDfileformat.BENTHIC_DATA <- ".csv"

SEL.STDfield.BENTHIC_METADATA <- c("STATION_ID","Year","longitude","latitude",TYPE) # ID	replicate	poslon	poslat	SUR_SWEPT_KM2	SUBSUR_SWEPT_KM2	MSFDhab
SEL.STDfilename.BENTHIC_METADATA <- c("Benthic_Data")
SEL.STDfileformat.BENTHIC_METADATA <- ".csv"



# CREATING REFERENCE FILES FOR FBIT workflow
for(GEAR in unique(BIOLOGICAL_METADATA[["ALL"]]$SamplingGearCAT))
{

TEMP.METADATA <- BIOLOGICAL_METADATA[["ALL"]][BIOLOGICAL_METADATA[["ALL"]]$SamplingGearCAT==GEAR,]
TEMP.BIOMASS <- BIOLOGICAL_BIOMASS[["ALLbta"]][is.na(match(BIOLOGICAL_BIOMASS[["ALLbta"]]$STATION_ID,TEMP.METADATA$STATION_ID))==F,]

TEMP.BIOMASS <- TEMP.BIOMASS[,SEL.STDfield.BENTHIC_DATA]
TEMP.BIOMASS[,"sample_m2"] <- 1  # TO ADD TO PREVIOUS FORMATTING STEP !!!!!
colnames(TEMP.BIOMASS)[match(BTA.NAMEstd[,2],colnames(TEMP.BIOMASS))] <- BTA.NAMEstd[,4] 

colnames(TEMP.BIOMASS)[colnames(TEMP.BIOMASS)=="STATION_ID"] <- "sample_ID"
names(TEMP.METADATA)[names(TEMP.METADATA)=="STATION_ID"] <- "ID"

TEMP.BIOMASS <- TEMP.BIOMASS[TEMP.BIOMASS$Biomass>0,]

setwd(DIR.SOURCE[[paste(7)]])
write.table(TEMP.BIOMASS,file=paste("Benthic_Data_",GEAR,"_",TYPE,".csv",sep=""),sep=";",row.names = F)
write.table(TEMP.METADATA,file=paste("Env_Data_",GEAR,"_",TYPE,".csv",sep=""),sep=";",row.names = F)
}


###







 if(SEL.PDF)
{for(i in 1:length(dev.list())){dev.off()
                                rm(i)}}



## DATASET from EVHOE IBTS
if(F)
{setwd(file.path(DIR.SOURCE[[paste(1)]],Type,"METADATA"))
write.table(DATASET$Station,file="METADATA_Megafauna_PL.csv",sep=";",row.names = F)

DATA_MegafaunaEVHOE <- as.data.frame(MATRIXspDW)
DATA_MegafaunaEVHOE$Taxa <- dataset.taxref1$L_VALIDE[match(DATA_MegafaunaEVHOE$Espece,dataset.taxref1$C_VALIDE)]
write.table(DATA_MegafaunaEVHOE,file="DATA_MegafaunaEVHOE.csv",sep=";",row.names = F)
      
}









#### TRASH #####
#### TRASH #####
#### TRASH #####
#### TRASH #####
#### TRASH #####



BIOLOGICAL_METADATA <- list()
Type <- "BIOMASS"

BiodatasetLIST <- paste(BiodatasetMetadata[,"OwnerCountry"],BiodatasetMetadata[,"SOURCEname"],sep="_")

for(name in list.files(file.path(DIR.SOURCE[[paste(3)]],Type,"METADATA")))
{BIOLOGICAL_METADATA[[name]] <- read.csv2(file=file.path(DIR.SOURCE[[paste(2)]],Type,"METADATA",name),header=T,sep=";",dec = ".",colClasses = "character")
rm(name)}

for(name in names(BIOLOGICAL_METADATA))
{
if(length(BIOLOGICAL_METADATA[["ALL"]])==0)
{BIOLOGICAL_METADATA[["ALL"]] <- BIOLOGICAL_METADATA[[name]]}else{
TEMP <- BIOLOGICAL_METADATA[[name]][,intersect(colnames(BIOLOGICAL_METADATA[["ALL"]]),colnames(BIOLOGICAL_METADATA[[name]]))]
BIOLOGICAL_METADATA[["ALL"]] <- BIOLOGICAL_METADATA[["ALL"]][,colnames(TEMP)]

BIOLOGICAL_METADATA[["ALL"]] <- rbind(BIOLOGICAL_METADATA[["ALL"]],TEMP[,match(colnames(BIOLOGICAL_METADATA[["ALL"]]),colnames(TEMP))])
}}

#
TEMP.MISSING <- NULL
for(VAR in c("longitude","latitude"))
{BIOLOGICAL_METADATA[["ALL"]][,VAR] <- as.numeric(BIOLOGICAL_METADATA[["ALL"]][,VAR])
TEMP.MISSING <- rbind(TEMP.MISSING,BIOLOGICAL_METADATA[["ALL"]][is.na(BIOLOGICAL_METADATA[["ALL"]][,VAR])==T,])
BIOLOGICAL_METADATA[["ALL"]] <- BIOLOGICAL_METADATA[["ALL"]][is.na(BIOLOGICAL_METADATA[["ALL"]][,VAR])==F,]}
#

# Sampling gear list
BIOLOGICAL_METADATA[["ALL"]]$SamplingGearCAT <- NA 
BIOLOGICAL_METADATA[["ALL"]]$SamplingGearCAT[BIOLOGICAL_METADATA[["ALL"]]$SamplingGear=="GOV"] <- "TRAWL"
BIOLOGICAL_METADATA[["ALL"]]$SamplingGearCAT[BIOLOGICAL_METADATA[["ALL"]]$SamplingGear!="GOV"] <- "GRAB"
##


## BIOMASS dataset
BIOLOGICAL_BIOMASS <- list()
Type <- "DATASET_Biological"

for(name in list.files(file.path(DIR.SOURCE[[paste(1)]],Type,"BIOMASS")))
{BIOLOGICAL_BIOMASS[[name]] <- read.csv2(file=file.path(DIR.SOURCE[[paste(1)]],Type,"BIOMASS",name),header=T,sep=";",dec = ".",colClasses = "character")

# Add genus level
BIOLOGICAL_BIOMASS[[name]]$Genus <- sub(" .*","",BIOLOGICAL_BIOMASS[[name]]$Taxa)
##


rm(name)}

for(name in names(BIOLOGICAL_BIOMASS))
{
if(length(BIOLOGICAL_BIOMASS[["ALL"]])==0)
{BIOLOGICAL_BIOMASS[["ALL"]] <- BIOLOGICAL_BIOMASS[[name]]}else{
TEMP <- BIOLOGICAL_BIOMASS[[name]][,intersect(colnames(BIOLOGICAL_BIOMASS[["ALL"]]),colnames(BIOLOGICAL_BIOMASS[[name]]))]
BIOLOGICAL_BIOMASS[["ALL"]] <- BIOLOGICAL_BIOMASS[["ALL"]][,colnames(TEMP)]

BIOLOGICAL_BIOMASS[["ALL"]] <- rbind(BIOLOGICAL_BIOMASS[["ALL"]],TEMP[,match(colnames(BIOLOGICAL_BIOMASS[["ALL"]]),colnames(TEMP))])
}}


       
##











