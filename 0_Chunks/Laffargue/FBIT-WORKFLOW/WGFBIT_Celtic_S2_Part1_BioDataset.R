#####################################
######### WGFBIT - October 2019 ##########
##################################### >)))°>
######## SCRIPT dataset compilation - Celtic/BoB Eco-region (WGFBIT)
##
### R version 3.6.2

#################################
##
### PART 1 - BIODATASET
##
#################################

 ### 2.2 - BIOLOGICAL DATASET
  
  ## 2.2.1 - METADATA

BiodatasetMetadata <- read.csv2(file.path(DIR.SOURCE[[paste(3)]],FILE.NAME[["BioMetadata"]]),stringsAsFactors=F)

# NEW files source

SEL.NEWfileUPDATE <- T

if(SEL.NEWfileUPDATE){

## Load previously created file
if(file.exists(file.path(file.path(DIR.OUPUT[[paste(1)]]),paste(SEL.ScriptName,"BIOdataset.RData",sep="")))){
load(file=file.path(file.path(DIR.OUPUT[[paste(1)]]),paste(SEL.ScriptName,"BIOdataset.RData",sep="")))}else{BIOdataset <- list()}


for(SURVEY in BiodatasetMetadata$SOURCEname){

print(paste("Processing data for survey:", SURVEY))

## 2.1 - Check for standard field names
# CONVERT / FORMAT biological dataset 
TEMP.LIST <- list()
for(FIELD in c(colnames(BiodatasetMetadata)[grep("FIELDname",colnames(BiodatasetMetadata))])){
TEMP.LIST <- rbind(TEMP.LIST,c(unlist(strsplit(FIELD,"_"))[2],BiodatasetMetadata[BiodatasetMetadata$SOURCEname==SURVEY,][[FIELD]]))}
BIOdataset[["SOURCE"]][[SURVEY]][["StandardName"]] <- TEMP.LIST
rm(TEMP.LIST)
#

## 2.2 - Load SOURCE files
TEMP.SEP <- BiodatasetMetadata[BiodatasetMetadata$SOURCEname==SURVEY,]$FileSeparator
for(FILE in  colnames(BiodatasetMetadata)[grep("FileName",colnames(BiodatasetMetadata))]){
TEMP.metadata <- BiodatasetMetadata[[FILE]][BiodatasetMetadata$SOURCEname==SURVEY]
TEMP.PATH <- file.path(DIR.SOURCE[[paste(3)]],BiodatasetMetadata$TYPEbiocomponent[BiodatasetMetadata$SOURCEname==SURVEY],SURVEY,TEMP.metadata)
if(file.exists(TEMP.PATH)){
BIOdataset[["SOURCE"]][[SURVEY]][[unlist(strsplit(FILE,"_"))[2]]] <-  read.csv2(TEMP.PATH,sep=TEMP.SEP,stringsAsFactors=F)}else{
BIOdataset[["SOURCE"]][[SURVEY]][[unlist(strsplit(FILE,"_"))[2]]] <- NULL}

rm(TEMP.metadata,TEMP.PATH)}
rm(TEMP.SEP)                                                                                





## 2.2.2 - Specific DATASET adds/modifications
 
   ## 2.2.2.1 - EVHOE
   
   # adapt bio variable (number, biomass)
   if(SURVEY=="EVHOE")
   {
   
   # SELECTION YEARS
   # SEL.YEARS <- c(2008:max(BIOdataset[["SOURCE"]][[SURVEY]][["Station"]]$Annee))
   SEL.YEARS <- c(2012:max(BIOdataset[["SOURCE"]][[SURVEY]][["Station"]]$Annee))
   SEL.YEARSexlusion <- c(2017)
   BIOdataset[["SOURCE"]][[SURVEY]][["Station"]] <- BIOdataset[["SOURCE"]][[SURVEY]][["Station"]][is.na(match(BIOdataset[["SOURCE"]][[SURVEY]][["Station"]]$Annee,SEL.YEARS))==F,]
   BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]] <-   BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]][is.na(match(BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]]$Trait,unique(BIOdataset[["SOURCE"]][[SURVEY]][["Station"]]$Trait)))==F,]
   #
  
   
   TEMP.VAR <- merge(BIOdataset[["SOURCE"]][[SURVEY]][["Station"]],BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]],by=c("Campagne", "Annee", "Trait"))
   
   # Compute density variable from hauling surface 
  for(VAR in c("Nombre","Poids")){TEMP.VAR[[paste(VAR,"1000m2",sep="")]] <- as.numeric(TEMP.VAR[[VAR]])/(as.numeric(TEMP.VAR$SurfaceBalayee)*1e+03)}
  BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]] <- TEMP.VAR
  rm(TEMP.VAR)
   #
   
   #BIOdataset[["SOURCE"]][[SURVEY]][["Station"]]
   
   # Taxa
   BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]]$TaxaSOURCE <- as.vector(BIOdataset[["SOURCE"]][[SURVEY]][["SpeciesList"]]$Nom_Scientifique[match(BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]]$Espece  ,BIOdataset[["SOURCE"]][[SURVEY]][["SpeciesList"]]$Code_Rubin)])
   
   # Controle absence REFTAX
   TEMP.MISSING <- unique(BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]]$Espece[is.na(BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]]$TaxaSOURCE)])

if(length(TEMP.MISSING)>0){
print(paste("WARNING !!!", "missing sp in species list",paste(TEMP.MISSING,collapse=",")))
}

  # List of species to correct manually
  SEL.SPcorr <- rbind(c("Aporrhais serresianus", "Aporrhais serresiana"))

  for(i in 1:nrow(SEL.SPcorr)){
  TEMP0 <- which(BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]]$TaxaSOURCE==SEL.SPcorr[,i][1])
  if(length(TEMP0)>0)
  {BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]]$TaxaSOURCE[TEMP0] <- SEL.SPcorr[,i][2]}
  rm(TEMP0)}

  } # END SPECIFIC CORRECTIONS EVHOE
##


 
   ## 2.2.2.2 - JNCC
   if(SURVEY=="JNCC")
   {
   
   ## Corrections species name (Catch table)
           ### VERIFIER /CORRIGER CI-DESSOUS !!!
                      ### VERIFIER /CORRIGER CI-DESSOUS !!!
                                 ### VERIFIER /CORRIGER CI-DESSOUS !!!
                                            ### VERIFIER /CORRIGER CI-DESSOUS !!!
                                                       ### VERIFIER /CORRIGER CI-DESSOUS !!!
                                                       
   TEMP <- grep("\\(",BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]]$Taxa)
      for(CORR in paste("\\(",unique(sub(".*\\(","",BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]]$Taxa[TEMP])),sep=""))
   {BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]]$Taxa <- sub(CORR, "",BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]]$Taxa)}
   rm(TEMP)

   } # END SPECIFIC CORRECTIONS JNCC


      # adapt bio variable (number, biomass)
   if(SURVEY=="DEMERSALES"){
  
  
  # Demersales adaptation script from raw dataset (Jose)
  if(F){}
  
  # "Raw dataset" from Jose adaptation script: DataDEM_2013_2020.csv 
   
   }
   ######
   ####
      
   
## 2.2.3 - CONVERT / FORMAT biological dataset with standard field names
# Station, Lon, Lat, Species, Biomasse (standard/area), Abundance (standard/area)


                                # 2.2.3.1 - Convert field names
# Standard field names definition
NAMEStation.Standard <- c("Survey", "Year", "Station", "Latitude", "Longitude")
                                
for(FILE in names(BIOdataset[["SOURCE"]][[SURVEY]]))
{TEMP.COLmatch <- match(colnames(BIOdataset[["SOURCE"]][[SURVEY]][[FILE]]),BIOdataset[["SOURCE"]][[SURVEY]]$StandardName[,2])
if(length(na.omit(TEMP.COLmatch))>0)
{colnames(BIOdataset[["SOURCE"]][[SURVEY]][[FILE]])[is.na(TEMP.COLmatch)==F] <-  BIOdataset[["SOURCE"]][[SURVEY]]$StandardName[,1][na.omit(TEMP.COLmatch)]}
rm(FILE,TEMP.COLmatch)}
                                                   
                                # 2.2.3.2 - Aggregate / final formatted dataset 
                                TEMP.NAME <- unlist(BIOdataset[["SOURCE"]][[SURVEY]]$StandardName[,1])
                                TEMP1 <- BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]][,is.na(match(colnames(BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]]),TEMP.NAME))==F]
                                TEMP2 <- BIOdataset[["SOURCE"]][[SURVEY]][["Station"]][,is.na(match(colnames(BIOdataset[["SOURCE"]][[SURVEY]][["Station"]]),TEMP.NAME))==F]
                                
                                if(is.null(TEMP1)==F&is.null(TEMP2)==F){BIOdataset[["FORMATED"]][[SURVEY]][["DATA"]] <- merge(TEMP1,TEMP2,by=intersect(colnames(TEMP1),colnames(TEMP2)))
                                TEMP2[,NAMEStation.Standard[is.na(match(NAMEStation.Standard,colnames(TEMP2)))]] <- NA
                                BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]] <- TEMP2[,NAMEStation.Standard]}
                                
                                if(is.null(TEMP1)&is.null(TEMP2)==F){
                                print(paste("!!! WARNING !!! No biological/observation data for dataset", SURVEY))
                                TEMP2[,NAMEStation.Standard[is.na(match(NAMEStation.Standard,colnames(TEMP2)))]] <- NA
                                BIOdataset[["FORMATED"]][[SURVEY]][["DATA"]] <- TEMP2
                                BIOdataset[["FORMATED"]][[SURVEY]][["STATION"]] <- TEMP2[,c("Survey", "Year", "Station", "Latitude", "Longitude")]}
                                
                                if(is.null(TEMP2)&is.null(TEMP1)==F){
                                print(paste("!!! WARNING !!! No station data for dataset", SURVEY))
                                nrow(TEMP1)
                                BIOdataset[["FORMATED"]][[SURVEY]][["DATA"]] <- TEMP1}
                                                                     
                                
                                BIOdataset[["FORMATED"]][[SURVEY]][["WORMSSpeciesList"]] <- BIOdataset[["SOURCE"]][[SURVEY]][["SpeciesList"]]
                                rm(TEMP1, TEMP2)
                                #



 
 
 ## 2.2.4 - HARMONIZATION/CORRECTION OF SPECIES NAME 

    # 2.2.4.1 - CORRECTION SPECIES WORMS

if(length(BIOdataset[["SOURCE"]][[SURVEY]][["Catch"]])>0){

TEMP.PATH <- file.path(DIR.SOURCE[[paste(3)]],BiodatasetMetadata$TYPEbiocomponent[BiodatasetMetadata$SOURCEname==SURVEY],SURVEY,paste("SpeciesList_",SURVEY,".csv",sep=""))

if(SEL.WORMSnew==T|file.exists(TEMP.PATH)==F){

TEMP.SP1 <- unique(as.vector(BIOdataset[["SOURCE"]][[SURVEY]]$TaxaSOURCE))

TEMP.SP2 <- list()
TEMP.SP3 <- NULL

for(i in c(1:length(TEMP.SP1)))
{print(TEMP.SP1[i])
 skip_to_next <- T
tryCatch({TEMP.SP2[[i]] <- WORMSextract(ListSP=TEMP.SP1[i],SEL.marineOnly=F,verbose=F)
TEMP.SP3 <- unique(c(TEMP.SP3 ,colnames(TEMP.SP2[[i]])))}
,error = function(e) { skip_to_next <<- T})
if(skip_to_next){ next }
}

TEMP.SP6 <- NULL
for(i in c(1:length(TEMP.SP2)))
{
TEMP.SP4 <- TEMP.SP3[is.na(match(TEMP.SP3,names(TEMP.SP2[[i]])))]
if(length(TEMP.SP4)>0)
{TEMP.SP5 <- t(data.frame(rep(NA,length(TEMP.SP4)))) 
colnames(TEMP.SP5) <- TEMP.SP4
TEMP.SP2[[i]] <- as.data.frame(cbind(TEMP.SP2[[i]],TEMP.SP5))
rm(TEMP.SP5)}
TEMP.SP2[[i]] <- TEMP.SP2[[i]][,order(names(TEMP.SP2[[i]]))]
TEMP.SP6 <- rbind(TEMP.SP6,TEMP.SP2[[i]])

rm(TEMP.SP4)}

 
TEMP.SP6 <- merge(data.frame("TaxaSOURCE"=TEMP.SP1),TEMP.SP6, by.x="TaxaSOURCE",by.y="name",all.x=T)
colnames(TEMP.SP6)[colnames(TEMP.SP6)=="Taxa"] <- "TaxaSOURCE"
                                                
#BIOdataset[["SOURCE"]][[SURVEY]][["SpeciesList"]] <- cbind(TAXREF[["TUTTItaxa"]][match(TEMP.SP6$TaxaSOURCE,as.vector(TAXREF[["TUTTItaxa"]]$Nom_Scientifique)),],TEMP.SP6)
BIOdataset[["SOURCE"]][[SURVEY]][["SpeciesList"]] <- TEMP.SP6
BIOdataset[["FORMATED"]][[SURVEY]][["WORMSSpeciesList"]] <- BIOdataset[["SOURCE"]][[SURVEY]][["SpeciesList"]]
write.csv2(BIOdataset[["SOURCE"]][[SURVEY]][["SpeciesList"]],file.path(DIR.SOURCE[[paste(3)]],BiodatasetMetadata$TYPEbiocomponent[BiodatasetMetadata$SOURCEname==SURVEY],SURVEY,paste("SpeciesList_",SURVEY,".csv",sep="")))



rm(TEMP.SP6,TEMP.SP1,TEMP.SP2,TEMP.SP3)}else{BIOdataset[["FORMATED"]][[SURVEY]][["WORMSSpeciesList"]] <- read.csv2(TEMP.PATH)}

rm(TEMP.PATH)

         # 2.5.2 - Check for missing Worms valid_name

TEMP.MISSING <- BIOdataset[["FORMATED"]][[SURVEY]][["WORMSSpeciesList"]][is.na(BIOdataset[["FORMATED"]][[SURVEY]][["WORMSSpeciesList"]]$valid_name),]

if(nrow(TEMP.MISSING)>0){print(paste("WARNING MISSING Worms valid names for following taxa"))
                         print(paste(TEMP.MISSING$TaxaSOURCE,collapse=","))}

#
rm(TEMP.MISSING)}




# Add valid name & AphiaID to dataset
for(VAR in c("valid_name","valid_AphiaID")){
BIOdataset[["FORMATED"]][[SURVEY]][["DATA"]][[VAR]] <- BIOdataset[["FORMATED"]][[SURVEY]][["WORMSSpeciesList"]][[VAR]][match(BIOdataset[["FORMATED"]][[SURVEY]][["DATA"]]$TaxaSOURCE, BIOdataset[["FORMATED"]][[SURVEY]][["WORMSSpeciesList"]]$TaxaSOURCE)] 
rm(VAR)}
#

rm(SURVEY)} # END LOOP by SURVEY

 
   
 
         # 2.2.5 -  Species list filtering  
         # Indicate unwanted species taxo class, order or taxa names (valid_names)

  # 2.2.5.1 - configure the species filter
   # Fish
   # "swimming" Cephalopods
   # Gelatinous
   
  SELECT.SpeciesFilter <- list()
  SELECT.SpeciesFilter[["class"]] <- c("Chlorophyta incertae sedis","Actinopterygii","Elasmobranchii","Holocephali","Petromyzonti","Myxini","Scyphozoa","Cephalopoda")
  SELECT.SpeciesFilter[["order"]] <- c("Myopsida","Oegopsida")
  SELECT.SpeciesFilter[["valid_name"]]  <- c("")
  #

  ## Specific filter (e.g. excluding highly motile / swimming species)

  #### TO COMPLETE
    #### TO COMPLETE
      #### TO COMPLETE
        #### TO COMPLETE

  ##


  # 2.2.5.2 - Apply filter by SURVEY
for(SURVEY in names(BIOdataset[["FORMATED"]])){
  BIOdataset[["FORMATED"]][[SURVEY]][["FilteredSpeciesList"]] <- BIOdataset[["FORMATED"]][[SURVEY]][["WORMSSpeciesList"]]
  for(FILTER in names(SELECT.SpeciesFilter)){
  TEMP.FILTER <- match(BIOdataset[["FORMATED"]][[SURVEY]][["FilteredSpeciesList"]][[FILTER]],SELECT.SpeciesFilter[[FILTER]])
  if(length(na.omit(TEMP.FILTER))>0){
  BIOdataset[["FORMATED"]][[SURVEY]][["FilteredSpeciesList"]] <- BIOdataset[["FORMATED"]][[SURVEY]][["FilteredSpeciesList"]][is.na(TEMP.FILTER),]}
  rm(TEMP.FILTER)}
  
  BIOdataset[["FORMATED"]][[SURVEY]][["DATA"]] <- BIOdataset[["FORMATED"]][[SURVEY]][["DATA"]][is.na(match(BIOdataset[["FORMATED"]][[SURVEY]][["DATA"]]$TaxaSOURCE,BIOdataset[["FORMATED"]][[SURVEY]][["FilteredSpeciesList"]]$TaxaSOURCE))==F,]
  
  # 2.2.5.3 - Filter redundant station data
  BIOdataset$FORMATED[[SURVEY]]$STATION <- unique(BIOdataset$FORMATED[[SURVEY]]$STATION)
  
  # Add survey name
  BIOdataset$FORMATED[[SURVEY]]$STATION$Survey <- SURVEY 
  
  
  
rm(SURVEY)} # End loop by SURVEY
  #
 
 # Save R object
TEMP.DIR <- file.path(DIR.OUPUT[[paste(1)]])
if(file.exists(TEMP.DIR)==F){dir.create(TEMP.DIR)}
save(BIOdataset,file=file.path(TEMP.DIR,paste(SEL.ScriptName,"BIOdataset.RData",sep="")))
rm(TEMP.DIR)
  #
}else{load(file=file.path(file.path(DIR.OUPUT[[paste(1)]]),paste(SEL.ScriptName,"BIOdataset.RData",sep="")))}
##################################
##################################
##################################

                               
                   
  
