##############################################################################
##############################################################################
############ Function to compute SE*RE indicator ("SEREfishPres")
############ adapted from Script WGBIODIV
############ ABA, PLAFF Fevrier-juin 2018
############ MODIF PLAFF Août 2019
############ MODIF PLAFF Février 2021
############ R 3.6.2
##############################################################################
##############################################################################

## dataset :  matrice station/espèces avec 
    # Nom scientifique de l'Espece 
    # variable numérique (abondance, biomasse, densité .... )


#####  FUNCTION

SEREfishPres <- function(dataset,SELtrans="log1p",TYPEoccu=F,SEL.CompleteTraits=T,SEL.WormsCorr=F,FILTERfish=T,FILTERcephalo=T,ScriptDIR="E:/PlaffIFREMER/R-SCRIPT/FUNCTIONS") {

SEL.StdNames <- c("Station","Taxa","Freq") # Standard names for data matrix and outputs
 
## 1 - SPECIFIC FUNCTIONS & DATASETS
    
   ## Fonctions
ScriptDIR <- file.path("E:/PlaffIFREMER/R-SCRIPT/FUNCTIONS")
source(file.path(ScriptDIR,"RFunction_WORMSextract.R"))
source(file.path(ScriptDIR,"RFunction_SpeciesLIST.R"))

   ## 1.2 - TRAITS standard table

    # Loading ### CORRIGER / COMPLETE to be included into function !!!!!
    ## Source traits table from Beauchard et al. xxx ("Tables_tableS1.csv")
    
    traits <- read.csv2(file=file.path(ScriptDIR,"RFunction_SEREfishPres.csv"),stringsAsFactors=F)
    
    # Mise à jour taxo table traits
    if(F){
    TEMP <- WORMSextract(as.vector(traits$Taxa))
    traits <- merge(traits,TEMP,by.x=colnames(traits)[1],by.y=colnames(TEMP)[1],all.x=T)
    write.csv2(traits, file=file.path(ScriptDIR,"SEREfishPres.csv"))
    rm(TEMP)}
    #

    # Standard traits code & score
    LIST.TRAITS <- data.frame("NAME"=c("BL","FR","BD","MO","AM","LS","OT","OS"),
                             "ScoreRangeMIN"=c(1,1,1,1,1,1,1,1),
                             "ScoreRangeMAX"=c(5,3,4,4,3,4,3,4))
                             rownames(LIST.TRAITS) <- LIST.TRAITS[,1]


   ## Apply standard fields names to bio dataset 
   TEMP.dimnames <- dimnames(dataset)
   names(TEMP.dimnames) <- SEL.StdNames[c(1,2)]
   dimnames(dataset) <- TEMP.dimnames
   ##
   
   
   ## Biological matrix transformation
   if(SELtrans=="log1p")
   {dataset <- log1p(dataset)}
   ##

   ## Filtre des espèces sans observation
   dataset <- dataset[,colSums(dataset)!=0]
   ##

   ## Filtre des stations sans observation
   LIST.STATIONempty <- unique(row.names(dataset)[rowSums(dataset)==0],row.names(dataset)[rowSums(dataset)==0])
   dataset <- dataset[is.na(match(row.names(dataset),LIST.STATIONempty)),]
   ##


   ## CREATION LISTE ESPECES
   ListSp <- as.data.frame(as.vector(colnames(dataset)))

   ## Ajout/correction information taxinomique
   if(SEL.WormsCorr|SEL.CompleteTraits){
   print("Complete taxonomic data from Worms")
   ListSp <- WORMSextract(as.vector(colnames(dataset)))}
   
   colnames(ListSp)[1] <- "Taxa"
   ##
   


   ## Filter Fish species and Cephalopods
   if(FILTERfish){
      LIST.FISH <- as.vector(na.omit(ListSp$Taxa[ListSp$Infraphylum=="Gnathostomata"]))
      dataset <- dataset[,is.na(match(colnames(dataset),LIST.FISH))]}
   
   
   if(FILTERcephalo){
      LIST.CEPHALO <- as.vector(na.omit(ListSp$Taxa[ListSp$Class=="Cephalopoda"]))
      dataset <- dataset[,is.na(match(colnames(dataset),LIST.CEPHALO))]}

   ##

 
   ## Complete traits matrix with selected level of aggregation
   
          TEMP.MATCHtraits <- match(ListSp$Taxa,traits[,"Taxa"])
          TEMP.TRAITS0 <- traits[,c("Taxa",as.vector(LIST.TRAITS[,1]))][TEMP.MATCHtraits,]
          colnames(TEMP.TRAITS0)[1] <- "TaxaMATCH"
          TEMP.TRAITS0 <- data.frame("Taxa"=ListSp$Taxa,TEMP.TRAITS0,stringsAsFactors =F)

if(SEL.CompleteTraits){
traitsSupLevel <- list()

   #if(length(traits$Genus)==0){traits$Genus <- sub(" .*","",as.vector(traits$Taxa))}
    
    # 3.2.1.2.2 - Complément au niveau genre
    # attribution du code Genre s'il existe dans la table traits, mediane des traits (ou valeur sup si non entier)
    
       #traits$CODE[is.na(traits$CODE)] <- spLIST$CODE[is.na(spLIST$Species)][match(traits$Genus[is.na(traits$CODE)],spLIST$Genus[is.na(spLIST$Species)])]           
       #TEMP.MATCHtraitsGenus <- match(spLIST$Species,traits$Taxa)

       # Generation table traits niveau du "valid_name", du genre et de la famille
          # Sélection niveau Taxa source
           #
          
       traitsSupLevel <- list()
       for(LEVEL in c("valid_name","Genus","Family")){
       traitsSupLevel[[LEVEL]] <- NULL
       for(TAXA in unique(traits[[LEVEL]])){
       TEMP0 <- which(traits[[LEVEL]]==TAXA)
       TEMP1 <- unique(traits[TEMP0,][,as.vector(LIST.TRAITS[,1])])
       TEMP1 <- ceiling(as.numeric(lapply(TEMP1,FUN="median")))
       traitsSupLevel[[LEVEL]] <- rbind(traitsSupLevel[[LEVEL]],c(TAXA,TEMP1))             
       rm(TAXA,TEMP0,TEMP1)}
       traitsSupLevel[[LEVEL]] <- traitsSupLevel[[LEVEL]][is.na(traitsSupLevel[[LEVEL]][,1])==F,]
       colnames(traitsSupLevel[[LEVEL]]) <- c("TaxaMATCH",as.vector(LIST.TRAITS[,1]))
  
   
   
   # Complément au niveau "LEVEL"
   TEMP.MATCHtraits2 <- match(ListSp[[LEVEL]],traitsSupLevel[[LEVEL]])
   TEMP.NAtest <- which(is.na(TEMP.TRAITS0[,"TaxaMATCH"]))
   if(length(TEMP.NAtest)>0){
   TEMP.TRAITS0[TEMP.NAtest,][,c(2:ncol(TEMP.TRAITS0))] <- traitsSupLevel[[LEVEL]][TEMP.MATCHtraits2[TEMP.NAtest],]}

   rm(LEVEL)} 
       
   
       # Nouvelle table de traits (sélection + aggrégation)
   ListSp <- data.frame(ListSp,TEMP.TRAITS0,stringsAsFactors=F)
   traitsSEL <- ListSp[,c("Species",c("valid_name","Genus","Family"),"Taxa","TaxaMATCH",as.vector(LIST.TRAITS[,1]))]
   traitsSEL <- traitsSEL[!is.na(ListSp[,as.vector(LIST.TRAITS[,1])][,1]),]
   rownames(traitsSEL) <- traitsSEL$Taxa  #
   #traitsSEL <- traitsSEL[!is.na(match(traitsSEL$Taxa,colnames(dataset))),]
   }else{
   ListSp <- data.frame(ListSp,TEMP.TRAITS0,stringsAsFactors=F)
   traitsSEL <- ListSp[,c("Taxa","TaxaMATCH",as.vector(LIST.TRAITS[,1]))]
   traitsSEL <- traitsSEL[!is.na(ListSp[,as.vector(LIST.TRAITS[,1])][,1]),]
   }

   
   ## FILTER dataset from available species in trait matrix
   MATRIXdatasetSEL <- dataset[,is.na(match(colnames(dataset), traitsSEL$Taxa))==F]
   
      # WARNING  !!!! NOUVELLES STATIONS SANS OBSERVATION APRES SELECTION !!!!
   LIST.STATIONempty2 <- unique(row.names( MATRIXdatasetSEL )[rowSums( MATRIXdatasetSEL )==0],row.names( MATRIXdatasetSEL )[rowSums( MATRIXdatasetSEL )==0])

   ##


 
   

   # CREATION Occurrence Matrix
   MATRIXspOCCU <- MATRIXdatasetSEL
   MATRIXspOCCU[MATRIXspOCCU>0] <- 1
   
   TEMP.occ <- aggregate(Freq~Taxa, data=data.frame(MATRIXspOCCU),FUN="sum")
   TEMP.occ[,"prop"] <- TEMP.occ[,"Freq"]/nrow(MATRIXspOCCU)
   #


   ### Test for species coverage from traits matrix
     ##
   LIST.SPECIEScovered <- colnames(MATRIXdatasetSEL[,na.omit(match(unique(traitsSEL$Taxa),colnames(MATRIXdatasetSEL)))])
   RES.SPcoverage <- rowSums(MATRIXdatasetSEL[,LIST.SPECIEScovered])/rowSums(dataset)
   TEMP.HISTcoverage <- hist(RES.SPcoverage, plot=F)
   TEMP.HISTcoverage$counts <- round(TEMP.HISTcoverage$counts/length(RES.SPcoverage)*100,2)
     ##
  
     ## Distribution taxo des espèces couvertes


     ##

   ###





  ### INDICATOR COMPUTATION

    ## FORMULA

   # SE
LIST.SEcomp <- c("BL*FR","BL*BD","FR*BD",
                 "SE~BL*FR*BD")
   # RE
LIST.REcomp <- c("MO*OT","MO*OS","MO*AM","OT*OS","OT*AM","OS*AM",
                 "MO*OT*OS","MO*OT*AM","MO*OS*AM","OT*OS*AM",
                 "MO*OT*OS*AM",
                 "RLS~LS-AM+1",
                 "RM~(AM^2)/(LS-AM+1)",
                 "RM*MO", "RM*OT","RM*OS",
                 "RM*MO*OT","RM*MO*OS","RM*OS*OT",
                 "RE~RM*MO*OT*OS")

   # Total formula
LIST.TotFormula <- c("SE.RE.a~SE+RE",
                     "SE.RE.m~SE*RE")



# Normalisation des scores de la table traits
traitsNorm <- traitsSEL
for(Tr in rownames(LIST.TRAITS))
{traitsNorm[[Tr]] <- as.numeric(traitsNorm[[Tr]])
traitsNorm[[Tr]] <- (traitsNorm[[Tr]]-LIST.TRAITS["ScoreRangeMIN"][Tr,])/(LIST.TRAITS["ScoreRangeMAX"][Tr,]-LIST.TRAITS["ScoreRangeMIN"][Tr,])
rm(Tr)}
#
traitsNorm[["SE"]] <- traitsNorm[["FR"]]*traitsNorm[["BD"]]*traitsNorm[["BL"]]
traitsNorm[["RE"]] <- (traitsNorm[["AM"]]/((traitsNorm[["LS"]] - traitsNorm[["AM"]]) +1)) * traitsNorm[["AM"]] * traitsNorm[["MO"]] * traitsNorm[["OT"]] * traitsNorm[["OS"]]
traitsNorm[["SE*RE"]] <- traitsNorm[["SE"]]*traitsNorm[["RE"]]
traitsNorm[["SE+RE"]] <- traitsNorm[["SE"]]+traitsNorm[["RE"]]

rownames(traitsNorm) <- traitsNorm$Taxa




      







   ## Summary taxa frequency
TEMPmean <- aggregate(Freq~Taxa,data=data.frame(MATRIXdatasetSEL),FUN=mean)
traitsNorm$Focc    <- TEMP.occ[,"prop"][match(traitsNorm$Taxa, TEMP.occ[,"Taxa"])]
traitsNorm$Fmean  <- TEMPmean[,"Freq"][match(traitsNorm$Taxa, TEMPmean[,"Taxa"])]
traitsNorm <- traitsNorm[is.na(traitsNorm$Focc)==F,] # WARNING SUPPRESSIONS DE NOUVELLES ESPECES A CE NIVEAU !!!!!
traitsNorm <- traitsNorm[traitsNorm$Focc>0,]
traitsNorm <- traitsNorm[order(traitsNorm$Focc,decreasing=T),]

rm(TEMP.occ,TEMPmean)





  ## 4.5 -  Pondération matrice de traits par variable biologique (abondance, biomasse, occurrence)
##

     # Calcul matrices valeurs relatives
MATRIXdatasetSELRel <- MATRIXdatasetSEL/rowSums(MATRIXdatasetSEL)
     #


datasetSEL <- data.frame(MATRIXdatasetSEL)
datasetSELRel <- data.frame(MATRIXdatasetSELRel)
datasetOCCU <- data.frame(MATRIXspOCCU)


#### Creation de l'objet résultat: traits sources et étapes d'aggrégations intermédiaires 
###
#

df.tr <- traitsNorm[,as.vector(LIST.TRAITS[,1])]
for(EXP in c(LIST.SEcomp,LIST.REcomp))
{df.tr[[gsub("\\*",".",gsub("~.*","",EXP))]] <- eval(parse(text=gsub(".*~","",EXP)), df.tr)
rm(EXP)}

# Normalisation
for(EXP in c(LIST.SEcomp,LIST.REcomp))
{df.tr[[gsub("\\*",".",gsub("~.*","",EXP))]]=(df.tr[[gsub("\\*",".",gsub("~.*","",EXP))]]-min(df.tr[[gsub("\\*",".",gsub("~.*","",EXP))]]))/(diff(range(df.tr[[gsub("\\*",".",gsub("~.*","",EXP))]])))
rm(EXP)}
#
##

for(EXP in LIST.TotFormula)
{df.tr[[gsub("\\*",".",gsub("~.*","",EXP))]] <- eval(parse(text=gsub(".*~","",EXP)), df.tr)

# Normalisation formule complète
df.tr[[gsub("\\*",".",gsub("~.*","",EXP))]]=(df.tr[[gsub("\\*",".",gsub("~.*","",EXP))]]-min(df.tr[[gsub("\\*",".",gsub("~.*","",EXP))]]))/(diff(range(df.tr[[gsub("\\*",".",gsub("~.*","",EXP))]])))
#
rm(EXP)}
##

datasetSEL=merge(datasetSEL[datasetSEL$Taxa%in%rownames(df.tr)==T,],df.tr,by.x="Taxa",by.y="row.names",all.x=T)
datasetSELRel=merge(datasetSELRel[datasetSELRel$Taxa%in%rownames(df.tr)==T,],df.tr,by.x="Taxa",by.y="row.names",all.x=T)
datasetOCCU=merge(datasetOCCU[datasetOCCU$Taxa%in%rownames(df.tr)==T,],df.tr,by.x="Taxa",by.y="row.names",all.x=T)

dataset.SEL <- list()
dataset.SEL[["abs"]] <- datasetSEL
dataset.SEL[["rel"]] <- datasetSELRel
dataset.SEL[["occu"]]<- datasetOCCU

TEMPselCOL <- match(colnames(df.tr),colnames(datasetSEL))

RESULTS.INDICATOR <- list()



  # 4.5.2 - CALCUL PONDERATION MATRICES (absolue et relatives)
for(NAME in names(dataset.SEL)){
print(paste("Indicator computation for", NAME))
resp <-  cbind(dataset.SEL[[NAME]][,"Freq"]*dataset.SEL[[NAME]][,TEMPselCOL] , Station=dataset.SEL[[NAME]]$Station)
rest=NULL
nCOL <- na.omit(match(colnames(df.tr),colnames(resp)))
for(i in sort(nCOL,decreasing=T)){
      ty <- droplevels(aggregate(resp[,i] ~ resp$Station, FUN=sum))
      rest= c(ty,rest)
}

RESULTS.INDICATOR[["Values"]][[NAME]] = matrix(unlist(rest[seq(2,length(rest),by=2)]),ncol=max(nCOL), dimnames=list(unlist(rest[1]),c(names(resp[,1:max(nCOL)])))) ;  #

rm(NAME,resp,rest,nCOL)}

## Calcul pour matrice occurrence relative
# Division de chaque valeur d'occurence par la somme des valeurs par station
# 
NAME <- "occuRel"
print(paste("Indicator computation for", NAME))
RESULTS.INDICATOR[["Values"]][[NAME]] <- RESULTS.INDICATOR[["Values"]][["occu"]]/aggregate(Freq~Station,data=dataset.SEL[["occu"]],FUN="sum")$Freq
rm(NAME)
#

RESULTS.INDICATOR[["dataset"]][["dataSelected"]] <- MATRIXdatasetSEL
RESULTS.INDICATOR[["dataset"]][["dataRelative"]] <- MATRIXdatasetSELRel
RESULTS.INDICATOR[["dataset"]][["dataOccurrence"]] <-MATRIXspOCCU
RESULTS.INDICATOR[["dataset"]][["traits"]] <- traitsNorm


return(RESULTS.INDICATOR)

} # End function



   #
  ###

#####










## EXEMPLE A PARTIR DONNEES EVHOE
  if(F){
      ## Série totale format RSUFI en local (AVEC Benthos, SANS SELECTION COSER)
 SEL.YEARS <- c(2008:2020)
 DATASET <- list()
      DIR.SOURCEevhoe <- "E:/PlaffIFREMER/2-EvHOE/RESEARCH_IN_PROGRESS/DATA_PAPER"
      SELECT.FILEPATHdataset <- file.path(DIR.SOURCEevhoe,"DATASET/RSUFI/RSUFI_Updated20210421")
      SEL.RSUFIfilesSOURCEstd <- c("Captures","Tailles","Traits","Strates")
      for(FILE in SEL.RSUFIfilesSOURCEstd){
      DATASET[["RAW"]][[FILE]] <- read.csv2(file.path(paste(SELECT.FILEPATHdataset),paste(FILE,".csv",sep="")),stringsAsFactors = F)
      rm(FILE)}

      # CREATE FIELDS DN & DW
      TEMP.DATA<- DATASET[["RAW"]]$Captures
      TEMP.DATA$SurfaceBalayee <- as.numeric(DATASET[["RAW"]][["Traits"]]$SurfaceBalayee[match(TEMP.DATA$Trait, DATASET[["RAW"]][["Traits"]]$Trait)])

      TEMP.DATA$DN <- as.numeric(TEMP.DATA$Nombre)/TEMP.DATA$SurfaceBalayee
      TEMP.DATA$DW <- as.numeric(TEMP.DATA$Poids)/TEMP.DATA$SurfaceBalayee
      TEMP.DATA <- TEMP.DATA[is.na(match(TEMP.DATA$Annee,SEL.YEARS))==F,]

      DATASET[["RAW"]]$Captures <-  TEMP.DATA

      # CREATE contingence DATA MATRIX

      for(VAR in c("DN","DW")){
     DATASET[["MATRIX"]][[paste("MATRIXsp",VAR,sep="")]] <- xtabs(as.formula(paste(VAR,"~",paste("Trait","Espece",sep="+"))),data=DATASET[["RAW"]]$Captures)
      rm(VAR)}


      ##

      MATRIXspDN <- DATASET[["MATRIX"]][[paste("MATRIXsp","DN",sep="")]]
      MATRIXspDW <- DATASET[["MATRIX"]][[paste("MATRIXsp","DW",sep="")]]
      
#

     RESULTS <- list()
     for(NAME in names(DATASET[["MATRIX"]])){
     RESULTS[[paste(NAME)]] <- SEREfishPres(DATASET[["MATRIX"]][[NAME]])
     rm(NAME)}

     }












##

































