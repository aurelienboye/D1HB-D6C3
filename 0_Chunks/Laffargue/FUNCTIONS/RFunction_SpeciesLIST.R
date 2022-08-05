### SCRIPT DE CREATION LISTE D'ESPECE A PARTIR DES DONNEES DE CAMPAGNE
## Plaffarg 2018-2022
## MaJ 01/2022
## CONTROLE DES DOUBLONS
# ELIMINER les redondances (établissement de la liste "minimale" de taxa)
# OPTION POUR CHAQUE LISTE PAR ANNEE (niveau de détermination change annuellement !!!)
# EXTRACTION de la liste des Gelatineux
##

## !!! WARNING !!!
## NEED WORMSextract function to be loaded
#DIRECTORY.ScriptSource <- "E:/PlaffIFREMER/R-SCRIPT/FUNCTIONS"
#source(file.path(DIRECTORY.ScriptSource,"RFunction_WORMSextract.R"))
##


# Data : données brutes en entrée avec au moins une liste d'espèce ou d'AphiaID, class, family ...
# éventuellement un champ année
# Data <- DATASET[["RSUFI"]][["1987-2019"]][["SERIE_REF"]][["TOTAL"]][["Captures"]]
# SELyear # distinction par année ou non
# SpeciesFieldName # Indication du champ nom d'espèce pas nécessaire si data est un vecteur de taxa, "OBLIGATOIRE si data est un data frame avec plusieurs champs


# YearFieldName # Indication du champ Année OPTIONNEL
# SELtaxoWorms (T or F) # Execute contrôle taxo worms Exécuté si pas données taxo complète renseignée (ou champ taxo attendu différent du fichier en entrée)





SpeciesLIST <- function(Data=NA,
                        SpeciesFieldName=NA,
                        TaxoFieldName=c("kingdom"=NA, "phylum"=NA, "class"=NA, "order"=NA, "family"=NA, "genus"=NA, "Species"=NA),
                        YearFieldName=NA,SELyear=F,SELtaxoWorms=F,NAomit=F)
{

 if(is.data.frame(Data)){if(nrow(Data)==0|is.na(SpeciesFieldName)){
 print("Missing Data or SpeciesFieldName parameters to execute SpeciesLIST function")
 stopifnot(nrow(Data)>0&is.na(SpeciesFieldName)==F)}}
 
 if(is.vector(Data) & length(Data)==0)
 {print("Missing Data or SpeciesFieldName parameters to execute SpeciesLIST function")
 stopifnot(nrow(Data)>0&is.na(SpeciesFieldName)==F)
 }
 
 


  # TEST
if(F)
{
load(file=file.path("E:/PlaffIFREMER/R-SCRIPT/FUNCTIONS","TESTlistSP.RData"))
Data=TESTspList
SpeciesFieldName="Nom_Scientifique" # "valid_name"
YearFieldName="Annee"
SELyear <- F
SELtaxoWorms <- T
}
 #

## Taxonomic table reference names (1 taxa name and n taxonomic fields)
if(is.vector(Data)){Data <- as.data.frame(Data)
SpeciesFieldName <- "Nom_Scientifique"
colnames(Data) <- "Nom_Scientifique"}

if(NAomit){Data <- Data[is.na(Data[,SpeciesFieldName])==F,]}

colnames(Data)[colnames(Data)==SpeciesFieldName] <- paste("SOURCE",SpeciesFieldName,sep="_")
SpeciesFieldName <- paste("SOURCE",SpeciesFieldName,sep="_")
LISTrefTaxo <-  c(SpeciesFieldName, names(TaxoFieldName))
#

if(length(which(is.na(TaxoFieldName)))>0|SELtaxoWorms==T)
{

 if(exists('WORMSextract', mode='function')==F)
 {print("WORMSextract function is needed")}
 stopifnot(exists('WORMSextract', mode='function'))

LISTsp <- WORMSextract(unique(Data[,SpeciesFieldName]))
colnames(LISTsp)[colnames(LISTsp)=="name"] <- SpeciesFieldName
Data <- merge(Data,LISTsp,by=SpeciesFieldName,suffixes = c(".x",""))


}else{
colnames(Data)[is.na(match(colnames(Data),TaxoFieldName))==F] <- names(TaxoFieldName)

if(length(which(is.na(match(names(TaxoFieldName),colnames(Data)))))>0)
{

 if(exists('WORMSextract', mode='function')==F)
 {print("WORMSextract function is needed")}
 stopifnot(exists('WORMSextract', mode='function'))

LISTsp <- WORMSextract(unique(Data[,SpeciesFieldName]))
Data <- merge(Data,LISTsp,by=SpeciesFieldName,suffixes = c(".x",""))
}
}


# Which species name's group == espece and which taxa name's group >= genre
if(SELyear & length(which(colnames(Data)==YearFieldName))>0)
{
TEMP.YEARsp <- unique(Data[,c(YearFieldName,SpeciesFieldName)])
TEMP.YEARsp <- table(TEMP.YEARsp)
TEMP.YEARsp <- rbind(TEMP.YEARsp,"ALL"=colSums(TEMP.YEARsp))}else{
if(SELyear){print("Check the year field name")}
TEMP.YEARsp <- t(data.frame("ALL"=rep(1,length(unique(Data[,c(SpeciesFieldName)])))))
colnames(TEMP.YEARsp) <- unique(Data[,c(SpeciesFieldName)])
}

LIST.SP <- list()
#LIST.SPECIES[["REDUNDANT"]] <- list()

for(YEAR in rownames(TEMP.YEARsp))
{
print(YEAR)
TEMP.YEARspSel <- TEMP.YEARsp[which(rownames(TEMP.YEARsp)==YEAR),]

TEMP.YEARspSel <- TEMP.YEARspSel[TEMP.YEARspSel>0]

#TEMP.SPECIES <- dataset.taxref1$C_Perm[match(names(TEMP.YEARspSel),dataset.taxref1$C_VALIDE)]
#TEMP.TAXREF  <- LIST.SPECIES[["ALL"]][match(TEMP.SPECIES,LIST.SPECIES[["ALL"]]$CPERM),]
#TEMP.ERROR <- TEMP.TAXREF[is.na(TEMP.SPECIES),]

TEMP.TAXREF <- Data[match(names(TEMP.YEARspSel),Data[,SpeciesFieldName]),]

## Recherche du niveau taxonomique du taxa considéré
# Recherche des possibles redondances
#
TEMP.LEVELnum <- NULL

for(i in 1:nrow(TEMP.TAXREF))
{TEMP0 <- which(is.na(TEMP.TAXREF[,c(LISTrefTaxo)][i,])==F)
if(length(TEMP0)>0)
{TEMP.LEVELnum <- c(TEMP.LEVELnum, max(TEMP0)-1)}else{
TEMP.LEVELnum <- c(TEMP.LEVELnum, ncol(TEMP.TAXREF[,c(LISTrefTaxo)])-1)}
rm(TEMP0)}

#TEMP.LEVELnum <- c(TEMP.LEVELnum, min(which(is.na(TEMP.TAXREF[i,][,c(LISTrefTaxo)][-c(1:(which(colnames(TEMP.TAXREF)=="Species")-1))])==F)))}
#TEMP.LEVELnum <- TEMP.LEVELnum+which(colnames(TEMP.TAXREF)=="Species")-1

TEMP.LEVELname <- colnames(TEMP.TAXREF[,c(LISTrefTaxo)])[TEMP.LEVELnum+1]

TEMP.isSPECIES <- which(TEMP.LEVELname=="Species")
TEMP.isnotSPECIES <- which(TEMP.LEVELname!="Species")

## Liste totale de taxa potentiels (pour tous les niveaux taxonomiques) à partir de la liste d'espèce
TEMP.isSPECIEStotLIST <- NULL
for(i in 1:nrow(TEMP.TAXREF[,c(LISTrefTaxo)][TEMP.isSPECIES,]))
{TEMP.isSPECIEStotLIST <- c(TEMP.isSPECIEStotLIST,na.omit(as.vector(unlist(TEMP.TAXREF[,c(LISTrefTaxo)][TEMP.isSPECIES,][i,]))))}
TEMP.isSPECIEStotLIST <- sort(unique(TEMP.isSPECIEStotLIST))
##

## Extraction des taxa potentiellement redondants
LIST.TAXAredundant <- NULL
for(i in 1:length(TEMP.isnotSPECIES))
{TEMP2 <- TEMP.TAXREF[,c(LISTrefTaxo)][TEMP.isnotSPECIES[i],]
TEMP0 <- which(is.na(TEMP2)==F) 
if(length(TEMP0)>0)
{TEMP3 <- TEMP2[[max(TEMP0)]] 
if(length(which(TEMP.isSPECIEStotLIST==TEMP3))>0)
 {LIST.TAXAredundant <- c(LIST.TAXAredundant,TEMP.isnotSPECIES[i])}
 rm(TEMP3)}
 #TEMP3 <- TEMP2[[min(which(is.na(TEMP2[-c(1:(which(colnames(TEMP.TAXREF[,c(LISTrefTaxo)])=="Species")-1))])==F))+(which(colnames(TEMP.TAXREF[,c(LISTrefTaxo)])=="Species")-1)]]
 rm(TEMP2, TEMP0)}
##

TEMP.TAXREF <- TEMP.TAXREF[,c("valid_AphiaID","valid_name","valid_authority",LISTrefTaxo)]
TEMP.TAXREF2 <- TEMP.TAXREF

LIST.SP[[paste(YEAR)]][["COMPLETE"]] <- TEMP.TAXREF
LIST.SP[[paste(YEAR)]][["REDUNDANT"]] <- TEMP.TAXREF[LIST.TAXAredundant,]


if(length(LIST.TAXAredundant)>0)
{TEMP.TAXREF2 <- TEMP.TAXREF2[-LIST.TAXAredundant,]}


####
## EXTRACTION FISH ##
# Following classes are selected: "Actinopteri", "Elasmobranchii","Holocephali","Myxini","Petromyzonti","Cephalaspidomorphi"
SEL.FISH <- c("Actinopteri", "Elasmobranchii","Holocephali","Myxini","Petromyzonti","Cephalaspidomorphi")
MISSINGclass <- SEL.FISH[is.na(match(SEL.FISH,TEMP.TAXREF2$class))]
if(length(MISSINGclass)>0){print(paste("WARNING: CHECK MISSING/ERRONEOUS CLASS for FISH:",MISSINGclass))}
LIST.SP[[paste(YEAR)]][["FISH"]]   <- TEMP.TAXREF2[is.na(match(TEMP.TAXREF2$class,SEL.FISH))==F,]

####

####
## EXTRACTION INVERTEBRATES
LIST.SP[[paste(YEAR)]][["INVERTEBRATE"]] <- TEMP.TAXREF2[is.na(match(TEMP.TAXREF2$class,SEL.FISH)),]
##
####

####
## EXTRACTION CEPHALOPODS
SEL.CEPHALOPOD <- "Cephalopoda"
LIST.SP[[paste(YEAR)]][["CEPHALOPOD"]]   <- TEMP.TAXREF2[is.na(match(TEMP.TAXREF2$class,SEL.CEPHALOPOD))==F,]
##
####


####
## EXTRACTION GELATINOUS
if(T){
   TEMP.Gelatinous <-  LIST.SP[[paste(YEAR)]][["INVERTEBRATE"]]
   
   TEMP.Gelatinous <-  TEMP.Gelatinous[is.na(match(TEMP.Gelatinous$phylum,c("Cnidaria", "Ctenophora","Mollusca")))==F,]
   TEMP.Gelatinous <-  TEMP.Gelatinous[which(TEMP.Gelatinous$class!="Anthozoa"),]
   TEMP.Gelatinous <-  TEMP.Gelatinous[which(TEMP.Gelatinous$order!="Anthoathecata"),]

   TEMPgel1 <- unique(c(which(TEMP.Gelatinous$phylum=="Mollusca"&TEMP.Gelatinous$order!="Pteropoda"),
                 which(TEMP.Gelatinous$phylum=="Mollusca"&is.na(TEMP.Gelatinous$order)==T),
                 which(TEMP.Gelatinous$order=="Leptothecata"&TEMP.Gelatinous$family!="Aequoreidae")))
   TEMP.Gelatinous <- if(length(TEMPgel1)>0){TEMP.Gelatinous[-TEMPgel1,]}
   LIST.SP[[paste(YEAR)]][["GELATINOUS"]] <- TEMP.Gelatinous
   
   LIST.SP[[paste(YEAR)]][["BENTHOS"]] <- LIST.SP[[paste(YEAR)]][["INVERTEBRATE"]][is.na(match(LIST.SP[[paste(YEAR)]][["INVERTEBRATE"]]$valid_name,LIST.SP[[paste(YEAR)]][["GELATINOUS"]]$valid_name)),]
}
####


####



##
####

####
## EXTRACTION DATRAS ("Commercial species")

## EXTRACTION COMMERCIAL INVERTEBRATES
# Specific commercial invertebrates list of species 
SEL.COMMinvertSP <- c("Palinurus elephas", "Homarus gammarus", "Nephrops norvegicus", "Maja brachydactyla", "Pecten maximus", "Aequipecten opercularis", "Necora puber", "Crangon crangon", "Palaemon serratus", "Buccinum undatum")

LISTspCom <- WORMSextract(SEL.COMMinvertSP)

TEMP.MATCH <- match(LISTspCom$valid_name,LIST.SP[[paste(YEAR)]][["COMPLETE"]]$valid_name) 
LISTspCom <- LIST.SP[[paste(YEAR)]][["COMPLETE"]][na.omit(TEMP.MATCH),]

if(length(which(is.na(TEMP.MATCH)))>0){
MISSINGspCOMinvert <- SEL.COMMinvertSP[is.na(TEMP.MATCH)]
print(paste("WARNING following commercial invertebrates missing in dataset:",paste(MISSINGspCOMinvert,collapse="/")))}

# "Commercial" species "DATRAS"
LIST.SP[[paste(YEAR)]][["DATRAS"]] <- rbind(LIST.SP[[paste(YEAR)]][["FISH"]],LIST.SP[[paste(YEAR)]][["CEPHALOPOD"]],LISTspCom)

##
####

rm(TEMP.isnotSPECIES,TEMP.isSPECIEStotLIST,TEMP.isSPECIES,TEMP.LEVELname,TEMP.LEVELnum)
rm(TEMP.YEARspSel,YEAR,TEMP.TAXREF2,TEMP.TAXREF)} # end loop by year
rm(TEMP.YEARsp)
return(LIST.SP)}




