## FUNCTION
## TITLE:Generate species list from WORMS database
## P.Laffargue 01/2020
## R3.5.1
## V0
####################################################################################################################################
####################################################################################################################################
# REFERENTIELS TAXO

### Generer liste taxa à partir de liste de taxa et DE WORMS sous format noms scientifiques ou Aphia ID
    # ListSP: vecteurs de noms de taxa (espece, genre, famille, etc ...)
    # SEL.taxosup: Ajoute compléments niveau taxo (pas forcément nécessaire, time consuming et nécessite autre library worrms)
    # FromAphiaID: T,F or AUTO à partir de l'APhiaID (T) ou à partir de taxon name (F) ou recherche auto du type de données en entrée (AUTO)
    # FuzzyMiss: recherche fuzzy pour les taxa non trouvés (que pour liste espèce par noms scientifiques)
    # na.omit: suppression des valeurs NA dans la liste d'espèces
    
#ListSP <- TEMP.SPlist1


### CORRECTIF A AJOUTER POUR SKIP ERROR
#for(i in c(1:length(TEMP.SP1)))
#{print(TEMP.SP1[i])
# skip_to_next <- F
#tryCatch({TEMP.SP2[[i]] <- WORMSextract(ListSP=TEMP.SP1[i],SEL.marineOnly=F,verbose=F)}
#,error = function(e) { skip_to_next <<- T})
#if(skip_to_next){ next }
#}
#

#library(worms)
#library(worrms)
#REFTAXworms <- wormsbymatchnames(data_species_GG_87_18)

#ListSP <- c("123432","234456")
#ListSP <- c("Genocidaris", "Leptochelastra toureti")
                                                                                                                                                                                                           
WORMSextract <- function(ListSP=NA,FromAphiaID="AUTO",SEL.taxosup=F, FuzzyMiss=T, na.omit=T,SEL.marineOnly=T,verbose=T,DIR.RESULTS=NA,bySP=F)
{

if(na.omit)
{ListSP <- as.vector(na.omit(ListSP))}

if(SEL.marineOnly){marine_only <- "true"}else{marine_only <- "false"}

options(warn=-1)
if(FromAphiaID=="AUTO")
{
if(length(na.omit(as.numeric(ListSP)))>0)
{FromAphiaID=T}else{FromAphiaID=F}
}
options(warn=0)

if(FromAphiaID==T){
REFTAXworms <- worms::wormsbyid(as.numeric(ListSP), verbose = TRUE, ids = T, sleep_btw_chunks_in_sec = 0.01)}else{

if(bySP){
for(SP in ListSP){
REFTAXworms <- NULL
tryCatch({REFTAXworms <- rbind(REFTAXworms, worms::wormsbynames(SP,ids = T, match = F, marine_only = marine_only,verbose=verbose))},error=function(e){print(paste("WORMS search tool error for taxa =",SP))})}}else{

# Retrieve AphiaID from taxon names (exact matching)
REFTAXworms <- worms::wormsbynames(as.character(ListSP), ids = T, match = F, marine_only = marine_only,verbose=verbose)}




# Ajout for nomatched species in the previous function from fuzzy match
if(FuzzyMiss)
{
print("Application of Fuzzy research for missing species")
TEMP0.nomatch <- which(is.na(REFTAXworms$match_type))
if(length(TEMP0.nomatch)>0)
{TEMP1.nomatch <- NULL
for(k in TEMP0.nomatch)
{
TEMP.k <- as.character(ListSP)[k]
print(TEMP.k)

tryCatch({REFTAXworms[k,] <- worms::wormsbymatchnames(TEMP.k,ids = T)},error=function(e){print(paste("WORMS search tool error for taxa =",TEMP.k))})

rm(TEMP.k)}
}#
}
##
}

## Warning missing species match
TEMP <- REFTAXworms[is.na(REFTAXworms$AphiaID),]
if(nrow(TEMP)>0)
{
print(paste("Warning no match found for the following species names"))
print(paste(TEMP$name,collapse=","))
}
rm(TEMP)
##


## Conserve seulement espèces trouvées sous WORMS
REFTAXworms <- REFTAXworms[is.na(REFTAXworms$AphiaID)==F,]

if(SEL.taxosup)
{
LIST.AphiaID <- REFTAXworms$AphiaID
REFTAXwormsCOMP <- list()

for(i in 1:length(LIST.AphiaID))
{REFTAXwormsCOMP[[paste(LIST.AphiaID[i])]] <- as.data.frame(worrms::wm_classification(LIST.AphiaID[i]))}

TEMP.NAMEStaxo <- NULL
for(j in 1:length(REFTAXwormsCOMP))
{TEMP.NAMEStaxo <- unique(c(TEMP.NAMEStaxo,REFTAXwormsCOMP[[j]]$rank))}

TEMP2 <- NULL
for(j in 1:length(REFTAXwormsCOMP))
{
TEMP0 <- NULL
for(k in TEMP.NAMEStaxo){
TEMP1 <- REFTAXwormsCOMP[[j]]$scientificname[REFTAXwormsCOMP[[j]]$rank==k][1]
if(length(TEMP1)==0)
{TEMP0 <- c(TEMP0,NA)}else{
TEMP0 <- c(TEMP0,TEMP1)}
rm(TEMP1)}
TEMP2 <- rbind(TEMP2,TEMP0)

rm(TEMP0)}

colnames(TEMP2) <- paste(TEMP.NAMEStaxo)
REFTAXworms <- data.frame("Nom_Scientifique"=REFTAXworms$name,REFTAXworms,TEMP2,stringsAsFactors=F)
}





# Extraction complète taxo worms
LIST.AphiaID <- na.omit(REFTAXworms$AphiaID)
REFTAXwormsCOMP <- list()
for(i in 1:length(LIST.AphiaID))
{REFTAXwormsCOMP[[paste(LIST.AphiaID[i])]] <- as.data.frame(worrms::wm_classification(LIST.AphiaID[i]))}


TEMP.NAMEStaxo <- NULL
for(j in 1:length(REFTAXwormsCOMP))
{TEMP.NAMEStaxo <- unique(c(TEMP.NAMEStaxo,REFTAXwormsCOMP[[j]]$rank))}

TEMP2 <- NULL
for(j in 1:length(REFTAXwormsCOMP))
{
TEMP0 <- NULL
for(k in TEMP.NAMEStaxo){
TEMP1 <- REFTAXwormsCOMP[[j]]$scientificname[REFTAXwormsCOMP[[j]]$rank==k][1]
if(length(TEMP1)==0)
{TEMP0 <- c(TEMP0,NA)}else{
TEMP0 <- c(TEMP0,TEMP1)}
rm(TEMP1)}
TEMP2 <- rbind(TEMP2,TEMP0)

rm(TEMP0)} 

colnames(TEMP2) <- paste(TEMP.NAMEStaxo)

TEMP2 <- TEMP2[match(as.character(LIST.AphiaID),names(REFTAXwormsCOMP)),]

REFTAXworms <- data.frame("Nom_Scientifique"=REFTAXworms$name,REFTAXworms,TEMP2,stringsAsFactors=F)

            print("Worms missing taxa")
            print(setdiff(ListSP,REFTAXworms$scientificname))
            print("Replaced by following taxa")
            print(setdiff(REFTAXworms$scientificname,ListSP))
if(F)
{write.csv2(REFTAXworms,file.path(DIR.RESULTS,"REFTAXworms.csv"),row.names=F)}

rm(TEMP2)

return(REFTAXworms)

}

## END FUNCTION