#####################################
######### WGFBIT - October 2019 ##########
##################################### >)))°>
######## SCRIPT dataset compilation - Celtic/BoB Eco-region (WGFBIT)
##
### R version 3.5.1
###############
####



#### 2 - DATASET download and formatting

 ### 2.1 - BTA DATASET

SEL.NewBTAfileUPDATE <- T

 # 1.2.5 - BTA table for each kind of biological dataset 
 FILE.NAME[["BTA"]] <- rbind(c("MACROFAUNA", "GRAB","BTA_EMODNET_LifeSpan-OBeauchard",".csv"),
                             c("MEGAFAUNA","TRAWL","BTA_BENTHIS_MEGAandMACROFAUNE_TABLE_VS4",".csv"))

 FILE.BTAsel <- FILE.NAME[["BTA"]][FILE.NAME[["BTA"]][,1]==SEL.COMPONENT]


# Standardize BTA class name
BTA.NAMEstd <- NULL



if(SEL.NewBTAfileUPDATE){


 

## Load previously created file
if(file.exists(file.path(DIR.OUPUT[[paste(1)]],paste(SEL.ScriptName,"BIOLOGICAL_BTA.RData",sep="")))){
load(file=file.path(DIR.OUPUT[[paste(1)]],paste(SEL.ScriptName,"BIOLOGICAL_BTA.RData",sep="")))}else{BIOLOGICAL_BTA <- list()}



   ## 2.1.0 - load BTA data

name <- FILE.BTAsel[1]

BIOLOGICAL_BTA[[name]][["MATRIX.FUNC"]] <- read.csv2(file=file.path(DIR.SOURCE[[paste(3)]],"BTA",paste(FILE.BTAsel[3],FILE.BTAsel[4],sep="")),header=T,sep=";",dec = ".",colClasses = "character",stringsAsFactors=F)



   ## 2.1.1 - Correct/Adapt traits matrix

TEMP.SELmatrix <- BIOLOGICAL_BTA[[name]][["MATRIX.FUNC"]]

    # 2.1.1.1 - Check / Suppress EMPTY Row
TEMP.SELmatrix <- data.frame("Taxa"=TEMP.SELmatrix[,1],as.data.frame(lapply(TEMP.SELmatrix[,-1],"as.numeric")))
TEMP.SELmatrix[is.na(TEMP.SELmatrix)] <- 0
TEMP.SELmatrix <- TEMP.SELmatrix[rowSums(TEMP.SELmatrix[,-1])>0,]

BIOLOGICAL_BTA[[name]][["MATRIX.FUNC"]] <- TEMP.SELmatrix

    # 2.1.1.2 - Update Species list
     # reuses by default the file already existing in the source directory (SEL.WORMSnew=F)

SEL.WORMSnew <- T # Select (T) to create the taxonomic table from Worms (if T !!! WARNING !!! Time consuming !!!)
               # OR reuses by default the file already existing in the source directory (SEL.WORMSnew=F)
               
# Select and test for already existing taxonomic data table
WormsFile <- file.path(DIR.SOURCE[[paste(3)]],"BTA",paste(FILE.BTAsel[3],"_LISTTaxa_WORMS",".csv",sep=""))

if(SEL.WORMSnew==T|file.exists(WormsFile)==F)
{

TEMP.SP1 <- unique(as.vector(TEMP.SELmatrix$Taxa))

TEMP.SP2 <- list()
TEMP.SP3 <- NULL

# Extract Worms information for each taxa
for(i in c(1:length(TEMP.SP1)))
{print(TEMP.SP1[i])
 skip_to_next <- F
tryCatch({TEMP.SP2[[i]] <- WORMSextract(ListSP=TEMP.SP1[i],SEL.marineOnly=F,verbose=F)
TEMP.SP3 <- unique(c(TEMP.SP3 ,colnames(TEMP.SP2[[i]])))}
,error = function(e) { skip_to_next <<- T})
if(skip_to_next){ next }
}

# Format taxonomic dataset
TEMP.SP6 <- NULL
TEMP.SP2b <- TEMP.SP2 
for(i in c(1:length(TEMP.SP2b)))
{
TEMP.SP4 <- TEMP.SP3[is.na(match(TEMP.SP3,names(TEMP.SP2b[[i]])))]
if(length(TEMP.SP4)>0)
{TEMP.SP5 <- t(data.frame(rep(NA,length(TEMP.SP4)))) 
colnames(TEMP.SP5) <- TEMP.SP4
TEMP.SP2b[[i]] <- as.data.frame(cbind(TEMP.SP2b[[i]],TEMP.SP5))
rm(TEMP.SP5)}

TEMP.SP2b[[i]] <- TEMP.SP2b[[i]][,order(names(TEMP.SP2b[[i]]))]
TEMP.SP6 <- rbind(TEMP.SP6,TEMP.SP2b[[i]])

rm(TEMP.SP4)}

## Suppress colnames TEMP2
## TO CHECK ERROR from previous part of the script
if(length(which(colnames(TEMP.SP6)=="TEMP2"))>0){
TEMP.SP6 <- TEMP.SP6[,-which(colnames(TEMP.SP6)=="TEMP2")]}

if(length(which(colnames(TEMP.SP6)=="Nom_Scientifique"))==1){
colnames(TEMP.SP6)[which(colnames(TEMP.SP6)=="Nom_Scientifique")] <- "TaxaSOURCE"}
                     
##
###

# Keep names not identified from WORMS dataset (na's in valid_name)
TEMP.TaxaCORR <- as.vector(TEMP.SP6$valid_name)
if(length(which(is.na(TEMP.SP6$valid_name)))>0)
{TEMP.TaxaCORR[is.na(TEMP.SP6$valid_name)] <- as.vector(TEMP.SP6$TaxaSOURCE[is.na(TEMP.SP6$valid_name)])}
TEMP.SP6$TaxaCORR <- TEMP.TaxaCORR
#

BIOLOGICAL_BTA[[name]][["MATRIX.FUNCtaxa"]] <- TEMP.SP6
write.csv2(TEMP.SP6, file.path(DIR.SOURCE[[paste(3)]],"BTA",paste(FILE.BTAsel[3],"_LISTTaxa_WORMS",".csv",sep="")))

rm(TEMP.SP1,TEMP.SP2,TEMP.SP3,TEMP.SP6)}else{
BIOLOGICAL_BTA[[name]][["MATRIX.FUNCtaxa"]] <- read.csv2(file.path(DIR.SOURCE[[paste(3)]],"BTA",paste(FILE.BTAsel[3],"_LISTTaxa_WORMS",".csv",sep="")))
}

     
    
    ## ## 2.1.2 - Redundancies test and suppression within functional matrix
    # Creation new corrected funtional matrix (MATRIX.FUNCnew)
    # & NAME CORRECTION from worms "valid_name" 
    
TEMP.MODIF <- NULL

TEMP.NAME0 <- as.vector(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCtaxa"]]$TaxaCORR[match(BIOLOGICAL_BTA[[name]][["MATRIX.FUNC"]]$Taxa,BIOLOGICAL_BTA[[name]][["MATRIX.FUNCtaxa"]]$TaxaSOURCE)])
TEMP.NAME0[is.na(TEMP.NAME0)] <- as.vector(rownames(BIOLOGICAL_BTA[[name]][["MATRIX.FUNC"]]))[is.na(TEMP.NAME0)]
TEMP.TEST0 <- table(TEMP.NAME0)[which(table(TEMP.NAME0)>1)]


TESTredundDIFF <- NULL # test / modifications info storage
ListredundTAXAsource <- NULL

for(SP in names(TEMP.TEST0)){
TEMP.TEST1 <- which(TEMP.NAME0==SP)

TEMP.TEST2 <- BIOLOGICAL_BTA[[name]][["MATRIX.FUNC"]][,-1][TEMP.TEST1,]

TEMP.TEST3 <- combn(1:length(TEMP.TEST1),2)

for(i in 1:ncol(TEMP.TEST3)){

TraitsDiff <- paste(colnames(setdiff(TEMP.TEST2[TEMP.TEST3[,i][1],],TEMP.TEST2[TEMP.TEST3[,i][2],])))
if(length(TraitsDiff)==0){TraitsDiff <- NA} 

ListredundTAXAsource <- c(ListredundTAXAsource, TEMP.TEST1)
TESTredundDIFF <- rbind(TESTredundDIFF,data.frame("validName"=SP,"redundantNames"=paste(BIOLOGICAL_BTA[[name]][["MATRIX.FUNC"]][,1][TEMP.TEST1],collapse="+"),"TraitsDiff"=TraitsDiff))

# Rules for redundant species: mean values / non empty traits

TEMP.TRAITS0 <- gsub('[[:digit:]]+', '', colnames(TEMP.TEST2))

TEMP.TRAITS5 <- NULL

for(TRAIT in unique(TEMP.TRAITS0))
{

TEMP.TRAITS1 <- which(is.na(match(gsub("([^A-Za])+", "", x = colnames(TEMP.TEST2)),TRAIT))==F)
TEMP.TRAITS2 <- TEMP.TEST2[,TEMP.TRAITS1]
TEMP.TRAITS2[is.na(TEMP.TRAITS2)] <- 0
TEMP.TRAITS2 <- t(apply(TEMP.TRAITS2, 1, function(i) i/sum(i))) # Harmonization / score 1

TEMP.TRAITS3 <- TEMP.TRAITS2[rowSums(TEMP.TRAITS2,na.rm=T)>0,]

if(length(nrow(TEMP.TRAITS3))==0&length(TEMP.TRAITS3)>0){
TEMP.TRAITS3 <- t(data.frame(TEMP.TRAITS3))}

if(nrow(TEMP.TRAITS3)>0){
#TEMP.TRAITS2[TEMP.TRAITS2==0] <- NA
TEMP.TRAITS4 <- colSums(TEMP.TRAITS3,na.rm=T)
TEMP.TRAITS5 <- c(TEMP.TRAITS5,TEMP.TRAITS4/sum(TEMP.TRAITS4))
}else{TEMP.TRAITS5 <-  c(TEMP.TRAITS5,colSums(TEMP.TRAITS2))}

}
TEMP.TRAITS5[is.na(TEMP.TRAITS5)] <- 0
TEMP.MODIF <- rbind(TEMP.MODIF,data.frame("Taxa"=SP,t(as.data.frame(TEMP.TRAITS5))))

}
}

if(length(ListredundTAXAsource)>0)
{BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]] <- rbind(BIOLOGICAL_BTA[[name]][["MATRIX.FUNC"]][-ListredundTAXAsource,],TEMP.MODIF)}

# Construct Traits Matrix from AphiaID
TEMP.MATCH <- match(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]]$Taxa,BIOLOGICAL_BTA[[name]][["MATRIX.FUNCtaxa"]]$valid_name)

if(length(which(is.na(TEMP.MATCH)==T))>0){
TEMP.MATCH[is.na(TEMP.MATCH)] <- match(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]]$Taxa[is.na(TEMP.MATCH)],BIOLOGICAL_BTA[[name]][["MATRIX.FUNCtaxa"]]$TaxaSOURCE)
TEMP.AphiaID <- BIOLOGICAL_BTA[[name]][["MATRIX.FUNCtaxa"]]$AphiaID[TEMP.MATCH]
MISSING_AphiaID <- BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]]$Taxa[is.na(TEMP.AphiaID)]
BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnewAphiaID"]] <- BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]][is.na(TEMP.AphiaID)==F,]
BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnewAphiaID"]]$Taxa <- as.vector(na.omit(TEMP.AphiaID))
rm(TEMP.MATCH)
}
#

if(length(which(is.na(TESTredundDIFF$TraitsDiff)==F))>0){
print("!!! Warning !!!")
print("Check discrepancies in trait matrix for redundant species ")
print("listed in TESTredundDIFF.csv")
print("!!! Warning !!!")} 

write.csv2(TESTredundDIFF, file.path(DIR.SOURCE[[paste(3)]],"BTA",paste(FILE.BTAsel[3],"_TESTredundDIFF_",".csv",sep="")))

BIOLOGICAL_BTA[[name]][["TESTredundDIFF"]] <- TESTredundDIFF





   ## 2.1.2 - Convert BTA dataset to matrix with species valid_name as rownames
   TEMP.ROWNAMES <- as.vector(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCtaxa"]]$TaxaCORR[match(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]]$Taxa,BIOLOGICAL_BTA[[name]][["MATRIX.FUNCtaxa"]]$TaxaSOURCE)])
   
   if(length(which(is.na(TEMP.ROWNAMES)))){
       TEMP.NAs <- BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]]$Taxa[is.na(TEMP.ROWNAMES)]
       TEMP.ROWNAMES[is.na(TEMP.ROWNAMES)] <- TEMP.NAs
   rm(TEMP.NAs)}
    
   rownames(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]]) <- TEMP.ROWNAMES
   BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]] <- BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]][,-which(colnames(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]])=="Taxa")]
   
   ## 2.1.2 - Convert BTA dataset to matrix with species AphiaID as rownames 
      ### TO COMPLETE
            ### TO COMPLETE
                  ### TO COMPLETE
                        ### TO COMPLETE
   
   ## 
   ##
   #

    
  
  
    

    ## 2.1.3 - Traits structures and modalities

# LISTE Fonctional groups of traits GROUPES F
TEMP.TRAITS <- gsub('[[:digit:]]+', '', colnames(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]]))

LIST.Fgroups <- unique(TEMP.TRAITS)
COL.BLOCKS <- table(TEMP.TRAITS)
COL.BLOCKS <- COL.BLOCKS[match(LIST.Fgroups, names(COL.BLOCKS))]

BIOLOGICAL_BTA[[name]][["BTAgroups"]] <- LIST.Fgroups
BIOLOGICAL_BTA[[name]][["BTAblocks"]] <- COL.BLOCKS

###


 

    ## 2.1.3.1 - Score harmonization

for(TRAIT in unique(LIST.Fgroups))
{

TEMPtrait <- which(is.na(match(gsub("([^A-Za])+", "", x = colnames(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]])),TRAIT))==F)

TEMP0 <- BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]][,TEMPtrait]
TEMP0[is.na(TEMP0)] <- 0
BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]][,TEMPtrait] <- t(apply(TEMP0, 1, function(i) i/sum(i)))     

t(apply(TEMP0, 1, function(i) i/sum(i)))

rm(TEMP0,TRAIT)}
BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]][is.na(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]])] <- 0
###



#TEMP$Taxa <- sub(" .*", "", TEMP$Taxa)
#BIOLOGICAL_BTA[[name]] <- unique(rbind(BIOLOGICAL_BTA[[name]][["MATRIX.FUNC"]],TEMP))
#BIOLOGICAL_BTA[[name]] <- data.frame("Taxa"=BIOLOGICAL_BTA[[name]][,1],as.data.frame(lapply(BIOLOGICAL_BTA[[name]][,-1],"as.numeric")))
#BIOLOGICAL_BTA[[name]] <- aggregate(.~Taxa,BIOLOGICAL_BTA[[name]],FUN="mean")
#BIOLOGICAL_BTA[[name]] <- cbind("Taxa"=BIOLOGICAL_BTA[[name]][,1],BIOLOGICAL_BTA[[name]][,-1]/rowSums(BIOLOGICAL_BTA[[name]][,-1])) # harmonization to 1
# aggregation of upper BTA level from mean score
## Add genus level
#BIOLOGICAL_BTA[[name]]$Genus <- sub(" .*","",BIOLOGICAL_BTA[[name]]$Taxa)
##

 
    ## 2.1.4 - Identification / automatic completion of taxa with missing information
##
#

# Select which levels to choose to complete traits matrix
SEL.LEVELS2COMPLETE <- c("genus","family","order")

# Nouvelle matrice avec compléments

if(T){
print("Complete missing traits by taxa")


# 2.1.4.1 - Search for missing traits information
for(TRAIT in unique(BIOLOGICAL_BTA[[name]][["BTAgroups"]]))
{print(TRAIT)
BIOLOGICAL_BTA[[name]][["MATRIX.FUNCcompleted"]][[TRAIT]] <- NULL
TEMPtrait <- which(is.na(match(gsub("([^A-Za])+", "", x = colnames(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]])),TRAIT))==F)

TEMP0 <- BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]][,TEMPtrait]
TEMP0[is.na(TEMP0)] <- 0
#LIST.MISSING[[TRAIT]] <- rownames(TEMP0)[rowSums(TEMP0)==0]
BIOLOGICAL_BTA[[name]][["MATRIX.FUNCmissing"]][[TRAIT]] <- rownames(TEMP0)[rowSums(TEMP0)==0]
BIOLOGICAL_BTA[[name]][["MATRIX.FUNCfilled"]][[TRAIT]] <- rownames(TEMP0)[rowSums(TEMP0)!=0]

BIOLOGICAL_BTA[[name]][["MATRIX.FUNCmissing"]][[TRAIT]] <- data.frame("taxa"=BIOLOGICAL_BTA[[name]][["MATRIX.FUNCmissing"]][[TRAIT]],BIOLOGICAL_BTA[[name]][["MATRIX.FUNCtaxa"]][,c("valid_name","genus","family","order")][match(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCmissing"]][[TRAIT]],BIOLOGICAL_BTA[[name]][["MATRIX.FUNCtaxa"]]$TaxaSOURCE),])

BIOLOGICAL_BTA[[name]][["MATRIX.FUNCfilled"]][[TRAIT]] <- data.frame("taxa"=BIOLOGICAL_BTA[[name]][["MATRIX.FUNCfilled"]][[TRAIT]],BIOLOGICAL_BTA[[name]][["MATRIX.FUNCtaxa"]][,c("valid_name","genus","family","order")][match(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCfilled"]][[TRAIT]],BIOLOGICAL_BTA[[name]][["MATRIX.FUNCtaxa"]]$TaxaSOURCE),])

    ## Complément auto missing data à partir espèces proches
if(nrow(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCmissing"]][[TRAIT]])>0){
for(TAXA in BIOLOGICAL_BTA[[name]][["MATRIX.FUNCmissing"]][[TRAIT]]$taxa){
TEMP4 <- NULL

TEMP1 <- BIOLOGICAL_BTA[[name]][["MATRIX.FUNCmissing"]][[TRAIT]][BIOLOGICAL_BTA[[name]][["MATRIX.FUNCmissing"]][[TRAIT]]$taxa==TAXA,]

for(LEVEL in SEL.LEVELS2COMPLETE)
{
if(length(TEMP4)==0){

if(is.na(TEMP1[[LEVEL]])==F){
TEMP2 <- BIOLOGICAL_BTA[[name]][["MATRIX.FUNCfilled"]][[TRAIT]][is.na(match(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCfilled"]][[TRAIT]][[LEVEL]],TEMP1[[LEVEL]]))==F,]
TEMP2 <- TEMP2[which(TEMP2$taxa!=TAXA),]

if(length(TEMP2)>1){

TEMP3 <- TEMP0[match(TEMP2$taxa, rownames(TEMP0)),]
TEMP4 <- colSums(t(apply(TEMP3, 1, function(i) i/sum(i))))
#TEMP4 <- t(as.data.frame((TEMP4/sum(TEMP4))*SELECT.POINTS))
TEMP4 <- t(as.data.frame((TEMP4/sum(TEMP4))))
rownames(TEMP4) <- TAXA
}}

if(length(TEMP4)>0)
{TEMP4 <- data.frame(TEMP4,"Level"=LEVEL,"Aggregate"=paste(rownames(TEMP3),collapse="+"))

# AJOUT nouvelle matrice fonctionnelle
BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]][,match(colnames(TEMP3),colnames(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]]))][match(rownames(TEMP4),rownames(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]])),] <- TEMP4[,match(colnames(TEMP3),colnames(TEMP4))]
}

}

rm(LEVEL)}



if(length(TEMP4)>0)
{BIOLOGICAL_BTA[[name]][["MATRIX.FUNCcompleted"]][[TRAIT]] <- rbind(TEMP4,BIOLOGICAL_BTA[[name]][["MATRIX.FUNCcompleted"]][[TRAIT]])}

rm(TAXA)}

rm(TEMP1)}


rm(TEMP0,TRAIT)} # END LOOP by TRAIT
} # SEL T/F

BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]][is.na(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]])] <- 0




## 2.1.5 - Complete new traits/modalitities from similar traits

    # 2.1.5.1 - NEW LONGEVITY CLASSIFICATION (6 classes "LG" from 4 classes "LO")

##
### A CORRIGER POUR CLASSIF A 5 MODALITES !!!
##

SEL.REF <- "LO"
SEL.NEW <- "LG"

TEMP.REF  <- BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]][,TEMP.TRAITS==SEL.REF]
TEMP.NEW <- BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]][,TEMP.TRAITS==SEL.NEW]

TEMP.MISSING <- which(rowSums(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]][,TEMP.TRAITS==SEL.NEW])==0)

if(length(TEMP.MISSING)>0){
TEMP.REF[is.na(TEMP.REF)] <- 0
TEMP.NEW[is.na(TEMP.NEW)] <- 0

    ## Attribution of LG modalities from initial LO values
    # Attribution rules
     # LG1   LO1+LO2/TOTAL
     TEMP.NEW[,1] <- TEMP.REF[,1]+TEMP.REF[,2]*(2/3)
      # LG2   LO3
     TEMP.NEW[,2] <- TEMP.REF[,2]*(1/3)+TEMP.REF[,3]*(2/3)
       # LG3   LO3
     TEMP.NEW[,3] <- TEMP.REF[,3]
        # LG4   L04 
     TEMP.NEW[,4] <- TEMP.REF[,4]
         # LG5   L04 2/6
     TEMP.NEW[,5] <- TEMP.REF[,4]*(1/3)
          # LG6   L04 1/6
     TEMP.NEW[,6] <- TEMP.REF[,4]*(1/6)

     # Score Harmonisation
     TEMP.NEW <- t(apply(TEMP.NEW, 1, function(i) i/sum(i)))
     #

## Ajout à la matrice complète
BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]][,which(is.na(match(gsub("([^A-Za])+", "", x = colnames(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]])),"LG"))==F)][match(rownames(TEMP.NEW),rownames(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]])),][TEMP.MISSING,] <- TEMP.NEW[TEMP.MISSING,]
BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]][is.na(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]])] <- 0
}

     ### Species to test/complete for higher longevity classes (>10)
SEL.CLASSsup10 <- c("LG4","LG5","LG6")

TEMP.TEST0 <- BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]][,grep("LG",colnames(BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]]))]
TEMP.TEST1 <- rownames(TEMP.TEST0)[which(rowSums(TEMP.TEST0[,match(SEL.CLASSsup10, colnames(TEMP.TEST0))])>0)]
##

         #
         par(mfrow=c(1,1))
         plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
         text(5.5, 8, paste("Test for species within higher longevity classes", sep=""), cex=1.2)
         text(5.5, 7, paste("considered higher classes are",paste(SEL.CLASSsup10,collapse=","), sep=" "), cex=1.1)
         text(5.5, 6, paste("From ", nrow( BIOLOGICAL_BTA[[name]][["MATRIX.FUNCnew"]]),"taxa total list", sep=" "), cex=1.1)
         text(5.5, 5, paste(length(TEMP.TEST1),"taxa have to be checked/completed for these classes"), cex=1.1)
        #


        # Save R object
if(file.exists(DIR.OUPUT[[paste(1)]])==F){dir.create(TEMP.DIR)}
save(BIOLOGICAL_BTA,file=file.path(DIR.OUPUT[[paste(1)]],paste(SEL.ScriptName,"BIOLOGICAL_BTA.RData",sep="")))
        #

}else{load(file=file.path(DIR.OUPUT[[paste(1)]],paste(SEL.ScriptName,"BIOLOGICAL_BTA.RData",sep="")))}
        ################################                      
        ### MERGING selected traits & Biomass per station
load(file=file.path(DIR.OUPUT[[paste(1)]],paste(SEL.ScriptName,"BIOLOGICAL_BTA.RData",sep="")))        
load(file=file.path(DIR.OUPUT[[paste(1)]],paste(SEL.ScriptName,"BIOdataset.RData",sep="")))        

for(SURVEY in names(BIOdataset$FORMATED)){

TEMP.BIO <- BIOdataset$FORMATED[[SURVEY]]$DATA
TEMP.BTAmatrix <- BIOLOGICAL_BTA[[SEL.COMPONENT]][["MATRIX.FUNCnew"]]

BIOdataset$FORMATED[[SURVEY]]$DATA_BTA <- cbind(TEMP.BIO,TEMP.BTAmatrix[,c("LO1","LO2","LO3","LO4")][match(TEMP.BIO$valid_name,rownames(TEMP.BTAmatrix)),])
rm(SURVEY,TEMP.BIO,TEMP.BTAmatrix)}
save(BIOdataset,file=file.path(DIR.OUPUT[[paste(1)]],paste(SEL.ScriptName,"BIOdataset.RData",sep="")))



##
rm(TEMP.NEW,TEMP.REF,TEMP.MISSING)

######################################################## 
############## END PART X BTA ##################
########################################################





