######################################
##### FUNCTION TO COMPUTE TDIs indices
##  Adapted from Jac et al. 2020
## update 2021/05
# derniere modification 05/07/2021
######################################

#########################################################
##############     SUMMARY            ###################
#########################################################

#### WARNING
# charger fonction WormsExtract

####


#### OPTIONS

# SELprotecStatus : "VME" or "OSPAR"


#### 1 - SOURCE DIRECTORIES and FILES
 # 1.1 - SOURCE DIRECTORY
 # 1.2 - SOURCE FILES
 # 1.3 - Harmonisation des champs de donn?es
 # 1.4 - Chargement des tableaux de donn?es ###
 
#### 2- TDI: Trawling Disturbance Impact (Foveau et al., 2017) ####
 # 2.1 - Correspondance entre le fichier final et la base de donn?es des traits des esp?ces


#### 3 - TDI: Trawling Disturbance Impact (Foveau et al., 2017) ####


#### 4 - Calcul du TDI selon la m?thode de De Juan #####

### 5 - pTDI (Jac et al. 2020) ###

#### 6 - mTDI: Indice de vuln?rabilit? modifi? d'apr?s Certain et al, 2015  (Jac et al., 2020) #####




# SELtaxaLevel = niveau taxo maximal utilis? pour compl?ter les don?nes de traits pour les taxa manquants
IndicTDI <- function(DataBio=NA,DataStation=NA,DataVul=NA,DataTaxref=NA,SELtaxaLevel="family",
FIELDtaxa=NA,FIELDstationID=NA,FIELDbioquanti=NA,VARprop=25,WormsCheck=T,SELprotecStatus="VME")
{


# Traits field names standard
TraitsFields <-  c("Mobility","Fragility","Position","Size","Feeding.mode","Protection.status.VME.","Protection.status..OSPAR.")
#


# Taxonomic standard levels
TaxoLevels <- c("genus","family","order","class","phylum")
TEMP <- which(TaxoLevels==SELtaxaLevel)
if(length(TEMP)>0){SELtaxaLevel2 <- TaxoLevels[c(1:TEMP)]}else{SEL.TaxoLevels <- "Taxa"}
#


#### 1 - PREPARATION DIRECTORIES and FILES
##

# 1.3 - Harmonisation des champs de donn?es

LIST.FIELDnames <- rbind(
cbind("Taxa",FIELDtaxa),
cbind("StationID",FIELDstationID),
cbind("VARbio",FIELDbioquanti))
colnames(LIST.FIELDnames) <- c("StandardFieldName", "DatasetFieldName")
#

for(OBJ in c("DataBio","DataStation","DataVul")){
TEMP.DATA <- get(OBJ)
if(length(na.omit(match(colnames(TEMP.DATA),LIST.FIELDnames[,2])))>0){
colnames(TEMP.DATA)[is.na(match(colnames(TEMP.DATA),LIST.FIELDnames[,2]))==F] <- LIST.FIELDnames[,1][na.omit(match(colnames(TEMP.DATA),LIST.FIELDnames[,2]))]}
assign(OBJ, TEMP.DATA)
rm(TEMP.DATA)}

#




#### 2- Correspondance entre le fichier final et la base de donn?es des traits des esp?ces

### merge catch and tdi
# DataVul$Taxa <- as.character(DataVul$Nom_Scientifique)
tab <- merge(DataBio,DataVul, by.x="Taxa", by.y="Taxa", all.x=T, all.y=F)


#### quels sont les taxons o? il manque le tdi
tab$Taxa <- as.character(tab$Taxa)
tabx <- tab[which(is.na(tab$Mobility)),]
MISStaxa <- unique(tabx$Taxa)  #### nom des taxons n'ayant pas de traits
COVEREDtaxa <- sort(setdiff(DataBio$Taxa,MISStaxa))


### PROBLEME EXEMPLE "Dichelopandalus bonnieri"

pb2 <- c()  #### cr?ation d'un vecteur vide pb
for(i in 1:length(MISStaxa)){
  taxon=MISStaxa[i]     #### taxon = nom des taxons n'ayant pas de traits
  names(taxon) <- "Taxa"
  taxon2 <- NULL
  if(WormsCheck&exists('WORMSextract', mode='function')){
  # options(warn=-1)
  tryCatch(Wormstaxon <- WORMSextract(taxon),error=function(e){print(paste("WORMS search tool error for taxa =",taxon))})
  taxon2 <- unique(Wormstaxon[,SELtaxaLevel2])
  # options(warn=0)
  }else{if(length(taxon2)==0){taxon2 <- taxon}}

  if(WormsCheck==F&nrow(DataTaxref)>0){taxon2 <- NULL
  taxon2 <- unique(nivsys[,SELtaxaLevel2])
  taxon2 <- as.data.frame(taxon2)
  names(taxon2) <- SELtaxaLevel2
  
  }else{if(length(taxon2)==0){taxon2 <- taxon}}

  Trait=NULL
  for(level in SELtaxaLevel2){
    if(length(Trait)==0){
      Trait = DataVul[which(DataVul[,level] %in% taxon2[[level]]),c(TraitsFields)] #### cr?ation du data frame Trait qui prends les valeurs de traits si la famille est retrouv?e dans la cat?gorie famille du tbaleau TDI
      Trait=unique(Trait)  ### permet de garder "qu'un exemplaire" des diff?rences
      
      if(nrow(Trait)==1){tab[which(tab$Taxa==taxon), TraitsFields] <- Trait}  ### Trait comporte qu'une ligne alors on peut lui attribuer les valeurs de Trait
      
      if(nrow(Trait)>1){  #### Trait comporte plusieurs lignes
        if(as.numeric(apply(Trait[1:5],2, sd))[1]<1.5 & as.numeric(apply(Trait,2, sd))[2]<1.5 & as.numeric(apply(Trait,2, sd))[3]<1.5 & as.numeric(apply(Trait,2, sd))[4]<1.5 & as.numeric(apply(Trait,2, sd))[5]<1.5)   {  ### l'?cart type de chaque trait doit ?tre inf?rieur ? 1.5
          if(sd(as.numeric(apply(Trait,1, sum)))<2.5){   #### l'?cart type de la somme des traits (donc du TDI) doit ?tre inf?rieur ? 2.5
            pb2 <- rbind(pb2, taxon)  #### tout les taxons qui r?pondent ? ces conditions vont se retrouver dans pb
            Trait[nrow(Trait)+1,] <- apply(Trait,2,max)  #### si toutes ces conditions r?unient alors on peut prendre la valeur maximale pour chaque trait et l'attribuer ? la famille
            tab[which(tab$Taxa==taxon), TraitsFields]<- Trait[nrow(Trait),]  #### et on peut le remettre dans le tableau du d?part
          }
          
        }else{Trait=NULL}
        }
      
      } # Test traits d?j? trouv? au niveau inf?rieur
} # End loop by levels

rm(Trait,taxon,taxon2)} # End loop by taxon

pb <- tab[which(is.na(tab$Mobility)),] #### permet de savoir quels taxons ont ?t? enlev?
MISStaxa2 <- unique(pb$Taxa)
if(F){table(pb[,c("Taxa","StationID")])}

### v?rification de la proportion de sp enlev?e par station
w<- aggregate(as.formula(paste("VARbio","~","StationID")),data=tab,FUN="sum")
colnames(w)[which(colnames(w)=="VARbio")] <- "VARbio_tot"
q <- merge(pb, w, by="StationID")


p <- aggregate(as.formula(paste("VARbio","~","StationID")),data=q,FUN="sum")
colnames(p)[which(colnames(p)=="VARbio")] <- "VARbio_enl"
q <- merge(p, w, by="StationID")

q$prop <- (q$VARbio_enl/q$VARbio_tot)*100

if(is.null(VARprop)|is.na(VARprop)){}else{
q <- subset(q, prop >= VARprop)} ### on conserve que les taxons o? moins de X% (d?fini par VARprop) de la biomasse a ?t? supprim?e par manque d'info
  
sta <- unique(q$StationID)
tab <- subset(tab, !StationID %in% sta)

MISSINGStation <- setdiff(DataStation$StationID,unique(tab$StationID))


#### suppression des taxons o? les traits sont encore vides
if (length(pb)>0) tab <- tab[-which(is.na(tab$Mobility)),] else tab <- tab


#### 3 - TDI: Trawling Disturbance Impact (Foveau et al., 2017) ####

# cr?ation de la colonne tdi
tab$TDI<-tab$Mobility + tab$Fragility + tab$Position + tab$Size + tab$Feeding.mode


# aggregate tab by trait, ann?e et taxon
Newtab<- aggregate(as.formula(paste("VARbio", "~", paste("StationID", "Taxa",sep="+"))) , data = tab, sum)

#compute total catch by trait et ann?e
Tot<- aggregate(as.formula(paste("VARbio", "~", paste("StationID",sep="+"))) , data = Newtab, sum)


# compute relative VARbio
Newtab<-merge(Newtab, Tot, by=c("StationID"), all=F)
Newtab$relwei<-Newtab$VARbio.x/Newtab$VARbio.y #division du poid du taxon par poid total de l'annee

Test<-aggregate(relwei ~ StationID , data = Newtab, sum)
summary(Test) ## permet de v?rifier que l'abondance relative totale est bien de 1


Newtab<-Newtab[order(Newtab[,4], decreasing=F),]  ####D ordonner Newtab pour que ca soit dans le m?me ordre que tab

Newtab<-merge(Newtab,tab) # fusion des tables avec le tdi et avec les biomasses


# Weighting tdi by relative biomass
#On applique la biomasse sur chaque indice pour avoir le nouveau TDI
Newtab$mobi<-as.numeric(Newtab$Mobility)*Newtab$relwei
Newtab$frag<-as.numeric(Newtab$Fragility)*Newtab$relwei
Newtab$posi<-as.numeric(Newtab$Position)*Newtab$relwei
Newtab$siz<-as.numeric(Newtab$Size)*Newtab$relwei
Newtab$feed<-as.numeric(Newtab$Feeding)*Newtab$relwei
Newtab$tdi2<-as.numeric(Newtab$TDI)*Newtab$relwei





### TDI par station
#aggregate per trawl
## jeu de donn?e avec les valeurs de TDI par station
RES_TDI <- aggregate(cbind(mobi,frag,posi,siz,feed,tdi2) ~ StationID, data = Newtab, sum)
names(RES_TDI) <- c("StationID", "Mobility", "Fragility", "Position", "Size", "Feeding", "TDIm" )




#### 4 - Calcul du TDI selon la methode de De Juan #####
# djTDI
## formation de groupe en fonction des niveaux de TDI des esp?ces 
t <- unique(Newtab$Taxa)

# Definitions groupes de sensibilit?
g5 <- subset(Newtab, TDI>13)  #### formation du groupe g5
g4 <- subset(Newtab, TDI>10 & TDI<14) #### groupe g5 + g4
g3 <- subset(Newtab, TDI<11 & TDI>7)
g2 <- subset(Newtab, TDI<8 & TDI>4)
g1 <- subset(Newtab, TDI<5)

#### calcul de l'abondance de chaque groupe ? chaque station 
g5_ab<-plyr::ddply(g5, .variables = c("StationID"), summarize, biomg5=sum(relwei))
g4_ab<-plyr::ddply(g4, .variables = c("StationID"), summarize, biomg4=sum(relwei))
g3_ab<-plyr::ddply(g3, .variables = c("StationID"), summarize, biomg3=sum(relwei))
g2_ab<-plyr::ddply(g2, .variables = c("StationID"), summarize, biomg2=sum(relwei))
g1_ab<-plyr::ddply(g1, .variables = c("StationID"), summarize, biomg1=sum(relwei))

#### calcul de l'abondance totale de chaque station 
ab_tot <- ddply(Newtab, .variables = c("StationID"), summarize, abontot=sum(VARbio), biomtot=sum(relwei))

#### reunion des diff?rents tableaux

ab_tot<-merge(ab_tot, g5_ab, by="StationID", all.x = T , all.y=F)
ab_tot<-merge(ab_tot, g4_ab, by="StationID", all.x = T , all.y=F)
ab_tot<-merge(ab_tot, g3_ab, by="StationID", all.x = T , all.y=F)
ab_tot<-merge(ab_tot, g2_ab, by="StationID", all.x = T , all.y=F)
ab_tot<-merge(ab_tot, g1_ab, by="StationID", all.x = T , all.y=F)

ab_tot <-replace(ab_tot, is.na(ab_tot),0)  ### remplacement des NA par 0


#### calcul du TDI #####

ab_tot$djTDI = (log(1)*log(ab_tot$biomg1 + 1) +log(2) *log(ab_tot$biomg2  + 1) + log(4) * log(ab_tot$biomg3  +1) +
                log(8)*log(ab_tot$biomg4  +1) + log(16)*log(ab_tot$biomg5  +1 ))/log(ab_tot$biomtot +1) ## jeu de donn?e avec les valeurs de TDI de de Juan par station

RES_djTDI <- ab_tot

if(F){RES_djTDI <- merge(RES_djTDI, DataStation, by="StationID")}




### 5 - pTDI (Jac et al. 2020) ###

## cr?ation d'une table avec uniquement les esp?ces des groupes 3, 4 et 5
tab2 <- subset (tab, TDI >7)

# aggregate tab by taxon,haul, year
Newtab2<- aggregate(as.formula(paste("VARbio","~",paste("StationID","Taxa",sep="+"))), data = tab2, sum)


#compute total catch by haul, year 
Tot2<- aggregate(as.formula(paste("VARbio","~",paste("StationID",sep="+"))), data = Newtab2, sum)

#compute relative biological variable (e.g. biomass)
Newtab3 <- merge(Newtab2, Tot2, by=c("StationID"), all=F)
Newtab3$relwei <- Newtab3$VARbio.x/Newtab3$VARbio.y #division du poid du taxon par poid total de l'annee

Test <- aggregate(as.formula(paste("relwei", "~", "StationID")), data = Newtab3, sum)
summary(Test)


Newtab3<-Newtab3[order(Newtab3[,4], decreasing=F),]  ####D ordonner Newtab3 pour que ca soit dans le m?me ordre que tab

Newtab3<-merge(Newtab3,tab2) # fusion des tables avec le tdi et avec les biomasses


#weighting tdi by relative biological variable
#On applique la pond?ration par la variable biologique sur chaque indice pour avoir le nouveau TDI
Newtab3$mobi<-as.numeric(Newtab3$Mobility)*Newtab3$relwei
Newtab3$frag<-as.numeric(Newtab3$Fragility)*Newtab3$relwei
Newtab3$posi<-as.numeric(Newtab3$Position)*Newtab3$relwei
Newtab3$siz<-as.numeric(Newtab3$Size)*Newtab3$relwei
Newtab3$feed<-as.numeric(Newtab3$Feeding)*Newtab3$relwei
Newtab3$tdi2<-as.numeric(Newtab3$TDI)*Newtab3$relwei

#test tdi indices range 0-3 or 0-15
#summary(Newtab3)# le tdi est influenc? maintenant par la biomasse

### TDI par station
## jeu de donn?e avec les valeurs de pTDI pour chaque station
#aggregate per trawl
RES_pTDI <- aggregate(as.formula(paste("cbind(mobi,frag,posi,siz,feed,tdi2)", "~", paste("StationID",sep="+"))), data = Newtab3, sum)
names(RES_pTDI) <- c("StationID","Mobility", "Fragility", "Position", "Size", "Feeding", "TDI_p" )


 


#### 6 - mTDI: Indice de vuln?rabilit? modifi? d'apr?s Certain et al, 2015  (Jac et al., 2020) #####

data <- tab

### remise des valeurs de traits entre 0 et 1
data$Mobility <- (data$Mobility+1)/4
data$Fragility <- as.numeric(data$Fragility)
data$Fragility <- (data$Fragility+1)/4
data$Size <- (data$Size+1)/4
data$Position<- (data$Position+1)/4
data$Feeding.mode <- (data$Feeding.mode+1)/4

## Ajout statut de protection (VME pour med et liste ospar pour Atlantique cf. base de traits)
# Protection.status.VME. Protection.status..OSPAR.
data[,SELprotecStatus] <- data[,grep(SELprotecStatus,colnames(data))]/4


a= data$Mobility*data$Position*data$Size  ### on a choisit que la mobilit?, la position et la taille ?tait des facteurs de sensibilit? direct primaire                                                  
g=data$Fragility  ### et la fragilit? un facteur direct aggravant

Vul.fct<-function(a,g,gamma=0.5){return(a^(1-(g/(g+gamma))))}
data$ti<- Vul.fct(a,g,gamma = 0.5)

b=(data[, SELprotecStatus]+data$Feeding.mode)/2   ### le statut de protection et le mode d'alimentation sont des facteurs indirects primaires
data$si<- Vul.fct(b, 0, gamma=0.5)

ab_rel <- aggregate(as.formula(paste("VARbio","~","StationID")),data=data,FUN="sum")
colnames(ab_rel)[which(colnames(ab_rel)=="VARbio")] <- "ab_tot"
data <- merge(data, ab_rel, by="StationID")
data$ab_relat <- data$VARbio/data$ab_tot  ## permet d'avoir l'abondance relative

data$tisi <- data$si*data$ti
data$vulnerabilite_sp <- data$ab_relat/(data$tisi)

### tableaux de donn?es avec les valeurs de mT pour chacune des stations ?chantillonn?es
## Calcul du Dj
RES_mTDI <- aggregate(as.formula(paste("vulnerabilite_sp","~","StationID")),data=data,FUN="sum")
colnames(RES_mTDI)[which(colnames(RES_mTDI)=="vulnerabilite_sp")] <- "mTDI"

RESULTS <- list(RES_TDI,RES_djTDI,RES_pTDI,RES_mTDI,MISStaxa2,MISSINGStation)
names(RESULTS) <- c("RES_TDI","RES_djTDI","RES_pTDI","RES_mTDI","MISSINGtaxa","MISSINGStation")
return(RESULTS)  
 } # END FUNCTION
#######


