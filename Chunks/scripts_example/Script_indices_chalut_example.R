setwd("") ## charger le dossier dans lequel on travaille

### packages ###
library("plyr")



##### Chargement des tableaux de données ###
donnee_medits200 <- read.table("donnee_final_gdl200.csv", header=T, sep=";", dec=".") ## fichier espèces
param <- read.table("emplacement_abrasion_station_gdl.csv", sep=";", dec=".", header=T)  ### fichier position des traits de chalut
TDI<- read.table("TDI_all.csv", sep=";", dec=",", header=T)  ## fichier traits biologiques/espèces



#### Trawling Disturbance Impact (Foveau et al., 2017) ####

###### Correspondance entre le fichier final et la base de données des traits des espèces
### merge catch and tdi
TDI$Taxon <- as.character(TDI$Taxon)
tab<- merge(donnee_medits200,TDI, by.x="taxon", by.y="Taxon", all.x=T, all.y=F)


#### quels sont les taxons où il manque le tdi
tab$taxon<- as.character(tab$taxon)
tabx<-tab[which(is.na(tab$Mobility)),]
b<-unique(tabx$taxon)  #### nom des taxons n'ayant pas de traits


a<-merge(tabx, TDI, by.x="taxon", by.y="Family", all.x=T, all.y=F)
match(b, TDI$Family) ### ici on fait au niveau de la famille mais on peut le faire au niveau taxonomique que l'on veut


pb <- c()  #### création d'un vecteur vide pb
for(i in 1:length(b)){
  taxon=b[i]     #### taxon = nom des taxons n'ayant pas de traits
  Trait = TDI[which(TDI$Family==taxon),c(7:11)] #### création du data frame Trait qui prends les valeurs de traits si la famille est retrouvée dans la catégorie famille du tbaleau TDI
  Trait=unique(Trait)  ### permet de garder "qu'un exemplaire" des diférences
  if(nrow(Trait)==1){tab[which(tab$taxon==taxon), c(18:23)]<- Trait}  ### Trait comporte qu'une ligne alors on peut lui attribuer les valeurs de Trait
  if(nrow(Trait)>1){  #### Trait comporte plusieurs lignes
    if(as.numeric(apply(Trait[1:5],2, sd))[1]<1.5 & as.numeric(apply(Trait,2, sd))[2]<1.5 & as.numeric(apply(Trait,2, sd))[3]<1.5 & as.numeric(apply(Trait,2, sd))[4]<1.5 & as.numeric(apply(Trait,2, sd))[5]<1.5)   {  ### l'écart type de chaque trait doit être inférieur à 1.5
      if(sd(as.numeric(apply(Trait,1, sum)))<2.5){   #### l'écart type de la somme des traits (donc du TDI) doit être inférieur à 2.5
        pb <- rbind(pb, taxon)  #### tout les taxons qui répondent à ces conditions vont se retrouver dans pb
        Trait[nrow(Trait)+1,] <-apply(Trait,2,max)  #### si toutes ces conditions réunient alors on peut prendre la valeur maximale pour chaque trait et l'attribuer à la famille
        tab[which(tab$taxon==taxon), c(18:24)]<- Trait[nrow(Trait),]  #### et on peut le remettre dans le tableau du départ
      }
    }
  }
}

pb<- tab[which(is.na(tab$Mobility)),] #### permet de savoir quels taxons ont été enlevé

### vérification de la proportion de sp enlever apr station
w<- ddply(tab, .variables="code_trait", summarize, biom_tot= sum(biom1km))
q<- merge(pb, w, by="code_trait")
p <- ddply(q, .variables = "code_trait", summarize, biom_enl=sum(biom1km))
q<- merge(p, w, by="code_trait")
q$prop <- (q$biom_enl/q$biom_tot)*100
q<- subset(q, prop >= 25) ### on conserve que les taxons où moins de 25% de la biomasse a été supprimée par manque d'info
i<-unique(q$code_trait)
tab<- subset(tab, !code_trait %in% i)

#### suppression des taxons où les traits sont encore vides
if (length(pb)>0) tab <- tab[-which(is.na(tab$Mobility)),] else tab <- tab



# création de la colonne tdi
tab$TDI<-tab$Mobility + tab$Fragility + tab$Position + tab$Size + tab$Feeding.mode


# aggregate tab by trait, année et taxon
Newtab<- aggregate(biom1km ~ code_trait + Annee + taxon , data = tab, sum)

#compute total catch by trait et année
Tot<- aggregate(biom1km ~ code_trait + Annee   , data = Newtab, sum)


#compute relative biomass
Newtab<-merge(Newtab, Tot, by=c("code_trait"), all=F)
Newtab$relwei<-Newtab$biom1km.x/Newtab$biom1km.y #division du poid du taxon par poid total de l'annee

Test<-aggregate(relwei ~ code_trait , data = Newtab, sum)
summary(Test) ## permet de vérifier que l'abondance relative totale est bien de 1


Newtab<-Newtab[order(Newtab[,4], decreasing=F),]  ####D ordonner Newtab pour que ca soit dans le même ordre que tab

Newtab<-merge(Newtab,tab) # fusion des tables avec le tdi et avec les biomasses


#weighting tdi by relative biomass
#On applique la biomasse sur chaque indice pour avoir le nouveau TDI
Newtab$mobi<-as.numeric(Newtab$Mobility)*Newtab$relwei
Newtab$frag<-as.numeric(Newtab$Fragility)*Newtab$relwei
Newtab$posi<-as.numeric(Newtab$Position)*Newtab$relwei
Newtab$siz<-as.numeric(Newtab$Size)*Newtab$relwei
Newtab$feed<-as.numeric(Newtab$Feeding)*Newtab$relwei
Newtab$tdi2<-as.numeric(Newtab$TDI)*Newtab$relwei

#test tdi indices range 0-3 or 0-15
summary(Newtab)# le tdi est influencé maintenant par la biomasse

#aggregate per trawl
Trawltdi<- aggregate(cbind(mobi,frag,posi,siz,feed,tdi2) ~ code_trait + Annee+ Trait, data = Newtab, sum)
names(Trawltdi)<-c("code_trait", "Annee", "Trait","Mobility", "Fragility", "Position", "Size", "Feeding", "TDIm" )

### TDI par station
Trawltdi  ## jeu de donnée avec les valeurs de TDI par station


##### Calcul du TDI selon la méthode de De Juan #####
## formation de groupe en fonction des niveaux de TDI des espèces 

t<- unique(Newtab$taxon)

g5<-subset(Newtab, TDI>13)  #### formation du groupe g5
g4 <- subset(Newtab, TDI>10 & TDI<14) #### groupe g5 + g4
g3 <- subset(Newtab, TDI<11 & TDI>7)
g2<- subset(Newtab, TDI<8 & TDI>4)
g1 <- subset(Newtab, TDI<5)

#### calcul de l'abondance de chaque groupe à chaque station 
g5_ab<-ddply(g5, .variables = c("code_trait"), summarize, biomg5=sum(relwei))
g4_ab<-ddply(g4, .variables = c("code_trait"), summarize, biomg4=sum(relwei))
g3_ab<-ddply(g3, .variables = c("code_trait"), summarize, biomg3=sum(relwei))
g2_ab<-ddply(g2, .variables = c("code_trait"), summarize, biomg2=sum(relwei))
g1_ab<-ddply(g1, .variables = c("code_trait"), summarize, biomg1=sum(relwei))

#### calcul de l'abondance totale de chaque station 
ab_tot <- ddply(Newtab, .variables = c("code_trait"), summarize, abontot=sum(densite1km), biomtot=sum(relwei))

#### reunion des différents tableaux

ab_tot<-merge(ab_tot, g5_ab, by="code_trait", all.x = T , all.y=F)
ab_tot<-merge(ab_tot, g4_ab, by="code_trait", all.x = T , all.y=F)
ab_tot<-merge(ab_tot, g3_ab, by="code_trait", all.x = T , all.y=F)
ab_tot<-merge(ab_tot, g2_ab, by="code_trait", all.x = T , all.y=F)
ab_tot<-merge(ab_tot, g1_ab, by="code_trait", all.x = T , all.y=F)

ab_tot <-replace(ab_tot, is.na(ab_tot),0)  ### remplacement des NA par 0


#### calcul du TDI #####

ab_tot$TDI = (log(1)*log(ab_tot$biomg1 + 1) +log(2) *log(ab_tot$biomg2  + 1) + log(4) * log(ab_tot$biomg3  +1) +
                log(8)*log(ab_tot$biomg4  +1) + log(16)*log(ab_tot$biomg5  +1 ))/log(ab_tot$biomtot +1) ## jeu de donnée avec les valeurs de TDI de de Juan par station


ab_tot<- merge(ab_tot, param, by="code_trait")




### pTDI (Jac et al. 2020) ###

## création d'une table avec uniquement les espèces des groupes 3, 4 et 5
tab2<- subset (tab, TDI >7)

# aggregate tab by taxon,haul, year
Newtab2<- aggregate(biom1km ~ code_trait+ Annee + taxon , data = tab2, sum)


#compute total catch by haul, year 
Tot2<- aggregate(biom1km ~ code_trait + Annee  , data = Newtab2, sum)

#compute relative biomass
Newtab3<-merge(Newtab2, Tot2, by=c("code_trait"), all=F)
Newtab3$relwei<-Newtab3$biom1km.x/Newtab3$biom1km.y #division du poid du taxon par poid total de l'annee

Test<-aggregate(relwei ~ code_trait , data = Newtab3, sum)
summary(Test)


Newtab3<-Newtab3[order(Newtab3[,4], decreasing=F),]  ####D ordonner Newtab3 pour que ca soit dans le même ordre que tab

Newtab3<-merge(Newtab3,tab2) # fusion des tables avec le tdi et avec les biomasses


#weighting tdi by relative biomass
#On applique la biomasse sur chaque indice pour avoir le nouveau TDI
Newtab3$mobi<-as.numeric(Newtab3$Mobility)*Newtab3$relwei
Newtab3$frag<-as.numeric(Newtab3$Fragility)*Newtab3$relwei
Newtab3$posi<-as.numeric(Newtab3$Position)*Newtab3$relwei
Newtab3$siz<-as.numeric(Newtab3$Size)*Newtab3$relwei
Newtab3$feed<-as.numeric(Newtab3$Feeding)*Newtab3$relwei
Newtab3$tdi2<-as.numeric(Newtab3$TDI)*Newtab3$relwei

#test tdi indices range 0-3 or 0-15
summary(Newtab3)# le tdi est influencé maintenant par la biomasse

#aggregate per trawl
Trawltdi_A_3_4_5<- aggregate(cbind(mobi,frag,posi,siz,feed,tdi2) ~ code_trait + Annee + Trait, data = Newtab3, sum)
names(Trawltdi_A_3_4_5)<-c("code_trait", "Annee","Trait", "Mobility", "Fragility", "Position", "Size", "Feeding", "TDI_p" )

### TDI par station
Trawltdi_A_3_4_5 ## jeu de donnée avec les valeurs de pTDI pour chaque station


#### Indice de vulnérabilité modifié d'après Certain et al, 2015  (Jac et al., 2020) #####

data <- tab [,c(1,2,10, 18:23)]  ## on conserve que les colonnes du fichier tab qui nous intéressent

### remise des valeurs de traits entre 0 et 1
data$Mobility <- (data$Mobility+1)/4
data$Fragility <- as.numeric(data$Fragility)
data$Fragility <- (data$Fragility+1)/4
data$Size <- (data$Size+1)/4
data$Position<- (data$Position+1)/4
data$Feeding.mode <- (data$Feeding.mode+1)/4
data$VME <- data$VME/4  ## statut de protection (VME pour med et liste ospar pour Atlantique cf. base de traits)

a= data$Mobility*data$Position*data$Size  ### on a choisit que la mobilité, la position et la taille était des facteurs de sensibilité direct primaire
g=data$Fragility  ### et la fragilité un facteur direct aggravant

Vul.fct<-function(a,g,gamma=0.5){return(a^(1-(g/(g+gamma))))}
data$ti<- Vul.fct(a,g,gamma = 0.5)


b=(data$VME+data$Feeding.mode)/2   ### le statut de protection et le mode d'alimentation sont des facteurs indirects primaires
data$si<- Vul.fct(b, 0, gamma=0.5)

ab_rel<- ddply(data, .variables=c("code_trait"), summarize, ab_tot = sum(biom1km))
data<- merge(data, ab_rel, by="code_trait")
data$ab_relat<- data$biom1km/data$ab_tot  ## permet d'avoir l'abondance relative


data$tisi <- data$si*data$ti
data$vulnerabilite_sp <- data$ab_relat/(data$tisi)


## Calcul du Dj
vulne_sp_positive<- ddply(data, .variables=c("code_trait"), summarize, Dj_1=-sum(vulnerabilite_sp))
vulne_sp_positive  ### tableaux de données avec les valeurs de mT pour chacune des stations échantillonnées
