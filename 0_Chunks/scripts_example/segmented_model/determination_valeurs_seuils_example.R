setwd("C://Users/cjac/Desktop/scripts_example/segmented_model")

TDI_cgfs <- read.table("TDI_cgfs.csv", dec=".", sep=";", header=T)
TDI_cgfs <- TDI_cgfs[, c(1,8)]

### rajouter la corse

library(e1071)
library(MASS)
library(gam)
library(segmented)
library(plyr)
library(robustHD)
library(optimx)


################### Exemple où le modèle linéaire négatif est sélectionné  ###################
TDI_habitat<- merge(Mdn, TDI_cgfs, by="code", all.x=F)
summary(TDI_habitat$HAB_TYPE)
summary(TDI_habitat$MFSD_BH17)

c<-subset(TDI_habitat, HAB_TYPE=="A5.14")
skewness(sqrt(c$risk))
kurtosis(sqrt(c$risk))

c$trans_risk <- sqrt(c$risk)
c$TDI_stand <- standardize(c$TDI)

plot(TDI_stand ~ trans_risk, data=c) ## représentation des données

#on va du plus simple au plus complexe

### modèle null (ou totalement plat)
y <- c$TDI_stand
x <- c$trans_risk
o<-lm(y~1)
summary(o)  ###   R²adj : NULL
lines(x, fitted(o), col=5)

## modèle lineaire simple
a<- lm(y ~ x)
summary(a) ###  verif de la pente et R²ad
### corrélation négative et significative
#adjr² = 0.024

lines(x, fitted(a), col=4) 

#modèle à un break point

g<- segmented(a, seg.Z = ~ x,control=seg.control(display=FALSE)) 
g



## 2 breakpoints
g2<- segmented(a, seg.Z = ~ x,control=seg.control(display=FALSE), psi=c(2, 4)) ## pas de breakpoint  
g2  

#la valeur du breakpoint estimée
g$psi[2]

# plat avant le breakpoint

fit.glm<-update(a,.~. -x)
fit.seg1 <- segmented.lm(fit.glm, seg.Z = ~x)


summary(fit.seg1)
# C'est ici que tu vas trouver le bon adjR2. adjR²= 0.021

# si la pente est positive: le modèle est disqualifié
#si la pente est negative il faut tester la significativité de la pente...
#pour ca, il va falloir ajuster un lm sur les données prises sur la pente (en supprimant les pallier) et ainsi tester la significativité de la pente.

yl<-y[x>=fit.seg1$psi[2]]
xl<-x[x>=fit.seg1$psi[2]]


lmf<-lm(yl~xl)
summary(lmf) ### négative et significative

#test breakpoint le pscore.test ne marche pas sur ces modèles contraints
davies.test(fit.seg1,seg.Z = ~ x,k=100) ## pas de breakpoint


#plot(fit.seg1)
lines(x[order(x, na.last = NA)],fitted(fit.seg1)[order(x, na.last = NA)], col=2)

## plat après le breakpoint
xx<- -x
fit.seg2<-segmented.lm(o,seg.Z=~xx)  ## pas de breakpoint


######################  Exemple où le modèle avec 1 breakpoint est sélectionné #####################
c<-subset(TDI_habitat, HAB_TYPE=="A5.15")

c$trans_risk <- sqrt(c$risk)
c$TDI_stand <- standardize(c$TDI)

plot(TDI_stand ~ trans_risk, data=c)

#on va du plus simple au plus complexe

### modèle null (ou totalement plat)
y <- c$TDI_stand
x <- c$trans_risk
o<-lm(y~1)
summary(o)  ###   R²adj : NULL
lines(x, fitted(o), col=5)

## modèle lineaire simple
a<- lm(y ~ x)
summary(a) ###  verif de la pente et R²ad
### corrélation négative et significative
# adjr² = 0.222

lines(x, fitted(a), col=4) 

#modèle à un break point

g<- segmented(a, seg.Z = ~ x,control=seg.control(display=FALSE)) 
g
#la valeur du breakpoint estimée
g$psi[2]

# plat avant le breakpoint

fit.glm<-update(a,.~. -x) 
fit.seg1 <- segmented.lm(fit.glm, seg.Z = ~x)
## pas ce modèle


## plat après le breakpoint
xx<- -x
fit.seg2<-segmented.lm(o,seg.Z=~xx)

# on verifie la pente et le ajusted R2 avec
summary(fit.seg2) ## adjr² = 0.259

# si la pente est positive: le modèle est disqualifié
#si la pente est negative il faut tester la significativité de la pente...
#pour ca, il va falloir ajuster un lm sur les données prises sur la pente (en supprimant les pallier) et ainsi tester la significativité de la pente.

yl<-y[x<=-fit.seg2$psi[2]]
xl<-x[x<=-fit.seg2$psi[2]]


lmf<-lm(yl~xl)
summary(lmf) ### pente négative et significative

#test breakpoint le pscore.test ne marche pas sur ces modèles contraints

davies.test(fit.seg2,seg.Z = ~ xx,k=100) ## ok


#Puis on fait le plot 

lines(x[order(x, na.last = NA)],fitted(fit.seg2)[order(x, na.last = NA)], col=3)



## 2 breakpoints
g2<- segmented(a, seg.Z = ~ x,control=seg.control(display=FALSE), psi=c(2, 6))   
g2  
summary(g2)


# ajuster un modèle à deux palliers complet avec optimx
#https://stats.stackexchange.com/questions/149627/piecewise-regression-with-constraints


#we need four parameters: the two breakpoints and the starting and ending intercepts
fun <- function(par, x) {
  #set all y values to starting intercept
  y1 <- x^0 * par["i1"]
  #set values after second breakpoint to ending intercept
  y1[x >= par["x2"]] <- par["i2"]
  #which values are between breakpoints?
  r <- x > par["x1"] & x < par["x2"]
  #interpolate between breakpoints
  y1[r] <- par["i1"] + (par["i2"] - par["i1"]) / (par["x2"] - par["x1"]) * (x[r] - par["x1"])
  y1
}

#sum of squared residuals
SSR <- function(par) {
  sum((y - fun(par, x))^2)
}


library(optimx)

#on fixe les paramètres initiaux x1, x2 (les breakpoints) sur la base du modèle g2
#on calcul des valeurs de palier de debut et de fin sur la base du modèle O (le modèle null)
#ie le pallier de depart = l'intercept du model null + la moyenne des residus du modele 
#et le pallier de fin = l'intercept du model null - la moyenne des residus du modele

fit.mod5<-optimx(par = c(x1 = g2$psi[1,2], x2 = g2$psi[2,2], i1 = mean(o$coefficients)+mean(o$residuals), i2 = mean(o$coefficients)-mean(o$residuals)), 
                 fn = SSR, 
                 method = "Nelder-Mead")

summary(fit.mod5)

#ici le modèle n'estime que les paramètres pas de test de significativité

yl<-y[x>=fit.mod5$x1]
xl<-x[x>=fit.mod5$x1]
yl<-yl[xl<=fit.mod5$x2]
xl<-xl[xl<=fit.mod5$x2]

lmf<-lm(yl~xl)
summary(lmf)

#pas de test de significativité des breakpoints :-(
#A defaut on peut verifier si la correlation obs/predicted est significative

cor.test(y,fun(c(x1 = fit.mod5$x1, x2 = fit.mod5$x2, i1 = fit.mod5$i1, i2 = fit.mod5$i2),x),method="pearson")

#Enfin, il faut se calculer le ajdR2 à la main avec des modif sur ci-dessous.


ybar<-mean(y)
totdev<-sum((y-ybar)^2, na.rm = TRUE)
fitted<-fun(c(x1 = fit.mod5$x1, x2 = fit.mod5$x2, i1 = fit.mod5$i1, i2 = fit.mod5$i2), x)
errdev<-sum((fitted-y)^2, na.rm = TRUE)

dftotal<-length(y)-1
dfres <-dftotal-4 #c'est un modèle à 4 paramètres

Rsq<-  1-(errdev/totdev)
adjRsq<- 1-(1-Rsq)*(dftotal/dfres)

adjRsq  ## 0.255 donc Adjr² plus faible que pour le modèle avec le plat après le breakpoint --> modèle avec 2 bkpt non sélectionné


################ Exemple où le modèle null a été sélectionné ###################
c<-subset(TDI_habitat, HAB_TYPE=="A5.27")
skewness((c$risk))
kurtosis((c$risk))

c$trans_risk <- (c$risk)
c$TDI_stand <- standardize(c$TDI)

plot(TDI_stand ~ trans_risk, data=c)

#on va du plus simple au plus complexe

### modèle null (ou totalement plat)
y <- c$TDI_stand
x <- c$trans_risk
o<-lm(y~1)
summary(o)  ###   R²adj : NULL
lines(x, fitted(o), col=5)

## modèle lineaire simple
a<- lm(y ~ x)
summary(a) ###  verif de la pente et R²ad
### corrélation négative mais pas significative
# adjr² = 0.029

