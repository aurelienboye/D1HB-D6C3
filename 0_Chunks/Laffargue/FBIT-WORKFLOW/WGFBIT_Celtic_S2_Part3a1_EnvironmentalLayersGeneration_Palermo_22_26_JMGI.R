#Script developed by González-Irusta to produce Longevity predictions across BoB and Celtic seas region
#Packages
library(ncdf4) 
library(raster) 
library(rgdal) 
library(ggplot2)
library(ClusterR)
library(rgeos)
library(sf)

#We load bathymetry

#Bathymetry was downloaded from EMODNET.

#Emodnet_c1 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/C/C1_2020.asc")
#crs(Emodnet_c1) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#Emodnet_c2 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/C/C2_2020.asc")
#crs(Emodnet_c2) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#Emodnet_c3 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/C/C3_2020.asc")
#crs(Emodnet_c3) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#Emodnet_c4 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/C/C4_2020.asc")
#crs(Emodnet_c4) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

#Emodnet_C <- raster::merge(Emodnet_c1, Emodnet_c2, Emodnet_c3, Emodnet_c4)

#Emodnet_d1 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/D/D1_2020.asc")
#crs(Emodnet_d1) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#Emodnet_d2 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/D/D2_2020.asc")
#crs(Emodnet_d2) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#Emodnet_d3 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/D/D3_2020.asc")
#crs(Emodnet_d3) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#Emodnet_d4 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/D/D4_2020.asc")
#crs(Emodnet_d4) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

#Emodnet_d <- raster::merge(Emodnet_d1, Emodnet_d2, Emodnet_d3, Emodnet_d4)

#Emodnet_e1 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/E/E1_2020.asc")
#crs(Emodnet_e1) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#Emodnet_e2 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/E/E2_2020.asc")
#crs(Emodnet_e2) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#Emodnet_e3 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/E/E3_2020.asc")
#crs(Emodnet_e3) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#Emodnet_e4 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/E/E4_2020.asc")
#crs(Emodnet_e4) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

#Emodnet_e <- raster::merge(Emodnet_e1, Emodnet_e2, Emodnet_e3, Emodnet_e4)

#Emodnet_f1 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/F/F1_2020.asc")
#crs(Emodnet_f1) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#Emodnet_f2 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/F/F2_2020.asc")
#crs(Emodnet_f2) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#Emodnet_f3 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/F/F3_2020.asc")
#crs(Emodnet_f3) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#Emodnet_f4 <- raster("E:/GIS/Bathymetry/Emodnet_Bathy_2020/F/F4_2020.asc")
#crs(Emodnet_f4) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

#Emodnet_e <- raster::merge(Emodnet_f1, Emodnet_f2, Emodnet_f3, Emodnet_f4)

#writeRaster(Emodnet_C, filename="E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/Depth/C.tif", format="GTiff", overwrite=TRUE)
#writeRaster(Emodnet_d, filename="E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/Depth/D.tif", format="GTiff", overwrite=TRUE)
#writeRaster(Emodnet_e, filename="E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/Depth/E.tif", format="GTiff", overwrite=TRUE)
#writeRaster(Emodnet_f, filename="E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/Depth/F.tif", format="GTiff", overwrite=TRUE)
#rm(list=ls())

#We load bathymetry after processing (the first time is needed to run previous lines and download the data from EMODNET bathymetry)

#C <- raster("E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/Depth/C.tif")
#D <- raster("E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/Depth/D.tif")
#E <- raster("E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/Depth/E.tif")
#F <- raster("E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/Depth/F.tif")

#FinalDepth <- merge(C,D,E,F)

#C <- NULL
#D <- NULL
#E <- NULL
#F <- NULL

#plot(FinalDepth)
#crs(FinalDepth)

#We crop the layer to our intrest area (Bob and Celtic Seas). First, we download the ICES regions from ICES website
#https://gis.ices.dk/sf/index.html

#ICESReg <- readOGR(dsn="E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers", layer="ICES_ecoregions_20171207_erase_ESRI")
#ICESReg_sel <- ICESReg[ICESReg$OBJECTID%in%c(2, 9),]
  
#plot(ICESReg_sel)

#DepthForArea <- crop(FinalDepth, ICESReg_sel)

#ICESReg_sel_LAEA <- spTransform(ICESReg_sel, LAEA)

#LAEA <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "

#beginCluster(n=4)
#Depth <- projectRaster(FinalDepth, crs=LAEA, res=3000, method="bilinear")
#endCluster()

#DepthForArea <- mask(Depth, ICESReg_sel_LAEA)
#plot(DepthForArea)

#writeRaster(DepthForArea, filename="E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/Results/Depth.tif", format="GTiff", overwrite=TRUE)


#In the first 100 lines we have downloaded several depth rasters from EMODNET, merge them, 
#project them to a final res of 3km and mask them using ICES regions
#If you want to repeat the process run previous lines after download the data. Otherwise, start here

Depth <- raster("E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/Results/Depth.tif", format="GTiff", overwrite=TRUE)
plot(Depth)
Depth 

#SEDIMENT TYPE
#The sediment type was also downloaded from EMODNET, in this case EMODNET habitats and rasterize in ARCGIS to a final resolution of 1 km
Habs_rast <- raster("E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/EMODNET/SedsAll.tif")
Habs_rast

Habs_rast_prj <- projectRaster(Habs_rast, Depth, method="ngb")
Habs_rast_prjMask <- crop(Habs_rast_prj, Depth)
Habs_rast_prjMask2 <- mask(Habs_rast_prjMask , Depth)

plot(Habs_rast_prjMask2)

Habs_rast <- Habs_rast_prjMask2 
writeRaster(Habs_rast, filename="E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/Results/SedsEmod.tif", overwrite=TRUE)

#Meaning of the current values in the raster
#0 NA value
#1 Fine mud or Sandy mud or Muddy sand
#2 Sand
#3 Coarse substrate
#4 Mixed sediment
#5 Rock
#6 Fine mud
#7 Muddy sand
#8 Biogenic
#9 Sandy mud
#10 Seabed
#11 Lophelia
#12 Sandy mud or muddy sand
#13 Sediment
#>13 Other biogenics

Habs <- as(Habs_rast, "SpatialPixelsDataFrame")
table(Habs$SedsAll)
summary(Habs$SedsAll)
#First, we pool all substrate (including all biogenic) that we don´t properly sample in just one unique category
Habs$SedsAll<- ifelse(Habs$SedsAll>13,0,Habs$SedsAll)
Habs$SedsAll <- ifelse(Habs$SedsAll==8,0,Habs$SedsAll)
Habs$SedsAll<- ifelse(Habs$SedsAll==11,0,Habs$SedsAll)
Habs$SedsAll<- ifelse(Habs$SedsAll==5,0,Habs$SedsAll)

#Second, we merge unknown substrate categories in just one unique category
Habs$SedsAll <- ifelse(Habs$SedsAll==13,10,Habs$SedsAll)
#Third, we group similar categories (e.g. Muddy sand with Fine mud or Sandy mud or Muddy sand)
Habs$SedsAll <- ifelse(Habs$SedsAll==6,1,Habs$SedsAll)
Habs$SedsAll<- ifelse(Habs$SedsAll==7,1,Habs$SedsAll)
Habs$SedsAll <- ifelse(Habs$SedsAll==9,1,Habs$SedsAll)
Habs$SedsAll <- ifelse(Habs$SedsAll==12,1,Habs$SedsAll)

#Finally, we remove ceros and no sampled substrates
Habs$SedsAll <- ifelse(Habs$SedsAll==0,NA, Habs$SedsAll)
Habs$SedsAll <- ifelse(Habs$SedsAll==10,5, Habs$SedsAll)


table(Habs$SedsAll)

#Final Class
#1 Fine mud or Sandy mud or Muddy sand
#2 Sand
#3 Coarse substrate
#4 Mixed sediment
#5 Sediment or Seabed



plot(raster(Habs[,"SedsAll"]))

FinalSeds <- raster(Habs[,"SedsAll"])


writeRaster(FinalSeds, filename="E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/Results/FinalSeds.tif", overwrite=TRUE)

###################################################################################################################################################################################
###################################################################################################################################################################################
###################################################################################################################################################################################
#We start with the oceanographic models
#CHOLOROPHYL
#The model is: "North Altantic Chlorophyll Concentration from Satellite observations (daily average) Reprocessed L4 (ESA-CCI): monthly"
#with Chl data downloaded from: htt#############################ps://resources.marine.copernicus.eu/product-detail/OCEANCOLOUR_ATL_CHL_L4_REP_OBSERVATIONS_009_091/INFORMATION
###################################################################################################################################################################################
######################################################################################################################################################
###################################################################################################################################################################################

PathToNC_file <-"E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/Chl/Chl_2015_2020_1km.nc"
data_Chl <- nc_open(PathToNC_file)
print(data_Chl)
head(data_Chl)
RasterStack <-brick(PathToNC_file)# Con level 1 elegimos el primer nivel
RasterStack

#We select data only from spring bloom (April and May)

AprMay <- c(4,5,16,17,28,29,40,41,52,53)

RasterStackFilt <- raster::subset(x=RasterStack ,subset=AprMay)
names(RasterStackFilt)

Chl_mean <- stackApply(RasterStackFilt, indices =  rep(1,nlayers(RasterStack)), fun = "mean", na.rm = T)
summary(Chl_mean)

FinalChl <- projectRaster(Chl_mean, Depth, method="bilinear")
plot(FinalChl)
FinalChl

#TEMPERATURE NEAR BOTTON. No regional models covering the whole area. We are using two different models

#Model 1. Atlantic-Iberian Biscay Irish- Ocean Physics Analysis and Forecast
#https://resources.marine.copernicus.eu/product-detail/IBI_ANALYSISFORECAST_PHY_005_001/INFORMATION

PathToNC_file_TempNB_IB <-("E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/TempNB/cmems_mod_ibi_phy_my_0.083deg-3D_P1M-m_1636978149762.nc")
data_TempNB_IB <- nc_open(PathToNC_file_TempNB_IB)
print(data_TempNB_IB)

#Temp near bottom is masked to certain depths, so I am creating it by myself

MeanByDeph <- list()

for (i in 1:50){ #There are 50 depths, each is a level
  
  RasterByDepth <- brick(PathToNC_file_TempNB_IB, var="thetao", level=i)# Con level 1 elegimos el primer nivel
  MeanByDeph[[51-i]] <- stackApply(RasterByDepth, indices =  rep(1,nlayers(RasterByDepth)), fun = "mean", na.rm = T)
  } #We inverse the order for the do.call. It will keep always the values in the first raster, because of this the deeper had to go first

FinalTempNB_IBI <- do.call(merge, MeanByDeph)

FinalTempNB_IBI 
plot(FinalTempNB_IBI)

FinalTempNB_Ibi_Prj <- projectRaster(FinalTempNB_IBI, Depth, method="bilinear")
plot(FinalTempNB_Ibi_Prj)

#Model 2. Atlantic - European North West Shelf - Ocean Physics Analysis and Forecast
#https://resources.marine.copernicus.eu/product-detail/NWSHELF_MULTIYEAR_PHY_004_009/INFORMATION

PathToNC_file_TempNB_CS <-("E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/TempNB/cmems_mod_nws_phy-bottomt_my_7km-2D_P1M-m_1636977209627.nc")
data_TempNB_CS <- nc_open(PathToNC_file_TempNB_CS)
print(data_TempNB_CS)
head(data_TempNB_CS)
RasterStack_TempNB_CS <-brick(PathToNC_file_TempNB_CS, var="bottomT")

TempNB_CSi_mean <- stackApply(RasterStack_TempNB_CS, indices =  rep(1,nlayers(RasterStack_TempNB_CS)), fun = "mean", na.rm = T)
summary(TempNB_CSi_mean)

FinalTempNB_CSi <- projectRaster(TempNB_CSi_mean, Depth, method="bilinear")
plot(FinalTempNB_CSi)

AllTemp <- merge(FinalTempNB_CSi, FinalTempNB_Ibi_Prj)
plot(AllTemp)

#layers from EMODNET habitats (https://www.emodnet-seabedhabitats.eu/access-data/launch-map-viewer/)
#Currents
Currents_Celtic <- raster("E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/EMODNET/ke_currents_atlantic.tif") 
crs(Currents_Celtic) <- "+proj=longlat +datum=WGS84 +no_defs"
Currents_BoB_N <- raster("E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/EMODNET/ke_currents_manga2500_2010_2015.tif")
crs(Currents_BoB_N) <- "+proj=longlat +datum=WGS84 +no_defs"
Currents_Iberia <- raster("E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/EMODNET/ke_currents_iberian_be_mesh_atlantic_maretec.tif")
crs(Currents_Iberia) <- "+proj=longlat +datum=WGS84 +no_defs"
Currents_Celtic_Proj <- projectRaster(Currents_Celtic, Depth, method="bilinear")
Currents_BoB_N_Proj <- projectRaster(Currents_BoB_N, Depth, method="bilinear")
Currents_Iberia_Proj <- projectRaster(Currents_Iberia, Depth, method="bilinear")
plot(stack(Currents_Celtic_Proj, Currents_BoB_N_Proj, Currents_Iberia_Proj))

Currents_Final <- merge(Currents_Celtic_Proj, Currents_BoB_N_Proj, Currents_Iberia_Proj)#For pixels share by more than 1 raster, merge follow the order in the code


plot(Currents_Final)

#Waves
Waves_Celtic <- raster("E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/EMODNET/ke_waves_atlantic.tif")
crs(Waves_Celtic) <- "+proj=longlat +datum=WGS84 +no_defs"
Waves_BoB <- raster("E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/EMODNET/ke_wave_homere_ubr_2010_2015.tif")
crs(Waves_BoB) <- "+proj=longlat +datum=WGS84 +no_defs"
Waves_Iberia <- raster("E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/OriginalLayers/EMODNET/iberian_ke_wave.tif")
crs(Waves_Iberia) <- "+proj=longlat +datum=WGS84 +no_defs"

Waves_Celtic_Proj <- projectRaster(Waves_Celtic, Depth, method="bilinear")
Waves_BoB_Proj <- projectRaster(Waves_BoB, Depth, method="bilinear")
Waves_Iberia_Proj <- projectRaster(Waves_Iberia, Depth, method="bilinear")
Waves_Final <- merge(Waves_Celtic_Proj, Waves_BoB_Proj, Waves_Iberia_Proj)
plot(Waves_Final)

#FinalEnergy. We add energy from Currents and Energy from Waves

EnergyFinal <- Currents_Final + Waves_Final

#AllLayersStack
FinalStack <- stack(Depth, FinalChl, AllTemp, Currents_Final, Waves_Final, EnergyFinal, FinalSeds)
names(FinalStack) <- c("Depth", "Chl", "Temp", "Currents", "Waves", "AllEnergy", "Substrate")

#We maskes using Depth for areas shallower than 900 m

Depth_SPDF <- as(Depth, "SpatialPixelsDataFrame")
Depth_SPDF$Depth <- ifelse(Depth_SPDF$Depth>= (-900), Depth_SPDF$Depth, NA)
MyMask <- raster(Depth_SPDF[,"Depth"])

FinalStack <- mask(FinalStack, MyMask)

windows()
plot(FinalStack)

writeRaster(FinalStack, filename="E:/OneDrive/IEO/Congresos&WK/WGFBIT/Palermo_2021/TraitAnalysis/Results/AllEnvLayers.tif", options="INTERLEAVE=BAND", overwrite=TRUE)


