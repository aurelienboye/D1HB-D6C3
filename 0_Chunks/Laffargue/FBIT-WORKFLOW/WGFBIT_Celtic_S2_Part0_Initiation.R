#####################################
######### WGFBIT - October 2019 ##########
##################################### >)))°>
######## SCRIPT dataset compilation - Celtic/BoB Eco-region (WGFBIT)
### R version 3.6.2

#################################
##
### PART 0 - INITIATION
##
#################################

##

###############
####


#################################
### SUMMARY
#################################


#### 0 - Libraries & Functions

 ### 0.1 - R Libraries
 ### 0.2 - Local Fonctions
  
#### 1 - Data selection sources
                      
 ### 1.1 - Data selection sources: File paths
 
 ### 1.2 - Data selection sources: Files names
 
  # 1.2.1 - EUNIS/MSFD
  # 1.2.2 - Depth
  # 1.2.3 - Fishing
  # 1.2.4 - Medata for biological datasets
  # 1.2.5 - BTA table for each kind of biological dataset (MEGAFAUNA or MACROFAUNA)

 
 ### 1.3 - R Libraries
 
 ### 1.4 - Local Fonctions


#### 2 - DATASET download and formatting
# (WGFBIT - Step 0.xx)
 
 ### 2.1 - BTA DATASET (STEP 3)
  
  ## 2.1.0 - load BTA data
  
  ## 2.1.1 - Correct/Adapt traits matrix
 
    # 2.1.1.1 - Check / Suppress EMPTY Row
  
    # 2.1.1.2 - Update Species list
     # reuses by default the file already existing in the source directory (SEL.WORMSnew=F)
 
  ## 2.1.2 - Redundancies test and suppression within functional matrix
                                             
  ## 2.1.4 - Identification / automatic completion of taxa with missing information
  
  ## 2.1.5 - Complete new traits/modalitities from similar traits
    
    # 2.1.5.1 - NEW LONGEVITY CLASSIFICATION (6 classes "LG" from 4 classes "LO")

    # => OUTPUT FINAL TRAIT MATRIX (links to FBIT STEP 2, line 7)
         # Creation new corrected funtional matrix (MATRIX.FUNCnew)
         # & NAME CORRECTION from worms "valid_name"

          
 ### 2.2 - BIOLOGICAL DATASET
# (WGFBIT - STEP 3)
  
  ## 2.2.1 - METADATA
  
  ## 2.2.2 - Specific DATASET adds/modifications

  ## 2.2.2.1 - EVHOE
  
  
  ## 2.2.3 - CONVERT / FORMAT biological dataset with standard field names
  
  ## 2.2.5 -  Species list filtering
   # 2.2.5.1 - configure the species filter
   # 2.2.5.2 - Apply filter by SURVEY
   
# => OUTPUT FINAL DATASET BIOLOGY (links to FBIT STEP 2, line 7)
# check for standard database format/fields (FBIT script)

# + mean biomass by C_SQUARE
# Species/Biomass/traits1/Traits2/traits3/traits4


 
 ### 2.3 - HABITAT DATASET "EUNIS" (MSFD broad habitat types)
  
  ## 2.3.1 - FILTER / MODIFY EUNIS HABITATS TYPES
  
  ## 2.3.2 - Region selection
   # 2.3.2.1 - Define polygon boundaries for subregion
   # 2.3.2.2 - Select the subregion to be treated
  
  ## 2.3.3 -  Surface by subregion and habitat type
  
  ## 2.3.4 - Sampling stations OVER EUNIS or MSFD habitats 
   # 2.3.4.1 - Format survey coordinates
   # 2.3.4.2 -  Attribute EUNIS/MSFD habitats to stations
   # 2.3.4.3 -  Attribute subregions to stations

 
 ### 2.4 - Depth dataset
 
  ## 2.4.1 - Attribute depths to stations
   # 2.1.4.1 - Search for missing traits information
    
  ## 2.4.2 - Habitat mean depth
    
        
    ## TO COMPLETE 
        ## TO COMPLETE
            ## TO COMPLETE
                ## TO COMPLETE
 
 
 
 ### 2.5 - Fishing pressure dataset
 
     ## TO COMPLETE
         ## TO COMPLETE
             ## TO COMPLETE
                 ## TO COMPLETE

 ### 2.6 - C-Square dataset

     ## TO COMPLETE
         ## TO COMPLETE
             ## TO COMPLETE
                 ## TO COMPLETE




#### 3 - PRELIMINARY tests

 ### 3.1 - Number of stations for each habitat TYPES
       ## 3.1.1 - covered Habitats
       ## 3.1.2 - Map / sampling coverage

 ### 3.2 - BTA matrix coverage for species lists by main habitats
       ## 3.2.1 -  Species richness coverage by habitat type
       ## 3.2.2 - Biomass coverage by habitat type

 ### 3.3 - Biomass Distribution / Longevity classes & Habitat
       ## 3.3.1 - For all stations with biological observations
       ## 3.3.2 - fuzzy pca to check for species/traits distribution for each habitat type ...
       ## 3.3.3 - with NULL (or low fishing pressure

#### 4 - Formatting and writing files for each FBIT step
  
       
###################################################################################################       
###################################################################################################
###################################################################################################
###

rm(list = ls(all=TRUE))


#### 0 - Libraries & Functions
 ### 0.1 - R Libraries

library(PBSmapping)
library(sp)
library(maptools)
library(rgeos)
library(raster)
library(ade4)
library(lme4)
library(remotes)

## Package VMS from ICES
install_github("ices-tools-prod/icesConnect")
install_github("ices-tools-prod/icesVMS")

 ### 0.2 - Local Fonctions

 # 0.2.1 - loads an RData file, and returns it
loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])}



## Options

SEL.ScriptName <- "WGFBIT_Celtic_S2_"  # useful to indicate script version into reports

SEL.PDF <- T # Select to create (T) or not (F) the pdfs reports

SEL.WORMSnew <- F # Utilize to update species lists and taxonomy from Worms online database. If T !!! WARNING !!! Time consuming !!!

# Select BTA dataset
SEL.TRAIT <- c("LO") # Select list of traits to considere (!!!! TO CHECK : LO or LG !!!!!)


# SEL.COMPONENT <- "MEGAFAUNA" # MEGAFAUNA or MACROFAUNA 
# selection of component to consider, utilized to select relevant BTA matrix

# SEL.WORMSnew <- F # if T !!! WARNING !!! Time consuming !!!
# Execution of taxa name update/correction from worms external sources

#
#

##


#### 1 - Data selection sources

### 1.1 - Data selection sources: File paths

DIR.SOURCE <- list()
DIR.SOURCE[[paste(0)]] <-  file.path("E:") # Local, network, github ...
DIR.SOURCE[[paste(1)]] <-  file.path(DIR.SOURCE[[paste(0)]],"PlaffIFREMER/4.1-ICES/ICES-WGFBIT/WorkInProgress/WesternWaters") # ROOT DIRECTORY

DIR.SOURCE[[paste(2)]] <-  file.path(DIR.SOURCE[[paste(1)]],"DATAlayers") # DATA
DIR.SOURCE[[paste(3)]] <-  file.path(DIR.SOURCE[[paste(2)]],"DATASET_Biological") # Biological Dataset
DIR.SOURCE[[paste(4)]] <-  file.path(DIR.SOURCE[[paste(2)]],"DATASET_FishingPressure_OSPAR") # Fishing pressure layer
DIR.SOURCE[[paste(5)]] <-  file.path(DIR.SOURCE[[paste(2)]],"DATASET_EnvironmentPressure") # Environment pressure layer
DIR.SOURCE[[paste("5a")]] <-  file.path(DIR.SOURCE[[paste(5)]],"Depth") # Depth

DIR.SOURCE[[paste("5b")]] <-  file.path(DIR.SOURCE[[paste(2)]],"DATASET_GRID") # 

DIR.SOURCE[[paste(6)]] <-  file.path(DIR.SOURCE[[paste(2)]],"DATASET_HABITAT","EUSEAMAP") # EUNIS HABITATS
DIR.SOURCE[[paste(7)]] <-  file.path(DIR.SOURCE[[paste(2)]],"outputs")
DIR.SOURCE[[paste(8)]] <-  file.path(DIR.SOURCE[[paste(7)]],"FBIT_WORKFLOW_formattedFILES/STEP2/Benthic data")
DIR.SOURCE[[paste(9)]] <-  file.path(DIR.SOURCE[[paste(1)]],"SCRIPTS/FUNCTIONS")

DIR.OUPUT <- list()                                                                                                             
DIR.OUPUT[[paste(1)]] <- file.path(DIR.SOURCE[[paste(1)]],"OUTPUTS") 


### 1.2 - Data selection sources: Files names

# 1.2.1 - EUNIS/MSFD
SOURCEFILE.EUNIS <-  "EUSeaMap2016_WesternWaters.shp"

# 1.2.2 - Depth
SOURCEFILE.DEPTH <- "GebcoDATA_2014_RASTER_TOTAL.RData"

# 1.2.3 - Fishing
REFname.FISHING <-  "OSPAR_intensity_total"
SOURCEFILE.FISHING <- list.files(DIR.SOURCE[[paste(4)]])[grep(REFname.FISHING,list.files(DIR.SOURCE[[paste(4)]]))]

# 1.2.4 - Medata for biological datasets
FILE.NAME <- list()
FILE.NAME[["DATASET_Biological"]][["METADATA"]] <- c("METADATA_UK_SB.csv")
FILE.NAME[["BioMetadata"]] <- "BioDataset_Metadata.csv"

# 1.2.5x - Selection of fauna types (MEGAFAUNA / MACROFAUNA)
 SEL.COMPONENT <- "MEGAFAUNA" # MEGAFAUNA or MACROFAUNA






 ## LOAD SPECIFIC FBIT FUNCTIONS
source(file.path(DIR.SOURCE[[paste(9)]],"Impact_continuous_longevity.R"))
source(file.path(DIR.SOURCE[[paste(9)]],"Habitatstatefishing.R"))

 ### 0.3 - External Fonctions
   # Source directory for functions
DIRECTORY.ScriptSource <- "E:/PlaffIFREMER/R-SCRIPT/FUNCTIONS"

   # Update species list from worms
source(file.path(DIR.SOURCE[[paste(9)]],"RFunction_WORMSextract.R"))


#


###########################################################
###########################################################
###########################################################
###########################################################