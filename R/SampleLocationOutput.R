library(raster);library(rgdal);library(tripack);library(SDMTools); library(manipulate);library(clhs);library(entropy);library(ggplot2)

InPSL <- c()
for (BGZ in c(7:14)) { #this will be 14 when we run this properly
  BGZ.St <- Sys.time()
  print(paste0("Starting Analysis of BGZ ",BGZ))
  #BGZ <- 2 #usually comment this out!
  
  #generate Dataframe of covariates
  ####create stacks of Covariates & groups of Covariates####
  CoV_SAR <- stack(lapply(list.files(path = "./AlignedNationWide",
                                     pattern = paste0("SAR_", BGZ,".tif$"),
                                     full.names = T, recursive = T), brick))
  
  CoV_DTM <- stack(lapply(list.files(path = "./AlignedNationWide",
                                     pattern = paste0("DTM_", BGZ,".tif$"), 
                                     full.names = T, recursive = T), raster))
  
  CoV_Slope <- stack(lapply(list.files(path = "./AlignedNationWide",
                                       pattern = paste0("Slope_", BGZ,".tif$"), 
                                       full.names = T, recursive = T), raster))
  
  CoV_PSL <- stack(lapply(list.files(path = "./AlignedNationWide",
                                     pattern = paste0("PSL_", BGZ,".tif$"),
                                     full.names = T, recursive = T), brick))
  
  CoV_SupGeo <- stack(lapply(list.files(path = "./AlignedNationWide",
                                        pattern = paste0("SupGeo_", BGZ,".tif$"),
                                        full.names = T, recursive = T), raster))
  
  CoV_BedGeo <- stack(lapply(list.files(path = "./AlignedNationWide",
                                        pattern = paste0("BedGeo_", BGZ,".tif$"), 
                                        full.names = T, recursive = T), raster))
  
  CoV_All <- stack(CoV_SAR, CoV_DTM, CoV_Slope, CoV_PSL, CoV_SupGeo, CoV_BedGeo) #RasterStack of All of the CoVariates
  
  #########################################
  #Check all covariates are loaded
  names(CoV_All)
  
  #creating a DataFrame of Covariates####
  CoV_DF <- as.data.frame(rasterToPoints(CoV_All)) #change the Cov_x layer as required
  CoV_DF[, grep("SupGeo", colnames(CoV_DF))] <- as.factor(CoV_DF[, grep("SupGeo", colnames(CoV_DF))])
  CoV_DF[, grep("Bed", colnames(CoV_DF))] <- as.factor(CoV_DF[, grep("Bed", colnames(CoV_DF))])
  
  #check for and remove columns that have only 1 level e.g all 0 (no coverage) in one of the PSLs
  no_uniques <- function(vec) length(unique(na.omit(vec))) > 1
  CoV_DF <- CoV_DF[,sapply(CoV_DF,no_uniques)]

#####
  #list of sample sizes by BGZ
  SampSizes <- c(252, 255, 257, 301, 330, 316, 295, 293, 277, 301, 350, 296, 299, 220)
   
  #PERFORM CLHS SIZE = NUMBER OF SAMPLES ITER = NUMBER OF ATTEMPTS TO COVER COVARIATE SPACE
  res <- clhs(CoV_DF[3:11], size = (SampSizes[BGZ]), iter = 100, progress = TRUE) #change size to var.
  
  #ADD CLHS RESULTS TO DATAFRAME CONTAINING XY VALUES & COVARIATES USED IN MODEL
  samplepoints <- CoV_DF[c(res),]
  
  #DF TO POINT SHP
  SPsamplepoints <- SpatialPoints(samplepoints[1:2], CRS("+init=epsg:27700"))
  
  shapefile(SPsamplepoints, file = paste0("F:/Projects/EPM/data/ModelOutputs/",BGZ,"/BGZ",BGZ,"_SamplePoints.shp"))

  #select PSL, add columns count non 0
  InPSL <- c(InPSL, sum(rowSums(samplepoints[,(grep("PSL_...$", names(samplepoints)))] > 0, na.rm = TRUE)))

  print(paste("Time to Complete BGZ",BGZ," = ", lubridate::as.duration(Sys.time() - BGZ.St)))
  
  #Clear object to make room in memory
  res <- c()
}

sum(InPSL)
  