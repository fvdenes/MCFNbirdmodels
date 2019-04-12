library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)
library(sp)
library(maptools)
library(rgeos)
library(tidyr)
library(beepr)

# Load BAM avian dataset ####

# including raw data for coni, hogr and seow
load("G:/My Drive/BAM.SharedDrive/DataStuff/Avian.Data/Processed/BAM-BBS-tables_20170630.Rdata")

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
coordinates(SS) <- c("X", "Y") 
proj4string(SS) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SSLCC <- spTransform(SS, LCC)

offl <- data.table(melt(OFF))
names(offl) <- c("PKEY","SPECIES","logoffset")
offl$SPECIES <- as.character(offl$SPECIES)
offl$PKEY <-as.character(offl$PKEY)
rm(OFF)


# Read MCFN homeland shp, convert to raster ####
MCFNhome_shp<-readOGR(dsn="D:/MCFN",layer="Homelands_Boundary_10252016")
MCFNhome_shp<-spTransform(MCFNhome_shp,LCC)
#writeOGR(LSFNhome_shp,dsn="D:/LSFN/Traditional_Territory_LacSeul_FN",layer="Traditional_Territory_Boundary_LacSeulFN_LCC", driver="ESRI Shapefile")

# Load Beaudoin layers (2011 and 2001), crop and mask for study region, save as rasters ####
b2011 <- list.files("D:/Beaudoin/2011/",pattern="tif$")
setwd("D:/Beaudoin/2011/")
bs2011 <- stack(raster(b2011[1]))
for (i in 2:length(b2011)) { bs2011 <- addLayer(bs2011, raster(b2011[i]))}
names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))

# Create buffers from MCFN homeland shp, convert to rasters, resample to Beaudoin, crop to max buffer extent, and create raster stack with all buffers #### 
incBuff<-c("0","50000","100000","150000","200000","250000","300000","350000","400000","450000","500000")

list1 <- lapply(incBuff, function(x) {
  srBuff <-gBuffer(MCFNhome_shp, width=as.numeric(x),quadsegs=35)
  r <- raster(ncol=500, nrow=500)
  extent(r)<-extent(srBuff)
  srBuff_rast <- rasterize(srBuff, r, 1)
  srBuff_rast_re <- resample(srBuff_rast,bs2011[[1]])
})

srBuff500000 <-gBuffer(MCFNhome_shp, width=as.numeric(300000),quadsegs=35) # re-create 5000000m buffer to use as cropping extent
srBuff_rast_crop<-crop(list1[[1]],srBuff500000)
MCFNbuffers<-stack(srBuff_rast_crop)
for(i in 2: length(list1)){MCFNbuffers <-addLayer(MCFNbuffers,crop(list1[[i]],srBuff500000))}

names(MCFNbuffers)<-paste0("MCFNBuff",incBuff)
writeRaster(MCFNbuffers, filename="D:/MCFN/rasters/MCFNbuffers.grd", format="raster",overwrite=TRUE)
#MCFNbuffers<-brick("D:/MCFN/rasters/MCFNbuffers.grd")

# New MCFN avian data ####
load("D:/MCFN/Additional bird data/MCFN_data_processed-20190325.RData")

#create spatial dataframe to project coords into LCC
SS_MCFNsp<-SpatialPointsDataFrame(cbind(SS_MCFN$X,SS_MCFN$Y),SS_MCFN, proj4string = CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"))
SS_MCFNsp<-spTransform(SS_MCFNsp,LCC)

names(SSLCC)

names(SS_MCFNsp)
SS_MCFNsp$Xcl<-SS_MCFNsp@coords[,1]
SS_MCFNsp$Ycl<-SS_MCFNsp@coords[,2]

SSMCFN_all<-rbind(SSLCC[,1:4],SS_MCFNsp[,c(2,1,8,9)])

#ADD NEW DATA TO PCTBL, OFFSET and PKEY DATAFRAMES
PCTBL_all<-PCTBL[,c(1:3,5,10)]
names(PCTBL_all)[5]<-"SPECIES"
PCTBL_all<-rbind(PCTBL_all,PCTBL_MCFN[,c(1,3,2,4:5)])


OFF_MCFN<-as.data.frame(OFF_MCFN)
OFF_MCFN$PKEY<-rownames(OFF_MCFN)


OFF_MCFN_long<-gather(OFF_MCFN,"SPECIES","logoffset",ALFL:YTVI,factor_key = F) # convert dataframe from wide to long format
OFF_ALL<-rbind(offl,OFF_MCFN_long)

names(PKEY_MCFN)[4:8]<-c("YEAR","MONTH","DAY","HOUR","MIN")
PKEY_all<-rbind(PKEY[,c(1:3,8)],PKEY_MCFN[,1:4])

# Extract point count SS for each buffer layer #### 
list2<- lapply(incBuff, function(x){
  srBuff <-gBuffer(MCFNhome_shp, width=as.numeric(x),quadsegs=35)
  srBuff<-SpatialPolygonsDataFrame(srBuff,data=data.frame(ID=1),match.ID=F)
  writeOGR(srBuff,dsn="D:/MCFN/buffers",layer=paste0("MCFNBuff",x), driver="ESRI Shapefile")
  MCFNSS <- crop(SSMCFN_all,srBuff)
  writeOGR(MCFNSS, dsn="D:/MCFN/buffers", layer=paste0("BAM_ptswBuff",x), driver="ESRI Shapefile")
  MCFNSS <- as.data.frame(MCFNSS)
})

names(list2)<-paste0("MCFNSS",incBuff)

beep(sound = 4, expr = NULL) # plays a sound to indicate command above finished running

# Finish importing Beaudoin layers 2011 and 2001 ####
#2011
MCFNbs2011 <- crop(bs2011,MCFNbuffers[[1]])
#abs2011<-mask(bs2011,MCFNbuffers[[1]])


#2001
b2001 <- list.files("D:/Beaudoin/2001/",pattern="tif$")
setwd("D:/Beaudoin/2001/")
bs2001 <- stack(raster(b2001[1]))
for (i in 2:length(b2001)) { bs2001 <- addLayer(bs2001, raster(b2001[i]))}
names(bs2001) <- gsub("NFI_MODIS250m_2001_kNN_","",names(bs2001))
MCFNbs2001 <- crop(bs2001,MCFNbuffers[[1]])
#bs2001<-mask(bs2001,MCFNbuffers[[1]])


# Calculate standard deviation of Beaudoin stacks to identify layers that do not occur in the study region, remove them
sd.test<-data.frame(var=names(MCFNbs2001[[5:79]]),sd=NA)
for (i in 1:nrow(sd.test)){
  sd.test[i,2]<-sd(getValues(MCFNbs2001[[i+5]]),na.rm = T)
}
drop<-which(names(MCFNbs2001)%in% sd.test$var[sd.test$sd==0])

MCFNbs2001<-dropLayer(MCFNbs2001,drop)
names(MCFNbs2001)
MCFNbs2011<-dropLayer(MCFNbs2011,drop)
names(MCFNbs2011)

writeRaster(MCFNbs2001, filename="D:/Beaudoin/2001/Processed/MCFN/bs2001_250m.grd", format="raster",overwrite=TRUE)
writeRaster(MCFNbs2011, filename="D:/Beaudoin/2011/Processed/MCFN/bs2011_250m.grd", format="raster",overwrite=TRUE)



# water layer - load NA land cover raster, reclassify to flag only water bodies, then crop to MCFN  ####
NALCC2005<-raster("M:/DataStuff/SpatialData/LCC05_NALCMS2010/Land_Cover_MXD/NA_LandCover_2005/data/NA_LandCover_2005/NA_LandCover_2005_LCC.img")
plot(NALCC2005)
m <- c(1:19,rep(0,17),1,0)
rclmat <- matrix(m, ncol=2, byrow=F)
wat<-reclassify(NALCC2005,rclmat)
watMCFN<- crop(wat,MCFNbuffers[[1]])
watMCFN<-resample(watMCFN,MCFNbuffers[[1]], method="ngb")
writeRaster(watMCFN, filename="D:/wat2011_lcc1/MCFN/watMCFN.tif", format="GTiff",overwrite=TRUE)
names(watMCFN)<-"water"

# CTI data #### NEED to ASSEMBLE MOSAIC IN ARCGIS
CTI<-raster("D:/CTI/CTI_MCFN.tif")
CTI<-projectRaster(CTI,MCFNbuffers[[1]])
CTI250<-resample(CTI,MCFNbuffers[[1]])

# Climate data- upload and resample to 250m resolution to match other layers ####
climate2010 <- list.files("D:/ClimateAdaptWest/baseline19812010/",pattern="asc$")
setwd("D:/ClimateAdaptWest/baseline19812010/")
clim2010 <- stack(raster(climate2010[1]))
for (i in 2:length(climate2010)) { clim2010 <- addLayer(clim2010, raster(climate2010[i]))}
proj4string(clim2010)<-LCC
MCFNclim2010 <- crop(clim2010,MCFNbuffers[[1]])
MCFNclim2010<-resample(MCFNclim2010,MCFNbuffers[[1]])
#aclim2010<-mask(aclim2010,abs2011$LandCover_NonVeg_v1)
plot(MCFNclim2010) 
writeRaster(MCFNclim2010, filename="D:/ClimateAdaptWest/baseline19812010/Processed/MCFN/MCFNclim2010.grd", format="raster",overwrite=TRUE)
MCFNclim2010<-stack("D:/ClimateAdaptWest/baseline19812010/Processed/MCFN/MCFNclim2010.grd")


beep(sound = 4, expr = NULL)

# Create Gaussian smoother with sigma = 750m for numeric predictors (Beaudoin, CTI)####
## 2001, sigma = 750m
fw750<-focalWeight(x=MCFNbs2001,d=750,type="Gauss")
MCFNbs2001_Gauss750<-brick(focal(MCFNbs2001[[1]],w=fw750,na.rm=TRUE))
names(MCFNbs2001_Gauss750)<-paste(names(MCFNbs2001)[[1]],"_Gauss750",sep="")  
for(i in 2:nlayers(MCFNbs2001)){
  MCFNbs2001_Gauss750<-addLayer(MCFNbs2001_Gauss750,focal(MCFNbs2001[[i]],w=fw750,na.rm=TRUE))
  names(MCFNbs2001_Gauss750)[i]<-paste(names(MCFNbs2001)[[i]],"_Gauss750",sep="")
}

writeRaster(MCFNbs2001_Gauss750, filename="D:/Beaudoin/2001/Processed/MCFN/MCFNbs2001_250_Gauss750m.grd", format="raster",overwrite=TRUE)
#MCFNbs2001_Gauss750<-brick("D:/Beaudoin/2001/Processed/MCFN/MCFNbs2001_250_Gauss750m.grd")
beep(sound = 4, expr = NULL)

## 2011, sigma = 750m
#fw750<-focalWeight(x=MCFNbs2011,d=750,type="Gauss")
MCFNbs2011_Gauss750<-brick(focal(MCFNbs2011[[1]],w=fw750,na.rm=TRUE))
names(MCFNbs2011_Gauss750)<-paste(names(MCFNbs2011)[[1]],"_Gauss750",sep="")  
for(i in 2:nlayers(MCFNbs2011)){
  MCFNbs2011_Gauss750<-addLayer(MCFNbs2011_Gauss750,focal(MCFNbs2011[[i]],w=fw750,na.rm=TRUE))
  names(MCFNbs2011_Gauss750)[i]<-paste(names(MCFNbs2011)[[i]],"_Gauss750",sep="")
}

writeRaster(MCFNbs2011_Gauss750, filename="D:/Beaudoin/2011/Processed/MCFN/MCFNbs2011_250_Gauss750m.grd", format="raster",overwrite=TRUE)
#MCFNbs2011_Gauss750<-brick("D:/Beaudoin/2011/Processed/MCFN/MCFNbs2011_250_Gauss750m.grd")
beep(sound = 4, expr = NULL)

## CTI
CTI250_Gauss750<-focal(CTI250,w=fw750,na.rm=TRUE)
names(CTI250_Gauss750)<-"CTI_MCFN_Gauss750"

# put together prediction rasterstack
pred_MCFN<-stack(MCFNbs2011,MCFNbs2011_Gauss750,CTI250,CTI250_Gauss750,watMCFN, MCFNclim2010, quick=TRUE)
#writeRaster(pred_MCFN,filename="D:/MCFN/rasters/Prediction dataset/pred_dat_MCFN", format="raster",overwrite=TRUE)

pred_MCFN_mask<-mask(pred_MCFN,MCFNbuffers[[1]])
pred_MCFN_mask<-crop(pred_MCFN_mask,MCFNhome_shp)

writeRaster(pred_MCFN_mask,filename="D:/MCFN/rasters/Prediction dataset/pred_MCFN_mask", format="raster",overwrite=TRUE)

save.image("D:/MCFN/ws1.RData")
#load("D:/MCFN/ws1.RData")
beep(sound = 4, expr = NULL) # plays a sound to indicate command above finished running
# Extracting data from rasters for surveyed locations for each buffer extent ####

# 
# # Reload rasters if connection is not working
 # MCFNbs2011<-brick("D:/Beaudoin/2011/Processed/MCFN/bs2011_250m.grd")
 # MCFNbs2011_Gauss750<-brick("D:/Beaudoin/2011/Processed/MCFN/MCFNbs2011_250_Gauss750m.grd")
 # MCFNclim2010<-brick("D:/ClimateAdaptWest/baseline19812010/Processed/MCFN/MCFNclim2010.grd")
 # MCFNbs2001<-brick("D:/Beaudoin/2001/Processed/MCFN/bs2001_250m.grd")
 # MCFNbs2001_Gauss750<-brick("D:/Beaudoin/2001/Processed/MCFN/MCFNbs2001_250_Gauss750m.grd")

list3<-list(NA)

for (j in 1:length(list2)){
  x<-list2[[j]]
  
  #2011
  dat2011 <- cbind(x, raster::extract(MCFNbs2011,as.matrix(cbind(x$X,x$Y)))) 
  for(i in 1:nlayers(MCFNbs2011_Gauss750)){
    dat2011 <- cbind(dat2011, raster::extract(MCFNbs2011_Gauss750[[i]],as.matrix(cbind(x$X,x$Y)))) # includes Beaudoin layers with Gaussian filter sigma=750m
    names(dat2011)[ncol(dat2011)] <- names(MCFNbs2011_Gauss750)[i]
  }
  dat2011<-cbind(dat2011,raster::extract(CTI250,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes CTI 250m resolution data
  names(dat2011)[ncol(dat2011)] <- "CTI_MCFN"
  dat2011<-cbind(dat2011,raster::extract(CTI250_Gauss750,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes CTI 250m resolution data, Gaussian filter sigma=750m
  names(dat2011)[ncol(dat2011)] <- "CTI_MCFN_Gauss750"
  dat2011<-cbind(dat2011,raster::extract(watMCFN,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes wat 250m resolution data
  names(dat2011)[ncol(dat2011)] <- "water"
  
  
  for(i in 1:nlayers(MCFNclim2010)){
    dat2011 <- cbind(dat2011, raster::extract(MCFNclim2010[[i]],as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes climate layers 
    names(dat2011)[ncol(dat2011)] <- names(MCFNclim2010)[i]
  }
  
  ### set up weight matrix for SS, and calculate weight values for each row in dat2011
  r2 <- MCFNbs2011[[1]]
  samprast2011 <- rasterize(cbind(dat2011$X,dat2011$Y), r2, field=1)
  sampsum25 <- focal(samprast2011, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE)
  
  dat2011 <- cbind(dat2011,raster::extract(sampsum25,as.matrix(cbind(dat2011$X,dat2011$Y))))
  names(dat2011)[ncol(dat2011)] <- "sampsum25"
  dat2011$wt <- 1/dat2011$sampsum25
  dat2011<-dat2011[which(!is.na(dat2011$wt)),]
  
  dat2011$SS <- as.character(dat2011$SS)
  dat2011$PCODE <- as.character(dat2011$PCODE)
  
  
  setwd("D:/MCFN/buffers")
  write.csv(dat2011,paste0("dat2011_buff",incBuff[j],".csv"))
  
  #2001
  dat2001 <- cbind(x, raster::extract(MCFNbs2001,as.matrix(cbind(x$X,x$Y)))) 
  for(i in 1:nlayers(MCFNbs2001_Gauss750)){
    dat2001 <- cbind(dat2001, raster::extract(MCFNbs2001_Gauss750[[i]],as.matrix(cbind(x$X,x$Y)))) # includes Beaudoin layers with Gaussian filter sigma=750m
    names(dat2001)[ncol(dat2001)] <- names(MCFNbs2001_Gauss750)[i]
  }
  dat2001<-cbind(dat2001,raster::extract(CTI250,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes CTI 250m resolution data
  names(dat2001)[ncol(dat2001)] <- "CTI_MCFN"
  dat2001<-cbind(dat2001,raster::extract(CTI250_Gauss750,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes CTI 250m resolution data, Gaussian filter sigma=750m
  names(dat2001)[ncol(dat2001)] <- "CTI_MCFN_Gauss750"
  dat2001<-cbind(dat2001,raster::extract(watMCFN,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes wat 250m resolution data
  names(dat2001)[ncol(dat2001)] <- "water"
  
  for(i in 1:nlayers(MCFNclim2010)){
    dat2001 <- cbind(dat2001, raster::extract(MCFNclim2010[[i]],as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes climate layers with Gaussian filter sigma=750m
    names(dat2001)[ncol(dat2001)] <- names(MCFNclim2010)[i]
  }
  
  ### set up weight matrix for SS, and calculate weight values for each row in dat2001
  r2 <- MCFNbs2001[[1]]
  samprast2001 <- rasterize(cbind(dat2001$X,dat2001$Y), r2, field=1)
  sampsum25 <- focal(samprast2001, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE)
  
  dat2001 <- cbind(dat2001,raster::extract(sampsum25,as.matrix(cbind(dat2001$X,dat2001$Y))))
  names(dat2001)[ncol(dat2001)] <- "sampsum25"
  dat2001$wt <- 1/dat2001$sampsum25
  dat2001<-dat2001[which(!is.na(dat2001$wt)),]
  
  
  dat2001$SS <- as.character(dat2001$SS)
  dat2001$PCODE <- as.character(dat2001$PCODE)
  
  
  setwd("D:/MCFN/buffers")
  write.csv(dat2001,paste0("dat2001_buff",incBuff[j],".csv"))
  
  # Prepare point count data for each SS and aggregate for 2001 and 2011.####
  
  PC <- inner_join(PCTBL_all,PKEY_all,by=c("PKEY","SS"))[,-1]
  colnames(PC)[5]<-"PCODE"
  #PC <- inner_join(PC,SSMCFN_all@data[,2],by="SS")
  
  
  
  
  MCFNPC <- PC[which(PC$SS%in%list2[[j]]$SS),]
  
  MCFNPC$SS <- as.character(MCFNPC$SS)
  MCFNPC$PKEY <- as.character(MCFNPC$PKEY)
  MCFNPC$PCODE <- as.character(MCFNPC$PCODE)
  MCFNPC$SPECIES <- as.character(MCFNPC$SPECIES)
  MCFNPC2001 <- MCFNPC[MCFNPC$YEAR < 2006,] 
  MCFNPC2011 <- MCFNPC[MCFNPC$YEAR > 2005,] 
  survey2001 <- aggregate(MCFNPC2001$ABUND, by=list("PKEY"=MCFNPC2001$PKEY,"SS"=MCFNPC2001$SS,"PCODE"=MCFNPC2001$PCODE), FUN=sum) 
  survey2011 <- aggregate(MCFNPC2011$ABUND, by=list("PKEY"=MCFNPC2011$PKEY,"SS"=MCFNPC2011$SS,"PCODE"=MCFNPC2011$PCODE), FUN=sum) 
  
  speclist<-c("CAWA","OSFL","RUBL","CONI")
  
  list3[[j]]<-list(MCFNPC2001,MCFNPC2011,survey2001,survey2011,dat2001,dat2011,speclist)
  names(list3[[j]])<-c("MCFNPC2001","MCFNPC2011","survey2001","survey2011","dat2001","dat2011","speclist")
}

names(list3)<-names(list2)

beep(sound = 4, expr = NULL)



rm(list=setdiff(ls(),c("pred_MCFN","pred_MCFN_mask","LCC","OFF_ALL","list3")))
gc()
save.image("D:/MCFN/BRTdata_pack.RData")
beep(sound = 4, expr = NULL) # plays a sound to indicate command above finished running
#load("D:/MCFN/BRTdata_pack.RData")
