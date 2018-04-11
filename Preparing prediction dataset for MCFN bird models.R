### Preparing prediction dataset for MCFN bird models
## Ecoregion datasets
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-4.1.1.Rdata")
dat1<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-4.1.2.Rdata")
dat2<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-5.1.5.Rdata")
dat3<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-5.1.2.Rdata")
dat4<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-5.2.1.Rdata")
dat5<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-8.1.1.Rdata")
dat6<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-5.1.3.Rdata")
dat7<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-5.2.3.Rdata")
dat8<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-5.1.6.Rdata")
dat9<-dat
load("M:/DataStuff/cawa_ms_prediction_grid_data_by_ecoregion/pgdat-3.4.2.Rdata")
dat10<-dat


dat<-rbind(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10)

require(sp)
require(rgeos)
require(rgdal)
require(raster)

# Set wd where MCFN spatial data are
wd<-"c:/Users/voeroesd/Dropbox/BAM/MCFN"

studyRegion<-readOGR(dsn=wd,layer="Homelands_Boundary_10252016")
str(studyRegion)
studyRegion2<-spTransform(studyRegion,CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "))
studyRegion2


# Vector of incremental buffer, must follow unit of the proj4srting chosen. Here buffers are in meters
incBuff<-c("0","50000","100000","150000","200000","250000","300000","350000","400000","450000","500000")
srBuff <-gBuffer(studyRegion2, width=as.numeric(500000),quadsegs=35)

plot(studyRegion2)


bufferlims <- lapply(incBuff, function(x){gBuffer(studyRegion2,width=as.numeric(x),quadsegs=35)})



plot_points<-function(buffer){
  mmi<- mm[samp(mm)$SS%in% buffer$SS ,]
  sppoints<- SpatialPoints(coords=samp(mmi)[,5:6],proj4string = CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "))
  plot(sppoints, add=T,pch=1,col=2,cex=0.5)
}

can0<-getData('GADM', country="CAN", level=0)
can0<-spTransform(can0,CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "))
can1<-getData('GADM', country="CAN", level=1)
can1<-spTransform(can1,CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "))

plot(srBuff, axes=F)
plot(studyRegion2,add=T, col="lightgrey")
lapply(bufferlims,plot,add=T)
lapply(buffers,plot_points)
#coords_lists<-locator(11)
text(coords_lists$x,coords_lists$y,labels=unlist(lapply(buffers,function(x){length(x$SS)})),cex=0.7)
#coords_lists2<-locator(11)
text(coords_lists2$x,coords_lists2$y,labels=as.numeric(incBuff)/1000,cex=0.7,col=4)
box()

plot(can0, add=T)
plot(can1, add=T)
#### Intersect ecoregions with buffers (examples)
spdat1<-SpatialPointsDataFrame(coords=dat1[,2:3],data=dat1,proj4string = CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "))
results1 <-intersect(spdat1, srBuff)
plot(results1,add=T)

spdat2<-SpatialPointsDataFrame(coords=dat2[,2:3],data=dat2,proj4string = CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "))
results2 <-intersect(spdat2, srBuff)
plot(results2,add=T)

spdat3<-SpatialPointsDataFrame(coords=dat3[,2:3],data=dat3,proj4string = CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "))
results3 <-intersect(spdat3, srBuff)
plot(results3,add=T)

# Buffer / intersetc / save SS as list in csv file / save spatial R obj in shapefile
## All ecoregions, all buffers (this takes a VERY long time to run)
spdat<-SpatialPointsDataFrame(coords=dat[,2:3],data=dat,proj4string = CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "))

inters<-function(x) {
  srBuff <-gBuffer(studyRegion2, width=as.numeric(x),quadsegs=35)
  results <-intersect(spdat, srBuff)
  if(nrow(results)==0){
    
  }
  else{
    outf <- paste0("BAM_pred_data_wBuff",x)
    write.table(results, file.path(wd,"output",paste0(outf,".csv")), row.names=FALSE, col.names=TRUE)
    results
  }
}

l1<-lapply(incBuff,inters)
str(l1)


l1<-lapply(incBuff, function(x) {
  srBuff <-gBuffer(studyRegion2, width=as.numeric(x),quadsegs=35)
  results <-intersect(spdat, srBuff)
  outf <- paste0("BAM_pred_data_wBuff",x)
  write.table(results, file.path(wd,"output",paste0(outf,".csv")), row.names=FALSE, col.names=TRUE)
  results
  })
str(l1)


#### Maps of data points per buffers

