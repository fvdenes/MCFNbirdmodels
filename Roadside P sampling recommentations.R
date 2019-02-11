# Roadside point count sampling recommendationsc



library(RColorBrewer)
library(viridis)
library(sp)
library(rgeos)
library(rgdal)
library(raster)
library(latticeExtra)
library(scales)

studyRegion<-readOGR(dsn="c:/Users/voeroesd/Dropbox/BAM/MCFN",layer="Homelands_Boundary_10252016")
studyRegion<-spTransform(studyRegion,CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "))


pred_data <- read.csv("C:/Users/voeroesd/Dropbox/BAM/MCFN/output/BAM_pred_data_wBuff0.csv", sep="")
sp_pred_data<-SpatialPixelsDataFrame(points=pred_data[,c(2,3)],data=pred_data, proj4string = CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "), tolerance=0.9 )

roadshape<-readOGR(dsn="c:/Users/voeroesd/Dropbox/BAM/MCFN",layer="ORN_ROAD_NET_ELEMENT")
proj4string(roadshape)
roadshape<-spTransform(roadshape,CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "))


burnshape<-readOGR(dsn="c:/Users/voeroesd/Dropbox/BAM/MCFN",layer="FIRE_DISTURBANCE_AREA")
proj4string(burnshape)
burnshape<-spTransform(burnshape,CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "))


points<-readOGR(dsn="c:/Users/voeroesd/Dropbox/BAM/MCFN",layer="sampling_points")
points<-spTransform(points,CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "))


habraster<-raster(sp_pred_data["HAB_NALC1"])



plot(habraster,legend=FALSE)
plot(burnshape,add=T, col=alpha("red",0.6),border = 'transparent')
plot(roadshape,add=T)
plot(points,pch=16,cex=0.5,col="blue" , add=T)

