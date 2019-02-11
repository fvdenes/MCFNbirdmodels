library(sp)
library(maptools)
library(rgdal)
library(raster)
# North French River representation analysis

# load rasters with predicted bird densities from regional and national scale models
listfiles <- list.files("D:/MCFN/northfrench/DianaNationalModels",pattern="asc$")
setwd("D:/MCFN/northfrench/DianaNationalModels")
national <- stack(raster(listfiles[1]))
for (i in 2:length(listfiles)) { national <- addLayer(national, raster(listfiles[i]))}
plot(national)

listfiles2 <- list.files("D:/MCFN/",pattern="tif$")
setwd("D:/MCFN")
regional <- stack(raster(listfiles2[1]))
for (i in 2:length(listfiles2)) { regional <- addLayer(regional, raster(listfiles2[i]))}
regional<-projectRaster(regional,crs=proj4string(national))
plot(regional)

regional_density<-regional[[c(1,3,5,7)]]

#load MCFN tt
MCFNtt<- readOGR(dsn="D:/MCFN",layer="Homelands_Boundary_10252016")

# load watershed shapefile
watershed<-readOGR(dsn="D:/MCFN/northfrench/NorthFrenchWatershed",layer="NorthFrenchWatershed")
plot(watershed)
area(watershed)
watershed<-spTransform(watershed,CRS(proj4string(national)))

# load MCFN traditional territory shapefile
MCFNtt<-readOGR(dsn="D:/MCFN",layer="Homelands_Boundary_10252016")
plot(MCFNtt)
MCFNtt<-spTransform(MCFNtt,CRS(proj4string(national)))

# load Ontario (Brandt) boreal shapefile:
ONboreal<- readOGR(dsn="D:/MCFN/northfrench",layer="OntarioBoreal")
proj4string(ONboreal)
ONboreal<- spTransform(ONboreal,CRS=proj4string(national))

#mask national raster to MCFNtt and ON boreal extents
nat_mask<-rasterize(MCFNtt,national,mask=TRUE)
plot(nat_mask)
nat_mask<-crop(nat_mask,regional)
plot(nat_mask)
nat_mask_density<-nat_mask[[c(2,4,6,8,10,12)]]

ON_nat_mask<-rasterize(ONboreal,national,mask=T)
plot(ON_nat_mask)
ON_nat_mask<-crop(ON_nat_mask,ONboreal)
ON_nat_mask_density<-ON_nat_mask[[c(2,4,6,8,10,12)]]

#mask national and regional rasters to watershed extent
nat_mask_w<-rasterize(watershed,national,mask=TRUE)
nat_mask_w<-crop(nat_mask_w,watershed)

nat_mask_w_density<-nat_mask_w[[c(2,4,6,8,10,12)]]

plot(nat_mask_w_density)

reg_mask_w<-rasterize(watershed,regional,mask=TRUE)
reg_mask_w<-crop(reg_mask_w,watershed)


reg_mask_w_density<-reg_mask_w[[c(1,3,5,7)]]
plot(reg_mask_w_density)

# Get populaton size in MCFNtt for each spp

#nat_mask is 4000x4000m resolution

# a function to obtain population size estimates for a density raster based on Poisson draws (with sd) and also cell sums
getPopsize<-function(pred.raster){
  nvalues<-length(getValues(pred.raster[[1]])[!is.na(getValues(pred.raster[[1]]))])
  out<-matrix(NA,nlayers(pred.raster),3)
  preds<-matrix(NA,nvalues,nlayers(pred.raster))
  areaconv<-res(pred.raster)[1]^2/100^2
  
  for(i in 1:nlayers(pred.raster)){
    preds[,i]<-getValues(pred.raster[[i]])[!is.na(getValues(pred.raster[[i]]))]
    preds_4km<-preds[,i]*areaconv
    
    popsize<-matrix(NA,nvalues,100)
    
    for(j in 1:nvalues){
      popsize[j,]<- rpois(100,preds_4km[j])
    }
    out[i,1:2]<-c(mean(colSums(popsize)),sd(colSums(popsize)))
    
    out[i,3]<-sum(getValues(pred.raster[[i]]*areaconv), na.rm= TRUE)
  }
  row.names(out)<- gsub("_currmean","",names(pred.raster))
  colnames(out)<- c("popsize mean", "popsize sd","cellsum")
  out
}

natpopsize<-getPopsize(nat_mask_density)
regpopsize<-getPopsize(regional_density)

w.natpopsize<-getPopsize(nat_mask_w_density)
w.regpopsize<-getPopsize(reg_mask_w_density)
round(w.regpopsize,2)

ON.natpopsize<-getPopsize(ON_nat_mask_density)

props.watershed.ONboreal.nat<-w.natpopsize[,1]/ON.natpopsize[,1]
props.watershed.MCFN.nat<-w.natpopsize[,1]/natpopsize[,1]
props.watershed.MCFN.reg<-w.regpopsize[,1]/regpopsize[,1]

round(rbind(props.watershed.ONboreal.nat,props.watershed.MCFN.nat),3)

round(props.watershed.MCFN.reg,3)



area.watershed<-area(watershed)/1000000
area.MCFNtt<- area(MCFNtt)/1000000
area.ONboreal<-area(ONboreal)[1]/1000000

round(area.watershed/area.MCFNtt,3)
round(area.watershed/area.ONboreal,3)
