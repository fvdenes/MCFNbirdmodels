library(mefa4)
library(opticut)
library(MASS)
library(pbapply)
library(pROC)
library(sp)

# raw data for coni
load("G:/My Drive/BAM.SharedDrive/DataStuff/AvianData/Processed/BAM-BBS-tables_20170630.Rdata")
rm(OFF,PKEY,SS)


TAX <- nonDuplicated(TAX, Species_ID, TRUE)
spp<-c("CAWA","OSFL","RUBL","CONI")
TAX2 <- droplevels(TAX[spp,])
#load data (with offsets)
load("C:/Users/voeroesd/Dropbox/BAM/Critical Habitat/CHID subunit delineation/pack_2016-12-01.Rdata")


YY2 <- Xtab(ABUND ~ PKEY + SPECIES_ALL, PCTBL)


YY3<-YY2[,spp]

pk <- rownames(DAT)

YY3 <- YY3[pk,]
YY3 <- YY3[,colSums(YY3) >= 25]
apply(as.matrix(YY3), 2, max)
#apply(as.matrix(YY), 2, table)


# import buffers to crop data
BUFF0 <- read.csv("C:/Users/voeroesd/Dropbox/BAM/MCFN/BAM_pts_wBuff0.csv", sep="", stringsAsFactors=FALSE)
BUFF50 <- read.csv("C:/Users/voeroesd/Dropbox/BAM/MCFN/BAM_pts_wBuff50000.csv", sep="", stringsAsFactors=FALSE)
BUFF100 <- read.csv("C:/Users/voeroesd/Dropbox/BAM/MCFN/BAM_pts_wBuff100000.csv", sep="", stringsAsFactors=FALSE)
BUFF150 <- read.csv("C:/Users/voeroesd/Dropbox/BAM/MCFN/BAM_pts_wBuff150000.csv", sep="", stringsAsFactors=FALSE)
BUFF200 <- read.csv("C:/Users/voeroesd/Dropbox/BAM/MCFN/BAM_pts_wBuff200000.csv", sep="", stringsAsFactors=FALSE)
BUFF250 <- read.csv("C:/Users/voeroesd/Dropbox/BAM/MCFN/BAM_pts_wBuff250000.csv", sep="", stringsAsFactors=FALSE)
BUFF300 <- read.csv("C:/Users/voeroesd/Dropbox/BAM/MCFN/BAM_pts_wBuff300000.csv", sep="", stringsAsFactors=FALSE)
BUFF350 <- read.csv("C:/Users/voeroesd/Dropbox/BAM/MCFN/BAM_pts_wBuff350000.csv", sep="", stringsAsFactors=FALSE)
BUFF400 <- read.csv("C:/Users/voeroesd/Dropbox/BAM/MCFN/BAM_pts_wBuff400000.csv", sep="", stringsAsFactors=FALSE)
BUFF450 <- read.csv("C:/Users/voeroesd/Dropbox/BAM/MCFN/BAM_pts_wBuff450000.csv", sep="", stringsAsFactors=FALSE)
BUFF500 <- read.csv("C:/Users/voeroesd/Dropbox/BAM/MCFN/BAM_pts_wBuff500000.csv", sep="", stringsAsFactors=FALSE)

buffers<-list(BUFF0,BUFF50,BUFF100,BUFF150,BUFF200,BUFF250,BUFF300,BUFF350,BUFF400,BUFF450,BUFF500)

DAT$ROAD<-factor(DAT$ROAD)
DAT$HAB_NALC1 <- DAT$HABTR
DAT$HAB_NALC2 <- DAT$HAB
DAT$YEAR <- DAT$YR+2013

mm <- Mefa(YY3, DAT, TAX2, "inner")

mm <- mm[!is.na(samp(mm)$HAB_NALC1),]


# Model function

model_1 <- function(spp, buffer = NULL, road="no", maxit=25) {
  if(is.null(buffer)==T){
    mmi<-mm }  
  else{
    mmi<- mm[samp(mm)$SS%in% buffer$SS ,]
    }
  road_table <- table(mmi@samp$ROAD)
  
  set.seed(1)
  
  #helper functions:
  find_levels <- function(spp, m=1000) {
    j <- rep(FALSE, nrow(mm))
    for (k in levels(droplevels(samp(mm)$HAB_NALC1))) {
      w <- which(samp(mm)$HAB_NALC1 == k)
      if (length(w) < m) {
        j[w] <- TRUE
      } else {
        j[sample(w, m)] <- TRUE
      }
    }
    mm <- mm[j,]
    y <- as.numeric(xtab(mm)[,spp])
    x <- droplevels(samp(mm)$HAB_NALC1)
    if(spp%in%colnames(OFF)==T){
        ol <- optilevels(y=y, x=x, dist="poisson", offset=OFF[rownames(mm),spp])
    }
    else{
        ol <- optilevels(y=y, x=x, dist="poisson")
    }
    ol
  }
  'logLik.try-error' <- function (object, ...) {
    structure(-.Machine$double.xmax^(1/3), df = 1,
              nobs = 1, class = "logLik")
  }
  logLik.glm_skeleton <- function (object, ...) {
    structure(object$logLik, df = object$df,
              nobs = object$nobs, class = "logLik")
  }
  glm_skeleton <- function(object, ..., CAICalpha=0.5, keep_call=TRUE, vcov=FALSE) {
    if (inherits(object, "try-error"))
      return(structure(as.character(object), class="try-error"))
    out <- structure(list(
      call=object$call,
      formula=formula(object),
      coef=coef(object),
      vcov=NULL,
      converge=object$converge,
      logLik=as.numeric(logLik(object)),
      df=attr(logLik(object), "df"),
      nobs=nobs(object)), class="glm_skeleton")
    if (!out$converge)
      return(structure("glm did not converge", class="try-error"))
    if (!keep_call) {
      out$call <- out$formula <- NULL
    }
    if (vcov)
      out$vcov <- vcov(object)
    out$class0 <- class(object)[1L]
    out$aic <- -2*out$logLik + 2*out$df
    out$bic <- -2*out$logLik + log(out$nobs)*out$df
    out$caic <- CAICalpha * out$aic + (1-CAICalpha) * out$bic
    out
  }
  
  ol <- find_levels(spp, m=1000) # use subset of offroad data
  rc <- ol$levels[[length(ol$levels)]]
  
  if(road=="yes"){
    ff <- y ~ x + ROAD + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT # with road
  }
  
  if(road=="no"){
    ff <- y ~ x +  CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT # without road
  }
    
  y <- as.numeric(xtab(mmi)[,spp])
  x <- droplevels(samp(mmi)$HAB_NALC1)
  levels(x) <- rc[levels(x)]
  if(spp%in%colnames(OFF)==T){
    model<- glm_skeleton(try(glm(ff, data=samp(mmi),family="poisson", offset=OFF[rownames(mmi),spp],maxit=maxit)), keep_call=FALSE, vcov=TRUE)
    
    if(class(model)=="glm_skeleton"){
      fitted<-exp(model.matrix(ff,samp(mmi))%*%model$coef + OFF[rownames(mmi),spp])
    
      mvsamps<- exp(model.matrix(ff,samp(mmi)) %*% t(mvrnorm(n=1000,mu=model$coef,Sigma=model$vcov)) + OFF[rownames(mmi),spp]) 
 
      SEs<-apply(mvsamps,1,sd)
      fits<-data.frame(count=y,fitted=fitted,SEs=SEs)
      
      #fitted values for BUFF0
      fits0<-fits[dimnames(fitted)[[1]]%in%dimnames(samp(mm[samp(mm)$SS%in%BUFF0$SS,]))[[1]],]
      
      
      out<-list(spp=spp,
                levels=ol,
                model=model,
                fitted=fits,
                fittedBUFF0=fits0)
    }
    else{
      out<-list(spp=spp,
                levels=ol,
                model=model
                )
    }
  }
  else{
    model<- glm_skeleton(try(glm(ff, data=samp(mmi),family="poisson", maxit=maxit)), keep_call=FALSE, vcov=TRUE)
    
    if(class(model)=="glm_skeleton"){
      fitted<-exp(model.matrix(ff,samp(mmi))%*%model$coef)
      
      mvsamps<- exp(model.matrix(ff,samp(mmi)) %*% t(mvrnorm(n=1000,mu=model$coef,Sigma=model$vcov))) 

      SEs<-apply(mvsamps,1,sd)
      
      fits<-data.frame(count=y,fitted=fitted,SEs=SEs)
      
      #fitted values for BUFF0
      fits0<-fits[dimnames(fitted)[[1]]%in%dimnames(samp(mm[samp(mm)$SS%in%BUFF0$SS,]))[[1]],]
      
      
      out<-list(spp=spp,
                levels=ol,
                model=model,
                fitted=fits,
                fittedBUFF0=fits0)
    }
    else{
      out<-list(spp=spp,
                levels=ol,
                model=model)
    }
  }
  out
}


# Run model for all buffers


modelallbuffers <- function(spp,road="no",maxit=25){
  t0 <- proc.time()
  mbuf0<-model_1(spp = spp, buffer=BUFF0, road=road)
  mbuf50<-model_1(spp = spp, buffer=BUFF50, road=road)
  mbuf100<-model_1(spp = spp, buffer=BUFF100, road=road)
  mbuf150<-model_1(spp = spp, buffer=BUFF150, road=road)
  mbuf200<-model_1(spp = spp, buffer=BUFF200, road=road)
  mbuf250<-model_1(spp = spp, buffer=BUFF250, road=road)
  mbuf300<-model_1(spp = spp, buffer=BUFF300, road=road)
  mbuf350<-model_1(spp = spp, buffer=BUFF350, road=road)
  mbuf400<-model_1(spp = spp, buffer=BUFF400, road=road)
  mbuf450<-model_1(spp = spp, buffer=BUFF450, road=road)
  mbuf500<-model_1(spp = spp, buffer=BUFF500, road=road)
  
  out<- list(spp=spp,
             time=as.numeric(proc.time() - t0)[3L],
             buffer0=mbuf0,
             buffer50=mbuf50,
             buffer100=mbuf100,
             buffer150=mbuf150,
             buffer200=mbuf200,
             buffer250=mbuf250,
             buffer300=mbuf300,
             buffer350=mbuf350,
             buffer400=mbuf400,
             buffer450=mbuf450,
             buffer500=mbuf500)
  out
}


mods_CAWA <- modelallbuffers(spp="CAWA", road="no",maxit=100)


mods_OSFL <- modelallbuffers(spp="OSFL", road="no",maxit=100)


mods_RUBL <- modelallbuffers(spp="RUBL", road="no",maxit=100)


mods_CONI <- modelallbuffers(spp="CONI", road="no",maxit=100)



# Predict density

# Load prediction dataset: MCFN homelands
pred_data <- read.csv("C:/Users/voeroesd/Dropbox/BAM/MCFN/output/BAM_pred_data_wBuff0.csv", sep="")

str(pred_data)



get_preds <- function(modlist,preddata){
  
  # Reclass land cover variable in prediction dataset according to reclassing done for model
  rc_pred<-function(modlistbuffer,preddata){
    rc<-modlistbuffer$levels$levels[[length(modlistbuffer$levels$levels)]]
    x <- droplevels(preddata$HAB_NALC1)
    levels(x) <- rc[levels(x)]
    preddata$x<-x
    preddata$x<-relevel(preddata$x,rc[1])
    preddata
  }
  
  rc_preddata<- lapply(modlist[-c(1,2)],rc_pred,preddata)
  
  rc_preddata<- lapply(rc_preddata,data.frame,y=1) # add intercept (y=1) column
  
  ff <- y ~ x +  CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT # model formula without road
  
  mod_matrix_pred<-lapply(rc_preddata,model.matrix,object=ff) # generate model matrix
  
  #####
  getCoef<- function(modlistbuffer){
    coef<-modlistbuffer$model$coef
    coef
  }
  
  getvcov<-function(modlistbuffer){
    vcov<-modlistbuffer$model$vcov
    vcov
  }
  
  z<-which(lapply(modlist[-c(1,2)],function(x){class(x$model)})=="glm_skeleton")
  
  coefs <-lapply(modlist[z+2],getCoef)
  vcovs <-lapply(modlist[z+2],getvcov)
  
  
 
  
  SEs<-CVs<-preds<-as.data.frame(matrix(NA,nrow=nrow(rc_preddata$buffer0),ncol = length(coefs),dimnames = list(NULL,names(coefs))))
  

  for (i in 1:length(coefs)){
    if(ncol(mod_matrix_pred[[names(coefs[i])]])==length(coefs[[i]])){
      preds[,names(coefs)[i]]<- exp(mod_matrix_pred[[names(coefs[i])]]%*%coefs[[i]])
      
      mvsamps<- exp(mod_matrix_pred[[names(coefs[i])]] %*% t(mvrnorm(n=1000,mu=coefs[[i]],Sigma=vcovs[[i]])))
      
      is.na(mvsamps)<-is.infinite(mvsamps)
      
      CVs[,names(CVs)[i]]<-apply(mvsamps,1,sd,na.rm=T)/apply(mvsamps,1,mean,na.rm=T)
      SEs[,names(SEs)[i]]<-apply(mvsamps,1,sd,na.rm=T)
    }
    else{
      x<-mod_matrix_pred[[names(coefs[i])]]
      x<-x[,-which(colnames(mod_matrix_pred[[names(coefs[i])]])%in%attr(coefs[[i]], "names")==F)]
      preds[,names(coefs)[i]]<- exp(x%*%coefs[[i]])
      
      mvsamps<- exp(x %*% t(mvrnorm(n=1000,mu=coefs[[i]],Sigma=vcovs[[i]])))
      
      is.na(mvsamps)<-is.infinite(mvsamps)
      
      CVs[,names(CVs)[i]]<-apply(mvsamps,1,sd,na.rm=T)/apply(mvsamps,1,mean,na.rm=T)
      SEs[,names(SEs)[i]]<-apply(mvsamps,1,sd,na.rm=T)
    }
  }
  
  out<-list(preds=preds,
            SEs=SEs,
            CVs=CVs)
  out
}



library(RColorBrewer)
library(viridis)
library(sp)
library(rgeos)
library(rgdal)
library(raster)

preds_CAWA<-get_preds(mods_CAWA,pred_data)
points_pred_CAWA<-data.frame(pred_data,density=preds_CAWA$preds$buffer450,SEs=preds_CAWA$SEs$buffer500,CVs=preds_CAWA$CVs$buffer500)
density_pixels_CAWA<-SpatialPixelsDataFrame(points=points_pred_CAWA[,c(2,3)],data=points_pred_CAWA, proj4string = CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "), tolerance=0.9 )
plot(density_pixels_CAWA["density"], col=rev(viridis(250, end=0.96)), scale.shrink=2, zlim=c(0,12),axes=T)
map.axes()
#plot(density_pixels_CAWA["SEs"], col=rev(magma(250)))
plot(density_pixels_CAWA["CVs"], col=rev(inferno(250)), scale.shrink=2, zlim=c(0,0.5))

plot(density_pixels_CAWA["CVs"], col= colorRampPalette(brewer.pal(9,"PRGn"))(250)[76:250] , scale.shrink=2, zlim=c(0,3.5))

box()

layout(matrix(1:4, 2, 2), heights = c(8,1))
plot(density_pixels_CAWA["density"],  what="image",col=rev(viridis(250, end=0.96)), zlim=c(0,0.12))
text(x=1050000,y=1100000,"Density of singing males (per ha)")
plot(density_pixels_CAWA["density"], col=rev(viridis(250, end=0.96)), zlim=c(0,0.12), what = "scale", axis.pos = 1)
plot(density_pixels_CAWA["CVs"], col= colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(250)[76:250],zlim=c(0,3.5), what="image")
text(x=1050000,y=1100000," Coefficient of variation")
plot(density_pixels_CAWA["CVs"], col= colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(250)[76:250],zlim=c(0,3.5),  what = "scale", axis.pos = 1)



writeRaster(raster(density_pixels_CAWA["density"]),filename = "predsCAWA500km",prj=TRUE,  format = "GTiff",overwrite=T)
writeRaster(raster(density_pixels_CAWA["CVs"]),filename = "predsCAWA500km_CVs",prj=TRUE,  format = "GTiff",overwrite=T)
#summary(density_pixels_CAWA$density)
#plot(density_pixels_CAWA$density)

preds_OSFL<-get_preds(mods_OSFL,pred_data)
points_pred_OSFL<-data.frame(pred_data,density=preds_OSFL$preds$buffer200, SEs=preds_OSFL$SEs$buffer200, CVs=preds_OSFL$CVs$buffer200)
density_pixels_OSFL<-SpatialPixelsDataFrame(points=points_pred_OSFL[,c(2,3)],data=points_pred_OSFL, proj4string = CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "), tolerance=0.9 )
plot(density_pixels_OSFL["density"], col=rev(viridis(250, end=0.96)))
plot(density_pixels_OSFL["density"], col=rev(viridis(250, end=0.96)),zlim=c(0,10))
#plot(density_pixels_OSFL["SEs"], col=rev(magma(250)))
plot(density_pixels_OSFL["CVs"], col=rev(inferno(250)))

layout(matrix(1:4, 2, 2), heights = c(8,1))
plot(density_pixels_OSFL["density"],  what="image",col=rev(viridis(250, end=0.96)), zlim=c(0,0.07))
text(x=1050000,y=1100000,"Density of singing males (per ha)")
plot(density_pixels_OSFL["density"], col=rev(viridis(250, end=0.96)), zlim=c(0,0.07), what = "scale", axis.pos = 1)
plot(density_pixels_OSFL["CVs"], col= colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(250)[76:250],zlim=c(0,3.5),  what="image")
text(x=1050000,y=1100000," Coefficient of variation")
plot(density_pixels_OSFL["CVs"],col= colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(250)[76:250],zlim=c(0,3.5),  what = "scale", axis.pos = 1)



writeRaster(raster(density_pixels_OSFL["density"]),filename = "predsOSFL200km",prj=TRUE,  format = "GTiff",overwrite=T)
writeRaster(raster(density_pixels_OSFL["CVs"]),filename = "predsOSFL200km_CVs",prj=TRUE,  format = "GTiff",overwrite=T)
#summary(density_pixels_OSFL$density)
#plot(density_pixels_OSFL$density)

preds_RUBL<-get_preds(mods_RUBL,pred_data)
points_pred_RUBL<-data.frame(pred_data,density=preds_RUBL$preds$buffer300, SEs=preds_RUBL$SEs$buffer300, CVs=preds_RUBL$CVs$buffer300)
density_pixels_RUBL<-SpatialPixelsDataFrame(points=points_pred_RUBL[,c(2,3)],data=points_pred_RUBL, proj4string = CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "), tolerance=0.9 )
plot(density_pixels_RUBL["density"], col=rev(viridis(250, end=0.96)))
plot(density_pixels_RUBL["density"], col=rev(viridis(250, end=0.96)),zlim=c(0,10))
#plot(density_pixels_RUBL["SEs"], col=rev(magma(250)))
plot(density_pixels_RUBL["CVs"], col=rev(inferno(250)))

layout(matrix(1:4, 2, 2), heights = c(8,1))
plot(density_pixels_RUBL["density"],  what="image",col=rev(viridis(250, end=0.96)), zlim=c(0,0.02))
text(x=1050000,y=1100000,"Density of singing males (per ha)")
plot(density_pixels_RUBL["density"], col=rev(viridis(250, end=0.96)), zlim=c(0,0.02),  what = "scale", axis.pos = 1)
plot(density_pixels_RUBL["CVs"], col=colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(250)[76:250],zlim=c(0,3.5),  what="image")
text(x=1050000,y=1100000," Coefficient of variation")
plot(density_pixels_RUBL["CVs"], col=colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(250)[76:250],zlim=c(0,3.5),  what = "scale", axis.pos = 1)

writeRaster(raster(density_pixels_RUBL["density"]),filename = "predsRUBL300km",prj=TRUE,  format = "GTiff",overwrite=T)
writeRaster(raster(density_pixels_RUBL["CVs"]),filename = "predsRUBL300km_CVs",prj=TRUE,  format = "GTiff",overwrite=T)
#summary(density_pixels_RUBL$density)
#plot(density_pixels_RUBL$density)

preds_CONI<-get_preds(mods_CONI,pred_data)
points_pred_CONI<-data.frame(pred_data,density=preds_CONI$preds$buffer500, SEs=preds_CONI$SEs$buffer500, CVs=preds_CONI$CVs$buffer500)
density_pixels_CONI<-SpatialPixelsDataFrame(points=points_pred_CONI[,c(2,3)],data=points_pred_CONI, proj4string = CRS("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "), tolerance=0.9 )
plot(density_pixels_CONI["density"], col=rev(viridis(250, end=0.96)))
plot(density_pixels_CONI["density"], col=rev(magma(250, end=0.96)),zlim=c(0,10))
#plot(density_pixels_CONI["SEs"], col=rev(magma(250)))
plot(density_pixels_CONI["CVs"], col=rev(inferno(250)))


layout(matrix(1:4, 2, 2), heights = c(8,1))
plot(density_pixels_CONI["density"],  what="image",col=rev(viridis(250, end=0.96)), zlim=c(0,0.008))
text(x=1050000,y=1100000,"Relative density (individuals per ha)")
plot(density_pixels_CONI["density"], col=rev(viridis(250, end=0.96)), zlim=c(0,0.008), what = "scale", axis.pos = 1)
plot(density_pixels_CONI["CVs"], col=colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(250)[76:250],zlim=c(0,3.5),  what="image")
text(x=1050000,y=1100000," Coefficient of variation")
plot(density_pixels_CONI["CVs"],col=colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(250)[76:250], zlim=c(0,3.5),  what = "scale", axis.pos = 1)

writeRaster(raster(density_pixels_CONI["density"]),filename = "predsCONI500km",prj=TRUE,  format = "GTiff",overwrite=T)
writeRaster(raster(density_pixels_CONI["CVs"]),filename = "predsCONI500km_CVs",prj=TRUE,  format = "GTiff",overwrite=T)
#summary(density_pixels_CONI$density)
#plot(density_pixels_CONI$density)

# a summary figure with the 4 density maps:
layout(matrix(1:8, 4, 2), heights = c(7,1,7,1))
par(mar=c(1, 1, 0, 0),oma=c(0,0,0,0))
plot(density_pixels_CAWA["density"],  what="image",col=rev(viridis(250, end=0.96)), zlim=c(0,0.12))
#text(x=1050000,y=1117000,"Canada Warbler")
#text(x=1050000,y=1105000,"Density of singing males (per ha)")
plot(density_pixels_CAWA["density"], col=rev(viridis(250, end=0.96)), zlim=c(0,0.12), what = "scale", axis.pos = 1,scale.shrink=0.5)

plot(density_pixels_OSFL["density"],  what="image",col=rev(viridis(250, end=0.96)), zlim=c(0,0.07))
#text(x=1050000,y=1117000,"Olive-sided Flycatcher")
#text(x=1050000,y=1105000,"Density of singing males (per ha)")
plot(density_pixels_OSFL["density"], col=rev(viridis(250, end=0.96)), zlim=c(0,0.07), what = "scale", axis.pos = 1,scale.shrink=0.5)

plot(density_pixels_RUBL["density"],  what="image",col=rev(viridis(250, end=0.96)), zlim=c(0,0.02))
#text(x=1050000,y=1117000,"Rusty Blackbird")
#text(x=1050000,y=1105000,"Density of singing males (per ha)")
plot(density_pixels_RUBL["density"], col=rev(viridis(250, end=0.96)), zlim=c(0,0.02),  what = "scale", axis.pos = 1,scale.shrink=0.5)

plot(density_pixels_CONI["density"],  what="image",col=rev(viridis(250, end=0.96)), zlim=c(0,0.008))
#text(x=1050000,y=1117000,"Common Nighthawk")
#text(x=1050000,y=1105000,"Relative density (individuals per ha)")
plot(density_pixels_CONI["density"], col=rev(viridis(250, end=0.96)), zlim=c(0,0.008), what = "scale", axis.pos = 1,scale.shrink=0.5)


# what is the land cover of grid cells with highest CVs:
lc_highCV_CAWA<-points_pred_CAWA$HAB_NALC1[which(points_pred_CAWA$CVs>quantile(points_pred_CAWA$CVs, probs=0.90))]
sort(round(table(lc_highCV_CAWA)/length(lc_highCV_CAWA),2),decreasing = T)

lc_highCV_OSFL<-points_pred_OSFL$HAB_NALC1[which(points_pred_OSFL$CVs>quantile(points_pred_OSFL$CVs, probs=0.90))]
sort(round(table(lc_highCV_OSFL)/length(lc_highCV_OSFL),2),decreasing = T)

lc_highCV_RUBL<-points_pred_RUBL$HAB_NALC1[which(points_pred_RUBL$CVs>quantile(points_pred_RUBL$CVs, probs=0.90))]
sort(round(table(lc_highCV_RUBL)/length(lc_highCV_RUBL),2),decreasing = T)

lc_highCV_CONI<-points_pred_CONI$HAB_NALC1[which(points_pred_CONI$CVs>quantile(points_pred_CONI$CVs, probs=0.90))]
sort(round(table(lc_highCV_CONI)/length(lc_highCV_CONI),2),decreasing = T)


# Predictions by land cover type (for plotting):
library(ggplot2)

## Function to summarize dat for ggplots.
## Gives count, mean, standard deviation, standard error of the mean, and confidence 
## interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
  library(doBy)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # Collapse the data
  formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
  datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  
  # Rename columns
  names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
  names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
  names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



plot_lc_dens<-function(preddf,modlistbuffer,response.density=TRUE){
  library(ggplot2)
  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
    library(doBy)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }
    
    # Collapse the data
    formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
    datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
    
    # Rename columns
    names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
    names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
    names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
  }
  # reclass landcover according to model reclassing
  rc_pred<-function(modlistbuffer,preddata){
    rc<-modlistbuffer$levels$levels[[length(modlistbuffer$levels$levels)]]
    x <- droplevels(preddata$HAB_NALC1)
    levels(x) <- rc[levels(x)]
    preddata$x<-x
    preddata$x<-relevel(preddata$x,rc[1])
    preddata
  }
  sumgg <- summarySE(preddf,measurevar="density",groupvars = "HAB_NALC1")
  sumgg <-rc_pred(modlistbuffer,sumgg)
  sumgg$HAB_NALC1<-factor(sumgg$HAB_NALC1, levels = sumgg$HAB_NALC1[order(sumgg$x)])
  sumgg$z<-as.factor(as.numeric(sumgg$x))
  
  if(response.density==TRUE){
    ggplot(sumgg,aes(x=HAB_NALC1,y=density, fill=z))+
      geom_bar(position=position_dodge(), stat="identity") +
      geom_errorbar(aes(ymin=density-sd, ymax=density+sd),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15,face="bold", colour = cividis(n = length(levels(sumgg$z)), end=0.9)[sort(sumgg$z)]),panel.background = element_rect(fill = "white", colour = "grey50"),axis.text.y=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"))+
      xlab("Land cover class") + ylab("Density of singing males (per ha)")+
      labs(fill = "Class grouping") +
      scale_fill_viridis(discrete=TRUE,  end=0.9, guide=F,option = "cividis") 
  }
  else{
    ggplot(sumgg,aes(x=HAB_NALC1,y=density, fill=z))+
      geom_bar(position=position_dodge(), stat="identity") +
      geom_errorbar(aes(ymin=density-sd, ymax=density+sd),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15,face="bold", colour = cividis(n = length(levels(sumgg$z)), end=0.9)[sort(sumgg$z)]),panel.background = element_rect(fill = "white", colour = "grey50"),axis.text.y=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"))+
      xlab("Land cover class") + ylab(" Relative density (per ha)")+
      labs(fill = "Class grouping")+
      scale_fill_viridis(discrete=TRUE,  end=0.9, guide=F,option = "cividis") 
    
  }
}


plot_lc_dens(points_pred_CAWA,mods_CAWA$buffer500)
plot_lc_dens(points_pred_OSFL,mods_OSFL$buffer200)
plot_lc_dens(points_pred_RUBL,mods_RUBL$buffer300)
plot_lc_dens(points_pred_CONI,mods_CONI$buffer500, response.density = F)

# Estimates and CIs for each buffer
dev.off()
get_CIs<-function(modlist){
  ests_CIs<- function(model){
    mvsamps<- mvrnorm(n=1000,mu=model$coef,Sigma=model$vcov) # CIs estimated by drawing from a multivariate distribution with coefs and vcov, and quantiles
    CIs<-t(apply(mvsamps,2,quantile,probs=c(0,0.9)))
    ests<-cbind(coef=model$coef,CIs)
    ests
  }
  z<-which(lapply(modlist[-c(1,2)],function(x){class(x$model)})=="glm_skeleton")
  
  lapply(modlist[z+2],function(x){ests_CIs(x$model)})
}


get_CIs(mods_CAWA)


## ROC/AUC

get_ROC_AUC<-function(modlist){
  simple_roc <- function(labels, scores){
    labels <- labels[order(scores, decreasing=TRUE)]
    data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
  }
  simple_auc <- function(ROC) {
    ROC$inv_spec <- 1-ROC$FPR
    dx <- diff(ROC$inv_spec)
    sum(dx * ROC$TPR[-1]) / sum(dx)
  }
  get_vectors<-function(modbuffer){
    v1 <- data.frame(labels=modbuffer$fitted$count,scores=modbuffer$fitted$fitted)
    v1$labels[which(v1$labels>0)]<-1
    v1
    v0 <- data.frame(labels=modbuffer$fittedBUFF0$count,scores=modbuffer$fittedBUFF0$fitted)
    v0$labels[which(v0$labels>0)]<-1
    v0
    out<-list(v1=v1,
             v0=v0)
    out
  }
 z<-which(lapply(modlist[-c(1,2)],function(x){class(x$model)})=="glm_skeleton")  
 vlist<-lapply(modlist[z+2],get_vectors)
 roc_list<-lapply(vlist,function(x){simple_roc(x$v1$labels,x$v1$scores)})
 auc<-lapply(roc_list,simple_auc)
 roc_list0<-lapply(vlist,function(x){simple_roc(x$v0$labels,x$v0$scores)})
 auc0<-lapply(roc_list0,simple_auc)
 out<-list(
   roc=roc_list,
   auc=auc,
   rocbuff0=roc_list0,
   aucbuff0=auc0
   )
}

roc_CAWA<-get_ROC_AUC(mods_CAWA)
roc_OSFL<-get_ROC_AUC(mods_OSFL)
roc_RUBL<-get_ROC_AUC(mods_RUBL)
roc_CONI<-get_ROC_AUC(mods_CONI)

plot_AUC<-function(roc_spp){
  plot(y=unlist(roc_spp$auc),x=1:length(unlist(roc_spp$auc))-0.1,xaxt="n",ylab="AUC",xlab="buffer (km)",ylim=c(0.5,1))
  axis(1,at=1:length(unlist(roc_spp$auc)),labels=attr(unlist(roc_spp$auc),"names"))
  points(y=unlist(roc_spp$aucbuff0),x=1:length(unlist(roc_spp$auc))+0.1,xaxt="n",ylab="AUC",pch=16)
  legend(1,0.6,legend=c("Buffer area", "MCFN homelands"), pch=c(1,16))
}



plot_AUC(roc_CAWA)

plot_AUC(roc_OSFL)
plot_AUC(roc_RUBL)
plot_AUC(roc_CONI)


plot_SE<-function(modlist, log=TRUE){
  
  plot_AUC<-function(roc_spp){
    plot(y=unlist(roc_spp$auc),x=1:length(unlist(roc_spp$auc))-0.1,xaxt="n",ylab="AUC",xlab="",ylim=c(0.5,1))
    #axis(1,at=1:length(unlist(roc_spp$auc)),labels=attr(unlist(roc_spp$auc),"names"))
    points(y=unlist(roc_spp$aucbuff0),x=1:length(unlist(roc_spp$auc))+0.1,xaxt="n",ylab="AUC",pch=16)
    legend(1,0.75,legend=c("Buffer area", "MCFN homelands"), pch=c(1,16))
  }
  
  inf.replace<-function(x){
    x$fitted$SEs[which(is.infinite(x$fitted$SEs))]<-NA
    x$fittedBUFF0$SEs[which(is.infinite(x$fittedBUFF0$SEs))]<-NA
    x
  }
  
  if(log==TRUE){
    z<-which(lapply(modlist[-c(1,2)],function(x){class(x$model)})=="glm_skeleton")  
    

    modlist[z+2]<-lapply(modlist[z+2], inf.replace)
    
    SElist2<-SElist<-lapply(modlist[z+2],function(x){mean(x$fitted$SEs,na.rm=T)})
    SElist0<-lapply(modlist[z+2],function(x){mean(x$fittedBUFF0$SEs,na.rm=T)})
    
    SE_median<-lapply(modlist[z+2],function(x){median(x$fitted$SEs,na.rm=T)})
    SElist0_median<-lapply(modlist[z+2],function(x){median(x$fittedBUFF0$SEs,na.rm=T)})
    
    is.na(SElist2)<-is.infinite(unlist(SElist))
    
    
    par(mfrow=c(3,1),mar=c(1.1,4.1,2.1,2.1))
    
    roc<-get_ROC_AUC(modlist)
    plot_AUC(roc)
    
    par(mar=c(2.1,4.1,1.1,2.1))
    
    plot(y=log(unlist(SElist)),x=1:length(SElist)-0.1,xaxt="n",ylab="log (mean SE of fitted values)",xlab="",ylim=c(min(log(c(unlist(SElist2),unlist(SElist0))),na.rm=T),max(log(c(unlist(SElist2),unlist(SElist0))),na.rm=T)))
    #axis(1,at=1:length(SElist),labels=attr(unlist(SElist),"names"))
    points(y=log(unlist(SElist0)),x=1:length(SElist0)+0.1,pch=16)
    #legend(2,60,legend=c("Buffer area", "MCFN homelands"), pch=c(1,16)) 
    
    par(mar=c(3.1,4.1,0.1,2.1))
    
    plot(y=log(unlist(SE_median)),x=1:length(SE_median)-0.1,xaxt="n",ylab="log (median SE of fitted values)",xlab="",ylim=c(min(log(c(unlist(SE_median),unlist(SElist0_median))),na.rm=T),max(log(c(unlist(SE_median),unlist(SElist0_median))),na.rm=T)))
    axis(1,at=1:length(SElist),labels=attr(unlist(SElist),"names"))
    points(y=log(unlist(SElist0_median)),x=1:length(SElist0)+0.1,pch=16)
    #legend(1,-4,legend=c("Buffer area", "MCFN homelands"), pch=c(1,16)) 
    par(mfrow=c(1,1))
  }
  if(log==FALSE){
    z<-which(lapply(modlist[-c(1,2)],function(x){class(x$model)})=="glm_skeleton")  
    SElist<-lapply(modlist[z+2],function(x){mean(x$fitted$SEs,na.rm=T)})
    SElist0<-lapply(modlist[z+2],function(x){mean(x$fittedBUFF0$SEs,na.rm=T)})
    
    SE_median<-lapply(modlist[z+2],function(x){median(x$fitted$SEs,na.rm=T)})
    SElist0_median<-lapply(modlist[z+2],function(x){median(x$fittedBUFF0$SEs,na.rm=T)})
    
    is.na(SElist2)<-is.infinite(unlist(SElist))
    
    par(mfrow=c(3,1),mar=c(1.1,4.1,2.1,2.1))
    
    roc<-get_ROC_AUC(modlist)
    plot_AUC(roc)
    
    par(mar=c(2.1,4.1,2.1,2.1))
    
    plot(y=unlist(SElist),x=1:length(SElist)-0.1,xaxt="n",ylab="mean SE of fitted values",xlab="",ylim=c(min(c(unlist(SElist2),unlist(SElist0)),na.rm=T),max(c(unlist(SElist2),unlist(SElist0)),na.rm=T)))
    #axis(1,at=1:length(SElist),labels=attr(unlist(SElist),"names"))
    points(y=unlist(SElist0),x=1:length(SElist0)+0.1,pch=16)
    #legend(2,60,legend=c("Buffer area", "MCFN homelands"), pch=c(1,16)) 
    
    par(mar=c(3.1,4.1,2.1,2.1))
    
    plot(y=unlist(SE_median),x=1:length(SE_median)-0.1,xaxt="n",ylab="median SE of fitted values",xlab="",ylim=c(min(c(unlist(SE_median),unlist(SElist0_median)),na.rm=T),max(c(unlist(SE_median),unlist(SElist0_median)),na.rm=T)))
    axis(1,at=1:length(SE_median),labels=attr(unlist(SE_median),"names"))
    points(y=unlist(SElist0_median),x=1:length(SElist0)+0.1,pch=16)
    #legend(1,-4,legend=c("Buffer area", "MCFN homelands"), pch=c(1,16)) 
    par(mfrow=c(1,1))
  }
}

plot_SE(mods_CAWA)
#plot_SE(mods_CAWA,F)
plot_SE(mods_OSFL)
plot_SE(mods_RUBL)
plot_SE(mods_CONI)

