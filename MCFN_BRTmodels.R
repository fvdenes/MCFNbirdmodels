library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)
library(blockCV)
library(gbm)
library(viridis)
library(beepr)
load("D:/MCFN/BRT outputs/out1.RData")

#load data ####
load("D:/MCFN/BRTdata_pack.RData")

# prepare list with all data ####
list4 <- list(NA)
for (k in 1:length(list3)){
  list5<-list(NA)
  for(i in 1:length(list3[[k]]$speclist)){
    specoff <- filter(OFF_ALL, SPECIES==as.character(list3[[k]]$speclist[i]))
    specoff <- distinct(specoff)
    
    specdat2001 <- filter(list3[[k]]$MCFNPC2001, SPECIES == as.character(list3[[k]]$speclist[i]))
    x1 <- try(specdat2001x <- aggregate(specdat2001$ABUND,by=list("PKEY"=specdat2001$PKEY,"SS"=specdat2001$SS), FUN=sum))
    
    if(class(x1)!="try-error"){
      names(specdat2001x)[3] <- "ABUND"
      dat1 <- right_join(specdat2001x,list3[[k]]$survey2001[,1:3],by=c("SS","PKEY")) 
      dat1$SPECIES <- as.character(list3[[k]]$speclist[i])
      dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND)) 
      s2001 <- left_join(dat1,specoff, by=c("SPECIES","PKEY"))
      d2001 <- left_join(s2001, list3[[k]]$dat2001, by=c("SS"))
    }

    if(class(x1)=="try-error"){
      dat1<-cbind(list3[[k]]$survey2001[,1:2],data.frame(ABUND=rep(0,nrow(list3[[k]]$survey2001))),data.frame(PCODE=list3[[k]]$survey2001[,3]),data.frame(SPECIES=rep(as.character(list3[[k]]$speclist[i]),nrow(list3[[k]]$survey2001))))
      s2001 <- left_join(dat1,specoff, by=c("SPECIES","PKEY"))
      d2001 <- left_join(s2001, list3[[k]]$dat2001, by=c("SS"))
    }
    
    specdat2011 <- filter(list3[[k]]$MCFNPC2011, SPECIES == as.character(list3[[k]]$speclist[i]))
    x2 <- try(specdat2011x <- aggregate(specdat2011$ABUND,by=list("PKEY"=specdat2011$PKEY,"SS"=specdat2011$SS), FUN=sum))
    
    if(class(x2)!="try-error"){
      names(specdat2011x)[3] <- "ABUND"
      dat1 <- right_join(specdat2011x,list3[[k]]$survey2011[,1:3],by=c("SS","PKEY")) 
      dat1$SPECIES <- as.character(list3[[k]]$speclist[i])
      dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND)) 
      s2011 <- left_join(dat1,specoff, by=c("SPECIES","PKEY"))
      d2011 <- left_join(s2011, list3[[k]]$dat2011, by=c("SS"))
    }

    if(class(x2)=="try-error"){
      dat1<-cbind(list3[[k]]$survey2011[,1:2],data.frame(ABUND=rep(0,nrow(list3[[k]]$survey2011))),data.frame(PCODE=list3[[k]]$survey2011[,3]),data.frame(SPECIES=rep(as.character(list3[[k]]$speclist[i]),nrow(list3[[k]]$survey2011))))
      s2011 <- left_join(dat1,specoff, by=c("SPECIES","PKEY"))
      d2011 <- left_join(s2011, list3[[k]]$dat2011, by=c("SS"))
    }
    
    datcombo<-rbind(d2001,d2011)
    datcombo<-datcombo[which(!is.na(datcombo$wt)),]
    
    
    list5[[i]]<-datcombo
    names(list5)[[i]]<- paste0("datcombo.",as.character(list3[[k]]$speclist[i]))
    
  }
  list4[[k]]<-list5
  rm(list5)
  names(list4)[[k]]<-paste0("datcombo",names(list3)[[k]])
}
beep(sound = 4, expr = NULL)

rm(list=setdiff(ls(),c("LCC","list3","list4","pred_MCFN_mask")))
gc()

# check sites with >0 detections per spp and buffer
lapply(list4$datcomboMCFNSS0,function(x){nrow(subset(x,ABUND>0))})
lapply(list4$datcomboMCFNSS50000,function(x){nrow(subset(x,ABUND>0))})
lapply(list4$datcomboMCFNSS100000,function(x){nrow(subset(x,ABUND>0))})
lapply(list4$datcomboMCFNSS150000,function(x){nrow(subset(x,ABUND>0))})
lapply(list4$datcomboMCFNSS200000,function(x){nrow(subset(x,ABUND>0))})
lapply(list4$datcomboMCFNSS250000,function(x){nrow(subset(x,ABUND>0))})
lapply(list4$datcomboMCFNSS300000,function(x){nrow(subset(x,ABUND>0))})
lapply(list4$datcomboMCFNSS350000,function(x){nrow(subset(x,ABUND>0))})
lapply(list4$datcomboMCFNSS400000,function(x){nrow(subset(x,ABUND>0))})
lapply(list4$datcomboMCFNSS450000,function(x){nrow(subset(x,ABUND>0))})
lapply(list4$datcomboMCFNSS500000,function(x){nrow(subset(x,ABUND>0))})

# extend prediction stack for better plotting
pred_MCFN_mask.ext<-extend(pred_MCFN_mask,100)


# A function that fits the BRT model ('gbm.step' from dismo package) on pre-defined folds, and saves outputs ####
brt_MCFN <- function(data=list4,buffer,species, pred.stack, seed = 1222, pred.variables ,output.folder, blocks=NULL, keep.out = TRUE, tc=3,lr=0.001,bf=0.5, save.points.shp=FALSE){ 
  # Arguments for this function
  ## data: data.frame object containing data for model fitting
  ## pred.stack: the raster stack/brick used as prediction dataset
  ## pred.variables: a character vector giving the names of predictor variables that will be included in BRT models
  ## blocks: object resulting from 'spatialBlocks' function that contains classification of sample points into folds
  ## output.folder: path to output folder (if folder does not exist it is created)
  ## keep.out: logical, whether to keep the output in the workspace. both blocks object and the brt output are automatically saved in the output folder regardless
  ## tc: BRT tree complexity
  ## lr: BRT learning rate
  ## bf: BRT bag fraction
  ## save.points.shp: logical, whether to save survey points as a shapefile in output folder
  
  # fit BRT models using pre-determined folds for 
  buffer<-as.character(buffer)
  
  bufferlist<-c("0","50000","100000","150000","200000","250000","300000","350000","400000","450000","500000")
  speclist<-c("CAWA","OSFL","RUBL","CONI")
  
  dat<-data[[which(bufferlist==buffer)]]
  datcombo<-dat[[which(speclist==species)]]
  
  if(any(is.na(datcombo$logoffset))==T){
    offset<-NULL
  }
  else{
    offset<-datcombo$logoffset
  }
  
  if (is.null(blocks)){
    folds<-NULL
    n.folds<-10
  }
  else {
    folds<-blocks$foldID
    n.folds<-blocks$k
  }
  #fit models ####
  x1 <-
    try(brt <-
          gbm.step(
            datcombo,
            gbm.y = "ABUND",
            gbm.x = pred.variables,
            fold.vector = folds,
            n.folds = n.folds,
            family = "poisson",
            tree.complexity = tc,
            learning.rate = lr,
            bag.fraction = bf,
            offset = offset,
            site.weights = datcombo$wt,
            keep.fold.models = T,
            keep.fold.fit = T
          ))
  
  if(class(x1)=="NULL"){#retry models that didn't converge with smaller learning rate
    x1 <-
      try(brt <-
            gbm.step(
              datcombo,
              gbm.y = "ABUND",
              gbm.x = pred.variables,
              fold.vector = folds,
              n.folds = n.folds,
              family = "poisson",
              tree.complexity = tc,
              learning.rate = lr/10,
              bag.fraction = bf,
              offset = offset,
              site.weights = datcombo$wt,
              keep.fold.models = T,
              keep.fold.fit = T
            ))
  }
  
  if(class(x1)=="NULL"){#retry models that didn't converge with smaller learning rate
    x1 <-
      try(brt <-
            gbm.step(
              datcombo,
              gbm.y = "ABUND",
              gbm.x = pred.variables,
              fold.vector = folds,
              n.folds = n.folds,
              family = "poisson",
              tree.complexity = tc,
              learning.rate = lr/100,
              bag.fraction = bf,
              offset = offset,
              site.weights = datcombo$wt,
              keep.fold.models = T,
              keep.fold.fit = T
            ))
  }
  
  if(class(x1)=="NULL"){#retry models that didn't converge with smaller learning rate
    x1 <-
      try(brt <-
            gbm.step(
              datcombo,
              gbm.y = "ABUND",
              gbm.x = pred.variables,
              fold.vector = folds,
              n.folds = n.folds,
              family = "poisson",
              tree.complexity = tc,
              learning.rate = lr/1000,
              bag.fraction = bf,
              offset = offset,
              site.weights = datcombo$wt,
              keep.fold.models = T,
              keep.fold.fit = T
            ))
  }
  
  if(class(x1)=="NULL"){#retry models that didn't converge with smaller learning rate
    x1 <-
      try(brt <-
            gbm.step(
              datcombo,
              gbm.y = "ABUND",
              gbm.x = pred.variables,
              fold.vector = folds,
              n.folds = n.folds,
              family = "poisson",
              tree.complexity = tc,
              learning.rate = lr/10000,
              bag.fraction = bf,
              offset = offset,
              site.weights = datcombo$wt,
              keep.fold.models = T,
              keep.fold.fit = T
            ))
  }
  
  x1 <-
    try(brt <-
          gbm.step(
            datcombo,
            gbm.y = "ABUND",
            gbm.x = pred.variables,
            fold.vector = folds,
            n.folds = n.folds,
            family = "poisson",
            tree.complexity = tc-1,
            learning.rate = lr,
            bag.fraction = bf,
            offset = offset,
            site.weights = datcombo$wt,
            keep.fold.models = T,
            keep.fold.fit = T
          ))
  
  if(class(x1)=="NULL"){#retry models that didn't converge with smaller learning rate
    x1 <-
      try(brt <-
            gbm.step(
              datcombo,
              gbm.y = "ABUND",
              gbm.x = pred.variables,
              fold.vector = folds,
              n.folds = n.folds,
              family = "poisson",
              tree.complexity = tc-1,
              learning.rate = lr/10,
              bag.fraction = bf,
              offset = offset,
              site.weights = datcombo$wt,
              keep.fold.models = T,
              keep.fold.fit = T
            ))
  }
  
  if(class(x1)=="NULL"){#retry models that didn't converge with smaller learning rate
    x1 <-
      try(brt <-
            gbm.step(
              datcombo,
              gbm.y = "ABUND",
              gbm.x = pred.variables,
              fold.vector = folds,
              n.folds = n.folds,
              family = "poisson",
              tree.complexity = tc-1,
              learning.rate = lr/100,
              bag.fraction = bf,
              offset = offset,
              site.weights = datcombo$wt,
              keep.fold.models = T,
              keep.fold.fit = T
            ))
  }
  
  if(class(x1)=="NULL"){#retry models that didn't converge with smaller learning rate
    x1 <-
      try(brt <-
            gbm.step(
              datcombo,
              gbm.y = "ABUND",
              gbm.x = pred.variables,
              fold.vector = folds,
              n.folds = n.folds,
              family = "poisson",
              tree.complexity = tc-1,
              learning.rate = lr/1000,
              bag.fraction = bf,
              offset = offset,
              site.weights = datcombo$wt,
              keep.fold.models = T,
              keep.fold.fit = T
            ))
  }
  
  if(class(x1)=="NULL"){#retry models that didn't converge with smaller learning rate
    x1 <-
      try(brt <-
            gbm.step(
              datcombo,
              gbm.y = "ABUND",
              gbm.x = pred.variables,
              fold.vector = folds,
              n.folds = n.folds,
              family = "poisson",
              tree.complexity = tc-1,
              learning.rate = lr/10000,
              bag.fraction = bf,
              offset = offset,
              site.weights = datcombo$wt,
              keep.fold.models = T,
              keep.fold.fit = T
            ))
  }
  
  if(class(x1)=="NULL"){
    stop("Restart model with even smaller learning rate, or other predictors!")
  }
  
  # Define/create folders for storing outputs ####
  if (class(x1) != "try-error") {
    z <- output.folder
    
    if (file.exists(z) == FALSE) {
      dir.create(z, recursive = T)
    }
    
    if (is.null(blocks)){
      save(brt, file = paste(z,buffer, species, "brtMCFN.R", sep = ""))
    }
    
    else {
      save(blocks, file = paste(z,buffer, species, "blocks.R", sep = ""))
      save(brt, file = paste(z,buffer, species, "brtAB.R", sep = ""))
    }  
    
    
    ## Model evaluation
    varimp <- as.data.frame(brt$contributions)
    write.csv(varimp, file = paste(z,buffer, species, "varimp.csv", sep = ""))
    cvstats <- t(as.data.frame(brt$cv.statistics))
    write.csv(cvstats, file = paste(z,buffer, species, "cvstats.csv", sep = ""))
    pdf(paste(z,buffer, species, "_plot.pdf", sep = ""))
    gbm.plot(
      brt,
      n.plots = length(pred.variables),
      smooth = TRUE,
      plot.layout = c(3, 3),
      common.scale = T
    )
    dev.off()
    pdf(paste(z,buffer, species, "_plot.var-scale.pdf", sep = ""))
    gbm.plot(
      brt,
      n.plots = length(pred.variables),
      smooth = TRUE,
      plot.layout = c(3, 3),
      common.scale = F,
      write.title = F
    )
    dev.off()
    
    
    ## Model prediction
    
    rast <-
      predict(pred.stack,
              brt,
              type = "response",
              n.trees = brt$n.trees)
    writeRaster(
      rast,
      filename = paste(z,buffer, species, "_pred1km", sep = ""),
      format = "GTiff",
      overwrite = TRUE
    )
    
    
    data_sp <-SpatialPointsDataFrame(coords = datcombo[which(datcombo$ABUND>0),8:9], data = datcombo[which(datcombo$ABUND>0),], proj4string = LCC)
    
    prev <- cellStats(rast, 'mean')    
    max <- 3*prev
    
    png(paste(z,buffer, species, "_preds.png", sep = ""), height=600, width=850)
    par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
    par(mar=c(0,0,5,0))
    plot(rast, col=viridis(15)[15], axes=TRUE, legend=FALSE, main=paste0(species," ",buffer,"m buffer", " current prediction"))
    
    plot(rast, col=viridis(15), zlim=c(0,max), axes=T, main="CAWA", add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.65,0.90,0.83,0.87), axis.args=list(cex.axis=1.2))
    text(350000,6850000,"Potential density (males/ha)", cex=1.3)
    
    points(data_sp$X, data_sp$Y, cex = 1, col="red", pch=3)
    dev.off()
    
    if(save.points.shp==T){
      writeOGR(
        data_sp,
        dsn = paste(z,buffer, "surveypoints.shp", sep = ""),
        layer = "data_sp",
        driver = "ESRI Shapefile"
      )
    }
    
    if(keep.out==T) {return(brt)}
  }
}

# assess SD among model predictors ####

SDtest<-apply(list4$datcomboMCFNSS0$datcombo.CONI[,16:62],2,sd, na.rm=T)
nonzeroSD<-names(list4$datcomboMCFNSS0$datcombo.CAWA[,16:62])[which(SDtest>0)]

# define prediction variables ####
pred.variables<-c(nonzeroSD,
                  "Structure_Biomass_TotalLiveAboveGround_v1",
                  "Structure_Stand_Age_v1",
                  "CTI_MCFN",
                  "water",
                  "AHM",
                  "DD18",
                  "MAT",
                  "MSP",
                  "FFP",
                  "TD"
)


# BRT models ####

list5<-list(NA)
bufferlist<-c("0","50000","100000","150000","200000","250000","300000","350000","400000","450000","500000")
speclist<-c("CAWA","OSFL","RUBL","CONI")
for(s in 1:length(speclist)){
  list6<-list(NA)
  for(b in 1:length(bufferlist)){
    out.folder<-paste0("D:/MCFN/BRT outputs/",speclist[s],"/buffer_",bufferlist[b],"m/")
    x1<- try(brt<-brt_MCFN(buffer=bufferlist[b],species = speclist[s], pred.variables=pred.variables, pred.stack = pred_MCFN_mask.ext, output.folder = out.folder,save.points.shp = TRUE,tc=3, lr=0.01 ))
    if (class(x1)=="NULL"){
      list6[[b]]<-NA
    }
    if (class(x1)=="try-error"){
      list6[[b]]<-NA
    }
    if(class(x1)!="try-error"&&class(x1)!="NULL"){
      list6[[b]]<-brt
    }
    names(list6)[[b]]<-paste0(speclist[s],"_buffer_",bufferlist[b],"m")
  }
  list5[[s]]<-list6
  rm(list6)
  names(list5)[[s]]<-speclist[s]
}


rm(list=setdiff(ls(),c("LCC","list3","list4","list5","pred_MCFN_mask","pred_MCFN_mask.ext","brt_MCFN")))
gc()
save.image("D:/MCFN/BRT outputs/out1.RData")
beep(sound = 4, expr = NULL) 

# Compare model CV statistics to find best model (less bias in predictions) ####
compCVstats<- function(modlist, species){
  library(ggplot2)
  lsmod<-modlist[[species]]
  
  if(any(is.na(lsmod))){
    lsmod<-lsmod[-which(is.na(lsmod))]
  }
  
  # # Deviance 
  # df1<-as.data.frame(matrix(nrow=length(lsmod),ncol=3))
  # colnames(df1)<-c("buffer","deviance","se")
  # df1[,1]<-gsub("m","",gsub(paste0(species,"_buffer_"),"",names(lsmod)))
  # df1[,1]<-factor(df1[,1],levels=df1[,1])
  # df1[,2]<-unlist(lapply(lsmod,function(x){x$cv.statistics$deviance.mean}))
  # df1[,3]<-unlist(lapply(lsmod,function(x){x$cv.statistics$deviance.se}))
  # k1<-ggplot(df1, aes(buffer,deviance,ymin=deviance-se,ymax=deviance+se))
  # k1+geom_pointrange()
  # ggsave(paste0("D:/MCFN/BRT outputs/", species, "/",species,"_CV_deviance.pdf"))
  
  # Deviance explained
  df5<-as.data.frame(matrix(nrow=length(lsmod),ncol=3))
  colnames(df5)<-c("buffer","deviance.explained","se")
  df5[,1]<-gsub("m","",gsub(paste0(species,"_buffer_"),"",names(lsmod)))
  df5[,1]<-factor(df5[,1],levels=df5[,1])
  df5[,2]<-unlist(lapply(lsmod,function(x){(x$self.statistics$null - x$cv.statistics$deviance.mean)/x$self.statistics$null}))
  df5[,3]<-unlist(lapply(lsmod,function(x){x$cv.statistics$deviance.se/x$self.statistics$null}))
  k5<-ggplot(df5, aes(buffer,deviance.explained,ymin=deviance.explained-se,ymax=deviance.explained+se))
  k5+geom_pointrange()+xlab("Buffer distance (km)")+scale_x_discrete(labels=as.numeric(levels(df5[,1]))/1000)+ylab("Deviance explained")
  ggsave(paste0("D:/MCFN/BRT outputs/", species, "/",species,"_CV_deviance_explained.pdf"))
  
  
  # Correlation
  df2<-as.data.frame(matrix(nrow=length(lsmod),ncol=3))
  colnames(df2)<-c("buffer","correlation","se")
  df2[,1]<-gsub("m","",gsub(paste0(species,"_buffer_"),"",names(lsmod)))
  df2[,1]<-factor(df2[,1],levels=df2[,1])
  df2[,2]<-unlist(lapply(lsmod,function(x){x$cv.statistics$correlation.mean}))
  df2[,3]<-unlist(lapply(lsmod,function(x){x$cv.statistics$correlation.se}))
  k2<-ggplot(df2, aes(buffer,correlation,ymin=correlation-se,ymax=correlation+se))
  k2+geom_pointrange()+xlab("Buffer distance (km)")+scale_x_discrete(labels=as.numeric(levels(df2[,1]))/1000)+ylab("Correlation")
  ggsave(paste0("D:/MCFN/BRT outputs/", species, "/",species,"_CV_correlation.pdf"))
  
  # Calibration:  The first two statistics were the estimated intercepts and slopes of linear regression models of predictions against observations. The intercept measures the magnitude and direction of bias, with values close to 0 indicating low or no bias. The slope yields information about the consistency in the bias as a function of the mean, with a value of 1 indicating a consistent bias if the intercept is a nonzero value.
  ## intercept
  df3<-as.data.frame(matrix(nrow=length(lsmod),ncol=3))
  colnames(df3)<-c("buffer","calibration.intercept","se")
  df3[,1]<-gsub("m","",gsub(paste0(species,"_buffer_"),"",names(lsmod)))
  df3[,1]<-factor(df3[,1],levels=df3[,1])
  df3[,2]<-unlist(lapply(lsmod,function(x){x$cv.statistics$calibration.mean[1]}))
  df3[,3]<-unlist(lapply(lsmod,function(x){x$cv.statistics$calibration.se[1]}))
  k3<-ggplot(df3, aes(buffer,calibration.intercept,ymin=calibration.intercept-se,ymax=calibration.intercept+se))+ylab("Calibration intercept")
  k3+geom_pointrange()+xlab("Buffer distance (km)")+scale_x_discrete(labels=as.numeric(levels(df3[,1]))/1000)
  ggsave(paste0("D:/MCFN/BRT outputs/", species, "/",species,"_CV_calibration.intercept.pdf"))
  
  ## slope
  df4<-as.data.frame(matrix(nrow=length(lsmod),ncol=3))
  colnames(df4)<-c("buffer","calibration.slope","se")
  df4[,1]<-gsub("m","",gsub(paste0(species,"_buffer_"),"",names(lsmod)))
  df4[,1]<-factor(df4[,1],levels=df4[,1])
  df4[,2]<-unlist(lapply(lsmod,function(x){x$cv.statistics$calibration.mean[2]}))
  df4[,3]<-unlist(lapply(lsmod,function(x){x$cv.statistics$calibration.se[2]}))
  k4<-ggplot(df4, aes(buffer,calibration.slope,ymin=calibration.slope-se,ymax=calibration.slope+se))
  k4+geom_pointrange()+xlab("Buffer distance (km)")+scale_x_discrete(labels=as.numeric(levels(df4[,1]))/1000)+ylab("Calibration slope")
  ggsave(paste0("D:/MCFN/BRT outputs/", species, "/",species,"_CV_calibration.slope.pdf"))
  
}

compCVstats(list5,"CAWA")
compCVstats(list5,"OSFL")
compCVstats(list5,"CONI")

# Buffer0 model has very high calibration values , cannot distinguish values of other buffers in plot. re-run without buff0
list6<-list5
list6$CONI$CONI_buffer_0m<-NULL
compCVstats(list6,"CONI")
rm(list6)
gc()

compCVstats(list5,"RUBL")


# Bootstrap and CI for best model ####
boot_brt<-function(data,brtmodel,pred.data,nsamples=1000,output.folder){
  
  z <- output.folder
  
  if (file.exists(z) == FALSE) {
    dir.create(z,recursive = TRUE)
  }
  
  for(i in 1:nsamples){
    
    
    cat("loop",i,"\n") # this prints the loop number on console to track function progress
    
    
    sample<-sample(1:nrow(data),size=nrow(data),replace=T)
    data2<-data[sample,]
    
    
    brt<-
      gbm.fit(x = data2[,brtmodel$gbm.call$gbm.x], 
              y = data2[,brtmodel$gbm.call$gbm.y],
              offset<-brtmodel$gbm.call$offset[sample],
              distribution ="poisson",
              w = data2$wt,
              n.trees = brtmodel$n.trees,
              interaction.depth = brtmodel$interaction.depth,
              shrinkage = brtmodel$gbm.call$learning.rate,
              bag.fraction = brtmodel$gbm.call$bag.fraction,
              var.names = brtmodel$gbm.call$gbm.x,
              response.name = brtmodel$gbm.call$gbm.y
      )
    
    rast <- predict(pred.data,
                    brt,
                    type = "response",
                    n.trees = brt$n.trees)
    
    
    
    #stack <- addLayer(stack, rast)
    #names(stack)[i+1]<-paste0("sample",i) 
    
    writeRaster(rast, filename=paste0(z, "sample",i,".tif"), format="GTiff",overwrite=TRUE)
    rm(rast)
    gc()
    
  }
  
  #fun0.05 <- function(x) {quantile(x, probs = 0.05, na.rm = TRUE)}
  #lower<- calc(stack[[-1]],fun0.05)
  #fun0.95 <- function(x) {quantile(x, probs = 0.95, na.rm = TRUE)}
  #upper<- calc(stack[[-1]],fun0.95)
  
  #writeRaster(lower, filename=paste0(z, " confint_lower.tif"), format="GTiff",overwrite=TRUE)
  #writeRaster(upper, filename=paste0(z, " confint_upper.tif"), format="GTiff",overwrite=TRUE)
  
  #return(stack)
}

#CAWA: buffer = 200000
start_time <- Sys.time()
boot_brt(data=list4$datcomboMCFNSS200000$datcombo.CAWA,brtmodel=list5$CAWA$CAWA_buffer_200000m,pred.data=pred_MCFN_mask,nsamples=500,output.folder = "D:/MCFN/BRT outputs/CAWA/buffer_200000m/confint/")
end_time <- Sys.time()
end_time - start_time

#OSFL: buffer = 100000
boot_brt(data=list4$datcomboMCFNSS100000$datcombo.OSFL,brtmodel=list5$OSFL$OSFL_buffer_100000m,pred.data=pred_MCFN_mask,nsamples=500,output.folder = "D:/MCFN/BRT outputs/OSFL/buffer_100000m/confint/")

#CONI
boot_brt(data=list4$datcomboMCFNSS450000$datcombo.CONI,brtmodel=list5$CONI$CONI_buffer_450000m,pred.data=pred_MCFN_mask,nsamples=500,output.folder = "D:/MCFN/BRT outputs/CONI/buffer_450000m/confint/")

#RUBL
boot_brt(data=list4$datcomboMCFNSS500000$datcombo.RUBL,brtmodel=list5$RUBL$RUBL_buffer_500000m,pred.data=pred_MCFN_mask,nsamples=500,output.folder = "D:/MCFN/BRT outputs/RUBL/buffer_500000m/confint/")

save.image("D:/MCFN/BRT outputs/out1.RData")


# Calculate CI based on quantiles for each model, calculate SE and map them ####

CI_SE_CV_plotfun<-function(species,buffer){
  
  fun0.05 <- function(x) {quantile(x, probs = 0.05, na.rm = TRUE)}
  fun0.95 <- function(x) {quantile(x, probs = 0.95, na.rm = TRUE)}
  
  path<-paste0("D:/MCFN/BRT outputs/",species,"/buffer_",buffer,"m/confint/")
  
  confintBRT<-brick(paste0(path,"sample1.tif"))
  files <- list.files(path,pattern="tif$")
  setwd(path)
  for (i in 2:length(files)) { confintBRT <- addLayer(confintBRT, raster(files[i]))}
  names(confintBRT)[1:500] <- paste0("sample",1:500) 
  
  lower<- calc(confintBRT,fun0.05)
  upper<- calc(confintBRT,fun0.95)
  SE<- calc(confintBRT,sd)
  CV<-calc(confintBRT,function(x){sd(x)/mean(x)})
  
  #plot(SE)
  writeRaster(SE, filename=paste0(path, "SE.tif"), format="GTiff",overwrite=TRUE)
  writeRaster(CV, filename=paste0(path, "CV.tif"), format="GTiff",overwrite=TRUE)
  writeRaster(lower, filename=paste0(path, "90CIlower.tif"), format="GTiff",overwrite=TRUE)
  writeRaster(upper, filename=paste0(path, "90CIupper.tif"), format="GTiff",overwrite=TRUE)
  
}


start_time2 <- Sys.time()
CI_SE_CV_plotfun("CAWA","200000")
end_time2 <- Sys.time()
end_time2- start_time2
beep(sound = 4, expr = NULL)

start_time2 <- Sys.time()
CI_SE_CV_plotfun("OSFL","100000")
end_time2 <- Sys.time()
end_time2- start_time2
beep(sound = 4, expr = NULL)

start_time2 <- Sys.time()
CI_SE_CV_plotfun("CONI","450000")
end_time2 <- Sys.time()
end_time2- start_time2
beep(sound = 4, expr = NULL)

start_time2 <- Sys.time()
CI_SE_CV_plotfun("RUBL","500000")
end_time2 <- Sys.time()
end_time2- start_time2
beep(sound = 4, expr = NULL)

save.image("D:/MCFN/BRT outputs/out1.RData")
beep(sound = 4, expr = NULL) 


# Plot predictions from best model for each spp ####
pred_MCFN_mask.ext<-extend(brick("D:/MCFN/rasters/Prediction dataset/pred_MCFN_mask"),100)

plot_preds<-function(species,buffer,modlist=list5,addpoints=T,SD=FALSE, CV=TRUE){
  lsmod<-modlist[[species]]
  bufferlist<-c("0","50000","100000","150000","200000","250000","300000","350000","400000","450000","500000")
  mod<-lsmod[[which(bufferlist==buffer)]]
  
  rast<-predict(pred_MCFN_mask.ext,
                mod,
                type = "response",
                n.trees = mod$n.trees)
  
  prev <- cellStats(rast, 'mean')    
  max <- 3*prev
  
  png(file=paste0("D:/MCFN/BRT outputs/",species,"/preds_",species,"_",buffer,".png"), height=600, width=850)
  par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
  par(mar=c(0,0,5,0))
  plot(rast, col=viridis(15)[15], axes=TRUE, legend=FALSE, main=paste0(species," ",as.numeric(buffer)/1000," Km buffer"))
  
  plot(rast, col=viridis(15), zlim=c(0,max), axes=T, main=species, add=TRUE, legend.width=1.5, horizontal = TRUE, smallplot = c(0.65,0.90,0.83,0.87), axis.args=list(cex.axis=1.2))
  text(1170000,6980000,"Potential density (males/ha)", cex=1.3)

  if(addpoints==TRUE){
    
    speclist<-c("CAWA","OSFL","RUBL","CONI")
    
    dat<-list4$datcomboMCFNSS0
    datcombo<-dat[[which(speclist==species)]]
    if(mean(datcombo$ABUND)==0){
      datcombo<-dat[[which(speclist=="CAWA")]]# any species here will do
      data_sp <-SpatialPointsDataFrame(coords = datcombo[,c("X","Y")], data = datcombo, proj4string = LCC)
      points(data_sp$X, data_sp$Y, cex = 1, col=1, pch=4)
      
    }
    else{
      data_sp <-SpatialPointsDataFrame(coords = datcombo[,c("X","Y")], data = datcombo, proj4string = LCC)
      points(data_sp$X[data_sp$ABUND==0], data_sp$Y[data_sp$ABUND==0], cex = 1, col=1, pch=4)
      points(data_sp$X[data_sp$ABUND>0], data_sp$Y[data_sp$ABUND>0], cex = 1, col=2, pch=4)
    }
    
  }
  dev.off()
  
  magma2 <- colorRampPalette(magma(15), space="Lab", bias=2)
  
  if (SD==TRUE){  
    SD<-raster(paste0("D:/MCFN/BRT outputs/",species,"/buffer_",buffer,"m/confint/SE.tif"))
    SD_ext<-extend(SD,100)
    
    png(file=paste0("D:/MCFN/BRT outputs/",species,"/SE_",species,"_",buffer,".png"), height=600, width=850)
    par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
    par(mar=c(0,0,5,0))
    plot(SD_ext, col=magma(15)[15], axes=TRUE, legend=FALSE, main=paste0(species," ",as.numeric(buffer)/1000," Km buffer"))
    plot(SD_ext, col=magma(15), axes=T,  main=paste0(species," ",buffer,"m buffer SD"), legend.width=1.5, horizontal = TRUE, smallplot = c(0.65,0.90,0.83,0.87), axis.args=list(cex.axis=1.2),add=T)
    text(1170000,6980000,"Density SE", cex=1.3)
    dev.off()
    
    png(file=paste0("D:/MCFN/BRT outputs/",species,"/SE_",species,"_",buffer,"_2.png"), height=600, width=850)
    par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
    par(mar=c(0,0,5,0))
    plot(SD_ext, col=magma2(15)[15], axes=TRUE, legend=FALSE, main=paste0(species," ",as.numeric(buffer)/1000," Km buffer"))
    plot(SD_ext, col=magma2(15), axes=T,  main=paste0(species," ",buffer,"m buffer SD"), legend.width=1.5, horizontal = TRUE, smallplot = c(0.65,0.90,0.83,0.87), axis.args=list(cex.axis=1.2),add=T)
    text(1170000,6980000,"Density SE", cex=1.3)
    dev.off()
    
  }
  
  if (CV==TRUE){  
    CV<-raster(paste0("D:/MCFN/BRT outputs/",species,"/buffer_",buffer,"m/confint/CV.tif"))
    CV_ext<-extend(CV,100)
    
    png(file=paste0("D:/MCFN/BRT outputs/",species,"/CV_",species,"_",buffer,".png"), height=600, width=850)
    par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
    par(mar=c(0,0,5,0))
    plot(CV_ext, col=magma(15)[15], axes=TRUE, legend=FALSE, main=paste0(species," ",as.numeric(buffer)/1000," Km buffer"))
    plot(CV_ext, col=magma(15),  axes=T,  main=paste0(species," ",buffer,"m buffer CV"), legend.width=1.5, horizontal = TRUE, smallplot = c(0.65,0.90,0.83,0.87),axis.args=list(cex.axis=1.2),add=T)
    text(1170000,6980000,"Density CV", cex=1.3)
    dev.off()
    
    
    png(file=paste0("D:/MCFN/BRT outputs/",species,"/CV_",species,"_",buffer,"_2.png"), height=600, width=850)
    par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0))
    par(mar=c(0,0,5,0))
    plot(CV_ext, col=magma2(15)[15], axes=TRUE, legend=FALSE, main=paste0(species," ",as.numeric(buffer)/1000," Km buffer"))
    plot(CV_ext, col=magma2(15),  axes=T,  main=paste0(species," ",buffer,"m buffer CV"), legend.width=1.5, horizontal = TRUE,smallplot = c(0.65,0.90,0.83,0.87), axis.args=list(cex.axis=1.2),add=T)
    text(1170000,6980000,"Density CV", cex=1.3)
    dev.off()
  }
  
}

plot_preds("CAWA","200000",SD=TRUE,CV=TRUE)
beep(sound = 4, expr = NULL) 
plot_preds("OSFL","100000",SD=TRUE,CV=TRUE)
plot_preds("CONI","400000",SD=TRUE,CV=TRUE)
plot_preds("RUBL","450000",SD=TRUE,CV=TRUE)
beep(sound = 4, expr = NULL) 

save.image("D:/MCFN/BRT outputs/out1.RData")
beep(sound = 4, expr = NULL) 


# Summarize predictions per land cover type ####
NALCC2005<-raster("M:/DataStuff/SpatialData/LCC05_NALCMS2010/Land_Cover_MXD/NA_LandCover_2005/data/NA_LandCover_2005/NA_LandCover_2005_LCC.img")
NALCC2005<-projectRaster(NALCC2005,pred_MCFN_mask.ext[[1]], method="ngb")
MCFN_LCC<- crop(NALCC2005,pred_MCFN_mask.ext[[1]])
MCFN_LCC<- mask(MCFN_LCC,pred_MCFN_mask.ext[[1]])
plot(MCFN_LCC)
ltnalc <- read.csv("C:/Users/voeroesd/Documents/Repos/bamanalytics/lookup/nalcms.csv")




# CAWA
mod<-list5$CAWA$CAWA_buffer_200000m
rast<-predict(pred_MCFN_mask.ext,
              mod,
              type = "response",
              n.trees = mod$n.trees)

dens<-data.frame(Value=getValues(MCFN_LCC),Density=getValues(rast))
dens<-dens[!is.na(dens),]
dens$class<-ltnalc$Label[match(dens$Value, ltnalc$Value)]
medianDens<-aggregate(dens$Density,by=list(dens$class), FUN=median, na.rm=T)

medianDens<-aggregate(dens$Density,by=list(dens$class), FUN=function(x){quantile(x,0.75,na.rm=T)})

areaHab<-aggregate(dens$Density,by=list(dens$class), FUN=function(x){6.25* length(x)})
hab_dens<-data.frame(Habitat = medianDens$Group.1,Median.Density=medianDens$x,areaHab=areaHab$x/sum(areaHab$x))
hab_dens$area2<-hab_dens$areaHab
hab_dens$area2[which(hab_dens$areaHab<0.01)]<-0.01
hab_dens$labels<-as.character(round(hab_dens$areaHab,digits=2))
hab_dens$labels[which(hab_dens$areaHab<0.01)]<-"<0.01"

# agriculture has only 4 grid cells, in which density is overpredicted. will replace them with 75% quantile value
quantile(dens$Density,0.75,na.rm=T)
hab_dens$Median.Density[1]<-quantile(dens$Density,0.75,na.rm=T)

png(filename = paste0("D:/MCFN/BRT outputs/CAWA/lc_mean_denspreds.png"), height=600, width=850)
ggplot(hab_dens,aes(x=Habitat,y=Median.Density))+
  geom_bar(aes(fill=Habitat),width=hab_dens$area2, stat="identity")+
  xlab("Land cover type")+
  ylab("Mean density (males/ha)")+
  theme(legend.position = "none")+
  geom_text(aes(label=hab_dens$labels),nudge_y = 0.0005,size=3.5)
dev.off()

#OSFL
mod<-list5$OSFL$OSFL_buffer_100000m
rast<-predict(pred_MCFN_mask.ext,
              mod,
              type = "response",
              n.trees = mod$n.trees)
dens<-data.frame(Value=getValues(MCFN_LCC),Density=getValues(rast))
dens<-dens[!is.na(dens),]
dens$class<-ltnalc$Label[match(dens$Value, ltnalc$Value)]
medianDens<-aggregate(dens$Density,by=list(dens$class), FUN=median, na.rm=T)
areaHab<-aggregate(dens$Density,by=list(dens$class), FUN=function(x){6.25* length(x)})
hab_dens<-data.frame(Habitat = medianDens$Group.1,Median.Density=medianDens$x,areaHab=areaHab$x/sum(areaHab$x))
hab_dens$area2<-hab_dens$areaHab
hab_dens$area2[which(hab_dens$areaHab<0.01)]<-0.01
hab_dens$labels<-as.character(round(hab_dens$areaHab,digits=2))
hab_dens$labels[which(hab_dens$areaHab<0.01)]<-"<0.01"

png(filename = paste0("D:/MCFN/BRT outputs/OSFL/lc_mean_denspreds.png"), height=600, width=850)
ggplot(hab_dens,aes(x=Habitat,y=Median.Density))+
  geom_bar(aes(fill=Habitat),width=hab_dens$area2, stat="identity")+
  xlab("Land cover type")+
  ylab("Mean density (males/ha)")+
  theme(legend.position = "none")+
  geom_text(aes(label=hab_dens$labels),nudge_y = 0.0002,size=3.5)
dev.off()

# CONI
mod<-list5$OSFL$OSFL_buffer_400000m
rast<-predict(pred_MCFN_mask.ext,
              mod,
              type = "response",
              n.trees = mod$n.trees)

dens<-data.frame(Value=getValues(MCFN_LCC),Density=getValues(rast))
dens<-dens[!is.na(dens),]
dens$class<-ltnalc$Label[match(dens$Value, ltnalc$Value)]
medianDens<-aggregate(dens$Density,by=list(dens$class), FUN=median, na.rm=T)
areaHab<-aggregate(dens$Density,by=list(dens$class), FUN=function(x){6.25* length(x)})
hab_dens<-data.frame(Habitat = medianDens$Group.1,Median.Density=medianDens$x,areaHab=areaHab$x/sum(areaHab$x))
hab_dens$area2<-hab_dens$areaHab
hab_dens$area2[which(hab_dens$areaHab<0.01)]<-0.01
hab_dens$labels<-as.character(round(hab_dens$areaHab,digits=2))
hab_dens$labels[which(hab_dens$areaHab<0.01)]<-"<0.01"

png(filename = paste0("D:/MCFN/BRT outputs/CONI/lc_mean_denspreds.png"), height=600, width=850)
ggplot(hab_dens,aes(x=Habitat,y=Median.Density))+
  geom_bar(aes(fill=Habitat),width=hab_dens$area2, stat="identity")+
  xlab("Land cover type")+
  ylab("Mean density (males/ha)")+
  theme(legend.position = "none")+
  geom_text(aes(label=hab_dens$labels),nudge_y = 0.0002,size=3.5)
dev.off()

# RUBL
mod<-list5$RUBL$RUBL_buffer_450000m
rast<-predict(pred_MCFN_mask.ext,
              mod,
              type = "response",
              n.trees = mod$n.trees)

dens<-data.frame(Value=getValues(MCFN_LCC),Density=getValues(rast))
dens<-dens[!is.na(dens),]
dens$class<-ltnalc$Label[match(dens$Value, ltnalc$Value)]
medianDens<-aggregate(dens$Density,by=list(dens$class), FUN=median, na.rm=T)
areaHab<-aggregate(dens$Density,by=list(dens$class), FUN=function(x){6.25* length(x)})
hab_dens<-data.frame(Habitat = medianDens$Group.1,Median.Density=medianDens$x,areaHab=areaHab$x/sum(areaHab$x))
hab_dens$area2<-hab_dens$areaHab
hab_dens$area2[which(hab_dens$areaHab<0.01)]<-0.01
hab_dens$labels<-as.character(round(hab_dens$areaHab,digits=2))
hab_dens$labels[which(hab_dens$areaHab<0.01)]<-"<0.01"

png(filename = paste0("D:/MCFN/BRT outputs/RUBL/lc_mean_denspreds.png"), height=600, width=850)
ggplot(hab_dens,aes(x=Habitat,y=Median.Density))+
  geom_bar(aes(fill=Habitat),width=hab_dens$area2, stat="identity")+
  xlab("Land cover type")+
  ylab("Mean density (males/ha)")+
  theme(legend.position = "none")+
  geom_text(aes(label=hab_dens$labels),nudge_y = 0.0005,size=3.5)
dev.off()


save.image("D:/MCFN/BRT outputs/out1.RData")
beep(sound = 4, expr = NULL) 

# Cleanup ws ####
load("D:/MCFN/BRT outputs/out1.RData")

rm(list=setdiff(ls(),c("pred_MCFN_mask.ext","LCC","list4","list5","brt_MCFN","CI_SE_CV_plotfun","plot_preds")))
gc()

