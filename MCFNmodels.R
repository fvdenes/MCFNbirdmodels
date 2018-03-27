library(mefa4)
library(opticut)
library(MASS)
library(pbapply)
library(pROC)

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



# helper functions


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
  ol <- optilevels(y=y, x=x, dist="poisson", offset=OFF[rownames(mm),spp])
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


# Model function


model_1 <- function(spp, buffer = NULL, road="no", maxit=25) {
  if(is.null(buffer)==T){
    mmi<-mm }  
  else{
    mmi<- mm[samp(mm)$SS%in% buffer$SS ,]}
  road_table <- table(mmi@samp$ROAD)
  
  set.seed(1)
  
  find_levels <- function(spp, m=1000) {
    j <- rep(FALSE, nrow(mmi))
    for (k in levels(droplevels(samp(mmi)$HAB_NALC1))) {
      w <- which(samp(mmi)$HAB_NALC1 == k)
      if (length(w) < m) {
        j[w] <- TRUE
      } else {
        j[sample(w, m)] <- TRUE
      }
    }
    mmi <- mmi[j,]
    y <- as.numeric(xtab(mmi)[,spp])
    x <- droplevels(samp(mmi)$HAB_NALC1)
    if(spp%in%colnames(OFF)==T){
        ol <- optilevels(y=y, x=x, dist="poisson", offset=OFF[rownames(mmi),spp])
    }
    else{
        ol <- optilevels(y=y, x=x, dist="poisson")
    }
    ol
  }
  
  
  ol <- find_levels(spp, m=1000) # use subset of offroad data
  rc <- ol$levels[[length(ol$levels)]]
  
  if(road=="yes"){
    ff <- y ~ x + ROAD + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT # with road
    y <- as.numeric(xtab(mmi)[,spp])
    x <- droplevels(samp(mmi)$HAB_NALC1)
    levels(x) <- rc[levels(x)]
    
    if(spp%in%colnames(OFF)==T){
      model<- glm_skeleton(try(glm(ff, data=samp(mmi),family="poisson", offset=OFF[rownames(mmi),spp],maxit=maxit)), keep_call=FALSE, vcov=TRUE)
      #if running model with ROAD covariate returns NA road coef, check variable format (should be factor), then road_table object to see on and off road proportions
    }
    else{
      model<- glm_skeleton(try(glm(ff, data=samp(mmi),family="poisson", maxit=maxit)), keep_call=FALSE, vcov=TRUE)
      #if running model with ROAD covariate returns NA road coef, check variable format (should be factor), then road_table object to see on and off road proportions
    } 
  }
 
  if(road=="no"){
    ff <- y ~ x +  CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT # without road
    y <- as.numeric(xtab(mmi)[,spp])
    x <- droplevels(samp(mmi)$HAB_NALC1)
    levels(x) <- rc[levels(x)]
    if(spp%in%colnames(OFF)==T){
      model<- glm_skeleton(try(glm(ff, data=samp(mmi),family="poisson", offset=OFF[rownames(mmi),spp],maxit=maxit)), keep_call=FALSE, vcov=TRUE)
    }
    else{
      model<- glm_skeleton(try(glm(ff, data=samp(mmi),family="poisson", maxit=maxit)), keep_call=FALSE, vcov=TRUE)
    }

  }
  
  out<-list(spp=spp,
            levels=ol,
            model=model)
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
mods_CAWA$buffer100$model

mods_OSFL <- modelallbuffers(spp="OSFL", road="no",maxit=100)
mods_OSFL$buffer100$model

mods_RUBL <- modelallbuffers(spp="RUBL", road="no",maxit=100)
mods_RUBL$buffer100$model

mods_CONI <- modelallbuffers(spp="CONI", road="no",maxit=100)
mods_CONI$buffer300$model




# Load prediction dataset: MCFN homelands
pred_data <- read.csv("C:/Users/voeroesd/Dropbox/BAM/MCFN/output/BAM_pred_data_wBuff0.csv", sep="")
str(pred_data)

modlist<-mods_CAWA
preddata<-pred_data


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
  coefs <-lapply(modlist[-c(1,2)],getCoef)
  
  ests<-data.frame(buffer0=rep(NA,nrow(rc_preddata$buffer0)),
                   buffer50=rep(NA,nrow(rc_preddata$buffer50)),
                   buffer100=rep(NA,nrow(rc_preddata$buffer100)),
                   buffer150=rep(NA,nrow(rc_preddata$buffer150)),
                   buffer200=rep(NA,nrow(rc_preddata$buffer200)),
                   buffer250=rep(NA,nrow(rc_preddata$buffer250)),
                   buffer300=rep(NA,nrow(rc_preddata$buffer300)),
                   buffer350=rep(NA,nrow(rc_preddata$buffer350)),
                   buffer400=rep(NA,nrow(rc_preddata$buffer400)),
                   buffer450=rep(NA,nrow(rc_preddata$buffer450)),
                   buffer500=rep(NA,nrow(rc_preddata$buffer500)))
  for (i in 1:length(coefs)){
    ests[,i]<- mod_matrix_pred[[i]]%*%coefs[[i]]
  }
 
}



model.matrix(ff,rc_preddata$buffer0)

modlist<-mods_CAWA
preddata<-rc_pred_CAWA_250

mods_CAWA$buffer0$levels$levels
str(rc_preddata)

levels(rc_preddata$x)


rc_pred_RUBL_250<-rc_pred(mods_RUBL$buffer250,pred_data)
# Predict density
ff <- y ~ x +  CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT # model formula without road
x<- model.matrix(ff,data.frame(y=1,rc_pred_CAWA_250))




x%*%mods_RUBL$buffer250$model$coef # this is in ha. need to convert to km2


cbind(attr(mod_matrix_pred, "dimnames")[[2]],attr(mods_CAWA$buffer250$model$coef,"names"))

mods_OSFL <- modelallbuffers(spp="OSFL", road="no",maxit=100)
mods_OSFL$buffer100$model

mods_RUBL <- modelallbuffers(spp="RUBL", road="no",maxit=100)
mods_RUBL$buffer100$model

mods_CONI <- modelallbuffers(spp="CONI", road="no",maxit=100)
mods_CONI$buffer300$model




# Estimates and CIs for each buffer
ests_CIs<- function(model){
  mvsamps<- mvrnorm(n=1000,mu=model$coef,Sigma=model$vcov) # CIs estimated by drawing from a multivariate distribution with coefs and vcov, and quantiles
  CIs<-t(apply(mvsamps,2,quantile,probs=c(0,0.9)))
  ests<-cbind(coef=model$coef,CIs)
  ests
}

ests_CIs(mods_CAWA$buffer0$model)
ests_CIs(mods_CAWA$buffer50$model)
ests_CIs(mods_CAWA$buffer100$model)
ests_CIs(mods_CAWA$buffer150$model)
ests_CIs(mods_CAWA$buffer200$model)
ests_CIs(mods_CAWA$buffer250$model)
ests_CIs(mods_CAWA$buffer300$model)
ests_CIs(mods_CAWA$buffer350$model)
ests_CIs(mods_CAWA$buffer400$model)
ests_CIs(mods_CAWA$buffer450$model)
ests_CIs(mods_CAWA$buffer500$model)




ff <- y ~ x +  CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT


predict(mods_CAWA$buffer0$model,newdata=pred_data)

mods_CAWA$buffer0$levels$levels



## ROC/AUC
rocAll1 <- pblapply(1:ncol(mn), function(i) {
  pp <- mn[ss1,i] * exp(off1[ss1])
  roc(Y1[ss1], pp)
})
names(rocAll1) <- c("NULL", names(mods))
auc <- sapply(rocAll1, function(z) as.numeric(z$auc))
barplot(auc, ylim=c(0,1), space=0, ylab="AUC", xlab="Stages")
lines(0:length(mods)+0.5, auc, col=2, lwd=2)


