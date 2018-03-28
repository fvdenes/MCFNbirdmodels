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
  preds<-data.frame(buffer0=rep(NA,nrow(rc_preddata$buffer0)),
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
  
  
  getCoef<- function(modlistbuffer){
    coef<-modlistbuffer$model$coef
    coef
  }
  

  z<-which(lapply(modlist[-c(1,2)],function(x){class(x$model)})=="glm_skeleton")
  
  coefs <-lapply(modlist[z+2],getCoef)
  
  
  for (i in 1:length(coefs)){
    if(ncol(mod_matrix_pred[[names(coefs[i])]])==length(coefs[[i]])){
      preds[,names(coefs)[i]]<- exp(mod_matrix_pred[[names(coefs[i])]]%*%coefs[[i]])*100 # km^2 vs ha diff is 100
    }
    else{
    x<-mod_matrix_pred[[names(coefs[i])]]
    x<-x[,-which(colnames(mod_matrix_pred[[names(coefs[i])]])%in%attr(coefs[[i]], "names")==F)]
        preds[,names(coefs)[i]]<- exp(x%*%coefs[[i]])*100 # km^2 vs ha diff is 100
    }
  }
  
 preds
}

preds_cawa<-get_preds(mods_CAWA,pred_data)


get_preds(mods_CAWA,pred_data)





# Estimates and CIs for each buffer

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
get_CIs(mods_RUBL)




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
  plot(y=unlist(roc_spp$auc),x=1:length(unlist(roc_spp$auc))-0.1,xaxt="n",ylab="AUC",xlab="buffer",ylim=c(0,1))
  axis(1,at=1:length(unlist(roc_spp$auc)),labels=attr(unlist(roc_spp$auc),"names"))
  points(y=unlist(roc_spp$aucbuff0),x=1:length(unlist(roc_spp$auc))+0.1,xaxt="n",ylab="AUC",pch=16)
  legend(1,0.4,legend=c("Buffer area", "MCFN homelands"), pch=c(1,16))
}

plot_AUC(roc_CAWA)

plot_AUC(roc_OSFL)
plot_AUC(roc_RUBL)
plot_AUC(roc_CONI)


plot_SE<-function(modlist, log=TRUE){
  if(log==TRUE){
    z<-which(lapply(modlist[-c(1,2)],function(x){class(x$model)})=="glm_skeleton")  
    SElist<-lapply(modlist[z+2],function(x){median(x$fitted$SEs,na.rm=T)})
    SElist0<-lapply(modlist[z+2],function(x){median(x$fittedBUFF0$SEs,na.rm=T)})
    
    plot(y=log(unlist(SElist)),x=1:length(SElist)-0.1,xaxt="n",ylab="log (median SE of fitted values)",xlab="buffer")
    axis(1,at=1:length(SElist),labels=attr(unlist(SElist),"names"))
    points(y=log(unlist(SElist0)),x=1:length(SElist0)+0.1,pch=16)
    legend(1,-2,legend=c("Buffer area", "MCFN homelands"), pch=c(1,16)) 
  }
  if(log==FALSE){
    z<-which(lapply(modlist[-c(1,2)],function(x){class(x$model)})=="glm_skeleton")  
    SElist<-lapply(modlist[z+2],function(x){median(x$fitted$SEs,na.rm=T)})
    SElist0<-lapply(modlist[z+2],function(x){median(x$fittedBUFF0$SEs,na.rm=T)})
    
    plot(y=unlist(SElist),x=1:length(SElist)-0.1,xaxt="n",ylab="median SE of fitted values",xlab="buffer")
    axis(1,at=1:length(SElist),labels=attr(unlist(SElist),"names"))
    points(y=unlist(SElist0),x=1:length(SElist0)+0.1,pch=16)
    legend(1,0.3,legend=c("Buffer area", "MCFN homelands"), pch=c(1,16))
  }
}
plot_SE(mods_CAWA)
plot_SE(mods_CAWA,F)
plot_SE(mods_OSFL)
plot_SE(mods_RUBL)
plot_SE(mods_CONI)
