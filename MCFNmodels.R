library(mefa4)
library(opticut)
library(MASS)

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
  
  
  ol <- find_levels(spp, m=1000) # use subset of offroad data
  rc <- ol$levels[[length(ol$levels)]]
  
  if(road=="yes"){
    ff <- y ~ x + ROAD + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT # with road
    y <- as.numeric(xtab(mmi)[,spp])
    x <- droplevels(samp(mmi)$HAB_NALC1)
    if (!is.null(reclass))
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
    if (!is.null(reclass))
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


## example: 
m1 <- model_1(spp = "CAWA", buffer=BUFF100)
m2 <- model_1(spp = "CONI", buffer=BUFF500)

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
mods_CONI$buffer100$model

"CONI"%in%colnames(OFF)




mods$buffer200$model
mods$buffer100$model
  
  
# Estimates and CIs for each buffer
ests_CIs<- function(model){
  mvsamps<- mvrnorm(n=1000,mu=model$coef,Sigma=model$vcov) # CIs estimated by drawing from a multivariate distribution with coefs and vcov, and quantiles
  CIs<-t(apply(mvsamps,2,quantile,probs=c(0,0.9)))
  ests<-cbind(coef=model$coef,CIs)
  ests
}

ests_CIs(mods$buffer0$model)
ests_CIs(mods$buffer50$model)
ests_CIs(mods$buffer100$model)
ests_CIs(mods$buffer150$model)
ests_CIs(mods$buffer200$model)
ests_CIs(mods$buffer250$model)
ests_CIs(mods$buffer300$model)
ests_CIs(mods$buffer350$model)
ests_CIs(mods$buffer400$model)
ests_CIs(mods$buffer450$model)
ests_CIs(mods$buffer500$model)

# the following functions are used to extract density predictions?
h <- function(x) {
    if (inherits(x, "try-error"))
        return(0)
    data.frame(D=round(sort(exp(c(x$coef[1], x$coef[1]+x$coef[-1]))), 4))
    }


h(mods$buffer500$model)



g <- function(x) {
    z <- x$model
    logD <- exp(c(z$coef[1], z$coef[1]+z$coef[-1]))
    names(logD) <- substr(names(logD), 2, nchar(names(logD)))
    rc <- x$levels$levels[[length(x$levels$levels)]]
    rc <- unique(unname(rc))
    names(logD)[1] <- rc[!(rc %in% names(logD))]
    data.frame(logD)
}


g(mods$buffer250)



LEV <- c("ConifDense", "Agr", "ConifOpen", "ConifSparse", "DecidDense",
    "DecidOpen", "DecidSparse", "Devel", "Grass", "MixedDense", "MixedOpen",
    "MixedSparse", "Shrub", "WetDense", "WetOpen", "WetSparse")
h2 <- function(x) {
    if (inherits(x, "try-error"))
        return(structure(rep(NA, length(LEV)), names=LEV))
    x <- data.frame(D=round(sort(exp(c(x$coef[1], x$coef[1]+x$coef[-1]))), 4))
    out <- list()
    for (i in 1:length(x$D)) {
        a <- rownames(x)[i]
        a <- substr(a, 2, nchar(a))
        a <- strsplit(a, "\\+")[[1]]
        out[[i]] <- structure(rep(x$D[i], length(a)), names=a)
    }
    out <- unlist(out)
    names(out)[names(out) == "Intercept)"] <- LEV[1]
    out <- out[match(LEV, names(out))]
    names(out) <- LEV
    out
}

h2(mods$buffer100$model)


