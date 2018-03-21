library(mefa4)
library(opticut)

#load data
load("C:/Users/voeroesd/Dropbox/BAM/Critical Habitat/CHID subunit delineation/pack_2016-12-01.Rdata")


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



#DAT$ROAD<-factor(DAT$ROAD) ##
DAT$HAB_NALC1 <- DAT$HABTR
DAT$HAB_NALC2 <- DAT$HAB
DAT$YEAR <- DAT$YR+2013

mm <- Mefa(YY, DAT, TAX, "inner")

mm <- mm[!is.na(samp(mm)$HAB_NALC1),]

mmi<- mm[samp(mm)$SS%in% BUFF100$SS ,]

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





model_1 <- function(spp, buffer = NULL) {
  if(is.null(buffer)==T){
    mmi<-mm }  
  else{
    mmi<- mm[samp(mm)$SS%in% buffer$SS ,]}
  road_table <- table(mmi@samp$ROAD)
  ol <- find_levels(spp, m=1000) # use subset of offroad data
  rc <- ol$levels[[length(ol$levels)]]
  #ff <- y ~ x + ROAD + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT
  ff <- y ~ x +  CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT
  y <- as.numeric(xtab(mmi)[,spp])
  x <- droplevels(samp(mmi)$HAB_NALC1)
  if (!is.null(reclass))
    levels(x) <- rc[levels(x)]
  model<- glm_skeleton(try(glm(ff, data=samp(mmi),family="poisson", offset=OFF[rownames(mmi),spp])), keep_call=FALSE, vcov=TRUE)
  out<-list(levels=ol,road_table=road_table,model=model)
  out
  }




m1 <- model_1(spp = "CAWA", buffer=BUFF100)

m1$road_table
m1$model$coef
m1$model$vcov

mvrnorm()

 # this is the command to generate the model. to obtain confidence intervals, need to use vcov matrix, draw samples from a multivariate distribution and estimate quantiles.


# the following functions are used to extract density predictions?
h <- function(x) {
    if (inherits(x, "try-error"))
        return(0)
    data.frame(D=round(sort(exp(c(x$coef[1], x$coef[1]+x$coef[-1]))), 4))
    }


h(m1$model)



g <- function(x) {
    z <- x$model
    logD <- exp(c(z$coef[1], z$coef[1]+z$coef[-1]))
    names(logD) <- substr(names(logD), 2, nchar(names(logD)))
    rc <- x$levels$levels[[length(x$levels$levels)]]
    rc <- unique(unname(rc))
    names(logD)[1] <- rc[!(rc %in% names(logD))]
    logD
}

g(m1)


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

h2(m1$model)


