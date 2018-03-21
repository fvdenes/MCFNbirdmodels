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

str(BUFF0)

for(i in 1:nrow(DAT)){
  
}

BUFF500

table(as.character(DAT$SS)%in%BUFF500$SS)

T <- DAT[which(DAT$SS %in% BUFF0$SS), ]
str(T)
table(as.character(T$SS))



# Add columns to DAT that identify data within each buffer
DAT$BUFF0<-NA
DAT$BUFF0

DAT$isBBS <- startsWith(rownames(DAT), "BBS")
table(DAT$isBBS, DAT$ROAD)
DAT$JBCR <- interaction(DAT$JURS, DAT$xBCR, sep="::", drop=TRUE)
DAT$RoadBBS <- interaction(DAT$ROAD, DAT$isBBS, sep="::", drop=TRUE)



DAT$HAB_NALC1 <- DAT$HABTR
DAT$HAB_NALC2 <- DAT$HAB
DAT$YEAR <- DAT$YR+2013

mm <- Mefa(YY, DAT, TAX, "inner")

mm <- mm[!is.na(samp(mm)$RoadBBS) & !is.na(samp(mm)$JBCR) & !is.na(samp(mm)$HAB_NALC1),]

mm<- mm[samp(mm)$SS%in% BUFF0$SS ,]

# helper functions
get_subset <- function(region, road=FALSE) {
  r <- if (road)
    samp(mm)$RoadBBS == "1::TRUE" else samp(mm)$RoadBBS == "0::FALSE"
  keep <- if (region == "all")
    rep(TRUE, nrow(mm)) else samp(mm)$JBCR == region
  mm[keep & r,]
}
tab_data <- function(spp, region) {
    r <- samp(mm)$RoadBBS %in% c("1::TRUE", "0::FALSE")
    keep <- if (region == "all")
      rep(TRUE, nrow(mm)) else samp(mm)$JBCR == region
    mmi <- mm[keep & r,]
    y <- factor(ifelse(as.numeric(xtab(mmi)[,spp]) > 0, 1, 0), 0:1)
    x <- samp(mmi)$HAB_NALC1
    list(off=table(NALC=samp(mmi)$HAB_NALC1[samp(mmi)$ROAD==0], Det=y[samp(mmi)$ROAD==0]),
        on=table(NALC=samp(mmi)$HAB_NALC1[samp(mmi)$ROAD==1], Det=y[samp(mmi)$ROAD==1]))
}
find_levels <- function(spp, region, road=FALSE, m=1000) {
    mmi <- get_subset(region, road)
    j <- rep(FALSE, nrow(mmi))
    for (k in levels(droplevels(samp(mm)$HAB_NALC1))) {
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
    ol <- optilevels(y=y, x=x, dist="poisson", offset=OFF[rownames(mmi),spp])
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

## model functions
model_null <- function(spp, region, road=FALSE, trend=FALSE) {
    ff0 <- y ~ 1
    ff1 <- y ~ yr
    ff <- if (trend)
        ff1 else ff0
    mmi <- get_subset(region, road)
    y <- as.numeric(xtab(mmi)[,spp])
    yr <- samp(mmi)$YEAR
    glm_skeleton(try(
        glm(ff, family="poisson", offset=OFF[rownames(mmi),spp])),
        keep_call=FALSE, vcov=TRUE)
}
model_lcc <- function(spp, region, road=FALSE, reclass=NULL, trend=FALSE) {
    ff0 <- y ~ x
    ff1 <- y ~ yr + x
    ff <- if (trend)
        ff1 else ff0
    mmi <- get_subset(region, road)
    y <- as.numeric(xtab(mmi)[,spp])
    x <- droplevels(samp(mmi)$HAB_NALC1)
    if (!is.null(reclass))
        levels(x) <- reclass[levels(x)]
    yr <- samp(mmi)$YEAR
    glm_skeleton(try(
        glm(ff, family="poisson", offset=OFF[rownames(mmi),spp])),
        keep_call=FALSE, vcov=TRUE)
}
model_clim <- function(spp, region, road=FALSE, trend=FALSE) {
    ff0 <- y ~ CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 +
        CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT
    ff1 <- y ~ yr + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 +
        CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT
    ff <- if (trend)
        ff1 else ff0
    mmi <- get_subset(region, road)
    y <- as.numeric(xtab(mmi)[,spp])
    yr <- samp(mmi)$YEAR
    glm_skeleton(try(glm(ff, data=samp(mmi),
        family="poisson", offset=OFF[rownames(mmi),spp])),
        keep_call=FALSE, vcov=TRUE)
}
model_lcclim <- function(spp, region, road=FALSE, reclass=NULL, trend=FALSE) {
    ff0 <- y ~ x + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 +
        CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT
    ff1 <- y ~ yr + x + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 +
        CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT
    ff <- if (trend)
        ff1 else ff0
    mmi <- get_subset(region, road)
    y <- as.numeric(xtab(mmi)[,spp])
    x <- droplevels(samp(mmi)$HAB_NALC1)
    if (!is.null(reclass))
        levels(x) <- reclass[levels(x)]
    yr <- samp(mmi)$YEAR
    glm_skeleton(try(glm(ff, data=samp(mmi),
        family="poisson", offset=OFF[rownames(mmi),spp])),
        keep_call=FALSE, vcov=TRUE)
}

model_all <- function(region, spp) {
    t0 <- proc.time()
    tab <- tab_data(spp, region)
    Err <- structure("0 detections", class="try-error")
    if (sum(tab$off[,"1"]) > 0) {
        set.seed(1)
        ol <- find_levels(spp, region, road=FALSE, m=1000) # use subset of offroad data
        rc <- ol$levels[[length(ol$levels)]]
        D00 <- model_null(spp, region, road=FALSE)
        D01 <- model_lcc(spp, region, road=FALSE, rc)
        D02 <- model_clim(spp, region, road=FALSE)
        D03 <- model_lcclim(spp, region, road=FALSE, rc)
        T00 <- model_null(spp, region, road=FALSE, trend=TRUE)
        T01 <- model_lcc(spp, region, road=FALSE, reclass=rc, trend=TRUE)
        T02 <- model_clim(spp, region, road=FALSE, trend=TRUE)
        T03 <- model_lcclim(spp, region, road=FALSE, reclass=rc, trend=TRUE)

        if (sum(tab$on[,"1"]) > 0) {
            D10 <- model_null(spp, region, road=TRUE)
            D11 <- model_lcc(spp, region, road=TRUE, reclass=rc)
            D12 <- model_clim(spp, region, road=TRUE)
            D13 <- model_lcclim(spp, region, road=TRUE, reclass=rc)
            T10 <- model_null(spp, region, road=TRUE, trend=TRUE)
            T11 <- model_lcc(spp, region, road=TRUE, reclass=rc, trend=TRUE)
            T12 <- model_clim(spp, region, road=TRUE, trend=TRUE)
            T13 <- model_lcclim(spp, region, road=TRUE, reclass=rc, trend=TRUE)
        } else {
            D10 <- D11 <- D12 <- D13 <- Err
            T10 <- T11 <- T12 <- T13 <- Err
        }
    } else {
        D00 <- D01 <- D02 <- D03 <- Err
        T00 <- T01 <- T02 <- T03 <- Err
        if (sum(tab$on[,"1"]) > 0) {
            set.seed(1)
            ol <- find_levels(spp, region, road=TRUE, m=1000) # use subset of offroad data
            rc <- ol$levels[[length(ol$levels)]]
            D10 <- model_null(spp, region, road=TRUE)
            D11 <- model_lcc(spp, region, road=TRUE, reclass=rc)
            D12 <- model_clim(spp, region, road=TRUE)
            D13 <- model_lcclim(spp, region, road=TRUE, reclass=rc)
            T10 <- model_null(spp, region, road=TRUE, trend=TRUE)
            T11 <- model_lcc(spp, region, road=TRUE, reclass=rc, trend=TRUE)
            T12 <- model_clim(spp, region, road=TRUE, trend=TRUE)
            T13 <- model_lcclim(spp, region, road=TRUE, reclass=rc, trend=TRUE)
        } else {
            ol <- NULL
            D10 <- D11 <- D12 <- D13 <- Err
            T10 <- T11 <- T12 <- T13 <- Err
        }
    }
    out <- list(
        species=spp,
        region=if (missing(region)) "AllRegions" else region,
        levels=ol,
        table=tab,
        time=as.numeric(proc.time() - t0)[3L],
        density=list(
            null_off=D00,
            lcc_off=D01,
            clim_off=D02,
            lcclim_off=D03,
            null_on=D10,
            lcc_on=D11,
            clim_on=D12,
            lcclim_on=D13),
        trend=list(
            null_off=T00,
            lcc_off=T01,
            clim_off=T02,
            lcclim_off=T03,
            null_on=T10,
            lcc_on=T11,
            clim_on=T12,
            lcclim_on=T13))
    out
}

spp <- "CAWA"
res0 <- model_all(spp=spp, region="all") # this is the command to generate the model. to obtain confidence intervals, need to use vcov matrix, draw samples from a multivariate distribution and estimate quantiles.

str(res0)


h <- function(x) {
    if (inherits(x, "try-error"))
        return(0)
    data.frame(D=round(sort(exp(c(x$coef[1], x$coef[1]+x$coef[-1]))), 4))
}

xx <- h(res0$density$lcc_on)

g <- function(x) {
    z <- x$density$lcc_on
    logD <- exp(c(z$coef[1], z$coef[1]+z$coef[-1]))
    names(logD) <- substr(names(logD), 2, nchar(names(logD)))
    rc <- x$levels$levels[[length(x$levels$levels)]]
    rc <- unique(unname(rc))
    names(logD)[1] <- rc[!(rc %in% names(logD))]
    logD
}

xx <- h(res0$density$lcc_on)

## todo:
## OK - add year effect
## - calculate geographic discrepancies
## - compare D_null and T_null with geographic sampling
## - calculate fit in other regions
## - process prediction grid
## - check spatial patterns and change climate if needed
## - pl/cl???

#fl <- list.files("e:/peter/bam/2017/foam")
#SPP <- substr(sapply(strsplit(fl, "_"), "[[", 2), 1, 4)

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

Den <- list()

spp <- "CAWA"
for (spp in SPP) {
    fn <- paste0("C:/Users/voeroesd/Dropbox/BAM/Critical Habitat/CHID subunit delineation/results/results_", spp, ".Rdata")
    load(fn)
    res4 <- res1[c("ON::8", "QC::8", "ON::12", "QC::12")]
    Den[[spp]] <- sapply(res4, function(z) h2(z$density$lcc_on))
}
save(Den, file="C:/Users/voeroesd/Dropbox/BAM/Critical Habitat/CHID subunit delineation/results_CAWA.Rdata/")

f <- function(x) {
    x <- t(x)
    x / apply(x, 1, max, na.rm=TRUE)
}
i <- "CAWA"
barplot(f(Den[[i]]), main=i, beside=TRUE, legend.text=TRUE,
    col=c("tomato", "gold", "grey", "turquoise"), ylab="density / max density")

barplot(t(Den[[i]]), main=names(Den)[i], beside=TRUE)





model_gw <- function(spp) {
  ff <- y ~ HAB + ROAD + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 +
    CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT
  y <- as.numeric(xtab(mm)[,spp])
  glm_skeleton(try(glm(ff, data=samp(mmi),
                       family="poisson", offset=OFF[rownames(mmi),spp])),
               keep_call=FALSE, vcov=TRUE)
}
