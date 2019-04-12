## - maxdur/maxdist
## - date
## - time
## - xy
## - Time zone
## - TREE NALC

library(maptools)
library(mefa4)

setwd("D:/MCFN/Additional bird data/")

# fl <- c("BCCA_PKEY_update2017.txt",
#         "BCCA_PTCOUNT_UPDATE2017.txt",
#         "BCCA_XYupdate2017.txt",
#         "MANBBA_PC_2016.txt",
#         "MANBBA_PKEY_2016.txt",
#         "MANBBA_XY_2016.txt",
#         "BCRoadside.txt",
#         "Pkey_QCAtlasv2.txt",
#         "PtCount_QcAtlasv2.txt",
#         "XY_QCAtlasv2.txt",
#         "OffsetIntersectionOct17_Atlasupdate.csv")
# 
# Ls <- lapply(fl, read.csv)
# names(Ls) <- sapply(strsplit(fl, "\\."), "[[", 1)
# 

fl<- c("MCFN_Pkey.csv",
       "MCFN_PtCount.csv",
       "MCFN_XY.csv",
       "MCFN_Tree_LC05_TZone.csv")
Ls <- lapply(fl, read.csv)
names(Ls) <- sapply(strsplit(fl, "\\."), "[[", 1)

# ## BC atlas has the incorrect SS -- it needs the pcode added
# ## should be pcode:site:station
# head(Ls$BCCA_XYupdate2017)
# head(Ls$BCCA_PKEY_update2017)
# head(Ls$BCCA_PTCOUNT_UPDATE2017)
# head(Ls$OffsetIntersectionOct17_Atlasupdate)
# head(Ls$OffsetIntersectionOct17_Atlasupdate[Ls$OffsetIntersectionOct17_Atlasupdate$PCODE=="BCCA",])
# 
# tmpfun <- function(x) {
#   x$SS <- as.character(x$SS)
#   PCODE <- paste0(as.character(x$PCODE), ":")
#   PCODE[PCODE != "BCCA:"] <- ""
#   x$SS <- paste0(PCODE, x$SS)
#   x$SS <- as.factor(x$SS)
#   x
# }
# Ls$BCCA_XYupdate2017 <- tmpfun(Ls$BCCA_XYupdate2017)
# head(Ls$BCCA_XYupdate2017)
# Ls$BCCA_PKEY_update2017 <- tmpfun(Ls$BCCA_PKEY_update2017)
# head(Ls$BCCA_PKEY_update2017)
# Ls$BCCA_PTCOUNT_UPDATE2017 <- tmpfun(Ls$BCCA_PTCOUNT_UPDATE2017)
# head(Ls$BCCA_PTCOUNT_UPDATE2017)
# Ls$OffsetIntersectionOct17_Atlasupdate <- tmpfun(Ls$OffsetIntersectionOct17_Atlasupdate)
# head(Ls$OffsetIntersectionOct17_Atlasupdate)
# head(Ls$OffsetIntersectionOct17_Atlasupdate[Ls$OffsetIntersectionOct17_Atlasupdate$PCODE=="BCCA",])

## SS

SS1 <- Ls$MCFN_XY
SS1$Y <- SS1$POINT_YNad83
SS1$X <- SS1$POINT_X_Nad83

cn <- c("SS", "PCODE", "X", "Y")
SS <-SS1[,cn]
rownames(SS) <- SS$SS
#write.csv(SS, row.names=FALSE, file="MCFN-update-xy.csv")

tmp <- Ls$MCFN_Tree_LC05_TZone
SS$TZ <- tmp$tzid[match(SS$SS, tmp$SS)]
SS$Tree <- tmp$TreeCov[match(SS$SS, tmp$SS)]
SS$NALCMS <- tmp$LCC05[match(SS$SS, tmp$SS)]

## PKEY

PK1 <- Ls$MCFN_Pkey
PK1$PCODE <- PK1$ProjectCode
PK1$YearCollected <- PK1$Year
PK1$MonthCollected <- PK1$Month
PK1$DayCollected <- PK1$Day
PK1$Min <- PK1$Minute

cn <- c("PKEY", "SS", "PCODE", "YearCollected", "MonthCollected", "DayCollected", "Hour", "Min")
PKEY <- PK1[,cn]
rownames(PKEY) <- PKEY$PKEY
PKEY <- data.frame(PKEY, SS[match(PKEY$SS, SS$SS),])
PKEY$SS.1 <- PKEY$PCODE.1 <- NULL

#levels(PKEY$TZ) <- c("", "America/Blanc-Sablon", "America/Dawson_Creek",
#    "America/Edmonton", "America/Goose_Bay", "America/Halifax", "America/Moncton",
#    "America/Montreal", "America/Rankin_Inlet", "America/Regina", "America/Toronto",
#    "America/Vancouver", "America/Winnipeg", "America/Yellowknife")
levels(PKEY$TZ) <- c("EDT")
# PKEY$TZ[PKEY$TZ==""] <- NA
# PKEY$TZ <- droplevels(PKEY$TZ)
# table(PKEY$TZ, useNA="a")
TZ <- PKEY$TZ
lttz <- read.csv("C:/Users/voeroesd/Documents/Repos/bamanalytics/lookup/tzone.csv")
lttz <- nonDuplicated(lttz, Timezone, TRUE)

YEAR <- PKEY$YearCollected
MM <- ifelse(PKEY$MonthCollected < 10,
             paste0("0", PKEY$MonthCollected), as.character(PKEY$MonthCollected))
DAY <- ifelse(PKEY$DayCollected < 10,
              paste0("0", PKEY$DayCollected), as.character(PKEY$DayCollected))
HH <- ifelse(PKEY$Hour < 10, paste0("0", PKEY$Hour), as.character(PKEY$Hour))
mm <- ifelse(PKEY$Min < 10, paste0("0", PKEY$Min), as.character(PKEY$Min))
DD <- with(PKEY, paste0(YEAR, "-", MM, "-", DAY, " ", HH, ":", mm, ":00"))
DD <- strptime(DD, "%Y-%m-%e %H:%M:%S", tz="EST5EDT") # EDt = UTC-4
PKEY$DATE <- DD
## Julian day
PKEY$JULIAN <- DD$yday # this is kept as original
PKEY$JDAY <- DD$yday / 365
summary(PKEY$JDAY)
## prevent too far extrapolation
#PKEY$JDAY[PKEY$JDAY < 0.35 | PKEY$JDAY > 0.55] <- NA
PKEY$JDAY[PKEY$JDAY < 0.35] <- 0.35
PKEY$JDAY[PKEY$JDAY > 0.55] <- 0.55
## TSSR = time since sunrise
Coor <- as.matrix(cbind(as.numeric(SS$X),as.numeric(SS$Y)))[match(PKEY$SS, rownames(SS)),]
JL <- as.POSIXct(DD)
subset <- rowSums(is.na(Coor))==0 & !is.na(JL) & !is.na(TZ)
sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", proj4string=CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs") ,POSIXct.out=FALSE) * 24
PKEY$srise <- NA
PKEY$srise[subset] <- sr
PKEY$start_time <- PKEY$Hour + PKEY$Min/60
PKEY$MDT_offset <- lttz$MDT_offset[match(TZ, rownames(lttz))]
table(TZ, PKEY$MDT_offset)
PKEY$TSSR <- (PKEY$start_time - PKEY$srise + PKEY$MDT_offset) / 24
PKEY$TSSR_orig <- PKEY$TSSR # keep a full copy
PKEY$TSSR[PKEY$start_time > 12] <- NA ## after noon
summary(PKEY$TSSR)
summary(PKEY$start_time)

PKEY$MAXDUR <- 5
PKEY$MAXDIS <- Inf

PKEY$TREE <- PKEY$Tree
PKEY$TREE[PKEY$TREE > 100] <- NA
PKEY$TREE[PKEY$TREE < 0] <- NA
PKEY$TREE <- PKEY$TREE / 100

ltnalc <- read.csv("C:/Users/voeroesd/Documents/Repos/bamanalytics/lookup/nalcms.csv")
PKEY$NALCMS2 <- ltnalc$Label[match(PKEY$NALCMS, ltnalc$Value)]

PKEY$WNALC <- PKEY$NALCMS2
levels(PKEY$WNALC)[levels(PKEY$WNALC) %in% c("Agr","Barren","Devel","Grass", "Shrub")] <- "Open"
PKEY$LCC2 <- as.factor(ifelse(PKEY$WNALC %in% c("Open", "Wet"), "OpenWet", "Forest"))
PKEY$LCC4 <- PKEY$WNALC
levels(PKEY$LCC4) <- c(levels(PKEY$LCC4), "DecidMixed")
PKEY$LCC4[PKEY$WNALC %in% c("Decid", "Mixed")] <- "DecidMixed"
PKEY$LCC4 <- droplevels(PKEY$LCC4)
PKEY$LCC4 <- relevel(PKEY$LCC4, "DecidMixed")

## offsets

ROOT <- "D:/MCFN/Additional bird data/"
library(mefa4)
library(QPAD)
source("C:/Users/voeroesd/Documents/Repos/bamanalytics/R/dataprocessing_functions.R")

load_BAM_QPAD(3)
getBAMversion()
sppp <- getBAMspecieslist()

offdat <- data.frame(PKEY[,c("PCODE","PKEY","SS","TSSR","JDAY","MAXDUR","MAXDIS","JULIAN",
                             "TREE","LCC2","LCC4")])
summary(offdat)

offdat$JDAY2 <- offdat$JDAY^2
offdat$TSSR2 <- offdat$TSSR^2
table(offdat$LCC4, offdat$LCC2)
offdat$MAXDIS <- offdat$MAXDIS / 100

Xp <- cbind("(Intercept)"=1, as.matrix(offdat[,c("TSSR","JDAY","TSSR2","JDAY2")]))
Xq <- cbind("(Intercept)"=1, TREE=offdat$TREE,
            LCC2OpenWet=ifelse(offdat$LCC2=="OpenWet", 1, 0),
            LCC4Conif=ifelse(offdat$LCC4=="Conif", 1, 0),
            LCC4Open=ifelse(offdat$LCC4=="Open", 1, 0),
            LCC4Wet=ifelse(offdat$LCC4=="Wet", 1, 0))

OFF <- matrix(NA, nrow(offdat), length(sppp))
rownames(OFF) <- offdat$PKEY
colnames(OFF) <- sppp

#spp <- "OVEN"
for (spp in sppp) {
  cat(spp, "\n");flush.console()
  p <- rep(NA, nrow(offdat))
  A <- q <- p
  
  ## constant for NA cases
  cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
  ## best model
  mi <- bestmodelBAMspecies(spp, model.sra=0:8, type="BIC")
  cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
  
  Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
  OKp <- rowSums(is.na(Xp2)) == 0
  Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
  OKq <- rowSums(is.na(Xq2)) == 0
  
  p[!OKp] <- sra_fun(offdat$MAXDUR[!OKp], cf0[1])
  unlim <- ifelse(offdat$MAXDIS[!OKq] == Inf, TRUE, FALSE)
  A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * offdat$MAXDIS[!OKq]^2)
  q[!OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[!OKq], cf0[2]))
  
  phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
  tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
  p[OKp] <- sra_fun(offdat$MAXDUR[OKp], phi1)
  unlim <- ifelse(offdat$MAXDIS[OKq] == Inf, TRUE, FALSE)
  A[OKq] <- ifelse(unlim, pi * tau1^2, pi * offdat$MAXDIS[OKq]^2)
  q[OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[OKq], tau1))
  
  ii <- which(p == 0)
  p[ii] <- sra_fun(offdat$MAXDUR[ii], cf0[1])
  
  OFF[,spp] <- log(p) + log(A) + log(q)
  
}

## checks
(Ra <- apply(OFF, 2, range))
summary(t(Ra))
which(!is.finite(Ra[1,]))
which(!is.finite(Ra[2,]))


PC1 <- Ls$MCFN_PtCount
colnames(PC1) <- toupper(colnames(PC1))
PC1$PCODE <- PC1$PROJECTCODE
PC1$SPECIES <- PC1$SPECIESCODE
PC1$DISTANCE <- PC1$BAMDISTANCE
PC1$DURATION <- PC1$PERIOD
PC1$ABUND <- PC1$SUMOFCOUNT

PC1$SS <- PKEY[match(PC1$PKEY, PKEY$PKEY),"SS"]
PC1$DURATION <- 1

speciesnames <- read.csv("C:/Users/voeroesd/Documents/Repos/bamanalytics/lookup/singing-species.csv")
PC1$SPECIES <-speciesnames[match(PC1$FULLNAME,speciesnames$English_Name),"Species_ID"]

cn <- c("PCODE", "PKEY", "SS", "SPECIES", "ABUND", "BEH", "DISTANCE","DURATION")
PCTBL <- PC1[,cn]

YY <- Xtab(ABUND ~ PKEY + SPECIES, PCTBL)

SS_MCFN<-SS
PCTBL_MCFN<-PCTBL
PKEY_MCFN<-PKEY
OFF_MCFN<-OFF
YY_MCFN<-YY

save(SS_MCFN, PCTBL_MCFN, PKEY_MCFN, OFF_MCFN, YY_MCFN,
     file="MCFN_data_processed-20190325.RData")

###
library(tidyr)
OFF_MCFN_long<-gather(OFF_MCFN,"SPECIES","logoffset",ALFL:YTVI,factor_key = F) # convert dataframe from wide to long format
OFF_ALL<-rbind(offl,OFF_MCFN_long)
