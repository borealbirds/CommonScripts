
rm(list=ls(all=TRUE))


library(mefa4)
library(QPAD)

workdir <- "M:\\Boreal\\DataStuff\\AvianData\\Processed"
setwd(workdir)

load("BAM-BBS-tables_20170630.Rdata")

load("new_offset_data_package_2017-03-01.Rdata")

nrow(nonDuplicated(dat, PKEY))
nrow(nonDuplicated(dat, SS))

nrow(nonDuplicated(pc, PKEY))
nrow(nonDuplicated(pc, SS))
nrow(nonDuplicated(PCTBL, PKEY))

# Specify QPAD version 3
load_BAM_QPAD(version = 3)    # 
getBAMversion()

# Get the species list from QPAD
sppp <- getBAMspecieslist()


# Create the dataframe for the offsets from the PKEY and SS table
offdat <- data.frame(PKEY[, c("PCODE", "PKEY", "SS", "TSSR", "JDAY", "MAXDUR", "MAXDIS", "JULIAN")], 
                     SS[match(PKEY$SS, rownames(SS)), c("TREE", "TREE3", "HAB_NALC1", "HAB_NALC2", "SPRNG")])
offdat$srise <- PKEY$srise + PKEY$MDT_offset           # Calculate sunrise time
offdat$DSLS <- (offdat$JULIAN - offdat$SPRNG) / 365    # Calculate days since local spring
summary(offdat)

# Create the quadratics of julian day, time since sunrise, and days since local spring
offdat$JDAY2 <- offdat$JDAY^2
offdat$TSSR2 <- offdat$TSSR^2
offdat$DSLS2 <- offdat$DSLS^2

# Reclassify HAB_NALC2 into 4 classes: DecidMixed, Conif, Open, and Wet
offdat$LCC4 <- as.character(offdat$HAB_NALC2)
offdat$LCC4[offdat$LCC4 %in% c("Decid", "Mixed")] <- "DecidMixed"
offdat$LCC4[offdat$LCC4 %in% c("Agr", "Barren", "Devel", "Grass", "Shrub")] <- "Open"
offdat$LCC4 <- factor(offdat$LCC4, c("DecidMixed", "Conif", "Open", "Wet"))    # Convert them to factors

# Reclassify LCC4 into 2 classes: Forest and OpenWet
offdat$LCC2 <- as.character(offdat$LCC4)
offdat$LCC2[offdat$LCC2 %in% c("DecidMixed", "Conif")] <- "Forest"
offdat$LCC2[offdat$LCC2 %in% c("Open", "Wet")] <- "OpenWet"
table(offdat$LCC4, offdat$LCC2)
offdat$MAXDIS <- offdat$MAXDIS / 100

# Create the variable matrix for time removal
Xp <- cbind("(Intercept)" = 1, as.matrix(offdat[, c("TSSR", "JDAY", "DSLS", "TSSR2", "JDAY2", "DSLS2")]))
dim(Xp)

# Create the variable matrix for distance sampling
Xq <- cbind("(Intercept)" = 1, TREE = offdat$TREE, LCC2OpenWet = ifelse(offdat$LCC2 == "OpenWet", 1, 0), 
            LCC4Conif = ifelse(offdat$LCC4 == "Conif", 1, 0), 
            LCC4Open = ifelse(offdat$LCC4 == "Open", 1, 0), 
            LCC4Wet = ifelse(offdat$LCC4 == "Wet", 1, 0))

# Create the empty offset database
OFF <- matrix(NA, nrow(offdat), length(sppp))
rownames(OFF) <- offdat$PKEY
colnames(OFF) <- sppp


# spp <- "ALFL" # If you want, try it with a single species first

# Loop over all the species. QPAD uses BIC to pick the best model and then predicts the offsets at all of the points
for(spp in sppp){
cat(spp, "\n"); flush.console()
p <- rep(NA, nrow(offdat))
A <- q <- p

# ????
cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))

# Find the best model
mi <- bestmodelBAMspecies(spp, type = "BIC")
cfi <- coefBAMspecies(spp, mi$sra, mi$edr)

# Get the variables for predicting the time removal offsets based on the estimated best model
Xp2 <- Xp[, names(cfi$sra), drop = FALSE]
OKp <- rowSums(is.na(Xp2)) == 0

# Get the variables for predicting the distance sampling offsets based on the estimated best model
Xq2 <- Xq[, names(cfi$edr), drop = FALSE]
OKq <- rowSums(is.na(Xq2)) == 0


p[!OKp] <- sra_fun(offdat$MAXDUR[!OKp], cf0[1])             # Add default coefficient of rowsum = 0
unlim <- ifelse(offdat$MAXDIS[!OKq] == Inf, TRUE, FALSE)    # Unlimited yes/no
A[!OKq] <- ifelse(unlim, pi*cf0[2]^2, pi*offdat$MAXDIS[!OKq]^2)
q[!OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[!OKq], cf0[2]))

phi1 <- exp(drop(Xp2[OKp, , drop=FALSE] %*% cfi$sra))
tau1 <- exp(drop(Xq2[OKq, , drop=FALSE] %*% cfi$edr))
p[OKp] <- sra_fun(offdat$MAXDUR[OKp], phi1)
unlim <- ifelse(offdat$MAXDIS[OKq] == Inf, TRUE, FALSE)
A[OKq] <- ifelse(unlim, pi*tau1^2, pi*offdat$MAXDIS[OKq]^2)
q[OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[OKq], tau1))


ii <- which(p == 0)
p[ii] <- sra_fun(offdat$MAXDUR[ii], cf0[1])

OFF[, spp] <- log(p) + log(A) + log(q)

}

head(OFF[, 1:10])

# Checks to make sure none of the estimates spun off into inifinity
(Ra <- apply(OFFY, 2, range))
summary(t(Ra))
which(!is.finite(Ra[1,])) # BARS GCSP
which(!is.finite(Ra[2,]))


