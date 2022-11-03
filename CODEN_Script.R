library(lme4)
library(merTools)
library(lmerTest)
library(sjstats)
library(car)
library(emmeans)
library(dplyr)
library(ggplot2)
library(interactions)
library(MuMIn)
library(assertive.types)
library(sjPlot)
library(Hmisc)
library(corrplot)
library(rcompanion)
library(maptools)
library(rgdal)
library(spdep)
library(ctv)
library(tmap)
library(sf)
library(ape)
library(ggplot2)
library(deldir)
library(spatialreg)
library(jtools)
library(mgcv)
library(MASS)
library(sandwich)
library(lmtest)
library(pscl)
library(sp)
library(broom)

options(scipen = 999)

###Import the correct geodatabase feature class
fgdb <- "FinalDataset.gdb"
subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(fgdb)
print(fc_list)
dat <- readOGR(dsn=fgdb,layer="FinalDataset_adj_na") #SPECIFIC FEATURE CLASS
summary(dat)
plot(dat)

#Assign coordinates
coords_C<-coordinates(dat)

#Create the nb object with four neighbors
fc_kn4_C <- knn2nb(knearneigh(coords_C, k = 4))

#Make the nb object symmetrical (https://r-spatial.github.io/spdep/reference/testnb.html)
fc_kn4_C_s <- make.sym.nb(fc_kn4_C)

#Create the listw object from the nb object (fc_kn4) for four closest neighbors
fc_kn4_C_s_w <- nb2listw(fc_kn4_C_s, zero.policy = TRUE)

###Declare some variables as factors and rescale some of them
#Variables as factors
dat$RE <- as.factor(dat$RE)
dat$Sex <- as.factor(dat$Sex)

#Re-scale variables
dat$pk_ovrlpK <- dat$pk_ovrlp/1000000 ##square km instead of square meters
dat$resDensK <- dat$resDens*1000000 ##Residential units per square kilometers instead of per square meters

#Recode park distance variable based on 1/2 mile distance
dat$park_dist_half[dat$park_dist< 804] <- 1
dat$park_dist_half[dat$park_dist> 804] <- 0

#Recode highway distance variable based on 1/4 mile distance
dat$hwy_dist_quarter[dat$hwy_dist<= 400] <- 1
dat$hwy_dist_quarter[dat$hwy_dist> 400] <- 0

###Standardize the continuous variables that go in the regression
dat$AgeAtDx_scaled <- scale(dat$AgeAtDx)
dat$BMI_scaled <-scale(dat$BMI)
dat$PM25_scaled <- scale(dat$PM25)
dat$ovrcrwd_scaled <- scale(dat$ovrcrwd)
dat$p_multif_scaled <- scale(dat$p_multif)
dat$per_hb_scaled <- scale(dat$per_hb)
dat$per_essen_scaled <- scale(dat$per_essen)
dat$per_tew_scaled <- scale(dat$per_tew)
dat$resDensK_scaled <- scale(dat$resDensK)
dat$bikeScore_int_scaled <- scale(dat$bikeScore_int)
dat$walkScore_int_scaled <- scale(dat$walkScore_int)
dat$tranScore_int_scaled <- scale(dat$tranScore_int)
dat$ndvi_scaled <- scale(dat$ndvi)
dat$pk_ovrlpK_scaled <- scale(dat$pk_ovrlpK)



##Run the ME model to calculate eigenvectors for all cases
model_splag_C <- ME(Hospitalized ~ AgeAtDx_scaled + Sex + RE + BMI_scaled + Tobacco + Diabetes + Hypertension + KidneyDisease + LungDisease + apartment + PM25_scaled + ovrcrwd_scaled + p_multif_scaled + per_hb_scaled + per_essen_scaled + per_tew_scaled + resDensK_scaled + bikeScore_int_scaled + walkScore_int_scaled + tranScore_int_scaled + hwy_dist_quarter + ndvi_scaled + pk_ovrlpK_scaled + park_dist_half, data = dat, family="poisson", listw=fc_kn4_C_s_w)

##Now calculate the model's coefficients for all cases - Poisson
glm_model_splag_C <- glm(Hospitalized ~ AgeAtDx_scaled + Sex + RE + BMI_scaled + Tobacco + Diabetes + Hypertension + KidneyDisease + LungDisease + apartment + PM25_scaled + ovrcrwd_scaled + p_multif_scaled + per_hb_scaled + per_essen_scaled + per_tew_scaled + resDensK_scaled + bikeScore_int_scaled + walkScore_int_scaled + tranScore_int_scaled + hwy_dist_quarter + ndvi_scaled + pk_ovrlpK_scaled + park_dist_half + fitted(model_splag_C), data = dat, family="poisson")
fsum <- summary(glm_model_splag_C)
sum <- tidy(glm_model_splag_C)
nk <- nagelkerke(glm_model_splag_C)
vif <- vif(glm_model_splag_C)
emmeans(glm_model_splag_C, list(pairwise ~ RE), adjust = "tukey")

##Now calculate the Poisson model with robust standard errors - all cases
cov.pois <- vcovHC(glm_model_splag_C, type="HC0")
std.err <- sqrt(diag(cov.pois))
r.est <- cbind(Estimate= coef(glm_model_splag_C), "Robust SE" = std.err,
"Pr(>|z|)" = 2 * pnorm(abs(coef(glm_model_splag_C)/std.err), lower.tail=FALSE),
LL = coef(glm_model_splag_C) - 1.96 * std.err,
UL = coef(glm_model_splag_C) + 1.96 * std.err)
r.est
ex_re <- exp(r.est)

save.image(file = "CODEN_Poisson.RData")
