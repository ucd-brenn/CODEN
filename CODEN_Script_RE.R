###Step 2: Run the above models for sub-samples by racial/ethnic group
load('CODEN_Poisson.RData')

##Start with HISPANIC

#Hispanic dataset
fc_H <- dat[ which(dat$RE=='Hispanic' ), ]

#Assign coordinates
coords_H<-coordinates(fc_H)

#Create the nb object with four neighbors
fc_kn4_H <- knn2nb(knearneigh(coords_H, k = 4))

#Make the nb object symmetrical
fc_kn4_H_s <- make.sym.nb(fc_kn4_H)

#Create the listw object from the nb object (fc_kn4) for four closest neighbors
fc_kn4_H_s_w <- nb2listw(fc_kn4_H_s, zero.policy = TRUE)

#Run the ME model to calculate eigenvectors for Hispanic
model_splag_H <- ME(Hospitalized ~ AgeAtDx_scaled + Sex + BMI_scaled + Tobacco + Diabetes + Hypertension + KidneyDisease + LungDisease + apartment + PM25_scaled + ovrcrwd_scaled + p_multif_scaled + per_hb_scaled + per_essen_scaled + per_tew_scaled + resDensK_scaled + bikeScore_int_scaled + walkScore_int_scaled + tranScore_int_scaled + hwy_dist_quarter + ndvi_scaled + pk_ovrlpK_scaled + park_dist_half, data = fc_H, family="poisson", listw=fc_kn4_H_s_w)

#Now calculate the model's coefficients for Hispanic - Poisson
glm_model_splag_H <- glm(Hospitalized ~ AgeAtDx_scaled + Sex + BMI_scaled + Tobacco + Diabetes + Hypertension + KidneyDisease + LungDisease + apartment + PM25_scaled + ovrcrwd_scaled + p_multif_scaled + per_hb_scaled + per_essen_scaled + per_tew_scaled + resDensK_scaled + bikeScore_int_scaled + walkScore_int_scaled + tranScore_int_scaled + hwy_dist_quarter + ndvi_scaled + pk_ovrlpK_scaled + park_dist_half + fitted(model_splag_H), data = fc_H, family="poisson")
fsum_H <- summary(glm_model_splag_H)
sum_H <- tidy(glm_model_splag_H)
nk_H <- nagelkerke(glm_model_splag_H)
vif_H <- vif(glm_model_splag_H)

##Now calculate the Poisson model with robust standard errors - Hispanic
cov.pois_H <- vcovHC(glm_model_splag_H, type="HC0")
std.err_H <- sqrt(diag(cov.pois_H))
r.est_H <- cbind(Estimate= coef(glm_model_splag_H), "Robust SE" = std.err_H,
"Pr(>|z|)" = 2 * pnorm(abs(coef(glm_model_splag_H)/std.err_H), lower.tail=FALSE),
LL = coef(glm_model_splag_H) - 1.96 * std.err_H,
UL = coef(glm_model_splag_H) + 1.96 * std.err_H)
r.est_H
ex_re_H <- exp(r.est_H)





###Step 3. Now do the same for NON HISPANIC WHITE

#non-Hispanic White dataset #THE SPECIFIC NAME OF THE DATASET BELOW MIGHT BE DIFFERENT FROM WHAT YOU HAVE
fc_W <- dat[ which(dat$RE=='White' ), ]

#Assign coordinates
coords_W<-coordinates(fc_W)

#Create the nb object with four neighbors
fc_kn4_W <- knn2nb(knearneigh(coords_W, k = 4))

#Make the nb object symmetrical
fc_kn4_W_s <- make.sym.nb(fc_kn4_W)

#Create the listw object from the nb object (fc_kn4) for four closest neighbors
fc_kn4_W_s_w <- nb2listw(fc_kn4_W_s, zero.policy = TRUE)

#Run the ME model to calculate eigenvectors for non-Hispanic White
model_splag_W <- ME(Hospitalized ~ AgeAtDx_scaled + Sex + BMI_scaled + Tobacco + Diabetes + Hypertension + KidneyDisease + LungDisease + apartment + PM25_scaled + ovrcrwd_scaled + p_multif_scaled + per_hb_scaled + per_essen_scaled + per_tew_scaled + resDensK_scaled + bikeScore_int_scaled + walkScore_int_scaled + tranScore_int_scaled + hwy_dist_quarter + ndvi_scaled + pk_ovrlpK_scaled + park_dist_half, data = fc_W, family="poisson", listw=fc_kn4_W_s_w)

#Now calculate the model's coefficients for non-Hispanic White - Poisson
glm_model_splag_W <- glm(Hospitalized ~ AgeAtDx_scaled + Sex + BMI_scaled + Tobacco + Diabetes + Hypertension + KidneyDisease + LungDisease + apartment + PM25_scaled + ovrcrwd_scaled + p_multif_scaled + per_hb_scaled + per_essen_scaled + per_tew_scaled + resDensK_scaled + bikeScore_int_scaled + walkScore_int_scaled + tranScore_int_scaled + hwy_dist_quarter + ndvi_scaled + pk_ovrlpK_scaled + park_dist_half + fitted(model_splag_W), data = fc_W, family="poisson")
fsum_W <- summary(glm_model_splag_W)
sum_W <- tidy(glm_model_splag_W)
nk_W <- nagelkerke(glm_model_splag_W)
vif_W <- vif(glm_model_splag_W)

##Now calculate the Poisson model with robust standard errors - non-Hispanic White
cov.pois_W <- vcovHC(glm_model_splag_W, type="HC0")
std.err_W <- sqrt(diag(cov.pois_W))
r.est_W <- cbind(Estimate= coef(glm_model_splag_W), "Robust SE" = std.err_W,
"Pr(>|z|)" = 2 * pnorm(abs(coef(glm_model_splag_W)/std.err_W), lower.tail=FALSE),
LL = coef(glm_model_splag_W) - 1.96 * std.err_W,
UL = coef(glm_model_splag_W) + 1.96 * std.err_W)
r.est_W
ex_re_W <- exp(r.est_W)

save.image(file = "CODEN_Poisson.RData")
