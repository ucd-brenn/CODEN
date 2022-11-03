###Step 4: Run the model without built environment characteristics
load('CODEN_Poisson.RData')

##Run the ME model to calculate eigenvectors for all cases
model_splag_C_nbe <- ME(Hospitalized ~ AgeAtDx_scaled + Sex + RE + BMI_scaled + Tobacco + Diabetes + Hypertension + KidneyDisease + LungDisease, data = dat, family="poisson", listw=fc_kn4_C_s_w)

##Now calculate the model's coefficients for all cases - Poisson
glm_model_splag_C_nbe <- glm(Hospitalized ~ AgeAtDx_scaled + Sex + RE + BMI_scaled + Tobacco + Diabetes + Hypertension + KidneyDisease + LungDisease + fitted(model_splag_C_nbe), data = dat, family="poisson")
fsum <- summary(glm_model_splag_C_nbe)
sum <- tidy(glm_model_splag_C_nbe)
nk <- nagelkerke(glm_model_splag_C_nbe)
vif <- vif(glm_model_splag_C_nbe)
emmeans(glm_model_splag_C_nbe, list(pairwise ~ RE), adjust = "tukey")

##Now calculate the Poisson model with robust standard errors - all cases
cov.pois_nbe <- vcovHC(glm_model_splag_C_nbe, type="HC0")
std.err_nbe <- sqrt(diag(cov.pois_nbe))
r.est_nbe <- cbind(Estimate= coef(glm_model_splag_C_nbe), "Robust SE" = std.err_nbe,
"Pr(>|z|)" = 2 * pnorm(abs(coef(glm_model_splag_C_nbe)/std.err_nbe), lower.tail=FALSE),
LL = coef(glm_model_splag_C_nbe) - 1.96 * std.err_nbe,
UL = coef(glm_model_splag_C_nbe) + 1.96 * std.err_nbe)
r.est_nbe
ex_re_nbe <- exp(r.est_nbe)

save.image(file = "CODEN_Poisson.RData")
