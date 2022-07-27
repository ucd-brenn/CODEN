# CODEN
The neighborhood built environment and COVID-19 hospitalizations - Statistical Analysis

## Scripts
1. CODEN_FULL.rmd - Full analysis of the entire dataset
    + _note individual RE models followed the same process, but with the dataset filtered using_
    + `w_fc_C = filter(fc_C, RE == 'White')` and
    + `h_fc_C = filter(fc_C, RE == 'Hispanic')` respectively
2. CODEN_DEM.rmd - Analysis of just the impact of individual factors on COVID-19 hospitalization, without neighborhood built environment effects
