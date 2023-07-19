# ====================================================================================================
# ================================= FIT FULL MODEL ===================================================
# ====================================================================================================


## ============= set up workflow and data ===================

# project root, dependencies, plot themes
setwd("/lustre/scratch/scratch/ucbtgib/dengue_viet_23")
#setwd("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/myriad_may2023/")

# inla install testing version and centos 7
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
# inla.binary.install()

library(dplyr); library(raster); library(rgdal); library(sf)
library(stringr); library(ggplot2); library(lubridate)
library(magrittr); library(INLA); library(spdep)
source("00_plot_themes.R")

# INLA modelling functions, priors and pardiso license
source("00_inla_setup_functions_r4.R")
inla.setOption(pardiso.license="pardiso.lic.txt")

# build dataframe calls 21_build_model_df.R, specifying 4 variables 
# projname: name of save directory for outputs
# region: either NA, north, south or central; whether to subset to specific region
# n_clim_bins: n bins for grouping climatic predictors for nonlinear effects
# plot_graph: visualise neighbourhood matrix?
projname = "fullmodel_northern"
region = "all"
region2 = c("Red River Delta", "North Central", "North/Northeast", "Northwest")
n_clim_bins = 40
plot_graph = FALSE
province_case_threshold = NA
source("00_build_model_df_clim.R")




# ================= customise dataframe for this modelling task ====================

# create dataframe for modelling, response and family
ddf = dd

# select only variables required for models and scale/log transform
ddf <- ddf %>%
  dplyr::select(
    province, 
    region2,
    region3,
    areaid, 
    district,
    year,
    month,
    cases, 
    popdens_census,
    logpop,
    logpopdens,
    polyid,
    yearx,
    provincex, 
    urban_pw,
    flux_grav1,
    urbanexp_3yr_d,
    urbanexp_10yr_d,
    urbanexp_meanannual,
    water_piped_well, 
    water_piped_year,
    sanitation_flushtoilet_indoor,
    flushtoilet_indoor_year,
    flushtoilet_any_year, 
    traffic_milperskm,
    traffic_thouskmperinhab,
    tmean_coolestmonth,
    tmean_1m, tmean_1m_g, 
    spei1_1m, spei1_1m_g,
    spei6_5m, spei6_5m_g
  ) %>%
  dplyr::mutate(
    urban_s = scale(urban_pw),
    gravityf_log = scale( log(flux_grav1 + 1) ),
    urbanexp10_log = scale( log(urbanexp_10yr_d + 1) ),
    traffic_kmperinhab_log = scale(log(traffic_thouskmperinhab * 1000)),
    water_g = inla.group(water_piped_year, n=40),
    flushany_g = inla.group(flushtoilet_any_year, n=40),
    tmean_coolestmonth_s = scale(tmean_coolestmonth)
  )

# extra region flags
ddf$regiony1 = ddf$regiony

# subset to required variables
ddf = ddf %>%
  dplyr::select(
    province, 
    region2,
    region3,
    areaid, 
    district,
    year,
    month,
    cases, 
    logpop,
    logpopdens,
    polyid,
    yearx,
    provincex, 
    tmean_coolestmonth_s,
    tmean_1m, tmean_1m_g, 
    spei1_1m, spei1_1m_g,
    spei6_5m, spei6_5m_g,
    urban_s,
    urbanexp3_log,
    urbanexp10_log,
    gravityf_log,
    traffic_kmperinhab_log,
    water_piped_year,
    flushtoilet_any_year,
    water_g,
    flushany_g
  )

# 
ddf$region3x = as.integer(as.factor(ddf$region3))



# =================== create dataframe of models to fit and compare =====================

# with a province level fixed effect to account for confounding effects of space
form_base = paste(
  c("y ~ 1",
    "offset(logpop)",
    "f(month, model='rw1', hyper=hyper2.rw, cyclic=TRUE, scale.model=TRUE, constr=TRUE, replicate=provincex)",
    "f(polyid, model='bym2', graph=nbmatrix_name, replicate=yearx, scale.model=TRUE, constr=TRUE, adjust.for.con.comp=TRUE, hyper=hyper.bym2)"),
  collapse = " + "
)

# full model 
effect_names = 
  c("tmean_coolestmonth_s", 
    "gravityf_log",
    "urban_s",
    "urbanexp10_log",
    "traffic_kmperinhab_log",
    "f(flushany_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
    "f(water_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
    "f(tmean_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
    "f(spei1_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
    "f(spei6_5m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)")
m1 = paste(effect_names, collapse=" + "); name1 = "full"
fx = m1

# create data frame
fx = data.frame(modid = 1:length(fx),
                fx = fx,
                candidate = c(name1),
                formula = paste(form_base, fx, sep=" + "))
bs = data.frame(modid = "baseline", fx = "baseline", candidate="baseline", formula=form_base)
fx = rbind(fx, bs)

# model name
fx$model_filename = paste("fullmod_nb_model_", fx$modid, ".R", sep="")



# ================== fit models in an iterative loop ======================

# run model selection loop
for(i in 1:nrow(fx)){

  # formula
  fx_i = fx[ i, ]
  form_i = formula(as.vector(fx_i$formula))
  
  # fit INLA model nested in tryCatch: time out after 2 hours (7200 secs)
  e = simpleError("error fitting")
  
  # data
  dd_i <- ddf %>% dplyr::mutate(y = cases)
  
  # fit model
  print("")
  print("==========================================")
  print(Sys.time())
  print("==========================================")
  print("")
  
  # set config to TRUE for posterior samples
  mod_i = tryCatch(
    fitINLAModel(form_i, dd_i, family="nbinomial", verbose=TRUE, config=TRUE, return.marginals=FALSE),
    error = function(e) return(e)
  )
  
  # # improve hyperpar estimates
  # mod_j = tryCatch(
  #   INLA::inla.hyperpar(mod_i, verbose=TRUE),
  #   error = function(e) return(e)
  #   )
  # if(class(mod_j)[1] != "simpleError"){
  #   mod_i = mod_j
  # }

  # write timeout to result
  if(class(mod_i)[1] == "simpleError"){
    
    ex = fx_i; ex$result = "error in fitting"
    ex_file_name = paste("fullmod_nb_err_", fx_i$modid, ".csv", sep="")
    write.csv(ex, paste(save_dir, "errors/", ex_file_name, sep=""), row.names=FALSE)
    
    # otherwise calculate and save fit metrics and model
  } else{
    
    fm = fitMetricsINLA(mod_i, data=dd_i, modname=fx_i$fx, inla.mode="experimental")
    res_i = cbind(fx_i, fm)
    fm_file_name = paste("fullmod_nb_fitmetrics_", fx_i$modid, ".csv", sep="")
    write.csv(res_i, paste(save_dir, "fitmetrics/", fm_file_name, sep=""), row.names=FALSE)
    
    # save model
    save(mod_i, file=paste(save_dir, "models/", fx_i$model_filename, sep=""))
  }
  
} # end of model fitting loop
