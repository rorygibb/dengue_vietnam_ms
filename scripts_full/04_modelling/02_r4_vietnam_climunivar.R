

# ====================================================================================================
# ====================== CLIMATIC VARIABLES: BASELINE + UNIVARIATE COVARIATES ========================
# ====================================================================================================

## Modelling script 1:
## Fits baseline model and baseline + single covariate ('univariate') models


## ============= set up workflow and data ===================

# project root, dependencies, plot themes
setwd("/lustre/scratch/scratch/ucbtgib/dengue_viet_23")
#setwd("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/myriad_sep2021/")

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
projname = "climunivar"
region = "all"
region2 = NA
n_clim_bins = 40
plot_graph = FALSE
province_case_threshold = NA
source("00_build_model_df_clim.R")

#dd = read.csv("./output/dataset_processed.csv")



# ================= customise dataframe for this modelling task ====================

# create dataframe for modelling, response and family
ddf = dd

# subset to key variables
# select only variables required for models and scale/log transform
ddf1 <- ddf %>%
  dplyr::select(
    province, 
    areaid, 
    district,
    year,
    cases, 
    logpop,
    logpopdens,
    polyid,
    month,
    yearx,
    provincex, 
    tmean_annualmean,
    tmin_annualmean,
    tmin_coolestmonth)
ddf2 <- ddf %>%
  dplyr::select(matches("tmean|tmin|tmax|tdrange|precip|spei")) %>%
  dplyr::select(matches("_g"))
ddf <- cbind(ddf1, ddf2)

# extra region flags
ddf$regiony1 = ddf$regiony

# other flags
ddf$areaidx = as.integer(as.factor(ddf$areaid))
ddf$polyidx = ddf$polyid

# add
ddf$tmin_coolestmonth_s = scale(ddf$tmin_coolestmonth)



# # =================== create dataframe of models to fit and compare =====================

# separate bym effects for each year (each year is independent replicate but there is spatial dependency within each year)
form_base = paste(
  c("y ~ 1",
    "offset(logpop)",
    "f(month, model='rw1', hyper=hyper2.rw, cyclic=TRUE, scale.model=TRUE, constr=TRUE, replicate=provincex)",
    "f(polyid, model='bym2', graph=nbmatrix_name, replicate=yearx, scale.model=TRUE, constr=TRUE, adjust.for.con.comp=TRUE, hyper=hyper.bym2)"), 
  collapse = " + "
)

# univariate climate fx
fx = c("f(tmean_0m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmean_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmean_2m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmean_3m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmean_4m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmean_5m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmean_6m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",

       "f(tmin_0m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmin_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmin_2m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmin_3m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmin_4m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmin_5m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmin_6m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       
       "f(tmax_0m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmax_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmax_2m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmax_3m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmax_4m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmax_5m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmax_6m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       
       "f(spei1_0m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei1_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei1_2m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei1_3m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei1_4m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei1_5m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei1_6m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",

       "f(spei6_0m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei6_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei6_2m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei6_3m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei6_4m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei6_5m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei6_6m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       
       "f(spei12_0m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei12_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei12_2m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei12_3m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei12_4m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei12_5m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(spei12_6m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       
       "f(precip_0m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(precip_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(precip_2m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(precip_3m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(precip_4m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(precip_5m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(precip_6m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)")

# covar names
covs = paste(rep(c("tmean", "tmin", "tmax", "spei1", "spei6", "spei12", "precip"), each=7), rep(0:6, 7), sep="_")

# create data frame
fx = data.frame(modid = 1:length(fx),
                covar = covs,
                fx = fx,
                formula = paste(form_base, fx, sep=" + "))
bs = data.frame(modid = "baseline", covar = "baseline", fx = "baseline", formula=form_base)
fx = rbind(fx, bs)

# model names
fx$model_filename = paste("climuni_nb_model_", fx$modid, ".R", sep="")



# ================== fit models in an iterative loop ======================
  
# run model selection loop
for(i in 1:nrow(fx)){

  # formula
  fx_i = fx[ i, ]
  form_i = formula(as.vector(fx_i$formula))

  # fit INLA model nested in tryCatch: time out after 2 hours (7200 secs)
  e = simpleError("error fitting")
  
  # save storage by cutting all unnecessary variables and only keep specified variable
  dd_i <- ddf %>%
    dplyr::select(
      province,
      areaid,
      areaidx,
      district,
      year,
      cases,
      logpop,
      logpopdens,
      polyid,
      polyidx,
      month,
      yearx,
      provincex,
      tmin_coolestmonth_s) %>%
    dplyr::mutate(y = cases) %>%
    cbind(ddf[ , grep(as.vector(fx_i$covar), names(ddf)), drop=FALSE])
  
  # fit model
  mod_i = tryCatch(
    fitINLAModel(form_i, dd_i, family="nbinomial", verbose=TRUE),
    error = function(e) return(e)
  )

  # write timeout to result
  if(class(mod_i)[1] == "simpleError"){

    ex = fx_i; ex$result = "error in fitting"
    ex_file_name = paste("climuni_nb_err_", fx_i$modid, ".csv", sep="")
    write.csv(ex, paste(save_dir, "errors/", ex_file_name, sep=""), row.names=FALSE)

    # otherwise calculate and save fit metrics and model
  } else{

    fm = fitMetricsINLA(mod_i, data=dd_i, modname=fx_i$fx, inla.mode="experimental")
    res_i = cbind(fx_i, fm)
    fm_file_name = paste("climuni_nb_fitmetrics_", fx_i$modid, ".csv", sep="")
    write.csv(res_i, paste(save_dir, "fitmetrics/", fm_file_name, sep=""), row.names=FALSE)

    # save model
    save(mod_i, file=paste(save_dir, "models/", fx_i$model_filename, sep=""))
  }

} # end of model fitting loop

