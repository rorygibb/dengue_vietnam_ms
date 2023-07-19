

# ====================================================================================================
# ===================== SOCIOECOLOGICAL VARIABLES: BASELINE + UNIVARIATE COVARIATES ==================
# ====================================================================================================

## Modelling script 1:
## Fits baseline model and baseline + single covariate ('univariate') models


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
projname = "sociounivar"
region = "all"
region2 = NA
n_clim_bins = 40
plot_graph = FALSE
province_case_threshold = NA
source("00_build_model_df_clim.R")

#write.csv(dd, "./output/dataset_processed.csv", row.names=FALSE)
#dd = read.csv("./output/dataset_processed.csv")



# ================= customise dataframe for this modelling task ====================

# create dataframe for modelling, response and family
ddf = dd

# tmean_g for all models
ddf$tmean_g = ddf$tmean_1m_g

# select only variables required for models and scale/log transform
ddf <- ddf %>%
  dplyr::select(
    province, 
    region3,
    areaid, 
    district,
    year,
    cases, 
    popdens_census,
    logpop,
    logpopdens,
    polyid,
    month,
    yearx,
    provincex, 
    tmean_g,
    tmean_annualmean,
    tmean_annualmean_g,
    tmin_annualmean, 
    tmin_annualmean_g,
    tmin_coolestmonth,
    tmin_coolestmonth_g,
    tmean_coolestmonth,
    tmean_coolestmonth_g,
    pop_propurban_census,
    urban_pw,
    urbanexp_3yr_d,
    urbanexp_10yr_d,
    urbanexp_meanannual,
    flux_grav1,
    flux_rad,
    water_piped_year,
    flushtoilet_indoor_year,
    flushtoilet_any_year,
    traffic_thouskmperinhab
  ) %>%
  dplyr::mutate(
    water_s = scale(water_piped_year),
    flushindoor_s = scale(flushtoilet_indoor_year),
    flushany_s = scale(flushtoilet_any_year),
    gravityf_s = scale(flux_grav1),
    radiationf_s = scale(flux_rad),
    urban_census_s = scale(pop_propurban_census), 
    urban_s = scale(urban_pw),
    urbanexp3_s = scale(urbanexp_3yr_d),
    urbanexp10_s = scale(urbanexp_10yr_d),
    traffic_thouskmperinhab_s = scale(traffic_thouskmperinhab),
    gravityf_log = log(flux_grav1 + 1),
    radiationf_log = log(flux_rad + 1),
    urbanexp3_log = log(urbanexp_3yr_d + 1),
    urbanexp10_log = log(urbanexp_10yr_d + 1),
    urbanexp_total = urbanexp_meanannual * 23,
    urbanexp_total_log = log(urbanexp_total + 1),
    urbanexp_meanannual_s = scale(urbanexp_meanannual),
    traffic_kmperinhab_log = log(traffic_thouskmperinhab * 1000),
    water_g = inla.group(water_piped_year, n=40),
    flushany_g = inla.group(flushtoilet_any_year, n=40),
    flushindoor_g = inla.group(flushtoilet_indoor_year, n=40),
    tmin_coolestmonth_s = scale(tmin_coolestmonth),
    tmin_annualmean_s = scale(tmin_annualmean),
    tmean_annualmean_s = scale(tmean_annualmean),
    tmean_coolestmonth_s = scale(tmean_coolestmonth)
  ) %>%
  dplyr::rename("popdens"=popdens_census) %>%
  dplyr::mutate(popdens_s = scale(popdens))

# # static mean variables for gravity/urban/pop/mobility
# static_vars = dd %>%
#   dplyr::group_by(areaid) %>%
#   dplyr::summarise(logpopdens_static = log(mean(popdens_census, na.rm=TRUE)),
#                    gravityf_log_static = log(mean(flux_grav1, na.rm=TRUE)+1),
#                    radiationf_log_static = log(mean(flux_rad, na.rm=TRUE)+1),
#                    urban_census_s_static = mean(pop_propurban_census),
#                    urban_s_static = mean(urban_pw, na.rm=TRUE),
#                    traffic_kmperinhab_log_static = log(mean(traffic_thouskmperinhab * 1000, na.rm=TRUE)),
#                    traffic_thouskmperinhab_static = mean(traffic_thouskmperinhab, na.rm=TRUE),
#                    urbanexp10_log_static = log(mean(urbanexp_10yr_d)+1),
#                    urbanexp3_log_static = log(mean(urbanexp_3yr_d)+1)) %>%
#   dplyr::mutate(urban_s_static = scale(urban_s_static),
#                 urban_census_s_static = scale(urban_census_s_static),
#                 traffic_thouskmperinhab_s_static = scale(traffic_thouskmperinhab_static))
# ddf = left_join(ddf, static_vars)

# grouped variables
ddf = ddf %>%
  dplyr::mutate(logpopdens_g = inla.group(logpopdens, n=n_clim_bins),
                urban_g = inla.group(urban_pw, n=n_clim_bins),
                urbancensus_g = inla.group(pop_propurban_census, n=n_clim_bins), 
                gravity_g = inla.group(flux_grav1, n=n_clim_bins),
                radiation_g = inla.group(flux_rad, n=n_clim_bins),
                traffic_perinhab_g = inla.group(traffic_thouskmperinhab, n=n_clim_bins))

# extra region flags
ddf$regiony1 = ddf$regiony

# other flags
ddf$areaidx = as.integer(as.factor(ddf$areaid))
ddf$polyidx = ddf$polyid



# # =================== create dataframe of models to fit and compare =====================

# separate bym effects for each year
form_base = paste(
  c("y ~ 1",
    "offset(logpop)",
    "f(month, model='rw1', hyper=hyper2.rw, cyclic=TRUE, scale.model=TRUE, constr=TRUE, replicate=provincex)",
    "f(polyid, model='bym2', graph=nbmatrix_name, replicate=yearx, scale.model=TRUE, constr=TRUE, adjust.for.con.comp=TRUE, hyper=hyper.bym2)"),
  collapse = " + "
)

# univariate fx
fx = c("popdens_s",
       "logpopdens",
       #"logpopdens_static",
       "urban_s",
       "urban_s_static",
       "urban_census_s",
       #"urban_census_s_static",
       "gravityf_s",
       "gravityf_log",
       #"gravityf_log_static",
       "radiationf_s",
       "radiationf_log",
       #"radiationf_log_static",
       "traffic_thouskmperinhab_s",
       "traffic_thouskmperinhab_s_static",
       #"traffic_kmperinhab_log",
       #"traffic_kmperinhab_log_static",
       "urbanexp3_s",
       "urbanexp10_s",
       "urbanexp3_log",
       "urbanexp10_log",
       #"urbanexp10_log_static",
       "water_s",
       #"water_s_static",
       "flushany_s", 
       #"flushany_s_static",
       "flushindoor_s", 
       #"flushindoor_s_static",
       "tmean_annualmean_s",
       "tmin_annualmean_s",
       "tmin_coolestmonth_s",
       "tmean_coolestmonth_s",
       "f(water_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(flushindoor_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(flushany_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(urbancensus_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(urban_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(logpopdens_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(gravity_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(radiation_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(traffic_perinhab_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmean_annualmean_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmin_annualmean_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmin_coolestmonth_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
       "f(tmean_coolestmonth_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)")

# create data frame
fx = data.frame(modid = 1:length(fx),
                fx = fx,
                effect_type = c(rep("linear", 22), rep("rw2", 13)),
                formula = paste(form_base, fx, sep=" + "))
bs = data.frame(modid = "baseline", fx = "baseline", effect_type = "baseline", formula=form_base)
fx = rbind(fx, bs)

# model name
fx$model_filename = paste("sociouni_nb_model_", fx$modid, ".R", sep="")



# ================== fit models in an interative loop ======================
  
# run model selection loop
for(i in 1:nrow(fx)){

  # formula
  fx_i = fx[ i, ]
  form_i = formula(as.vector(fx_i$formula))

  # fit INLA model nested in tryCatch
  e = simpleError("error fitting")
  
  # save storage by cutting all unnecessary variables and only keep specified variable
  if(fx_i$effect_type %in% c("rw2" , "baseline")){ 
    
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
        tmean_g,
        water_g,
        tmean_annualmean_g,
        tmin_annualmean_g,
        tmin_coolestmonth_g,
        tmean_coolestmonth_g,
        flushindoor_g, 
        flushany_g,
        urban_g,
        urbancensus_g,
        logpopdens_g,
        gravity_g,
        radiation_g,
        traffic_perinhab_g) %>%
      dplyr::mutate(y = cases)
    
  } else{
    
    dd_i <- ddf %>%
      dplyr::select(
        province, 
        areaid, 
        areaidx,
        district,
        year,
        cases, 
        tmin_coolestmonth_s,
        logpop,
        logpopdens,
        polyid,
        polyidx,
        month,
        yearx,
        provincex, 
        tmean_g) %>%
      dplyr::mutate(y = cases)
    dd_i <- cbind(dd_i, ddf[ , which(names(ddf) %in% unlist(strsplit(as.vector(fx_i$fx), split="[ + ]"))), drop=FALSE])
    
    }
  
  # fit model
  mod_i = tryCatch(
    fitINLAModel(form_i, dd_i, family="nbinomial", verbose=TRUE),
    error = function(e) return(e)
  )

  # write timeout to result
  if(class(mod_i)[1] == "simpleError"){

    ex = fx_i; ex$result = "error in fitting"
    ex_file_name = paste("sociouni_nb_err_", fx_i$modid, ".csv", sep="")
    write.csv(ex, paste(save_dir, "errors/", ex_file_name, sep=""), row.names=FALSE)

    # otherwise calculate and save fit metrics and model
  } else{

    fm = fitMetricsINLA(mod_i, data=dd_i, modname=fx_i$fx, inla.mode="experimental")
    res_i = cbind(fx_i, fm)
    fm_file_name = paste("sociouni_nb_fitmetrics_", fx_i$modid, ".csv", sep="")
    write.csv(res_i, paste(save_dir, "fitmetrics/", fm_file_name, sep=""), row.names=FALSE)

    # save model
    save(mod_i, file=paste(save_dir, "models/", fx_i$model_filename, sep=""))
  }

} # end of model fitting loop

