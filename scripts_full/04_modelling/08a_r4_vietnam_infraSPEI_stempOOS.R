
# ====================================================================================================
# ============================== OUT OF SAMPLE TEST FOR INTERACTIONS =================================
# ====================================================================================================


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
projname = "infraspei_stempOOS"
region = "all"
region2 = c("South", "South Central Coast", "Central Highlands")
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
    date,
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
    water_piped_year,
    flushtoilet_indoor_year,
    flushtoilet_any_year, 
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
    date,
    cases, 
    logpop,
    logpopdens,
    polyid,
    yearx,
    provincex, 
    tmean_coolestmonth_s,
    tmean_1m, tmean_1m_g, 
    urban_pw,
    urban_s,
    urbanexp10_log,
    gravityf_log,
    traffic_kmperinhab_log,
    water_piped_year,
    flushtoilet_any_year,
    water_g,
    flushany_g,
    spei1_1m, spei1_1m_g,
    spei6_5m, spei6_5m_g
  )

# 
ddf$region3x = as.integer(as.factor(ddf$region3))

# three categories (low, normal, high) of piped water and urbanicity
cats = ddf %>%
  dplyr::select(areaid, year, water_piped_year, urban_pw) %>%
  distinct()

# terciles of water
#tercs = quantile(aa$water_piped_year, c(0.33, 0.67))
tercs = c(0.25, 0.75)
cats$water_c = 1
cats$water_c[ cats$water_piped_year >= tercs[1] ] = 2
cats$water_c[ cats$water_piped_year >= tercs[2] ] = 3
table(cats$water_c)

# urban >= 0.25 or 0.75
cats$urban_c = 1
cats$urban_c[ cats$urban_pw >= 0.25 ] = 2
cats$urban_c[ cats$urban_pw >= 0.75 ] = 3

# combine
ddf = ddf %>%
  dplyr::left_join(cats[ , c(1, 2, 5, 6)])




# ================== set up 5fold exclusion via different regimes ===============

# specify type (random observations, district-by-year, district)
#type = "kfoldOOSq"
type = "stempOOS"

# unique identifier
unique_id = sample(1:10^6, 1)
mod_names = paste(unique_id, type, sep="_")

# by quarters
if(type == "kfoldOOSq"){
  
  # 3-month groupings
  my = ddf %>%
    dplyr::select(date, month) %>%
    distinct() %>%
    dplyr::arrange(date) %>%
    dplyr::left_join(
      data.frame(month = 1:12, subyear=rep(1:4, each=3))
    ) %>%
    dplyr::select(-month)
  ddf = left_join(ddf, my)
  
  # 5-fold quarterly groupings per district
  kfd = ddf %>%
    dplyr::select(areaid, year, subyear) %>%
    distinct() %>%
    dplyr::group_by(areaid) %>%
    dplyr::mutate(kfold = kfold_func(subyear, k=5))
  ddf = left_join(ddf, kfd)
  
  write.csv(
    ddf %>% dplyr::select(areaid, date, kfold),
    paste(save_dir, mod_names, "_folds.csv", sep=""),
    row.names=FALSE
  )
}

if(type == "stempOOS"){
  folds = ddf %>% dplyr::select(areaid, year) %>% distinct()
  folds$kfold = kfold_func(folds, k = 5)
  write.csv(folds, paste(save_dir, mod_names, "_folds.csv", sep=""), row.names=FALSE)
  ddf = left_join(ddf, folds)
}




# =================== create dataframe of models to fit and compare =====================

form_base = paste(
  c("y ~ 1",
    "offset(logpop)",
    "f(month, model='rw1', hyper=hyper2.rw, cyclic=TRUE, scale.model=TRUE, constr=TRUE, replicate=provincex)",
    "f(polyid, model='bym2', graph=nbmatrix_name, replicate=yearx, scale.model=TRUE, constr=TRUE, adjust.for.con.comp=TRUE, hyper=hyper.bym2)"),
  collapse = " + "
)

# full model effects 
effect_names = 
  c("tmean_coolestmonth_s", 
    "gravityf_log", 
    "traffic_kmperinhab_log", 
    "urbanexp10_log",
    "f(flushany_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
    "f(water_c, model='rw1', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
    "f(urban_c, model='rw1', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
    "f(tmean_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
    "f(spei1_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
    "f(spei6_5m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)")

# model list
fx = vector("list", length=length(effect_names)+1)

# models
m1 = paste(effect_names, collapse=" + "); name1 = "full"
m1.1 = paste(effect_names[ -grep("spei6", effect_names)], collapse=" + "); name1.1 = "no_spei6"
m1.2 = paste(effect_names[ -grep("spei1", effect_names)], collapse=" + "); name1.2 = "no_spei1"

# long-term drought
ef2 = replace(effect_names, grepl("spei6", effect_names), "f(spei6_5m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE, replicate=water_c)")
ef3 = replace(effect_names, grepl("spei6", effect_names), "f(spei6_5m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE, replicate=urban_c)")

# short-term drought
ef5 = replace(effect_names, grepl("spei1", effect_names), "f(spei1_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE, replicate=water_c)")
ef6 = replace(effect_names, grepl("spei1", effect_names), "f(spei1_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE, replicate=urban_c)")

# combine into models
m2 = paste(ef2, collapse=" + "); name2 = "spei6water"
m3 = paste(ef3, collapse=" + "); name3 = "spei6urban"
m5 = paste(ef5, collapse=" + "); name5 = "spei1water"
m6 = paste(ef6, collapse=" + "); name6 = "spei1urban"

# unlist and add into df
fx = c(m1, m1.1, m1.2, m2, m3, m5, m6)

# create data frame
fx = data.frame(modid = 1:length(fx),
                fx = fx,
                candidate = c(name1, name1.1, name1.2, name2, name3, name5, name6),
                formula = paste(form_base, fx, sep=" + "))

bs = data.frame(modid = "baseline", fx = "baseline", candidate="baseline", formula=form_base)
fx = rbind(fx, bs)

# model name
fx$model_filename = paste(mod_names, "_infraspei_nb_model_", fx$modid, sep="")



# ================== fit models in an iterative loop ======================

# run model selection loop
for(i in 1:nrow(fx)){
  
  # formula
  fx_i = fx[ i, ]
  form_i = formula(as.vector(fx_i$formula))
  
  # fit INLA model nested in tryCatch: time out after 2 hours (7200 secs)
  e = simpleError("error fitting")
  
  # for each fold
  for(k in 1:5){
    
    # data and set group k to NA
    dd_i <- ddf %>% dplyr::mutate(y = cases) %>%
      dplyr::mutate(y = replace(y, kfold == k, NA))
    
    print("")
    print("==========================================")
    print(Sys.time())
    print("==========================================")
    print("")
    
    # fit 
    mod_i = tryCatch(
      fitINLAModel(form_i, dd_i, family="nbinomial", verbose=TRUE),
      error = function(e) return(e)
    )
    
    # if at first you don't succeed
    if(class(mod_i)[1] == "simpleError"){
      mod_i = tryCatch(
        fitINLAModel(form_i, dd_i, family="nbinomial", verbose=TRUE),
        error = function(e) return(e)
      )
    }
    
    # write timeout to result
    if(class(mod_i)[1] == "simpleError"){
      
      ex = fx_i; ex$result = "error in fitting"
      ex_file_name = paste(k, "_", mod_names, "_nb_err_", fx_i$modid, ".csv", sep="")
      write.csv(ex, paste(save_dir, "errors/", ex_file_name, sep=""), row.names=FALSE)
      
      # otherwise calculate and save fit metrics and model
    } else{
      
      fm = fitMetricsINLA(mod_i, data=dd_i, modname=fx_i$fx, inla.mode="experimental")
      res_i = cbind(fx_i, fm)
      fm_file_name = paste(k, "_", mod_names, "_nb_fitmetrics_", fx_i$modid, ".csv", sep="")
      write.csv(res_i, paste(save_dir, "fitmetrics/", fm_file_name, sep=""), row.names=FALSE)
      
      # save model
      #save(mod_i, file=paste(save_dir, "models/", paste(k, fx_i$model_filename, sep="_"), sep=""))
    }
    
    # extract observed and fitted and save
    dd_o = ddf %>%
      dplyr::select(areaid, province, date, logpop, cases, kfold) %>%
      dplyr::mutate(oos = ifelse(kfold == k, TRUE, FALSE),
                    holdout_type = type,
                    holdout_id = unique_id,
                    model = fx_i$candidate,
                    mean = mod_i$summary.linear.predictor$mean,
                    lower = mod_i$summary.linear.predictor$`0.025quant`,
                    upper = mod_i$summary.linear.predictor$`0.975quant`) %>%
      dplyr::filter(oos == TRUE)
    write.csv(dd_o, paste(save_dir, "models/", paste(k, "output", fx_i$model_filename, ".csv", sep="_"), sep=""), row.names=FALSE)
    
  } # end of kfold loop
  
} # end of model fitting loop


