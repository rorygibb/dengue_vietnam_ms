# ====================================================================================================
# ============================ CROSS VALIDATION MODEL SCRIPT =========================================
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
projname = "tempOOS_1"
region = "all"
region2 = NA
n_clim_bins = 40
plot_graph = FALSE
province_case_threshold = NA
source("00_build_model_df_clim.R")

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
    spei1_1m, spei1_1m_g,
    spei6_5m, spei6_5m_g,
    urban_s,
    urbanexp10_log,
    gravityf_log,
    traffic_kmperinhab_log,
    water_piped_year,
    flushtoilet_any_year,
    water_g,
    flushany_g
  )




# ================= initialises 5fold exclusion =================

# n.b. the below only runs once, the first time this script is run for a specified save_dir
# creates candidate models list, tracker and completed files to track and store model runs

# parameters for setup

n_reps = 12
k_folds = 5 
#type = "spOOS" 
# type = "stempOOS"
type = "tempOOS"

# create initialise files have not all created already
ll = list.files(save_dir)
if(!"models_list_all.csv" %in% ll){
  
  # unique IDs vec
  kfold_ids = c()
  
  # ----- create n_reps k-fold sets ------
  
  for(nn in 1:n_reps){
    
    unique_id_n = sample(1:10^6, 1)
    mod_names_n = paste(unique_id_n, type, sep="_")
    kfold_ids = c(kfold_ids, unique_id_n)
    
    # partition dataset
    if(type == "stempOOS"){
      folds = ddf %>% dplyr::select(areaid, year) %>% distinct()
      folds$kfold = kfold_func(folds, k = k_folds)
      write.csv(folds, paste(save_dir, unique_id_n, "_folds.csv", sep=""), row.names=FALSE)
      #ddf = left_join(ddf, folds)
    }
    
    if(type == "spOOS"){
      folds = ddf %>% dplyr::select(areaid) %>% distinct() 
      folds$kfold = kfold_func(folds, k = k_folds)
      write.csv(folds, paste(save_dir, unique_id_n, "_folds.csv", sep=""), row.names=FALSE)
      #ddf = left_join(ddf, folds)
    }
    
    # by quarters
    if(type == "tempOOS"){
      
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
      
      # 5-fold quarterly groupings across the board
      # kfd = ddf %>%
      #   dplyr::select(year, subyear) %>%
      #   distinct() %>%
      #   dplyr::mutate(kfold = kfold_func(subyear, k=k_folds))
      
      ddf = left_join(ddf, kfd)
      
      write.csv(
        ddf %>% dplyr::select(areaid, date, kfold),
        paste(save_dir, unique_id_n, "_folds.csv", sep=""),
        row.names=FALSE
      )
      
      ddf = ddf %>% dplyr::select(-kfold)
    }
    
  } # end partitioning loop
  
  
  # ----- create models dataframe -----
  
  form_base = paste(
    c("y ~ 1",
      "offset(logpop)",
      "f(month, model='rw1', hyper=hyper2.rw, cyclic=TRUE, scale.model=TRUE, constr=TRUE, replicate=provincex)",
      "f(polyid, model='bym2', graph=nbmatrix_name, replicate=yearx, scale.model=TRUE, constr=TRUE, adjust.for.con.comp=TRUE, hyper=hyper.bym2)"),
    collapse = " + "
  )
  
  effect_names = 
    c("f(tmean_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
      "f(spei1_1m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
      "f(spei6_5m_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
      "tmean_coolestmonth_s", 
      "gravityf_log",
      "urban_s",
      "urbanexp10_log",
      "traffic_kmperinhab_log",
      "f(flushany_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)",
      "f(water_g, model='rw2', hyper=hyper2.rw, scale.model=TRUE, constr=TRUE)")
  
  # model list
  fx = vector("list", length=4)
  
  # models
  m1 = paste(effect_names, collapse=" + "); name1 = "full"
  
  # add models
  fx[1] = list(m1)
  
  # individual holdouts (n.b. quarterly model so only including monthly varying climate vars)
  for(i in 1:3){
    fx[[ i + 1 ]] = paste(effect_names[ -i ], collapse=" + ")
  }
  
  # unlist and add into df
  fx = unlist(fx)
  
  # create data frame including formulae
  fx = data.frame(modid = 1:length(fx),
                  fx = fx,
                  candidate = c(name1, "tmean", "spei1", "spei6"),
                  formula = paste(form_base, fx, sep=" + "))
  
  bs = data.frame(modid = "baseline", fx = "baseline", candidate="baseline", formula=form_base)
  fx = rbind(fx, bs)
  
  # repeat for each n_reps ID to create full set of models to fit, cross referenced to unique holdout set
  fx_mods = data.frame()
  for(ii in kfold_ids){
    fx_mods = rbind(fx_mods, 
                    fx %>% dplyr::mutate(unique_id = ii,
                                         model_filename = paste(ii, type, "nb_model", modid, sep="_")))
  }
  fx_mods$model_identifier = 1:nrow(fx_mods)
  
  
  # ------- save initialisation objects ---------
  
  # list of models to fit
  write.csv(fx_mods, paste(save_dir, "models_list_all.csv", sep=""), row.names=FALSE)
  
  # tracker (which models are currently running)
  tracker = data.frame(unique_id = "dummy", model_filename = "dummy", model_identifier = "dummy")
  write.csv(tracker, paste(save_dir, "models_tracker.csv", sep=""), row.names=FALSE)
  
  # completed
  completed = data.frame(unique_id = "dummy", model_filename = "dummy", model_identifier = "dummy", candidate = "dummy",
                         n_models_fitted = "dummy", mae_oos = "dummy", rmse_oos = "dummy")
  write.csv(completed, paste(save_dir, "models_completed.csv", sep=""), row.names=FALSE)
  
} # end init block





# ===================== chooses and fits model under k-fold CV ============================

# check currently running and completed models
tracker = read.csv(paste(save_dir, "models_tracker.csv", sep=""))
completed = read.csv(paste(save_dir, "models_completed.csv", sep=""))

# list of models to fit
fx = read.csv(paste(save_dir, "models_list_all.csv", sep="")) %>%
  dplyr::filter(!model_identifier %in% c(tracker$model_identifier, completed$model_identifier)) 

# select model and fit if > 0 models left in list
if(nrow(fx) > 0){
  
  # select 1 model
  fx = fx %>% sample_n(1)
  
  # append selected model to tracker (removed after model completed)
  to_append = fx[ , c("unique_id", "model_filename", "model_identifier")]
  write.table(to_append, file=paste(save_dir, "models_tracker.csv", sep=""), 
              append=TRUE, col.names=FALSE, row.names=FALSE, sep=",")
  
  # add correct folds information to dd
  folds = read.csv(paste(save_dir, fx$unique_id, "_folds.csv", sep=""))
  ddf = left_join(ddf, folds)
  
  # fit each of the models in dataframe fx (by default 1)
  for(i in 1:nrow(fx)){
    
    # formula
    fx_i = fx[ i, ]
    form_i = formula(as.vector(fx_i$formula))
    
    # fit INLA model nested in tryCatch
    e = simpleError("error fitting")
    
    # the folds to fit
    folds_seq = unique(ddf$kfold)[ order(unique(ddf$kfold)) ]
    
    # oos_results data frame
    oos_results = data.frame()
    
    # for each fold
    for(k in folds_seq){
      
      # data and set group k to NA
      dd_i <- ddf %>% dplyr::mutate(y = cases) %>%
        dplyr::mutate(y = replace(y, kfold == k, NA))
      
      print("")
      print("==========================================")
      print(Sys.time())
      print("==========================================")
      print("")
      
      # fit model
      mod_i = tryCatch(
        fitINLAModel(form_i, dd_i, family="nbinomial", verbose=TRUE),
        error = function(e) return(e)
      )
      
      # try again if failed
      if(class(mod_i)[1] == "simpleError"){
        mod_i = tryCatch(
          fitINLAModel(form_i, dd_i, family="nbinomial", verbose=TRUE),
          error = function(e) return(e)
        )
      }
      
      # if failed again write timeout to result
      if(class(mod_i)[1] == "simpleError"){
        
        ex = fx_i; ex$result = "error in fitting"
        ex_file_name = paste(k, "_", type, "_nb_err_", fx_i$model_identifier, ".csv", sep="")
        write.csv(ex, paste(save_dir, "errors/", ex_file_name, sep=""), row.names=FALSE)
        
        # otherwise calculate and save fit metrics and model
      } else{
        
        fm = fitMetricsINLA(mod_i, data=dd_i, modname=fx_i$fx, inla.mode="experimental")
        res_i = cbind(fx_i, fm)
        fm_file_name = paste(k, "_", type, "_nb_fitmetrics_", fx_i$model_identifier, ".csv", sep="")
        write.csv(res_i, paste(save_dir, "fitmetrics/", fm_file_name, sep=""), row.names=FALSE)
        
        # save model
        #save(mod_i, file=paste(save_dir, "models/", paste(k, fx_i$model_filename, sep="_"), sep=""))
        
        # extract observed and fitted and save
        dd_o = ddf %>%
          dplyr::select(areaid, date, logpop, cases, kfold) %>%
          dplyr::mutate(oos = ifelse(kfold == k, TRUE, FALSE),
                        model = fx_i$candidate,
                        holdout_id = fx_i$unique_id,
                        model_identifier = fx_i$model_identifier,
                        holdout_type = type,
                        mean = mod_i$summary.linear.predictor$mean,
                        lower = mod_i$summary.linear.predictor$`0.025quant`,
                        upper = mod_i$summary.linear.predictor$`0.975quant`) %>%
          dplyr::filter(oos == TRUE)
        write.csv(dd_o, paste(save_dir, "models/", paste(k, "output", fx_i$model_identifier, fx_i$model_filename, ".csv", sep="_"), sep=""), row.names=FALSE)
        
        # add to growing OOS results dataframe
        oos_results = rbind(oos_results, dd_o)
        
      }
      
    } # end of kfold loop
    
    
    # ============ final operations ==============
    
    print("Saving results")
    
    # create dataframe to add to "completed" csv
    completed_i = data.frame(unique_id = fx_i$unique_id, 
                             model_filename = fx_i$model_filename, 
                             model_identifier = fx_i$model_identifier, 
                             candidate = fx_i$candidate,
                             n_models_fitted = n_distinct(oos_results$kfold), 
                             mae_oos = NA, 
                             rmse_oos = NA)
    
    # calculate predicted and residual error (run through link function)
    oos_results$predicted = exp(oos_results$logpop + oos_results$mean)
    oos_results$resid = oos_results$predicted - oos_results$cases
    
    # calculate summary stats
    completed_i$mae_oos = mean(abs(oos_results$resid), na.rm=TRUE)
    completed_i$rmse_oos = sqrt(mean(oos_results$resid^2, na.rm=TRUE))
    
    # remove from 
    #if(completed_i$n_models_fitted == k_folds){
      
    # append to completed
    write.table(completed_i, file=paste(save_dir, "models_completed.csv", sep=""), 
                append=TRUE, col.names=FALSE, row.names=FALSE, sep=",")
      
    #}
    
    # remove model from tracker
    read.csv(paste(save_dir, "models_tracker.csv", sep="")) %>%
      dplyr::filter(model_identifier != fx_i$model_identifier) %>%
      write.csv(paste(save_dir, "models_tracker.csv", sep=""))
    
  } # end of model fitting loop
  
  
} # end of if statement



# ===================== ENDS =============================
