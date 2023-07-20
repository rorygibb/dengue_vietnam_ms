

# ================== Script containing setup priors and functions for INLA spatiotemporal models ====================


# ----------------- 1. Project hyperpriors -------------------

# iid model 
hyper.iid = list(theta = list(prior="pc.prec", param=c(1, 0.01)))

# ar1 model
hyper.ar1 = list(theta1 = list(prior='pc.prec', param=c(0.5, 0.01)),
                  rho = list(prior='pc.cor0', param = c(0.5, 0.75)))

# bym model
hyper.bym = list(theta1 = list(prior="pc.prec", param=c(1, 0.01)),
                 theta2 = list(prior="pc.prec", param=c(1, 0.01)))

# bym2 model
# probability of SD of theta1 > 1 = 0.01
hyper.bym2 = list(theta1 = list(prior="pc.prec", param=c(1, 0.01)),
                  theta2 = list(prior="pc", param=c(0.5, 0.5)))

# hyperpriors for model grouping (iid / ar1) if used
# group.control.iid = list(model='iid', hyper = list(prec = list(prior='pc.prec',param=c(1, 0.01))))
# group.control.ar1 = list(model='ar1', hyper = list(theta1 = list(prior='pc.prec', param=c(1, 0.01)), rho = list(prior='pc.cor0', param = c(0.5, 0.75))))

# rw1/rw2 model: three levels of constraint on precision parameter 
# (puts more or less prior probability density on more or less wiggly)
hyper1.rw = list(prec = list(prior='pc.prec', param=c(0.1, 0.01))) # strictest smoothing; sd constrained to be low
hyper2.rw = list(prec = list(prior='pc.prec', param=c(0.3, 0.01))) # medium
hyper3.rw = list(prec = list(prior='pc.prec', param=c(1, 0.01))) # weaker (suggested INLA default) 
hyper4.rw = list(prec = list(prior='pc.prec', param=c(2, 0.01))) # weakest; sd can be quite wide 




# ------------------- 2. INLA model fitting function --------------------

### fitINLAModel: fit and return INLA model with specified formula and family

#' @param formx inla formula object; i.e. created using formula(y ~ x + f())
#' @param data model data frame
#' @param family likelihood
#' @param config boolean (default FALSE): set config in compute to TRUE for inla.posterior.sample()
#' @param verbose verbose reporting on or off? default FALSE
#' @param return.marginals specify whether model should save and return marginals

fitINLAModel = function(formx, data, family, config=FALSE, verbose=FALSE, return.marginals=FALSE){
  
  return(
    
    INLA::inla(
      
      # formula, data and model family
      formx,
      data = data,
      family = family,
      
      # fixed effects calibration
      control.fixed = list(mean.intercept=0, 
                           prec.intercept=1, # precision 1
                           mean=0, 
                           prec=1), # weakly regularising on fixed effects (sd of 1)
      
      # save predicted values on response scale
      control.predictor = list(compute=FALSE, 
                               link=1),
      
      # items to compute
      control.compute = list(cpo=TRUE, 
                             waic=TRUE, 
                             dic=TRUE, 
                             config=config, # set config to TRUE if want to do INLA::inla.posterior.sample()
                             return.marginals=return.marginals), # do not return marginals unless specified (saves memory)
      
      # configure inla approx
      control.inla = list(strategy='adaptive', # adaptive gaussian
                          cmin=0), # fixing Q factorisation issue https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls)
      
      # verbose?
      inla.mode = "experimental", # new version of INLA algorithm (requires R 4.1 and INLA testing version)
      verbose = verbose
    )
  )
  
  # # quick viz of prior sd distribution with precision
  # precx = 1
  # hist(rnorm(100000, 0, sd=sqrt(1/precx)), 100)
}




# ------------------ 3. Wrapper functions to extract useful stuff from fitted INLA models -----------------


### fitMetricsINLA: extract and report on metrics of fit

#' @param mod fitted INLA model
#' @param data data used to fit the model; n.b. assumes that response variable column is called "y"
#' @param modname name to give model in resulting dataframe
#' @param inla.mode either "default" or "experimental" (if experimental, needs to add offset to fitted values and exponentiate)

fitMetricsINLA = function(mod, data, modname="mod", inla.mode="default"){
  
  # if(inla.mode == "default"){
  #   
  #   dx = data[ , c("y"), drop=FALSE] %>%
  #     dplyr::mutate(#fitted = mod$summary.fitted.values$mean,
  #                   #abs_err = abs(y - fitted),
  #                   cpo = mod$cpo$cpo,
  #                   cpo_fail = mod$cpo$failure)
  # }
  # 
  # if(inla.mode == "experimental"){
  #   
  #   dx = data[ , c("y", "logpop"), drop=FALSE] %>%
  #     dplyr::mutate(#fitted = exp(mod$summary.fitted.values$mean + logpop), # add offset and exp
  #                   #abs_err = abs(y - fitted),
  #                   cpo = mod$cpo$cpo,
  #                   cpo_fail = mod$cpo$failure)
  # }
  
  fit = data.frame(modname = modname,
                   dic = mod$dic$dic, 
                   waic = mod$waic$waic,
                   waic_neffp = mod$waic$p.eff,
                   #mae = mean(dx$abs_err, na.rm=TRUE),
                   logscore = -mean(log(mod$cpo$cpo), na.rm=TRUE),
                   cpo_fail = sum(mod$cpo$failure == 1 & !is.na(mod$cpo$failure)))
  return(fit)
}


### extractRandomINLA: extract random effect and rename columns

# if effect is grouped/replicated by a factor, automatically assign each subgroup to its grouping factor (labelled 1:n) 
# if BYM model, further partition into u and v components

#' @param summary_random points to model$summary.random$effect_of_interest
#' @param effect_name name to assign to fitted effect in dataframe (can be anything)
#' @param model_is_bym boolean; to specify if model is joint Besag-York-Mollie
#' @param transform specify whether to exponentiate coefficients (i.e. back transform to relative risk)
#' 
extractRandomINLA = function(summary_random, effect_name, model_is_bym=FALSE, transform=FALSE){
  
  # extract model effect
  rf = summary_random %>%
    dplyr::rename("value"=1, "lower"=4, "median"=5, "upper"=6)

  # label by grouping factor (if not replicated, group is 1 for all observations)
  rf$group = rep(1:as.vector(table(rf$value)[1]), each=n_distinct(rf$value))
  
  # partition BYM into u and v components
  if(model_is_bym){
    rf$component = rep(c("uv_joint", "u_besag"), each=n_distinct(rf$value)/2)
    rf$value = rep(1:(n_distinct(rf$value)/2), n_distinct(rf$group)*2)
  }
  
  # back transform if specified
  if(transform == TRUE){
    rf[ , 2:7 ] = exp(rf[ , 2:7])
  }
  
  # name and return
  rf$effect = effect_name
  return(rf)
}


### extractFixedINLA: extract fixed effects and rename columns

extractFixedINLA = function(model, model_name="mod", transform=FALSE){
  ff = model$summary.fixed
  ff$param = row.names(ff)
  ff$param[ ff$param == "(Intercept)" ] = "Intercept"
  names(ff)[3:5] = c("lower", "median", "upper")
  if(transform == TRUE){
    ff[ 1:5 ] = exp(ff[ 1:5 ])
  }
  ff
}


###################################################################################################

# kfold function from dismo package

kfold_func = function (x, k = 5, by = NULL) 
{
  singlefold <- function(obs, k) {
    if (k == 1) {
      return(rep(1, obs))
    }
    else {
      i <- obs/k
      if (i < 1) {
        stop("insufficient records:", obs, ", with k=", 
             k)
      }
      i <- round(c(0, i * 1:(k - 1), obs))
      times = i[-1] - i[-length(i)]
      group <- c()
      for (j in 1:(length(times))) {
        group <- c(group, rep(j, times = times[j]))
      }
      r <- order(runif(obs))
      return(group[r])
    }
  }
  if (is.vector(x)) {
    if (length(x) == 1) {
      if (x > 1) {
        x <- 1:x
      }
    }
    obs <- length(x)
  }
  else if (inherits(x, "Spatial")) {
    if (inherits(x, "SpatialPoints")) {
      obs <- nrow(coordinates(x))
    }
    else {
      obs <- nrow(x@data)
    }
  }
  else {
    obs <- nrow(x)
  }
  if (is.null(by)) {
    return(singlefold(obs, k))
  }
  by = as.vector(as.matrix(by))
  if (length(by) != obs) {
    stop("by should be a vector with the same number of records as x")
  }
  un <- unique(by)
  group <- vector(length = obs)
  for (u in un) {
    i = which(by == u)
    kk = min(length(i), k)
    if (kk < k) 
      warning("lowered k for by group: ", u, "  because the number of observations was  ", 
              length(i))
    group[i] <- singlefold(length(i), kk)
  }
  return(group)
}
