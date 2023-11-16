


# ====================== Analyse full model =======================

# 1. Calculate and tabulate model information criteria/summary metrics
# 2. Visualise fitted effects
# 4. Tabulate parameter estimates
# 5. Run posterior predictive checks for calibration


library(dplyr)
library(sf)
library(ggplot2)
library(vroom)
library(raster)

PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

source("./scripts/04_modelling/00_inla_setup_functions_r4.R")
source("./scripts/04_modelling/00_plot_themes.R")

# directory where model objects are stored (large files so external)
models_dir = "C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/output/model_outputs/results_23/"



# ================= shapefile and dengue data =================

# districts to be excluded (offshore)
# only ones with substantial dengue cases are in Kien Giang; could be worth exploring including them but thisis sufficient for now
offshore_areas = c(70154, 70339, 70273, 70355, 70698)

# districts shapefile for Vietnam
shp = st_read("./data/shapefiles/dengue_districts_shapefile.shp") %>%
  dplyr::filter(!areaid %in% offshore_areas)

# provinces shapefile with regions info, and vietnam outline for mapping
shp_prov = st_read("./data/shapefiles/provinces.shp")
st_crs(shp_prov) = st_crs(shp)
shp_prov = st_crop(shp_prov, shp)
rr = read.csv("./data/shapefiles/regions_lookup.csv")
shp_prov = left_join(shp_prov, rr[ , c("provincena", "region")])
shp_vt = st_read("./data/shapefiles/vt_national.shp") %>%
  st_crop(shp)

# dengue, regions, climate, landuse, connectivity data
dd = read.csv("./output/model_data/ModelData_Dengue_VietAll.csv") %>%
  dplyr::left_join(read.csv("./output/model_data/ModelData_ClimaticRegions.csv")) %>%
  dplyr::left_join(read.csv("./output/model_data/ModelData_SocioEcologicalCovar_VietAll.csv")) %>%
  dplyr::left_join(vroom::vroom("./output/model_data/ModelData_ClimateLags_VietAll.csv.gz", col_types = list(areaid=vroom::col_character(), date=vroom::col_character()))) %>%
  dplyr::mutate(date = as.Date(date))

# provincex
dd$provincex = as.integer(as.factor(dd$province))

# polyid indicator
shpf = shp[ shp$areaid %in% dd$areaid, ]
id_ref = data.frame(areaid = shpf$areaid, polyid = 1:nrow(shpf))
dd = left_join(dd, id_ref)

# add polygon area, lat lon, region information
#shp$area_km2 = as.vector(st_area(shp) / 10^6)
shp = cbind(shp, as.data.frame(st_coordinates(st_centroid(shp))) %>% dplyr::rename("longitude"=1, "latitude"=2))
dd = left_join(dd, shp[ , c("areaid", "latitude", "longitude")] %>% st_drop_geometry())
shp = left_join(shp, dd[ !duplicated(dd$areaid), c("areaid", "region1", "region2", "region3") ])





# ======================== Read full model ======================== 

# full model and baseline model
load(paste(models_dir, "fullmodel/models/fullmod_nb_model_1.R", sep=""))
mx = mod_i
rm(mod_i)

load(paste(models_dir, "fullmodel/models/fullmod_nb_model_baseline.R", sep=""))
mb = mod_i
rm(mod_i)



# ======================= Tablulate change in metrics from baseline =====================

# extract ranef mae
calc_MAERE = function(model){
  strf = extractRandomINLA(model$summary.random$polyid, effect_name ="", model_is_bym = TRUE) %>%
    dplyr::filter(component == "uv_joint") %>%
    dplyr::left_join(data.frame(group=1:23, year=1998:2020)) %>%
    dplyr::left_join(
      dd %>% dplyr::select(polyid, areaid, year_useable_from) %>% distinct(),
      by=c("value"="polyid")
    ) %>%
    dplyr::mutate(mean = replace(mean, year < year_useable_from, NA)) %>%
    dplyr::filter(group > 4)
  mae_re = mean(abs(strf$mean), na.rm=TRUE)
  return(mae_re)
}

# metrics
metrics = rbind(
  fitMetricsINLA(mb, dd, modname = "Baseline")[ , c(1, 2, 3, 5)],
  fitMetricsINLA(mx, dd, modname = "Full")[ , c(1, 2, 3, 5)]
) %>%
  cbind(
    data.frame(MAE_re = c(calc_MAERE(mb), calc_MAERE(mx)))
  ) %>%
  dplyr::rename(
    "Model"=1, "DIC"=2, "WAIC"=3, "log-score"=4
  )
metrics = rbind(metrics, as.numeric(unlist(c("Difference", metrics[2, 2:5]-metrics[1, 2:5]))))
metrics[ , 2:3 ] = round(metrics[ , 2:3 ], 1)
metrics[ , 4:5 ] = round(metrics[ , 4:5 ], 3)
metrics$Model[3] = "Difference"

write.csv(metrics, "./output/figures/SuppTable_FullModelSummaryMetrics.csv", row.names=FALSE)



# ====================== Extract and save fitted climate effects and parameters ========================

# ranefs
ranefs = mx$summary.random
rf = lapply(ranefs[3:length(ranefs)], extractRandomINLA, effect_name = "x", transform=FALSE)
for(i in 1:length(rf)){
  rf[[i]]$effect = names(ranefs[ 3:length(ranefs)][i] )
}
rf = do.call(rbind.data.frame, rf)
row.names(rf) = c()
write.csv(rf, "./output/model_outputs/fullmodel/fitted_params/fullmodel_fittedclimatefunctions_rw2.csv", row.names=FALSE)

# fixed effects
fixefs = extractFixedINLA(mx, model_name = "full_model") %>%
  dplyr::mutate(description = "slope parameter for scaled covariate")
row.names(fixefs) = c()
write.csv(fixefs, "./output/model_outputs/fullmodel/fitted_params/fullmodel_fixedeffects.csv", row.names=FALSE)



# ======================== Visualise full results ===========================

# colours for different variable types
col_clim = viridis::viridis(200)[40]
col_socio = viridis::viridis(200)[105]

# ranefs
ranefs = mx$summary.random
rf = lapply(ranefs[3:length(ranefs)], extractRandomINLA, effect_name = "x", transform=TRUE)
for(i in 1:length(rf)){
  rf[[i]]$effect = names(ranefs[ 3:length(ranefs)][i] )
}
rf = do.call(rbind.data.frame, rf)

# extract ranefs and create standardised plots
effs = unique(rf$effect)
plots = vector("list", length(effs))
plotRF = function(x){
  
  rr=rf[ rf$effect == effs[x], ]
  
  if(effs[x] %in% c("tmean_1m_g", "spei1_1m_g", "spei6_5m_g")){
    colx = col_clim
  } else{
    colx = col_socio
  }
  
  # if(effs[x] == "tmean_1m_g"){
  #   # rr$upper[ rr$value < 14 ] = NA
  #   # rr$lower[ rr$value < 14 ] = NA
  #   rr = rr %>% dplyr::filter(value >= quantile(ddf$tmean_1m, 0.001))
  # }
  
  rr = left_join(
    rr,
    data.frame(
      effect = c("flushany_g", "water_g", "tmean_1m_g", "spei1_1m_g", "spei6_5m_g"),
      effname = c("Hygienic toilet access\n(indoor/outdoor)", "Piped or borehole-derived\nwater access", "Tmean (1-month lag)", "SPEI-1 (1-month lag)", "SPEI-6 (5-month lag)"),
      type = c("Infrastructure", "Infrastructure", "Climate", "Climate", "Climate"),
      units = c("Prop. households", "Prop. households", "Tmean (°C)", "SPEI-1", "SPEI-6")
    )
  )
  
  # visualise plot
  px = ggplot(rr) + 
    geom_line(aes(value, median)) +
    geom_ribbon(aes(value, ymin=lower, ymax=upper), alpha=0.2, fill=colx) +
    geom_hline(yintercept=1, lty=2) +
    #facet_wrap(~effname, scales="free") +
    theme_classic() +
    ggtitle(rr$effname[1]) + 
    ylab("Relative risk") + xlab(rr$units[1]) +
    theme(plot.title = element_text(size=12.25, hjust=0.5),
          axis.text = element_text(size=10.25),
          axis.title = element_text(size=11))
  
  # add vline for 0 in SPEI demarcating wet/dry
  if(grepl("spei", effs[x])){
    px = px + geom_vline(xintercept=0, lty=1, size=0.25, col="grey70", alpha=0.4)
  }
  
  # add density of climate vars
  if(effs[x] == "tmean_1m_g"){
    densx = dd %>% 
      dplyr::filter(tmean_1m >= min(rr$value) & tmean_1m <= max(rr$value)) %>%
      ggplot() + 
      geom_density(aes(tmean_1m), fill=colx, alpha=0.2,  size=0.4, color="grey50") + 
      theme_classic() +
      theme(plot.title = element_text(size=12.25, hjust=0.5),
            axis.text.y = element_text(size=8),
            axis.text.x = element_text(size=10.25),
            axis.title = element_text(size=11), 
            axis.title.y = element_text(size=10, color="black")) + 
      xlab(rr$units[1]) + 
      ylab("Density") + 
      scale_y_continuous(n.breaks=3)
    
    px = gridExtra::grid.arrange(px + theme(axis.title.x = element_blank()), 
                                 densx, nrow=2, heights=c(0.8, 0.25))
  }
  
  if(effs[x] == "spei1_1m_g"){
    densx = dd %>% 
      dplyr::filter(spei1_1m >= min(rr$value) & spei1_1m <= max(rr$value)) %>%
      ggplot() + 
      geom_vline(xintercept=0, lty=1, size=0.25, col="grey70", alpha=0.4) +
      geom_density(aes(spei1_1m), fill=colx, alpha=0.2,  size=0.4, color="grey50") + 
      theme_classic() +
      theme(plot.title = element_text(size=12.25, hjust=0.5),
            axis.text.y = element_text(size=10.25),
            axis.text.x = element_text(size=10.25),
            axis.title = element_text(size=11), 
            axis.title.y = element_text(size=10, color="black")) + 
      xlab(rr$units[1]) + 
      ylab("Density") + 
      scale_y_continuous(n.breaks=2)
    
    px = gridExtra::grid.arrange(px + theme(axis.title.x = element_blank()), 
                                 densx, nrow=2, heights=c(0.8, 0.25))
  }
  
  if(effs[x] == "spei6_5m_g"){
    densx = dd %>% 
      dplyr::filter(spei6_5m >= min(rr$value) & spei6_5m <= max(rr$value)) %>%
      ggplot() + 
      geom_vline(xintercept=0, lty=1, size=0.25, col="grey70", alpha=0.4) +
      geom_density(aes(spei6_5m), fill=colx, alpha=0.2, size=0.4, color="grey50") + 
      theme_classic() +
      theme(plot.title = element_text(size=12.25, hjust=0.5),
            axis.text.y = element_text(size=10.25),
            axis.text.x = element_text(size=10.25),
            axis.title = element_text(size=11), 
            axis.title.y = element_text(size=10, color="black")) + 
      xlab(rr$units[1]) + 
      ylab("Density") + 
      scale_y_continuous(n.breaks=3)
    
    px = gridExtra::grid.arrange(px + theme(axis.title.x = element_blank()), 
                                 densx, nrow=2, heights=c(0.8, 0.25))
  }  
  
  # return
  plots[[x]] = px
}
prf = lapply(1:length(effs), plotRF)


library(ggh4x)

# fixed effects plot
pfx = extractFixedINLA(mx, model_name="mod", transform=TRUE) %>%
  dplyr::filter(param != "Intercept") %>%
  dplyr::filter(!grepl("province|areaid", param)) %>%
  dplyr::left_join(
    data.frame(
      param=c("tmean_coolestmonth_s", "gravityf_log", "urban_s", "urbanexp10_log", "traffic_kmperinhab_log"),
      paramname = c("Tmean\ncoolest\nmonth", "Gravity\n(log)", "Built-up\nland", "Urban\nexpansion\nrate (log)", "Traffic\nper inhab\n(log)"),
      type = c("Climate", "Socio-env.", "Socio-env.", "Socio-env.", "Socio-env."),
      facet = c("Slope (risk ratio) ", "Slope (risk ratio)", "Slope (risk ratio)", "Slope (risk ratio)", "Slope (risk ratio)")
    )
  ) %>%
  dplyr::mutate(paramname = factor(paramname, levels=c("Tmean\ncoolest\nmonth", "Gravity\n(log)", "Traffic\nper inhab\n(log)",  "Built-up\nland", "Urban\nexpansion\nrate (log)"), ordered=TRUE),
                facet=factor(facet, levels=c("Slope (risk ratio) ", "Slope (risk ratio)"), ordered=TRUE)) %>%
  ggplot() +
  geom_point(aes(paramname, mean, col=type), size=2) +
  geom_linerange(aes(paramname, ymin=lower, ymax=upper, col=type), show.legend=FALSE) +
  geom_hline(yintercept=1, lty=2) +
  theme_classic() +
  scale_color_manual(values=c("Climate"=col_clim, "Socio-env."=col_socio)) +
  #scale_y_continuous(limits=c(0.7, 1.8), breaks=c(0.75, 1, 1.25, 1.5, 1.75), labels=c(0.75, 1, 1.25, 1.5, 1.75)) +
  ylab("Slope (risk ratio)") + 
  xlab("") + 
  #scale_color_viridis_d(begin=0.2, end=0.6) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0, size=11),
        axis.title.y = element_blank(),
        #axis.title.y = element_text(size=11),
        axis.text.y = element_text(size=10)) +
  facet_wrap(~facet, ncol=2, scales="free", strip.position="left") +
  force_panelsizes(cols = c(0.2, 1)) +
  #ggtitle("Waffle\nwaffles") +
  theme(#strip.text = element_blank(),
        strip.text = element_text(size=11),
        #plot.title = element_text(size=12.15, color="white"),
        strip.background = element_blank(),
        strip.placement="outside",
        legend.text = element_text(size=10.5), legend.title = element_blank(), #legend.background=element_rect(fill=NA, color="grey70"),
        legend.position=c(0.83, 0.9)) + guides(colour = guide_legend(override.aes = list(size=5, alpha=0.8)))

# pfx2 = extractFixedINLA(mx, model_name="mod", transform=TRUE) %>%
#   dplyr::filter(param != "Intercept") %>%
#   dplyr::filter(grepl("tmean", param)) %>%
#   dplyr::left_join(
#     data.frame(
#       param=c("tmean_coolestmonth_s", "urban_s", "urbanexp10_log", "traffic_kmperinhab_log"),
#       paramname = c("Tmean\ncoolest\nmonth", "Built-up\nland", "Urban\nexpansion\nrate (log)", "Traffic\nper inhab\n(log)"),
#       type = c("Climate", "Urbanisation", "Urbanisation", "Mobility")
#     )
#   ) %>%
#   ggplot() +
#   geom_point(aes(paramname, mean), size=2, color=col_clim) +
#   geom_linerange(aes(paramname, ymin=lower, ymax=upper), color=col_clim) +
#   geom_hline(yintercept=1, lty=2) +
#   theme_classic() +
#   ylab("Slope (risk ratio)") + 
#   xlab("") + 
#   #scale_color_viridis_d(begin=0.2, end=0.6) + 
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(angle=0, size=11),
#         axis.text.y = element_text(size=10),
#         axis.title.y = element_text(size=11))

# fixed effects
#pfx = gridExtra::grid.arrange(pfx2, pfx1, widths=c(0.3, 1))

# row 1: fixed + socioenv
pc1 = gridExtra::grid.arrange(pfx, 
                              prf[[1]], 
                              prf[[2]], nrow=1, widths=c(1.2, 0.9, 0.9))

# row 2: climatic
pc2 = gridExtra::grid.arrange(prf[[3]], 
                              prf[[4]], 
                              prf[[5]], nrow=1)

# combine and save
pcomb = gridExtra::grid.arrange(pc1, pc2, nrow=2, heights=c(1, 1.15))
pcomb = ggpubr::as_ggplot(pcomb)  +
  cowplot::draw_plot_label(label = c("a", "b", "c", "d", "e", "f", "g"), 
                           fontface = "bold", size = 20, 
                           x = c(0.01, 0.14, 0.442, 0.75, 0.04, 0.39, 0.715), y = c(1, 1, 0.97, 0.97, 0.53, 0.53, 0.53))
ggsave(pcomb, file="./output/figures/Figure3_DengueDrivers.png", dpi=600, device="png", units="in", width=11.5, height=7.2, scale=0.95)
ggsave(pcomb, file="./output/figures/Figure3_DengueDrivers.pdf", device="pdf", units="in", width=11.5, height=7.2, scale=0.95)



# ================= Visualise ST random effects ====================

rf = extractRandomINLA(mx$summary.random$polyid, effect_name ="", model_is_bym = TRUE) %>%
  dplyr::filter(component == "uv_joint") %>%
  dplyr::left_join(data.frame(group=1:23, year=1998:2020)) %>%
  dplyr::left_join(
    dd %>% dplyr::select(polyid, areaid, year_useable_from) %>% distinct(),
    by=c("value"="polyid")
  ) %>%
  dplyr::mutate(mean = replace(mean, year < year_useable_from, NA)) 

#cs = colorRampPalette(RColorBrewer::brewer.pal(11, "BrBG"))(200)
lims = max(abs(rf$mean), na.rm=TRUE)
lims = c(-lims, lims)

p1 = shp %>%
  dplyr::select(areaid) %>%
  dplyr::full_join(rf) %>%
  sf::st_as_sf() %>%
  ggplot()+
  geom_sf(aes(fill=mean), color=NA) + 
  geom_sf(data = shp_vt, color="grey20", size=0.5, fill=NA) + 
  scale_fill_gradientn(colors=rev(colorRampPalette(MetBrewer::met.brewer("Benedictus", 11))(200)), na.value="white", limits=lims, name="Posterior\nmean") +
  theme_void() + 
  facet_wrap(~year, nrow=3) +
  theme(strip.text = element_text(size=13), 
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))
ggsave(p1, file="./output/figures/SuppFigure_STRanefs.jpg", device="jpg", units="in", dpi=300, width=12, height=9, scale=0.85)


# ================== Tabulate full model param results ===============

# fxx 
fxx = extractFixedINLA(mx) %>%
  dplyr::left_join(
    data.frame(
      param=c("gravityf_log", "urban_s", "urbanexp10_log", "traffic_kmperinhab_log", "tmean_coolestmonth_s", "Intercept"),
      paramname = c("Gravity (log)", "Built-up land", "Urban expansion rate (10 year, log)", "Traffic per inhab (log)",
                    "Tmean coolest month", "Intercept"))
  ) %>%
  dplyr::mutate(Type = "Fixed effect") %>%
  dplyr::select(paramname, Type, median, lower, upper) %>%
  dplyr::mutate(median = round(median, 3), lower=round(lower, 3), upper=round(upper, 3)) %>%
  dplyr::rename("Parameter" = 1, "Median"=3, "CI_0.025"=4, "CI_0.975"=5) 

# hyperparams
hyp = mx$summary.hyperpar
hyp = round(hyp, 3)
hyp$param = row.names(hyp)
hyp = hyp %>%
  dplyr::left_join(
    data.frame(
      param = c("size for the nbinomial observations (1/overdispersion)", 
                "Precision for monthdengue",
                "Precision for polyid",
                "Phi for polyid",
                "Precision for flushany_g", 
                "Precision for water_g",
                "Precision for tmean_1m_g",
                "Precision for spei1_1m_g",
                "Precision for spei6_5m_g"),
      Parameter = c("Size (1/overdisp) for neg. binom",
                    "Precision for Month (RW1)", 
                    "Precision for District (BYM2)",
                    "Phi for District (BYM2)", 
                    "Precision for Hygienic toilet (RW2)", 
                    "Precision for Piped water (RW2)",
                    "Precision for Tmean 1m (RW2)", 
                    "Precision for SPEI-1 1m (RW2)", 
                    "Precision for SPEI-6 5m (RW2)")
    )
  ) %>%
  dplyr::mutate(Type = "Hyperparameter") %>%
  dplyr::select(8, 9, 4, 3, 5) %>%
  dplyr::rename(
    "Median"=3, "CI_0.025"=4, "CI_0.975"=5
  )

# combine and save
tab = rbind(fxx, hyp)
write.csv(tab, "./output/figures/SuppTable_FullModelTabulation.csv", row.names=FALSE)






# ======================== Posterior predictive check for model calibration =====================

# extract 500 samples
# ss = INLA::inla.posterior.sample(n=500, result=mx)
# save(ss, file="./finalmodel_posteriorsamples_20220419.R")
load(file="./finalmodel_posteriorsamples_20220419.R")

# remove model
rm(mx); rm(mod_i); rm(mb)

# for each extract latent field and simulate 25 outcomes
#for(i in 1:length(ss)){
for(i in 1:100){

  #i = 1
  cat(paste(i, "...", sep=""))
  latent_i = ss[[i]]$latent[ 1:nrow(ddf) ]
  od_i = as.vector(ss[[i]]$hyperpar[1])
  
  # does not need pop adding
  # ddf_i = ddf %>%
  #   dplyr::select(areaid, date, logpop, cases) %>%
  #   dplyr::mutate(latent = latent_i + logpop) 
  # for_sim = exp(ddf_i$latent)
  
  # exponentiate to mean
  for_sim = exp(latent_i)
  
  # simulate 25 outcomes
  sim_nbinom = function(x){
    rnbinom(n=25, mu=for_sim[x], size=od_i) 
  }
  res_i = lapply(1:length(for_sim), sim_nbinom)
  res_i = do.call(rbind, res_i)
  
  # combine and cbind to growing matrix of distributions
  if(i == 1){
    result = res_i
  } else{
    result = cbind(result, res_i)
  }
}



# #example
# x = sample(which(ddf$cases > 30), 1)
# #x = 42312      
# hist(result[x, ], 50)
# abline(v= ddf$cases[x], col="red")
# quantile(result[ x, ], c(0.025, 0.16667, 0.5, 0.83333, 0.975))

# quantiles = 95% and 67% predictive interval
qq = t(apply(result, 1, quantile, c(0.025, 0.16667, 0.5, 0.83333, 0.975)))

# combine
ppd = ddf %>%
  dplyr::select(areaid, region, province, date, cases) %>%
  cbind(qq) %>%
  dplyr::rename(
    "lower_95"=6, "lower_67"=7, "median"=8, "upper_67"=9, "upper_95"=10
  ) %>%
  dplyr::mutate(
    within_95 = cases >= lower_95 & cases <= upper_95,
    within_67 = cases >= lower_67 & cases <= upper_67,
    region2 = replace(region, region %in% c("Southeast", "Mekong River Delta"), "South"),
    region2 = replace(region2, region2 %in% c("Northeast", "Northwest"), "Northeast/Northwest")
  ) 

# regional PPD 
cal_plot = ppd %>%
  dplyr::select(-province, -areaid, -within_95, -within_67, -region2) %>%
  dplyr::group_by(region, date) %>%
  dplyr::summarise_all("sum", na.rm=TRUE) %>%
  dplyr::mutate(date=as.Date(date)) %>%
  # dplyr::mutate(region2 = factor(region2, levels=c("Northeast/Northwest", "Red River Delta",
  #                                                  "Central Highlands", "North Central",
  #                                                  "South", "South Central Coast"))) %>%
  dplyr::mutate(region = factor(region, levels=c("Northwest", "Northeast", "North Central", "Red River Delta",
                                                   "Central Highlands", "South Central Coast",
                                                   "Mekong River Delta", "Southeast"))) %>%
  ggplot() + 
  theme_classic() + 
  geom_ribbon(aes(date, ymin=lower_95/1000, ymax=upper_95/1000), alpha=0.2, fill="skyblue4", col=NA) + 
  geom_ribbon(aes(date, ymin=lower_67/1000, ymax=upper_67/1000), alpha=0.4, fill="skyblue4", col=NA) + 
  geom_line(aes(date, cases/1000), col="red", size=0.35) + 
  facet_wrap(~region, scales="free_y", ncol=2) +
  ylab("Dengue cases (thousands)") + 
  xlab("Month") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=14),
        axis.title = element_text(size=13),
        axis.text = element_text(size=11))
ggsave(cal_plot, file="./output/figs_fin/SuppFigureX_ModelPPDCalibration.png", device="png", units="in", width=9, height=8, dpi=600, scale=0.95)

# overall summaries
ppd %>%
  dplyr::group_by(region2) %>%
  dplyr::summarise(
    within_95 = sum(within_95, na.rm=TRUE) / length(within_95) ,
    within_67 = sum(within_67, na.rm=TRUE) / length(within_67) 
  )
ppd %>%
  dplyr::summarise(
    within_95 = sum(within_95, na.rm=TRUE) / length(within_95) ,
    within_67 = sum(within_67, na.rm=TRUE) / length(within_67) 
  )



# ============== more granular examination of distribution across pred intervals

# what interval falls in?
qq = t(apply(result, 1, quantile, c(0.025, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.975)))
get_interval = function(x){ 
  obs =  ddf$cases[x]
  lt = obs < as.vector(qq[x, ])
  return(sum(lt==FALSE))
  }
ppd$interval = sapply(1:nrow(qq), get_interval)

table(ppd$interval[ ppd$cases > 0 ])
hist(ppd$interval[ ppd$cases > 0 ])

# int names
int_names = c("<2.5", "2.5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", ">97.5")

p0 = ppd %>%
  dplyr::filter(region == "South Central Coast") %>%
  dplyr::filter(lubridate::year(as.Date(date)) >= 2002) %>%
  dplyr::filter(cases > 0) %>%
  dplyr::group_by(date, interval) %>%
  dplyr::summarise(nobs = length(interval)) %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(nobs_total = sum(nobs)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(prop_obs = nobs / nobs_total) 
expand.grid(unique(p0$date), unique(p0$interval)) %>%
  dplyr::rename("date"=1, "interval"=2) %>%
  dplyr::left_join(p0) %>%
  dplyr::mutate(prop_obs = replace(prop_obs, is.na(prop_obs), 0),
                date = as.Date(date)) %>%
  ggplot() + 
  geom_tile(aes(date, interval, fill=prop_obs), width=31) + 
  scale_fill_gradientn(colors = viridisLite::turbo(200)) + 
  #scale_fill_gradient(low="white", high="black") + 
  theme_classic()
