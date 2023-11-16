

# ============================== UNIVARIATE MODEL ANALYSES =====================================

library(dplyr)
library(sf)
library(ggplot2)
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
  dplyr::left_join(vroom::vroom("./output/model_data/ModelData_ClimateLags_VietAll.csv.gz", col_types = list(areaid=col_character(), date=col_character()))) %>%
  dplyr::mutate(date = as.Date(date))

# provincex
dd$provincex = as.integer(as.factor(dd$province))

# polyid indicator
shpf = shp[ shp$areaid %in% dd$areaid, ]
id_ref = data.frame(areaid = shpf$areaid, polyid = 1:nrow(shpf))
dd = left_join(dd, id_ref)

# name fix
dd$district[ dd$areaid == 70266 ] = "Ninh Binh"

# add polygon area, lat lon, region information
#shp$area_km2 = as.vector(st_area(shp) / 10^6)
shp = cbind(shp, as.data.frame(st_coordinates(st_centroid(shp))) %>% dplyr::rename("longitude"=1, "latitude"=2))
dd = left_join(dd, shp[ , c("areaid", "latitude", "longitude")] %>% st_drop_geometry())
shp = left_join(shp, dd[ !duplicated(dd$areaid), c("areaid", "region1", "region2", "region3") ])




# ============ view goodness of fit (DIC) metrics for all univariate models ===============
  
# functions for reading in
countvars = function(x){ sapply(strsplit(as.vector(x), "[+]"), length) }
readfile = function(x){
  foo = read.csv(x)
  if(!"covar" %in% names(foo)){ foo$covar = "" }
  if(!"model_sub" %in% names(foo)){ foo$model_sub = "" }
  if(!"model_filename" %in% names(foo)){ foo$model_filename = "" }
  foo
}

# socio-environmental covariate models
ff = list.files("./output/model_outputs/sociounivar/fitmetrics/", full.names=TRUE, pattern=".csv")

fx1 = do.call(rbind.data.frame, lapply(ff, readfile)) %>%
  dplyr::mutate(num_predictors = countvars(modname),
                covar = modname) %>%
  dplyr::arrange(waic) %>%
  dplyr::select(-formula)

fx1$deltawaic = fx1$waic - fx1$waic[ fx1$modname == "baseline"]
fx1$deltadic = fx1$dic - fx1$dic[ fx1$modname == "baseline"]
fx1$deltamae = fx1$mae - fx1$mae[ fx1$modname == "baseline"]
fx1$deltaLS = fx1$logscore - fx1$logscore[ fx1$modname == "baseline"]

fx1$covar_grp = NA
fx1$covar_grp[ grep("radiation|gravity|flights|traffic", fx1$modname) ] = "Mobility"
fx1$covar_grp[ grep("popdens", fx1$modname) ] = "Demography"
fx1$covar_grp[ grep("urban", fx1$modname) ] = "Urbanisation"
fx1$covar_grp[ grep("infra|water|flush", fx1$modname) ] = "Infrastructure"
fx1$covar_grp[ grep("tmin|tmean", fx1$modname) ] = "Temperature"

# climate models
ff = list.files("./output/model_outputs/climunivar/fitmetrics/", full.names=TRUE, pattern=".csv")

fx2 = do.call(rbind.data.frame, lapply(ff, readfile)) %>%
  dplyr::mutate(num_predictors = countvars(modname), effect_type = "rw2") %>%
  dplyr::arrange(dic) %>%
  dplyr::select(-formula)

fx2$deltawaic = fx2$waic - fx2$waic[ fx2$modname == "baseline"]
fx2$deltadic = fx2$dic - fx2$dic[ fx2$modname == "baseline"]
fx2$deltamae = fx2$mae - fx2$mae[ fx2$modname == "baseline"]
fx2$deltaLS = fx2$logscore - fx2$logscore[ fx2$modname == "baseline"]

fx2$modname[ fx2$modname == "baseline"] = "climate_baseline"
fx2$covar_grp = "Temperature"
fx2$covar_grp[ grep("spei|precip", fx2$covar) ] = "Hydrometeorology"

# remove 6 month lag from temp
fx2 = fx2[ !fx2$covar %in% c("tmean_6", "tmax_6", "tmin_6", "tmean_5", "tmax_5", "tmin_5"), ]

# keep only the best lag for each
fx2 = fx2 %>% 
  dplyr::mutate(xx = unlist(lapply(strsplit(covar, "_"), "[", 1))) %>%
  dplyr::group_by(xx) %>%
  dplyr::filter(deltawaic == min(deltawaic)) %>%
  #dplyr::filter(deltadic == min(deltadic)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-xx)

# combine
fx = rbind(fx1, fx2)

# plot change in DIC
fx %>% 
  dplyr::arrange(deltawaic) %>%
  tibble::rowid_to_column("id") %>%
  ggplot() + 
  geom_point(aes(id, deltawaic, col=num_predictors), size=2) +
  geom_line(aes(id, deltawaic)) + 
  theme_minimal() +
  scale_color_viridis_c(option="magma", end=0.8)

# plot change in dic for different vars in univariate
dat_subset = fx %>%
  #dplyr::filter(type == "uni") %>%
  dplyr::filter(!modname %in% c("climate_baseline", "baseline")) %>%
  dplyr::filter(!grepl("spei2_", covar)) %>%
  dplyr::mutate(covar = replace(covar, grep("infra_g", covar), "infra_rw"), 
                covar = replace(covar, grep("flushindoor_g", covar), "flushindoor_rw"),
                covar = replace(covar, grep("flushany_g", covar), "flushany_rw"),
                covar = replace(covar, grep("traffic_perinhab_g", covar), "traffic_rw"),
                covar = replace(covar, grep("water_g", covar), "pipedwater_rw"),
                covar = replace(covar, grep("popdens_g", covar), "popdens_rw"),
                covar = replace(covar, grep("urbancensus_g", covar), "urbancensus_rw"),
                covar = replace(covar, grep("traffic_milperskm_g", covar), "traffic_milperskm_rw"),
                covar = replace(covar, grep("gravity_g", covar), "gravity_rw"),
                covar = replace(covar, grep("radiation_g", covar), "radiation_rw"),
                covar = replace(covar, grep("urban_g", covar), "urban_rw"),
                covar = replace(covar, grep("flushoutdoor_g", covar), "flushoutdoor_rw"),
                covar = replace(covar, grep("tmin_annualmean_g", covar), "tmin_annual_rw"),
                covar = replace(covar, grep("tmin_coolestmonth_g", covar), "tmin_coolest_rw"),
                covar = replace(covar, grep("tmean_coolestmonth_g", covar), "tmean_coolest_rw"),
                covar = replace(covar, grep("tmean_annualmean_g", covar), "tmean_annual_rw")) %>%
  dplyr::arrange(covar_grp, desc(deltawaic)) %>%
  dplyr::mutate(covar = factor(covar, levels=covar, ordered=TRUE)) 

dat_subset = dat_subset %>%
  dplyr::filter(!covar %in% c("urbanexp10_s",
                              "urbanexp3_s",
                              "radiationf_s",
                              "gravityf_s",
                              "popdens_s",
                              "infra_s",
                              "infra_rw",
                              "water_s",
                              "traffic_rw",
                              "flushany_s",
                              #"flushany_rw",
                              "flushindoor_s",
                              "infra_s_static",
                              "tmin_annual_rw",
                              "tmean_annual_rw",
                              "tmin_coolest_rw",
                              "tmin_annualmean_all",
                              "tmin_coolestmonth_all",
                              "flushoutdoor_rw",
                              "flushoutdoor_s",
                              "flushoutdoor_s_static",
                              "urban_rw",
                              "radiation_rw",
                              "gravity_rw",
                              "traffic_milperskm_rw",
                              "popdens_rw",
                              "urbancensus_rw",
                              "urban_census_s",
                              "urban_census_s_static",
                              "urbanexp10_dev",
                              "tmean_coolest_rw",
                              "traffic_kmperinhab_s",
                              "traffic_milperskm_log_static",
                              "traffic_milperskm_log",
                              "traffic_thouskmperinhab_s_static")) %>%
  dplyr::mutate(covar = as.vector(covar),
                covar = replace(covar, covar == "flushindoor_rw", "Indoor hygeinic toilet access (rw)"),
                covar = replace(covar, covar == "flushany_rw", "Hygeinic toilet access (rw)"),
                #covar = replace(covar, covar == "flushany_s", "Hygeinic toilet access"),
                covar = replace(covar, covar == "flushany_s_static", "Hygeinic toilet access (temporally-fixed)"),
                covar = replace(covar, covar == "flushindoor_s_static", "Indoor hygeinic toilet access (temporally-fixed)"),
                covar = replace(covar, covar == "pipedwater_rw", "Piped water access (rw)"),
                covar = replace(covar, covar == "water_s_static", "Piped water access (temporally-fixed)"),
                # covar = replace(covar, covar == "infra_rw", "Infrastructure index (rw)"),
                # covar = replace(covar, covar == "infra_s_static", "Infrastructure index (temporally-fixed)"),
                #covar = replace(covar, covar == "traffic_rw", "Traffic distance per inhab (rw)"),
                covar = replace(covar, covar == "urban_s", "Built-up land"),
                covar = replace(covar, covar == "urbanexp10_log", "Urban expansion rate (10yr) log"),
                covar = replace(covar, covar == "urbanexp10_log_static", "Urban expansion rate (temporally-fixed)"),
                covar = replace(covar, covar == "urban_s_static", "Built-up land (temporally-fixed)"),
                covar = replace(covar, covar == "urbanexp3_log", "Urban expansion rate (3yr) log"),
                covar = replace(covar, covar == "gravityf_log", "Gravity log"),
                covar = replace(covar, covar == "gravityf_log_static", "Gravity log (temporally-fixed)"),
                covar = replace(covar, covar == "radiationf_log", "Radiation log"),
                covar = replace(covar, covar == "radiationf_log_static", "Radiation log (temporally-fixed)"),
                covar = replace(covar, covar == "tmin_annualmean_s", "Tmin (annual mean)"),
                covar = replace(covar, covar == "tmean_annualmean_s", "Tmean (annual mean)"),
                covar = replace(covar, covar == "tmin_coolestmonth_s", "Tmin (coolest month)"),
                covar = replace(covar, covar == "tmean_coolestmonth_s", "Tmean (coolest month)"),
                covar = replace(covar, covar == "logpopdens_s", "Population density log"),
                covar = replace(covar, covar == "logpopdens_static", "Population density log (temporally-fixed)"),
                covar = replace(covar, covar == "traffic_milperskm_log_static", "Traffic total distance (temporally-fixed)"),
                covar = replace(covar, covar == "traffic_thouskmperinhab_s_static", "Road travel per inhab (temporally-fixed)"),
                covar = replace(covar, covar == "traffic_thouskmperinhab_s", "Road travel per inhab"),
                covar = replace(covar, covar == "traffic_milperskm_log", "Road travel total distance"),
                covar = replace(covar, covar == "traffic_kmperinhab_log", "Road travel per inhab log"),
                covar = replace(covar, covar == "traffic_kmperinhab_log_static", "Road travel per inhab log (temporally-fixed)"),
                covar = replace(covar, covar == "tmean_1", "Tmean (1m lag, rw)"),
                covar = replace(covar, covar == "tmax_1", "Tmax (1m lag, rw)"),
                covar = replace(covar, covar == "tdrange_6", "T diurnal range (6m lag, rw)"),
                covar = replace(covar, covar == "spei6_5", "SPEI-6 (5m lag, rw)"),
                covar = replace(covar, covar == "spei12_4", "SPEI-12 (4m lag, rw)"),
                covar = replace(covar, covar == "tmin_1", "Tmin (1m lag, rw)"),
                covar = replace(covar, covar == "spei3_6", "SPEI-3 (6m lag, rw)"),
                covar = replace(covar, covar == "spei1_1", "SPEI-1 (1m lag, rw)"),
                covar = replace(covar, covar == "precip_1", "Precipitation (1m lag, rw)"))

waic_plot = dat_subset %>%
  dplyr::filter(covar != "Multivariate") %>%
  dplyr::arrange(desc(covar_grp), desc(deltawaic)) %>%
  #dplyr::filter(!grepl("total distance", covar)) %>%
  dplyr::mutate(covar = factor(covar, levels=covar, ordered=TRUE)) %>%
  ggplot() +
  geom_point(aes(covar, deltawaic, fill=covar_grp), size=6, pch=21) +
  geom_hline(yintercept=0, lty=2, size=1) +
  coord_flip() +
  theme_minimal() +
  ylab(expression(paste(Delta, "WAIC"))) + xlab("") +
  theme(axis.text.y = element_text(size=13),
        panel.border = element_rect(color="grey20", fill=NA),
        axis.text.x = element_text(size=11),
        axis.title.x = element_text(size=14),
        plot.title = element_text(size=15, hjust=0.5),
        legend.title=element_blank(),
        legend.text = element_text(size=13),
        legend.position="top")



# ============= for each covariate, calculate the change in district-level random effects (i.e. unexplained variation) =================

# univar baseline
load(paste(models_dir, "sociounivar/models/sociouni_nb_model_baseline.R", sep=""))
b_bym = extractRandomINLA(mod_i$summary.random$polyid, effect_name = "baseline_uni", model_is_bym = TRUE) %>%
  dplyr::filter(component == "uv_joint") %>%
  dplyr::select(value, group, mean) %>%
  dplyr::rename("baseline_uni"=mean)

# climate baseline
load(paste(models_dir, "climunivar/models/climuni_nb_model_baseline.R", sep=""))
cb_bym = extractRandomINLA(mod_i$summary.random$polyid, effect_name = "baseline_clim", model_is_bym = TRUE) %>%
  dplyr::filter(component == "uv_joint") %>%
  dplyr::select(value, group, mean) %>%
  dplyr::rename("baseline_clim"=mean)

# combine
bym = left_join(b_bym, cb_bym)

# add in region info
#ddf = read.csv("./output/dataset_processed.csv")
rg = dd %>% 
  dplyr::select(areaid, province, district, polyid, region2, region3) %>%
  dplyr::filter(!duplicated(areaid)) %>%
  dplyr::rename("value"=polyid)
bym = left_join(bym, rg)

# result
result_df = data.frame()
ranefs_output = data.frame()

# create filepath for models
dat_subset$filepath=NA
dat_subset$filepath[ grep("sociouni", dat_subset$model_filename) ] = paste(models_dir, "sociounivar/models/", dat_subset$model_filename[ grep("sociouni", dat_subset$model_filename) ], sep="")
dat_subset$filepath[ grep("climuni", dat_subset$model_filename) ] = paste(models_dir, "climunivar/models/", dat_subset$model_filename[ grep("climuni", dat_subset$model_filename) ], sep="")

# for each model in dat_subset
for(i in 1:nrow(dat_subset)){
  
  # data subset
  print(i)
  ds = dat_subset[ i, ]
  
  # load model
  load(ds$filepath)
  
  # extract bym
  bym_i = extractRandomINLA(mod_i$summary.random$polyid, effect_name = "bym", model_is_bym = TRUE) %>%
    dplyr::filter(component == "uv_joint") %>%
    dplyr::select(value, group, mean) %>%
    dplyr::rename("modi"=mean)
  
  # combine with baseline byms
  bym_i = left_join(bym, bym_i)
  
  # define baseline
  if(grepl("climuni", ds$model_filename)){
    bym_i$bl = bym_i$baseline_clim
  } else{
    bym_i$bl = bym_i$baseline_uni
  }
    
  # exclude pre-2002 (because large areas of country have no data before this)
  bym_i = bym_i %>% dplyr::filter(group > 4) 

  # total unexplained spatial variation (proportion change in MAE)
  ds$MAEre = mean(abs(bym_i$modi))
  #ds$SDre = sd(bym_i$modi)
  ds$MAEre_bl = mean(abs(bym_i$bl))
  #ds$SDre_bl = sd(bym_i$bl)
  
  result_df = rbind(result_df, ds)
  ranefs_output = rbind(ranefs_output, bym_i %>% dplyr::select(areaid, district, group, modi, bl) %>% dplyr::mutate(model = ds$covar))
}

# calculate prop. change in spatiotemporal variation relative to baseline
result_df$ve_MAE = (result_df$MAEre - result_df$MAEre_bl) / result_df$MAEre_bl



# ======================= as above, but for the seasonal random effect ========================

# result
result_dft = data.frame()

# seasonal effect climate baseline
# climate baseline
load(paste(models_dir, "climunivar/models/climuni_nb_model_baseline.R", sep=""))
cb_seas = extractRandomINLA(mod_i$summary.random$month, effect_name = "baseline_clim") %>%
  dplyr::rename("provincex"=group) %>%
  dplyr::left_join(dd %>% dplyr::filter(!duplicated(province)) %>% dplyr::select(province, provincex, region, region3)) %>%
  dplyr::rename("bl"=mean)

ds_clim = dat_subset[ grep("Temperature|Hydro", dat_subset$covar_grp), ] 
ds_clim = ds_clim[ !ds_clim$covar %in% c("Tmin (coolest month)", "Tmean (annual mean)", "Tmin (annual mean)", "Tmean (coolest month)"), ]

# for each model in dat_subset
for(i in 1:nrow(ds_clim)){
  
  # load model
  ds = ds_clim[ i, ]
  load(ds$filepath)
  
  # seasoal effect
  seas_i = extractRandomINLA(mod_i$summary.random$month, effect_name = "seas_clim") %>%
    dplyr::rename("provincex"=group) %>%
    dplyr::left_join(dd %>% dplyr::filter(!duplicated(province)) %>% dplyr::select(province, provincex, region, region3)) %>%
    dplyr::select(value, province, mean) %>%
    dplyr::rename("modi"=mean)
  
  # combine with baseline byms
  seas_i = left_join(cb_seas, seas_i)
  
  # total unexplained spatial variation (proportion change in MAE)
  ds$seas_ve = (mean(abs(seas_i$modi)) - mean(abs(seas_i$bl))) / mean(abs(seas_i$bl))

  result_dft = rbind(result_dft, ds)
}


# ================== create visualisation ================

fac_levels = c("Demography", "Urbanisation", "Mobility", "Infrastructure", "Temperature", "Hydrometeorology")
plot = result_df %>%
  dplyr::filter(!grepl("Indoor", covar)) %>%
  dplyr::left_join(result_dft %>% dplyr::select(covar, seas_ve)) %>%
  dplyr::mutate(covar_grp = factor(covar_grp, levels=fac_levels, ordered=TRUE)) %>%
  dplyr::arrange(desc(covar_grp), desc(deltawaic)) %>%
  dplyr::mutate(covar = factor(covar, levels=covar, ordered=TRUE),
                ve_MAE = ve_MAE*100, 
                seas_ve = seas_ve*100) %>%
  dplyr::select(covar, covar_grp, deltawaic, ve_MAE, seas_ve) %>%
  dplyr::filter(covar != "Traffic distance per inhab (nonlinear)") %>%
  reshape2::melt(id.vars=1:2) %>%
  dplyr::mutate(variable = as.vector(variable), 
                variable = replace(variable, variable == "deltawaic", "ΔWAIC"),
                variable = replace(variable, variable == "ve_MAE", "% change in spatiotemporal\nrandom effects variation"),
                variable = replace(variable, variable == "seas_ve", "% change in seasonal\nrandom effects variation"),
                variable = factor(variable, levels=c("ΔWAIC", "% change in spatiotemporal\nrandom effects variation",
                                                     "% change in seasonal\nrandom effects variation"),
                                                     ordered=TRUE) ) %>%
  ggplot() + 
  geom_hline(yintercept=0, lty=2, size=0.75) +
  geom_point(aes(covar, value, fill=covar_grp), size=6, pch=21) +
  coord_flip() +
  facet_wrap(~variable, nrow=1, scales="free_x") +
  #scale_fill_viridis_d(option="viridis") + 
  scale_fill_met_d(name="Archambault") + 
  theme_minimal() +
  ylab("Change relative to baseline model") + xlab("") +
  theme(axis.text.y = element_text(size=14),
        panel.border = element_rect(color="grey20", fill=NA),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=15),
        plot.title = element_text(size=15, hjust=0.5),
        legend.title=element_blank(), 
        legend.text = element_text(size=15),
        legend.position="bottom", 
        strip.text = element_text(size=17)) 
plot

ggsave(plot, file="./output/figures/SuppFigure_UnivariateResults.jpg", device="jpg", dpi=600, units="in", width=16.5, height=9, scale=0.85)






# =================== Extract fitted climate covariate functions ======================

ff = list.files("./output/model_outputs/climunivar/fitmetrics/", full.names=TRUE, pattern=".csv")

fx2 = do.call(rbind.data.frame, lapply(ff, readfile)) %>%
  dplyr::mutate(num_predictors = countvars(modname), effect_type = "rw2") %>%
  dplyr::arrange(dic) %>%
  dplyr::select(-formula)

fx2$deltawaic = fx2$waic - fx2$waic[ fx2$modname == "baseline"]
fx2$deltadic = fx2$dic - fx2$dic[ fx2$modname == "baseline"]
fx2$deltamae = fx2$mae - fx2$mae[ fx2$modname == "baseline"]
fx2$deltaLS = fx2$logscore - fx2$logscore[ fx2$modname == "baseline"]

fx2$filepath = paste(models_dir, "climunivar/models/", fx2$model_filename, sep="")

clim_funcs = data.frame()

for(i in 1:nrow(fx2)){
  
  print(fx2$covar[i])
  load(fx2$filepath[i])
  func = extractRandomINLA(mod_i$summary.random[[3]], effect_name = fx2$covar[i])
  func$deltawaic = fx2$deltawaic[i]
  func$deltadic = fx2$deltadic[i]
  func$deltaLS = fx2$deltaLS[i]
  clim_funcs = rbind(clim_funcs, func)
  
}

clim_funcs$clim = unlist(lapply(strsplit(clim_funcs$effect, "_"), "[", 1))
clim_funcs$lag = unlist(lapply(strsplit(clim_funcs$effect, "_"), "[", 2))

write.csv(clim_funcs, "./output/model_outputs/climunivar/fitted_functions/clim_functions.csv", row.names=FALSE)

# plot
clim_funcs %>%
  dplyr::filter(clim == "spei1") %>%
  ggplot() + 
  geom_ribbon(aes(value, ymin=exp(lower), ymax=exp(upper), fill=deltawaic), alpha=0.4) +
  geom_line(aes(value, exp(mean))) + 
  geom_hline(yintercept=1, lty=2) +
  facet_wrap(~lag, nrow=1) + 
  theme_minimal() +
  scale_fill_gradientn(colors=viridisLite::mako(200))
  

