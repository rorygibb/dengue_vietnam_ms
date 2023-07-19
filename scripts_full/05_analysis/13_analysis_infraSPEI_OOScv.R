

# ============================== Out of sample tests for water/SPEI interaction =====================================

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
models_dir_seas = "C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/output/model_outputs/results_23/infraspei_tempOOS/"
models_dir_st = "C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/output/model_outputs/results_23/infraspei_stempOOS/"



# ================= shapefile and dengue data =================

# districts to be excluded (offshore)
# only ones with substantial dengue cases are in Kien Giang; could be worth exploring including them but thisis sufficient for now
offshore_areas = c(70154, 70339, 70273, 70355, 70698)

# dengue, regions, climate, landuse, connectivity data
dd = read.csv("./output/model_data/ModelData_Dengue_VietAll.csv") %>%
  dplyr::left_join(read.csv("./output/model_data/ModelData_ClimaticRegions.csv")) %>%
  dplyr::left_join(read.csv("./output/model_data/ModelData_SocioEcologicalCovar_VietAll.csv")) %>%
  dplyr::left_join(vroom::vroom("./output/model_data/ModelData_ClimateLags_VietAll.csv.gz", col_types = list(areaid=vroom::col_character(), date=vroom::col_character()))) %>%
  dplyr::mutate(date = as.Date(date)) %>%
  dplyr::filter(region2 %in% c("South", "South Central Coast", "Central Highlands")) %>%
  dplyr::filter(!areaid %in% offshore_areas)

# provincex
dd$provincex = as.integer(as.factor(dd$province))

# districts shapefile for Vietnam
shp = st_read("./data/shapefiles/dengue_districts_shapefile.shp") %>%
  dplyr::filter(areaid %in% dd$areaid)

# provinces shapefile with regions info, and vietnam outline for mapping
shp_prov = st_read("./data/shapefiles/provinces.shp") %>%
  dplyr::filter(provincena %in% dd$province)
st_crs(shp_prov) = st_crs(shp)
shp_prov = st_crop(shp_prov, shp)
rr = read.csv("./data/shapefiles/regions_lookup.csv")
shp_prov = left_join(shp_prov, rr[ , c("provincena", "region")])

# polyid indicator
shpf = shp[ shp$areaid %in% dd$areaid, ]
id_ref = data.frame(areaid = shpf$areaid, polyid = 1:nrow(shpf))
dd = left_join(dd, id_ref)

# add polygon area, lat lon, region information
#shp$area_km2 = as.vector(st_area(shp) / 10^6)
shp = cbind(shp, as.data.frame(st_coordinates(st_centroid(shp))) %>% dplyr::rename("longitude"=1, "latitude"=2))
dd = left_join(dd, shp[ , c("areaid", "latitude", "longitude")] %>% st_drop_geometry())
shp = left_join(shp, dd[ !duplicated(dd$areaid), c("areaid", "region1", "region2", "region3") ])

# three categories (low, normal, high) of piped water and urbanicity
cats = dd %>%
  dplyr::select(areaid, year, water_piped_year, urban_pw) %>%
  distinct()

# terciles of water
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
dd = dd %>%
  dplyr::left_join(cats[ , c(1, 2, 5, 6)])



# ============= read in stOOS and tempOOS files ==============

readOOSFiles = function(dir){
  
  # model outputs
  ff = list.files(paste(dir, "models/", sep=""), full.names=TRUE, pattern=".csv")
  
  # read preds and keep only OOS preds
  readfilex = function(x){ 
    
    print(x)
    
    # read in
    csvx = read.csv(x) %>% dplyr::filter(oos == TRUE) #%>% dplyr::mutate(modelname = mod_deets)
    
    if(!"mean" %in% names(csvx)){ csvx$mean = NA }
    if(!"lower" %in% names(csvx)){ csvx$lower = NA }
    if(!"upper" %in% names(csvx)){ csvx$upper = NA }
    
    return(csvx)
  }
  pp = lapply(ff, readfilex)
  pp = do.call(rbind.data.frame, pp)
  
  # calculated predicted cases
  pp$mean = exp(pp$mean + pp$logpop)
  pp$lower = exp(pp$lower + pp$logpop)
  pp$upper = exp(pp$upper + pp$logpop)
  
  return(pp)
  
}

# 
pp1 = readOOSFiles(models_dir_seas)
pp2 = readOOSFiles(models_dir_st)
pp = rbind(pp1, pp2)

# incompletes
pp = pp %>% dplyr::filter(!holdout_id %in% c(132165, 39574) )


# =============== Calculate summary metrics =================

# 1. Across all observations

sm = pp %>%
  dplyr::filter(!is.na(cases)) %>%
  dplyr::mutate(resid = mean - cases) %>%
  dplyr::group_by(holdout_id, model) %>%
  dplyr::summarise(holdout_type = head(holdout_type, 1), 
                   mae = mean(abs(resid), na.rm=TRUE),
                   rmse = sqrt(mean(resid^2, na.rm=TRUE))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(holdout_id) %>%
  dplyr::mutate(deltamae = mae - mae[ model == "full" ],
                deltarmse = rmse - rmse[ model == "full" ],
                deltamae_baseline = mae - mae[ model == "baseline" ],
                deltarmse_baseline = rmse - rmse[ model == "baseline" ])

fac_levs = data.frame(
  model_name = c("Full (no interactions)", "SPEI excluded", "SPEI-1 excluded", "SPEI-1 (1m) x Water", "SPEI-1 (1m) x Urban", "SPEI-6 excluded", "SPEI-6 (5m) x Water", "SPEI-6 (5m) x Urban"),
  model = c("full", "no_spei", "no_spei1", "spei1water", "spei1urban", "no_spei6", "spei6water", "spei6urban")
)

sm = left_join(sm, fac_levs)
sm$model_name[ sm$model == "baseline" ] = "Baseline"
write.csv(sm, file=paste(dir, "OOSCV_countrywide_summary.csv", sep=""), row.names=FALSE)



# 2. summary metrics at district-level

# mean error for full model
full_mods = pp %>% 
  dplyr::filter(model == "full") %>%
  dplyr::mutate(resid = mean - cases) %>% 
  dplyr::select(areaid, date, holdout_id, holdout_type, resid) %>%
  dplyr::group_by(areaid, holdout_id) %>%
  dplyr::summarise(holdout_type = head(holdout_type, 1),
                   mae_full = mean(abs(resid), na.rm=TRUE),
                   rmse_full = sqrt(mean(resid^2, na.rm=TRUE)))

# mean error for model with spei6_excluded
# mods2 = pp %>% 
#   dplyr::filter(model == "spei1_only") %>%
#   dplyr::mutate(resid = mean - cases) %>% 
#   dplyr::select(areaid, date, holdout_id, holdout_type, resid) %>%
#   dplyr::group_by(areaid, holdout_id) %>%
#   dplyr::summarise(holdout_type = head(holdout_type, 1),
#                    mae_spei6ex = mean(abs(resid), na.rm=TRUE),
#                    rmse_spei6ex = sqrt(mean(resid^2, na.rm=TRUE)))
mods2 = pp %>%
  dplyr::filter(model == "baseline") %>%
  dplyr::mutate(resid = mean - cases) %>%
  dplyr::select(areaid, date, holdout_id, holdout_type, resid) %>%
  dplyr::group_by(areaid, holdout_id) %>%
  dplyr::summarise(holdout_type = head(holdout_type, 1),
                   mae_base = mean(abs(resid), na.rm=TRUE),
                   rmse_base = sqrt(mean(resid^2, na.rm=TRUE)))
full_mods = left_join(full_mods, mods2)

# mean for all models per holdout
per_holdout = pp %>%
  dplyr::mutate(resid = mean - cases) %>%
  dplyr::select(areaid, date, holdout_id, holdout_type, model, resid) %>%
  dplyr::group_by(areaid, holdout_id, model) %>%
  dplyr::summarise(holdout_type = head(holdout_type, 1),
                   mae = mean(abs(resid), na.rm=TRUE),
                   rmse = sqrt(mean(resid^2, na.rm=TRUE)))

# combine
sm_district = full_join(full_mods, per_holdout)
#write.csv(sm_district, file=paste(dir, "OOSCV_districts_summary.csv", sep=""), row.names=FALSE)




# ================ Visualise mean error change across all areas ================

# read in sm
#sm = read.csv(paste(dir, "OOSCV_countrywide_summary.csv", sep=""))

# for water only
smx = sm %>%
  dplyr::filter(!grepl("urban|full", model))
flx = fac_levs %>%
  dplyr::filter(!grepl("urban|full", model))

# rename holdout vards
sm2 = smx %>%
  dplyr::filter(!model %in% c("baseline")) %>%
  dplyr::mutate(
    holdout_type = replace(holdout_type, holdout_type=="kfoldOOSq", "Seasonal"),
    holdout_type = replace(holdout_type, holdout_type=="stempOOS", "Spatiotemporal"),
    holdout_type = factor(holdout_type, 
                          levels=c("Spatiotemporal", "Seasonal"), 
                          ordered=TRUE)
  ) %>%
  dplyr::mutate(model_name = factor(model_name, levels=rev(flx$model_name), ordered=TRUE))

# mean and se
sm_summary = sm2 %>%
  dplyr::group_by(holdout_type, model_name) %>%
  dplyr::summarise(deltamae_mean = mean(deltamae),
                   deltamae_se = plotrix::std.error(deltamae),
                   deltarmse_mean = mean(deltarmse),
                   deltarmse_se = plotrix::std.error(deltarmse))

# delta-MAE
p1 = sm2 %>%
  ggplot() + 
  geom_point(aes(model_name, deltamae, group=holdout_id), size=3.5, pch=21, color="skyblue4", fill="skyblue4", alpha=0.2) +
  geom_point(data=sm_summary, aes(model_name, deltamae_mean), pch=18, color="red", alpha=1, size=3) +
  geom_linerange(data=sm_summary, aes(model_name, ymin=deltamae_mean-1.96*deltamae_se, ymax=deltamae_mean+1.96*deltamae_se), color="red", alpha=1, size=0.7) +
  #geom_rect(aes(xmin=xmin1, xmax=xmax1, ymin=ymin1, ymax=0), alpha=0.1, fill="green") + 
  geom_hline(yintercept=0, lty=2) + 
  theme_bw() + 
  #ylab(expression(paste(Delta, "MAE (out-of-sample) relative to full model"))) +
  ylab("Change in mean prediction error (ΔMAE) from full model") +
  #ylab("Improvement in MAE (out-of-sample)") +
  xlab("Model") + 
  theme(strip.background = element_blank(), 
        strip.text=element_text(size=12),
        axis.text = element_text(size=11), 
        axis.title = element_text(size=12)) + 
  facet_wrap(~holdout_type) +
  coord_flip()




# ======================= District-level analysis ============================

# sm_district
#sm_district = read.csv(paste(dir, "OOSCV_districts_summary.csv", sep=""))

# calculate summary pointwise per model
smx = sm_district %>%
  dplyr::mutate(mae_diff = mae - mae_full,
                rmse_diff = rmse - rmse_full) %>%
  dplyr::group_by(holdout_type, model, areaid) %>%
  dplyr::summarise(mae_diff = mean(mae_diff, na.rm=TRUE),
                   rmse_diff = mean(rmse_diff, na.rm=TRUE))


# 1.  % of districts where predictive error is reduced compared to full model
# could also be framed as relative to baseline

n_areas = n_distinct(smx$areaid)
p2 = smx %>% 
  dplyr::mutate(mae_imp = mae_diff<0, rmse_imp = rmse_diff<0) %>%
  dplyr::group_by(holdout_type, model) %>%
  dplyr::summarise(mae_imp_n = sum(mae_imp),
                   mae_imp = sum(mae_imp)/n_areas, 
                   rmse_imp_n = sum(rmse_imp),
                   rmse_imp = sum(rmse_imp)/n_areas) %>%
  dplyr::filter(!grepl("urban", model)) %>%
  dplyr::left_join(fac_levs) %>%
  dplyr::mutate(model_name = replace(model_name, model=="baseline", "Baseline")) %>%
  dplyr::select(-model) %>%
  reshape2::melt(id.vars = c("holdout_type", "model_name")) %>%
  dplyr::filter(variable == "mae_imp") %>%
  dplyr::mutate(
    holdout_type = replace(holdout_type, holdout_type=="kfoldOOSq", "Seasonal"),
    holdout_type = replace(holdout_type, holdout_type=="stempOOS", "Spatiotemporal"),
    holdout_type = factor(holdout_type, 
                          levels=c("Spatiotemporal", "Seasonal"), 
                          ordered=TRUE)
  ) %>%
  dplyr::filter(model_name != "Full (no interactions)") %>%
  dplyr::mutate(model_name = factor(model_name, levels=rev(c("Baseline",
                                                             "SPEI-1 excluded",
                                                             "SPEI-1 (1m) x Water",
                                                             "SPEI-6 excluded",
                                                             "SPEI-6 (5m) x Water")), ordered=TRUE)) %>%
  ggplot() + 
  geom_point(aes(model_name, value*100), size=3, pch=21, col="black", fill="skyblue4", alpha=0.9) + 
  facet_wrap(~holdout_type) + 
  #geom_hline(yintercept=50, lty=1, size=0.7, col="grey70") +
  theme_bw() + 
  theme(strip.background = element_blank(), 
        strip.text=element_text(size=12),
        axis.text = element_text(size=11), 
        axis.title = element_text(size=12)) +
  coord_flip() +
  ylab("% districts with decreased prediction error\ncompared to full model") + xlab("Model")




#=============== mapping improvements over a full model ===============

#
pw = sm_district %>%
  dplyr::mutate(mae_diff = mae - mae_full) %>%
  dplyr::filter(model=="spei6water") %>%
  dplyr::group_by(areaid, holdout_type) %>%
  dplyr::summarise(mean_diff = mean(mae_diff))

# categories
pw$md_cat = "< -0.25"
pw$md_cat[ pw$mean_diff > -0.25 ] = "-0.1 to -0.25"
pw$md_cat[ pw$mean_diff > -0.1 ] = "0 to -0.1"
pw$md_cat[ pw$mean_diff > 0 ] = "0 to +0.1"
pw$md_cat[ pw$mean_diff > 0.1 ] = "+0.1 to +0.25"
pw$md_cat[ pw$mean_diff > 0.25 ] = "> +0.25"

pw$md_cat = "< -0.4"
pw$md_cat[ pw$mean_diff > -0.4 ] = "-0.1 to -0.4"
pw$md_cat[ pw$mean_diff > -0.1 ] = "0 to -0.1"
pw$md_cat[ pw$mean_diff > 0 ] = "0 to +0.1"
pw$md_cat[ pw$mean_diff > 0.1 ] = "+0.1 to +0.4"
pw$md_cat[ pw$mean_diff > 0.4 ] = "> +0.4"


# 0.25 captures 10% of obs; 0.1 captures 68% of obs
# thresh = 0.25
# pw$test= pw$mean_diff < -thresh | pw$mean_diff > thresh
# table(pw$test) / nrow(pw)

# shapefile
shp = sf::st_read("./data/shapefiles/dengue_districts_shapefile.shp") %>%
  dplyr::filter(areaid %in% pw$areaid)
shp_prov = sf::st_read("./data/shapefiles/provinces.shp") %>%
  dplyr::filter(provincena %in% ddf$province) 
st_crs(shp_prov) = crs(shp)
shp_prov = st_crop(shp_prov, shp)

# viz
p5.2 = shp %>% 
  dplyr::filter(areaid %in% pw$areaid) %>%
  full_join(pw) %>%
  dplyr::filter(holdout_type == "kfoldOOSq") %>%
  dplyr::mutate(holdout_type = "SPEI-6 (5m) x Water\n(Seasonal)") %>%
  #dplyr::mutate(md_cat = factor(md_cat, levels=c("< -0.25", "-0.1 to -0.25", "0 to -0.1", "0 to +0.1", "+0.1 to +0.25", "> +0.25"), ordered=TRUE)) %>%
  dplyr::mutate(md_cat = factor(md_cat, levels=c("< -0.4", "-0.1 to -0.4", "0 to -0.1", "0 to +0.1", "+0.1 to +0.4", "> +0.4"), ordered=TRUE)) %>%
  ggplot() + 
  geom_sf(aes(fill = md_cat), col=NA) + 
  geom_sf(data=shp_prov, fill=NA, col="grey50", alpha=0.5, size=0.25) + 
  #scale_fill_discrete(type=colorRampPalette(c("red", "white", "blue"))(200)) + 
  scale_fill_brewer(palette="BrBG", direction=-1, name=expression(paste(Delta, "MAE"))) + 
  maptheme + 
  facet_wrap(~holdout_type) + 
  theme(legend.position="none", legend.title = element_text(size=13), legend.text = element_text(size=12),
        strip.text=element_text(size=13))

p5.1 = shp %>% 
  dplyr::filter(areaid %in% pw$areaid) %>%
  full_join(pw) %>%
  dplyr::filter(holdout_type == "stempOOS") %>%
  dplyr::mutate(holdout_type = "SPEI-6 (5m) x Water\n(Spatiotemporal)") %>%
  #dplyr::mutate(md_cat = factor(md_cat, levels=c("< -0.25", "-0.1 to -0.25", "0 to -0.1", "0 to +0.1", "+0.1 to +0.25", "> +0.25"), ordered=TRUE)) %>%
  dplyr::mutate(md_cat = factor(md_cat, levels=c("< -0.4", "-0.1 to -0.4", "0 to -0.1", "0 to +0.1", "+0.1 to +0.4", "> +0.4"), ordered=TRUE)) %>%
  ggplot() + 
  geom_sf(aes(fill = md_cat), col=NA) + 
  geom_sf(data=shp_prov, fill=NA, col="grey50", alpha=0.5, size=0.25) + 
  #scale_fill_discrete(type=colorRampPalette(c("red", "white", "blue"))(200)) + 
  scale_fill_brewer(palette="BrBG", direction=-1, name=expression(paste(Delta, "MAE"))) + 
  maptheme + 
  facet_wrap(~holdout_type) + 
  #theme(legend.position="none") + 
  theme(legend.position="left", legend.title = element_text(size=13), legend.text = element_text(size=12),
        strip.text=element_text(size=13))

# combine
#p3 = gridExtra::grid.arrange(p5.1, p5.2, ncol=2, widths=c(0.92, 1.5))
p3 = gridExtra::grid.arrange(p5.1, p5.2, ncol=2, widths=c(1.5, 0.92))




# =========== combine all out of sample tests into a full figure ==================

# combine into full figure
pc = gridExtra::grid.arrange(p1, p2, p3, nrow=3, heights=c(0.9, 1, 1.5))
pc = ggpubr::as_ggplot(pc)  +
  cowplot::draw_plot_label(label = c("a", "b", "c"), 
                           fontface = "bold", size = 22.5, 
                           x = c(0.015, 0.015, 0.015), y = c(1, 0.72, 0.42))
ggsave(pc, file="./output/figures/SuppFigure_InfraSPEI_OOS.jpg", device="jpg", 
       units="in", width=7*0.9, height=9*0.8, dpi=600)
