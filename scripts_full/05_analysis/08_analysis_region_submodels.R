
# ================== View combined results of region-wise sub-models (how robust are relationships in different areas) ================

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



# ============= Models ============

# n.b. models are large files so stored externally from GitHub
# model parameters are stored in the "./output/model_outputs/regional_fullmodel" folder

# southern
load(paste(models_dir, "fullmodel_southern/models/fullmod_nb_model_1.R", sep=""))
mod_s = mod_i

# northern
load(paste(models_dir, "fullmodel_northern/models/fullmod_nb_model_1.R", sep=""))
mod_n = mod_i



# =============== Extract ranefs and fixefs =================

ranefs = data.frame()
byms = data.frame()
fixefs = data.frame()

for(region in c("Southern", "Northern")){

  if(region == "Southern"){ mod_i = mod_s }
  if(region == "Northern"){ mod_i = mod_n }
  
  # ranefs
  rr = mod_i$summary.random
  rf = lapply(rr[3:length(rr)], extractRandomINLA, effect_name = "x", transform=TRUE)
  for(i in 1:length(rf)){
    rf[[i]]$effect = names(rr[ 3:length(rr)][i] )
  }
  rf = do.call(rbind.data.frame, rf) %>%
    dplyr::mutate(region = region)
  
  # byms
  byx = extractRandomINLA(mod_i$summary.random$polyid, effect_name="bym", model_is_bym = TRUE) %>%
    dplyr::mutate(region = region) 
  byms = rbind(byms, byx)
  
  # fixefs
  fx = extractFixedINLA(mod_i, model_name=region, transform=FALSE) %>%
    dplyr::mutate(region = region)
  
  # combine
  ranefs = rbind(ranefs, rf)
  fixefs = rbind(fixefs, fx)
  
}

# save model parameters
# write.csv(ranefs, "./output/model_outputs/regional_fullmodel/fitted_params/nonlinear_funcs.csv", row.names=FALSE)
# write.csv(fixefs, "./output/model_outputs/regional_fullmodel/fitted_params/fixed_effects.csv", row.names=FALSE)


# ================== visualise ====================

# read model parameters from storage on GitHib
ranefs = read.csv("./output/model_outputs/regional_fullmodel/fitted_params/nonlinear_funcs.csv")
fixefs = read.csv("./output/model_outputs/regional_fullmodel/fitted_params/fixed_effects.csv")

# extract ranefs and create standardised plots
rf = ranefs
effs = unique(rf$effect)
plots = vector("list", length(effs))
plotRF = function(x){
  
  rr=rf[ rf$effect == effs[x], ]
  
  rr = left_join(
    rr,
    data.frame(
      effect = effs,
      effname = c("Hygeinic toilet access\n(indoor/outdoor)", "Piped or drilled well\nwater access", "Tmean (1-month lag)", "SPEI-1 (1-month lag)", "SPEI-6 (5-month lag)"),
      xlab = c("Prop. households", "Prop. households", "Tmean (°C)", "SPEI-1", "SPEI-6"),
      type = c("Infrastructure", "Infrastructure", "Climate", "Climate", "Climate")
    )
  )
  
  px = rr %>%
    dplyr::mutate(region = factor(region, levels=c("Southern", "Northern"), ordered=TRUE)) %>%
    ggplot() + 
    geom_ribbon(aes(value, ymin=lower, ymax=upper, group=region), fill="grey80", alpha=0.3) +
    geom_line(aes(value, mean, group=region, col=region), size=0.9) +
    geom_hline(yintercept=1, lty=2) +
    facet_wrap(~effname) +
    theme_classic() +
    scale_color_viridis_d(option="magma", begin=0.15, end=0.7,name="Region") +
    ylab("Relative risk") + xlab(rr$xlab[1]) +
    theme(axis.text = element_text(size=11),
          axis.title = element_text(size=12), 
          strip.background = element_blank(),
          strip.text = element_text(size=13),
          legend.position="none") 
  
  # add vline for 0 in SPEI demarcating wet/dry
  if(grepl("spei", effs[x])){
    px = px + geom_vline(xintercept=0, lty=1, size=0.25, col="grey70", alpha=0.4)
  }
  
  plots[[x]] = px
}
prf = lapply(1:length(effs), plotRF)

# non intercept fixed effects
p1 = fixefs %>%
  dplyr::filter(param != "Intercept") %>%
  dplyr::left_join(
    data.frame(
      param=c("gravityf_log", "urban_s", "urbanexp10_log", "traffic_kmperinhab_log", "tmean_coolestmonth_s"),
      paramname = c("Gravity\n(log)", "Built-up\nland", "Urban\nexpansion\nrate (log)", "Traffic\nper inhab\n(log)", "Tmean\ncoolest\nmonth")
    )
  ) %>%
  dplyr::mutate(
    paramname = factor(paramname, levels=c("Tmean\ncoolest\nmonth", "Gravity\n(log)", "Traffic\nper inhab\n(log)", "Built-up\nland", "Urban\nexpansion\nrate (log)"), ordered=TRUE),
    #region = factor(region, levels=c("South", "Central", "North"), ordered=TRUE) 
    region = factor(region, levels=c("Southern", "Northern"), ordered=TRUE)
  ) %>%
  ggplot() +
  geom_point(aes(paramname, median, group=region, col=factor(region)), position=position_dodge(width=0.5), size=2.5) +
  geom_linerange(aes(paramname, ymin=lower, ymax=upper, group=region, col=factor(region)), position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0, lty=2) +
  theme_classic() +
  scale_color_viridis_d(option="magma", begin=0.15, end=0.7,name="Subregion") +
  ylab("Linear fixed effect\n(posterior median + 95% CI)") + xlab("") + 
  theme(legend.position=c(0.85, 0.8), 
        legend.background = element_rect(color="grey50", fill=NA),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text = element_text(angle=0, size=11.5),
        axis.title = element_text(size=12))
        
# intercept
p2 = fixefs %>%
  dplyr::filter(param == "Intercept") %>%
  dplyr::mutate(
    param = "Intercept\n \n ",
    region = factor(region, levels=c("Southern", "Northern"), ordered=TRUE)
  ) %>%
  ggplot() +
  geom_point(aes(param, median, group=region, col=factor(region)), position=position_dodge(width=0.5), size=2.5) +
  geom_linerange(aes(param, ymin=lower, ymax=upper, group=region, col=factor(region)), position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0, lty=2) +
  theme_classic() +
  scale_color_viridis_d(option="magma", begin=0.15, end=0.7, name="Subregion") +
  ylab("Linear fixed effect\n(posterior median + 95% CI)") + xlab("") + 
  theme(legend.position="none", 
        legend.background = element_rect(color="grey50", fill=NA),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text = element_text(angle=0, size=11.5),
        axis.title = element_text(size=12))

# visualise map
map = shp %>%
  dplyr::left_join(
    dd %>% dplyr::select(province, region) %>% distinct(), 
    by=c("areaprovin"="province")
  ) %>%
  dplyr::mutate(
    regionx = "Northern",
    regionx = replace(regionx, region %in% c("Southeast", "Mekong River Delta", "South Central Coast", "Central Highlands"), "Southern"),
    regionx = factor(regionx, levels=c("Southern","Northern"), ordered=TRUE) 
  ) %>%
  ggplot() + 
  geom_sf(aes(fill=regionx), col=NA) + 
  maptheme +
  scale_fill_viridis_d(option="magma", begin=0.15, end=0.75, name="Subregion", guide=guide_legend(reverse = TRUE)) +
  theme(legend.position=c(0.1, 0.5)) +
  geom_sf(data=shp_prov, fill=NA, color="grey60", alpha=0.5, size=0.25) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12))


pc1 = gridExtra::grid.arrange(p2,
                              p1 + theme(axis.title.y = element_blank()),
                              prf[[1]], 
                              prf[[2]] + ylab(""), nrow=1, 
                              widths=c(0.35, 1, 0.9, 0.9))

# row 2: climatic
pc2 = gridExtra::grid.arrange(prf[[3]], 
                              prf[[4]] + ylab(""), 
                              prf[[5]] + ylab(""), 
                              map, 
                              nrow=1)

# combine
pc = gridExtra::grid.arrange(pc1, pc2, nrow=2, heights=c(1.1, 1))

# combine 
ggsave(pc, file="./output/figures/SuppFigure_SubRegionModels.png", device="png", dpi=600, width=15, height=8.3, units="in", scale=0.85)

