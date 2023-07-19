

# ============================== ANNUAL TEMPERATURE ANALYSIS =====================================

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
shp_vt = st_read("./data/shapefiles/gadm36_VNM_0.shp") %>%
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


fx1 = fx1 %>%
  dplyr::filter(covar %in% c("baseline", "tmin_annualmean_s", "tmean_annualmean_s", "tmean_coolestmonth_s")) %>%
  dplyr::mutate(
    covar = replace(covar, covar == "baseline", "Baseline"),
    covar = replace(covar, covar == "tmin_annualmean_s", "Tmin (annual mean)"),
    covar = replace(covar, covar == "tmean_annualmean_s", "Tmean (annual mean)"),
    covar = replace(covar, covar == "tmean_coolestmonth_s", "Tmean (coolest month)")
    )

fx1$filepath = paste(models_dir, "sociounivar/models/", fx1$model_filename, sep="")




# ============= calculate ST ranefs for each model ================

# unexplained st var
extract_MAERE = function(model){
  strf = extractRandomINLA(model$summary.random$polyid, effect_name ="", model_is_bym = TRUE) %>%
    dplyr::filter(component=="uv_joint") %>%
    dplyr::left_join(data.frame(group=1:23, year=1998:2020)) %>%
    dplyr::left_join(
      dd %>% dplyr::select(polyid, areaid, year_useable_from) %>% distinct(),
      by=c("value"="polyid")
    ) %>%
    dplyr::mutate(mean = replace(mean, year < year_useable_from, NA)) %>%
    dplyr::filter(group > 4) 
  #mae_re = mean(abs(strf$mean), na.rm=TRUE)
  return(strf)
}

fx1$MAE_re = NA
ranefs = data.frame()

for(i in 1:nrow(fx1)){
  
  load(fx1$filepath[i])
  ranefs_x = extract_MAERE(mod_i) %>% dplyr::mutate(covar = fx1$covar[i])
  fx1$MAE_re[i] = mean(abs(ranefs_x$mean), na.rm=TRUE)
  ranefs = rbind(ranefs, ranefs_x)
  
}



# ==================== Visualise ======================

cov_names = data.frame(
  covar = c("Tmean (coolest month)", "Tmin (annual mean)", "Tmean (annual mean)", "Baseline"),
  covar2 = c("Tmean\n(coolest month)", "Tmin\n(annual)", "Tmean\n(annual)", "Baseline"),
  covar3 = c("Tmean (coolest month)", "Tmin (annual)", "Tmean (annual)", "Baseline")
)

# plot 1 = change in WAIC
p1 = fx1 %>%
  dplyr::arrange(deltawaic) %>%
  dplyr::mutate(covar = Hmisc::capitalize(covar)) %>%
  dplyr::left_join(
    cov_names
  ) %>%
  dplyr::mutate(covar2 = factor(covar2, levels=cov_names$covar2, ordered=TRUE)) %>%
  ggplot() + 
  geom_point(aes(covar2, deltawaic), size=7, pch=21, fill="skyblue4", alpha=0.9) + 
  theme_bw() +
  geom_hline(yintercept=0, lty=2) + 
  ylab(expression(paste(Delta, "WAIC", sep=""))) +
  xlab("Model") + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=13)) +
  coord_flip()

# random effects magnitude
p2 = fx1 %>%
  dplyr::arrange(deltawaic) %>%
  dplyr::mutate(covar = Hmisc::capitalize(covar)) %>%
  dplyr::left_join(
    cov_names
  ) %>%
  dplyr::mutate(covar2 = factor(covar2, levels=cov_names$covar2, ordered=TRUE)) %>%
  ggplot() + 
  geom_point(aes(covar2, MAE_re), size=7, pch=21, fill="skyblue4", alpha=0.9) + 
  theme_bw() +
  geom_hline(aes(yintercept=MAE_re[ covar == "Baseline"]), lty=2) + 
  ylab(expression('MAE'[RE])) +
  xlab("Model") + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=13)) +
  coord_flip()

# plot 3 = overlaid densities
p3 = ranefs %>%
  dplyr::left_join(
    ranefs %>% dplyr::filter(covar=="Baseline") %>% dplyr::select(areaid, year, mean) %>% dplyr::rename("bl"=mean)
  ) %>%
  dplyr::left_join(
    cov_names
  ) %>%
  dplyr::filter(covar != "Baseline") %>%
  ggplot() + 
  geom_vline(xintercept=0, lty=2) + 
  geom_density(aes(x=bl), alpha=0.2, fill="skyblue4", col="black") +
  geom_density(aes(x=mean), alpha=0.2, fill="red", col="black") +
  theme_minimal() + 
  facet_wrap(~covar3, ncol=1) +
  xlab("Random effect posterior mean") + ylab("Density") +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size=11.5),
        axis.title = element_text(size=13), 
        strip.text = element_text(size=13.5))

# plot 4: spatiotemporal change in random effects from including tcoolest
change = ranefs %>%
  dplyr::filter(!is.na(mean)) %>%
  dplyr::filter(covar == "Tmean (coolest month)") %>%
  dplyr::left_join(
    ranefs %>% dplyr::filter(!is.na(mean)) %>% dplyr::filter(covar=="Baseline") %>% dplyr::select(areaid, year, mean) %>% dplyr::rename("bl"=mean)
  ) %>%
  dplyr::group_by(areaid) %>%
  dplyr::summarise(meanval = mean(mean),
                   meanval_bl = mean(bl),
                   modi = mean(abs(mean)), bl=mean(abs(bl))) %>%
  dplyr::mutate(diff = (modi - bl)/bl)
# lims = c(-max(abs(change$diff)), max(abs(change$diff)))
# shp %>% 
#   left_join(change) %>%
#   dplyr::filter(!areanameen %in% c("Truong Sa", "Hoang Sa")) %>%
#   ggplot() + 
#   geom_sf(aes(fill=diff), col=NA) + 
#   scale_fill_gradientn(colors=rev(colorRampPalette(RColorBrewer::brewer.pal(11, "BrBG"))(200)), limits=lims, na.value="grey97") + 
#   maptheme

ch = change %>%
  dplyr::select(1:3) %>%
  dplyr::rename("Baseline"=3, "Tmean (coolest month)"=2) %>%
  reshape2::melt(id.vars=1) 
lims = c(-max(abs(ch$value)), max(abs(ch$value)))
map_plot = shp %>% 
  full_join(ch) %>%
  dplyr::filter(!areanameen %in% c("Truong Sa", "Hoang Sa")) %>%
  dplyr::filter(!is.na(variable)) %>%
  dplyr::mutate(variable = factor(variable, levels=c("Baseline", "Tmean (coolest month)", ordered=TRUE))) %>%
  ggplot() + 
  geom_sf(aes(fill=value), col=NA) +
  geom_sf(data=shp_vt, fill=NA, col="grey20") + 
  scale_fill_gradientn(colors=rev(colorRampPalette(RColorBrewer::brewer.pal(11, "BrBG"))(200)), limits=lims, na.value="grey97", name="u+v") + 
  maptheme +
  theme(legend.position = c(0.2, 0.4), legend.title=element_blank(), legend.text = element_text(size=11)) +
  facet_wrap(~variable)

pc1 = gridExtra::grid.arrange(p3, map_plot, ncol=2, widths=c(0.6, 1))
pc2 = gridExtra::grid.arrange(p1, p2, nrow=2)
pc = gridExtra::grid.arrange(pc2, pc1, widths=c(0.5, 1.2))
pc = ggpubr::as_ggplot(pc)  +
  cowplot::draw_plot_label(label = c("a", "b", "c", "d", "e"), 
                           fontface = "bold", size = 22.5, 
                           x = c(0.115, 0.115, 0.34, 0.56, 0.77), y = c(0.985, 0.48, 0.96, 0.93, 0.93))
ggsave(pc, file="./output/figures/SuppFigure_TemperatureAnnualResults.jpg", device="jpg", width=12,  height=6, units="in", dpi=300, scale=1.1)



