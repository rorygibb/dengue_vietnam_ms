

# ====================== Analyse full model =======================

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
models_dir = "C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/output/model_outputs/results_23/speilags_south/"



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




# ================= model fit metrics ================

# get information criteria for models with different lags of SPEI
ff = list.files(paste(models_dir, "fitmetrics", sep=""), full.names=TRUE, pattern=".csv")
fx = do.call(rbind.data.frame, lapply(ff, read.csv)) %>%
  dplyr::mutate(filename = ff) %>%
  dplyr::arrange(dic) 



# ================= hazard functions ==================

clim_effects = data.frame()
fx = fx %>% 
  dplyr::mutate(deltawaic = waic - waic[ modname  == "baseline" ],
                deltadic = dic - dic[ modname  == "baseline" ]) %>%
  dplyr::filter(candidate %in% c("spei1", "precip", "spei6"))

for(i in 1:nrow(fx)){
  
  print(paste("Processing", i, "of", nrow(fx), sep=" "))

  ith = fx[i, ]
  mm = list.files(paste(models_dir, "models", sep=""), pattern=".R", full.names=TRUE)
  mm = mm[ grep(ith$model_filename, mm) ]
  load(mm)

  nlf = extractRandomINLA(mod_i$summary.random[[ grep("spei|precip", names(mod_i$summary.random)) ]], effect_name = ith$fx) %>%
    dplyr::mutate(rank=i, 
                  deltawaic=ith$deltawaic, 
                  deltadic = ith$deltadic,
                  clim_var = ith$candidate,
                  lag = strsplit(ith$fx, "_")[[1]][2])
  clim_effects = rbind(clim_effects, nlf)

}

# viz
xlims = range(clim_effects$value[ grep("spei", clim_effects$clim_var) ])
ylims = exp(range(c(clim_effects$lower[ grep("spei", clim_effects$clim_var) ], clim_effects$upper[ grep("spei", clim_effects$clim_var) ])))

spei_plot <- clim_effects %>%
  dplyr::filter(grepl("spei", clim_var)) %>%
  dplyr::filter(lag != "0m") %>%
  #dplyr::filter(value > qq[1] & value < qq[2]) %>%
  dplyr::filter(clim_var != "spei12") %>%
  # dplyr::mutate(lag = paste(substr(lag, 1, nchar(lag)-1), "months", sep=" "),
  #               clim_var = toupper(clim_var),
  #               clim_var = paste(substr(clim_var, 1, 4), "-", substr(clim_var, 5, 5), sep=""),
  #               clim_var = factor(clim_var, levels=c("SPEI-1", "SPEI-3", "SPEI-6", "SPEI-12"), ordered=TRUE)) %>%
  dplyr::mutate(lag = paste(substr(lag, 1, nchar(lag)-1), "months", sep=" "),
                #clim_var = toupper(clim_var),
                #clim_var = paste(substr(clim_var, 1, 4), "-", substr(clim_var, 5, 5), sep=""),
                clim_var = replace(clim_var, clim_var == "spei1", "SPEI-1\n(short\ntimescale)"),
                clim_var = replace(clim_var, clim_var == "spei6", "SPEI-6\n(long\ntimescale)"),
                clim_var = replace(clim_var, clim_var == "spei12", "SPEI-12\n(long\ntimescale)"),
                clim_var = factor(clim_var, levels=c("SPEI-1\n(short\ntimescale)", "SPEI-6\n(long\ntimescale)"), ordered=TRUE)) %>%
  ggplot() + xlab("Standardised Precipitation Evapotranspiration Index (SPEI)") + ylab("Dengue relative risk") +
  # geom_rect(aes(xmin=-2, xmax=-1, ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey95", alpha=0.1) +
  # geom_rect(aes(xmin=xlims[1], xmax=-2, ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey90", alpha=0.1) +
  geom_rect(aes(xmin=-2, xmax=2, ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey92", alpha=0.1) +
  geom_rect(aes(xmin=-1, xmax=1, ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey86", alpha=0.1) +
  # geom_rect(aes(xmin=1, xmax=2, ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey95", alpha=0.1) +
  # geom_rect(aes(xmin=2, xmax=xlims[2], ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey90", alpha=0.1) +
  geom_hline(yintercept=c(2), col="grey92", size=0.3) +
  geom_hline(yintercept=1, lty=2, col="grey50", size=0.4) +
  geom_vline(xintercept=0, lty=1, col="grey80") +
  geom_ribbon(aes(value, ymin=exp(lower), ymax=exp(upper), fill=deltawaic), alpha=0.7) +
  geom_line(aes(value, exp(median)), col="black", size=0.7) +
  lemon::facet_rep_grid(clim_var~lag) +
  theme_classic() +
  #scale_y_continuous(position="right") +
  # scale_fill_viridis_c(end=0.8, name=expression(paste(Delta, "WAIC"))) +
  # scale_color_viridis_c(end=0.8, name=expression(paste(Delta, "WAIC"))) +
  scale_fill_gradientn(colors=viridisLite::mako(200)[1:190], name=expression(paste(Delta, "WAIC"))) +
  scale_color_gradientn(colors=viridisLite::mako(200)[1:190], name=expression(paste(Delta, "WAIC"))) +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.line = element_line(color="grey44", size=0.35),
        axis.ticks = element_line(color="grey44", size=0.35),
        #panel.grid.major.y = element_line(color="grey93", size=0.5),
        # panel.grid.minor.y = element_line(color="grey93", size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle=0, size=16),
        axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        strip.text.x = element_text(size=14),
        legend.title = element_text(size=16),
        legend.text = element_text(size=11),
        legend.position = "right")

ggsave(spei_plot, file="./output/figures/SuppFigure_SPEIlags.jpg", units="in", width=13, height=4.2, dpi=600, device="jpg")
