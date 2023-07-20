
# ==========================================================================================#
#                                                                                           #
#           EXAMPLE PIPELINE 4: VISUALISES DENGUE TRENDS FOR SUBSET OF PROVINCES            #
#                                                                                           #
# ==========================================================================================#


# project root and dependencies
PATH = dirname(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(PATH)
pacman::p_load("dplyr", "raster", "rgdal", "sf", "ecmwfr", "stringr", "ggplot2", "lubridate", "magrittr", "vroom")

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))

# theme for mapping
maptheme = theme_classic() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.ticks = element_blank(),
        plot.title = element_text(hjust=0.5, size=12),
        legend.title = element_text(size=10), 
        strip.background = element_blank())



# ================= key objects =================

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

# add polygon area, lat lon, region information
#shp$area_km2 = as.vector(st_area(shp) / 10^6)
shp = cbind(shp, as.data.frame(st_coordinates(st_centroid(shp))) %>% dplyr::rename("longitude"=1, "latitude"=2))
dd = left_join(dd, shp[ , c("areaid", "latitude", "longitude")] %>% st_drop_geometry())
shp = left_join(shp, dd[ !duplicated(dd$areaid), c("areaid", "region1", "region2", "region3") ])

# subset to provinces
shp = shp[ shp$areaid %in% dd$areaid, ]




# ============== HEATMAP ===============

# latitude order
xy = cbind(shp[ , c("areaid"), drop=F] %>% sf::st_drop_geometry(), st_coordinates(st_centroid(shp))) %>%
  dplyr::arrange(Y)

# tmean order
# xy = dd %>% 
#   dplyr::select(areaid, tmean_annualmean) %>%
#   dplyr::group_by(areaid) %>%
#   dplyr::summarise(tmean = mean(tmean_annualmean)) %>%
#   dplyr::arrange(desc(tmean))

heatmap = dd %>%
  dplyr::mutate(date=as.Date(date), 
                areaid = factor(areaid, levels=xy$areaid, ordered=TRUE),
                province = factor(province, levels=c("Ha Noi", "Dak Lak", "Khanh Hoa", "Dong Nai"), ordered=TRUE)
  ) %>% 
  ggplot() + 
  geom_tile(aes(x = date, y=district, fill=log(incidence+1)), width=31) + 
  scale_fill_gradientn(colors=viridisLite::turbo(200), na.value="grey60", name="Dengue\nincidence\n(log)") +
  theme_minimal() + 
  theme(axis.text.y = element_text(size=6), 
        axis.text.x = element_text(size=13), 
        panel.grid = element_blank(),
        legend.title = element_text(size=13), 
        legend.text = element_text(size=12.5),
        strip.text = element_text(size=14),
        axis.title = element_text(size=14), 
        axis.ticks.y = element_blank()) +
  xlab("Month") + ylab("District") +
  facet_wrap(~province, ncol=1, scales="free_y")

ggsave(heatmap, file="./output/figures_example/DengueHeatmap.jpg", device="jpg", units="in", width=7, height=8, dpi=600, scale=0.90)



# ========================= MAPS =================================

# calculate yearly incidence (excluding Cam My in 2019 because extremely high so skews viz)
dd1 = dd[ -which(dd$district == "Cam My" & dd$year == 2019), ]

dd1 = dd1 %>%
  dplyr::filter(!is.na(cases)) %>%
  dplyr::group_by(year, areaid) %>%
  dplyr::summarise(cases = sum(cases, na.rm=TRUE),
                   pop = head(population_census, 1))
dd1$incidence = ( dd1$cases / dd1$pop ) * 100000

# calculate average over study period
dd1 = dd1 %>% dplyr::group_by(areaid) %>% dplyr::summarise(incidence = mean(incidence, na.rm=TRUE))

shpx = shp %>% left_join(dd1)

# individual provinces
pp = "Ha Noi"
p1 = ggplot() + 
  geom_sf(data=shpx[ shpx$areaprovin == pp, ], aes(fill=log(incidence+1)), col=NA) + 
  #geom_sf(data=shp_prov[ shp_prov$provincena %in% pp, ], fill=NA, col="grey30", alpha=0.8, size=0.4) + 
  scale_fill_gradientn(colors = rev(viridisLite::mako(200)), name="Mean\nannual\nincidence\n(log)") +
  maptheme +
  facet_wrap(~areaprovin, nrow=1) + 
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        strip.text = element_text(size=13),
        legend.position="right")

pp = "Dak Lak"
p2 = ggplot() + 
  geom_sf(data=shpx[ shpx$areaprovin == pp, ], aes(fill=log(incidence+1)), col=NA) + 
  #geom_sf(data=shp_prov[ shp_prov$provincena %in% pp, ], fill=NA, col="grey30", alpha=0.8, size=0.4) + 
  scale_fill_gradientn(colors = rev(viridisLite::mako(200)), name="Mean\nannual\nincidence\n(log)") +
  maptheme +
  facet_wrap(~areaprovin, nrow=1) + 
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        strip.text = element_text(size=13),
        legend.position="right")

pp = "Khanh Hoa"
p3 = ggplot() + 
  geom_sf(data=shpx[ shpx$areaprovin == pp, ], aes(fill=log(incidence+1)), col=NA) + 
  #geom_sf(data=shp_prov[ shp_prov$provincena %in% pp, ], fill=NA, col="grey30", alpha=0.8, size=0.4) + 
  scale_fill_gradientn(colors = rev(viridisLite::mako(200)), name="Mean\nannual\nincidence\n(log)") +
  maptheme +
  facet_wrap(~areaprovin, nrow=1) + 
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        strip.text = element_text(size=13),
        legend.position="right")

pp = "Dong Nai"
p4 = ggplot() + 
  geom_sf(data=shpx[ shpx$areaprovin == pp, ], aes(fill=log(incidence+1)), col=NA) + 
  #geom_sf(data=shp_prov[ shp_prov$provincena %in% pp, ], fill=NA, col="grey30", alpha=0.8, size=0.4) + 
  scale_fill_gradientn(colors = rev(viridisLite::mako(200)), name="Mean\nannual\nincidence\n(log)") +
  maptheme +
  facet_wrap(~areaprovin, nrow=1) + 
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        strip.text = element_text(size=13),
        legend.position="right")

pc = gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)
ggsave(pc, file="./output/figures_example/ProvincesMap.jpg", device="jpg", units="in", width=8, height=8, dpi=600, scale=0.90)


