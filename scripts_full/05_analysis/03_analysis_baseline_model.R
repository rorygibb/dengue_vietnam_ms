

# ================= Visualise fitted spatial and seasonal effects of baseine model ===============

# Script produces Supp. Figure X: baseline model fitted ranefs

PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

source("./scripts/04_modelling/00_inla_setup_functions_r4.R")
source("./scripts/04_modelling/00_plot_themes.R")



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





# ================= Examine random effects from baseline model ===============

# 1. Offset
# 2. Monthly effect by province
# 3. District effect (BYM2) by year (dengue year)

# read model (stored externally as very large file)
load("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/output/model_outputs/results_23/sociounivar/models/sociouni_nb_model_baseline.R")
m1 = mod_i

# seasonal effect
s = extractRandomINLA(m1$summary.random$month, effect_name = "district", transform=TRUE) %>%
  dplyr::left_join(dd %>% dplyr::select(provincex, province, region) %>% distinct() %>% dplyr::rename("group"=provincex))

# spatiotemporal effects
bym = extractRandomINLA(m1$summary.random$polyid, effect_name = "bym", model_is_bym = TRUE) %>%
  dplyr::left_join(data.frame(group=1:23, year=1998:2020)) %>%
  dplyr::left_join(dd %>% dplyr::select(polyid, areaid, district, province, year_useable_from) %>% distinct(), by=c("value"="polyid")) 
#bym$mean[ bym$year_useable_from < bym$year ] = NA
bym = bym[ -which(bym$year < bym$year_useable_from), ]






# =============== Seasonal by geography =================

xy = cbind(shp[ , c("areaprovin"), drop=F] %>% sf::st_drop_geometry(), st_coordinates(st_centroid(shp))) %>%
  dplyr::group_by(areaprovin) %>%
  dplyr::summarise(y = mean(Y)) %>%
  dplyr::arrange(desc(y))

region_facs = 
  c("Mekong River Delta",
    "Southeast",
    "Central Highlands", 
    "South Central Coast", 
    "North Central",
    "Red River Delta", 
    "Northwest",
    "Northeast")

# set dengue month
mm = as.vector(lubridate::month(c(5:12, 1:4), label=TRUE, abbr=TRUE))
s = s %>%
  left_join(
    data.frame(value=c(1:12), month = mm)
  ) %>%
  dplyr::mutate(month = factor(month, levels=mm, ordered=TRUE))

p1 = s %>% 
  dplyr::mutate(province = factor(province, levels=rev(xy$areaprovin), ordered=TRUE),
                value = factor(value, levels=1:12, ordered=TRUE),
                region = factor(region, levels=rev(region_facs), ordered=TRUE)) %>%
  ggplot() + 
  ggridges::geom_density_ridges(aes(month, province, height=mean, group=province, fill=region), 
                                stat="identity", scale=3, alpha=0.5, size=0.15, col=NA) + 
  ggridges::geom_density_ridges(aes(month, province, height=mean, group=province, fill=region),
                                stat="identity", scale=3, alpha=0.4, size=0.2, col="grey10") +
  theme_classic() + 
  xlab("Month") + 
  ylab("Province") +
  #scale_fill_viridis_d(option="viridis", end=0.95, direction=-1) +
  scale_fill_discrete(type = rev(met.brewer(name="Archambault", n=8))) +
  theme(legend.position="none", 
        axis.title = element_text(size=13), 
        axis.text.x = element_text(size=12))



# ============== Visualise map of geography of regions ==============

p2 = shp %>%
  dplyr::left_join(
    dd %>% dplyr::select(province, region) %>% distinct(), by=c("areaprovin"="province")
  ) %>%
  dplyr::mutate(
    region = factor(region, levels=rev(region_facs), ordered=TRUE)
  ) %>%
  dplyr::filter(!is.na(region)) %>%
  ggplot() + 
  geom_sf(aes(fill=region), color=NA) + 
  geom_sf(data=shp_prov, fill=NA, color="grey60", alpha=0.5, size=0.25) +
  maptheme + 
  scale_fill_discrete(type = rev(met.brewer(name="Archambault", n=8))) +
  theme(legend.position=c(0.01, 0.5), legend.title=element_blank(), legend.text=element_text(size=9.5),
        legend.background = element_blank())

  


# ============ visualise map of ranefs =================

ch = bym %>%
  dplyr::filter(component == "uv_joint") %>%
  dplyr::group_by(areaid) %>%
  dplyr::summarise(mean = mean(mean, na.rm=TRUE))
lims = c(-max(abs(ch$mean)), max(abs(ch$mean)))

p3 = shp %>%
  dplyr::left_join(ch) %>%
  ggplot() + 
  geom_sf(aes(fill=mean), col=NA) + 
  geom_sf(data=shp_prov, fill=NA, color="grey60", alpha=0.5, size=0.25) +
  #scale_fill_gradientn(colors=rev(colorRampPalette(RColorBrewer::brewer.pal(11, "BrBG"))(200)), limits=lims, na.value="grey97", name="u+v") + 
  scale_fill_gradientn(colors=rev(colorRampPalette(MetBrewer::met.brewer("Benedictus", 11))(200)), na.value="white", limits=lims, name="u+v") +
  maptheme + 
  theme(legend.position=c(0.01, 0.5), legend.title=element_text(size=12), legend.text=element_text(size=9.5),
        legend.background = element_blank())




# ================ combine ================

# v2
pc1 = gridExtra::grid.arrange(p2, p3, ncol=1)
pc = gridExtra::grid.arrange(p1, pc1, ncol=2, widths=c(1.3, 1))
pc = ggpubr::as_ggplot(pc)  +
  cowplot::draw_plot_label(label = c("a", "b", "c"), 
                           fontface = "bold", size = 22.5, 
                           x = c(0.03, 0.585, 0.585), y = c(0.99, 0.99, 0.48))
ggsave(pc, file="./output/figures/SuppFigure_BaselineModel.jpg", device="jpg", units="in", width=10, height=7.6, dpi=300)




# # =============== visualise interannual ranefs variation ================
# 
# dy = data.frame(year = 1998:2020) %>%
#   dplyr::mutate(dengue_year = paste(year, year+1, sep="-"))
# 
# lims = c(-8, 8)
# ranefs_dyn = shp %>%
#   full_join(
#     bym %>%
#       dplyr::filter(component == "uv_joint") %>%
#       dplyr::select(areaid, mean, year)
#   ) %>%
#   dplyr::mutate(
#     mean = replace(mean, mean < lims[1], lims[1]),
#     mean = replace(mean, mean > lims[2], lims[2])
#   ) %>%
#   dplyr::filter(!is.na(year)) %>%
#   dplyr::left_join(dy) %>%
#   ggplot() +
#   geom_sf(aes(fill=mean), col=NA) + 
#   geom_sf(data=shp_vt, fill=NA, color="grey60", alpha=0.5, size=0.25) +
#   scale_fill_gradientn(colors=rev(colorRampPalette(RColorBrewer::brewer.pal(11, "BrBG"))(200)), limits=lims, na.value="grey97", name="u+v") + 
#   maptheme + 
#   theme(legend.title=element_text(size=12), legend.text=element_text(size=10)) +
#   facet_wrap(~dengue_year, nrow=3)
# ggsave(ranefs_dyn, file="./output/figures/SuppFigure_SpatiotemporalRanefs.jpg", device="jpg", units="in", width=14, height=9, dpi=600)
# 
# 
