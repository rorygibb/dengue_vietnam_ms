

# ======================= Visualise spatial, temporal and province-level dengue trends ====================

# Script produces MS Figure 1, and Supp. Figure 1 (dengue trends)

# project root and dependencies
PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
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



# ============ Calculate general annual estimates ==========

dd %>%
  dplyr::group_by(areaid, year) %>%
  dplyr::filter(!is.na(cases)) %>%
  dplyr::summarise(region3 = head(region3, 1), cases = sum(cases), population = head(population_census, 1)) %>%
  dplyr::group_by(region3, year) %>%
  dplyr::summarise(cases = sum(cases), population = sum(population)) %>%
  dplyr::filter(year >= 2002) %>%
  dplyr::mutate(incidence = (cases/population)*100000) %>%
  ggplot() + 
  geom_line(aes(year, cases)) +
  facet_wrap(~region3)


dd %>%
  dplyr::group_by(areaid, year) %>%
  dplyr::filter(!is.na(cases)) %>%
  dplyr::summarise(region3 = head(region3, 1), cases = sum(cases), population = head(population_census, 1)) %>%
  dplyr::group_by(region3, year) %>%
  dplyr::summarise(cases = sum(cases), population = sum(population)) %>%
  dplyr::filter(year >= 2002) %>%
  dplyr::mutate(incidence = (cases/population)*100000) %>%
  dplyr::group_by(region3) %>%
  dplyr::mutate(cases_rm = data.table::frollmean(cases, 5, align="left"),
                incidence_rm = data.table::frollmean(incidence, 5, align="left")) %>%
  ggplot() + 
  geom_line(aes(year, incidence_rm)) +
  facet_wrap(~region3)


ss = dd %>%
  dplyr::group_by(areaid, year) %>%
  dplyr::filter(!is.na(cases)) %>%
  dplyr::summarise(region3 = head(region3, 1), cases = sum(cases), population = head(population_census, 1)) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(cases = sum(cases), population = sum(population)) %>%
  dplyr::filter(year >= 2002) %>%
  dplyr::mutate(incidence = (cases/population)*100000) %>%
  dplyr::mutate(cases_rm = data.table::frollmean(cases, 5, align="right"),
                incidence_rm = data.table::frollmean(incidence, 5, align="right"))



# ============= 1. Time series visualisation of dengue within districts stratified by province =================

# For supplementary material
# n.b. Cam My excluded because of extremely high incidence affecting the ylims

# prep and plot
lat = as.data.frame(st_coordinates(st_centroid(shp)))
lat$areaid = shp$areaid; lat$province = shp$areaprovin
lat = lat %>% dplyr::group_by(province) %>% dplyr::summarise(lat=mean(Y)) %>% dplyr::arrange(desc(lat))
dd2 = dd; dd2$province = factor(dd2$province, levels=lat$province, ordered=TRUE)
p1 = ggplot(dd2[ dd2$district != "Cam My", ]) + 
  geom_line(aes(date, incidence, group=areaid), col="darkblue", alpha=0.4) + 
  facet_wrap(~province, scales="free_y", ncol=7) +
  theme_minimal() + 
  theme(legend.position="none",
        strip.text = element_text(size=15.5),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=10.5),
        axis.title = element_text(size=22)) + 
  xlab("Month") + 
  ylab("Dengue incidence (cases per 100,000 inhabitants)")
ggsave(p1, file="./output/figures/plots/DengueIncidence_PerProvince_Facet.jpg", device="jpg", width=16, height=13, units="in", dpi=300)




# =========== 2. Dengue heatmap showing how transmission settings vary by geography ===============

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
  dplyr::filter(district != "Cam My") %>%
  dplyr::mutate(date=as.Date(date), 
                areaid = factor(areaid, levels=xy$areaid, ordered=TRUE),
                # rm = (region3 == "South" & date < as.Date("2001-01-01")) | (region3 == "Central" & date < as.Date("2000-01-01")) |  (region3 == "North" & date < as.Date("1998-01-01")),
                # incidence = replace(incidence, rm==TRUE, NA),
                region = factor(region, 
                                levels=rev(c("Mekong River Delta", "Southeast", "South Central Coast",
                                         "Central Highlands", "North Central", "Red River Delta", "Northwest", "Northeast")))
                ) %>% 
  ggplot() + 
  geom_tile(aes(x = date, y=areaid, fill=log(incidence+1)), width=31) + 
  scale_fill_gradientn(colors=viridisLite::turbo(200), na.value="grey60", name="Dengue\nincidence\n(log)") +
  theme_minimal() + 
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(size=13), 
        panel.grid = element_blank(),
        legend.title = element_text(size=13), 
        legend.text = element_text(size=12.5),
        strip.text = element_text(size=14),
        axis.title = element_text(size=14), 
        axis.ticks.y = element_blank()) +
  xlab("Month") + ylab("District (ordered from low to high latitude)") +
  facet_wrap(~region, ncol=1, scales="free_y")

ggsave(heatmap, file="./output/figures/SuppFigure_DengueIncidenceHeatmap.jpg", device="jpg", units="in", width=9, height=13, dpi=600, scale=0.90)




# ================= 3. Plot the long-term geography of dengue incidence and geographical expansion =====================

# calculate yearly incidence
dd1 = dd %>%
  dplyr::filter(!is.na(cases)) %>%
  dplyr::group_by(year, areaid) %>%
  dplyr::summarise(cases = sum(cases, na.rm=TRUE),
                   pop = head(population_census, 1))
dd1$incidence = ( dd1$cases / dd1$pop ) * 100000



# ======= 2.1 Map incidence in 5-year epochs ========

# map incidence in 5 year epochs
epochs = data.frame(year = 1998:2020,
                    epoch=c(rep("1998-2003", length(1998:2003)), rep("2004-2009", length(2004:2009)), 
                            rep("2010-2014", length(2010:2014)), rep("2015-2020", length(2015:2020))))

# mean annual incidence across years
dx = left_join(dd1, epochs) %>%
  dplyr::filter(!is.na(incidence)) %>%
  dplyr::group_by(areaid, epoch) %>%
  dplyr::summarise(incidence = mean(incidence, na.rm=TRUE))
shpx = full_join(shp, dx, by="areaid")
shpx = shpx[ !is.na(shpx$epoch), ]

colScale = colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(200)
colScale = colorRampPalette(RColorBrewer::brewer.pal(9, "PuBuGn"))(200)

# map figure
p2 = ggplot() + 
  geom_sf(data=shpx, aes(fill=log(incidence+1)), col=NA) + 
  geom_sf(data=shp_prov, fill=NA, col="grey70", alpha=0.2, size=0.2) + 
  geom_sf(data = shp_vt, fill=NA, col="grey20", size=0.2) +
  scale_fill_gradientn(colors = colScale, name="Mean\nannual\nincidence\n(log)") +
  maptheme +
  facet_wrap(~epoch, nrow=1) + 
  theme(legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        strip.text = element_text(size=16),
        legend.position="right")



# ============== map geography of slopes =================

# use linear regression to estimate slope of change in log incidence over time

# calculate number of nonzero obs
dd2 = dd1 %>% 
  dplyr::arrange(areaid, year) %>%
  dplyr::group_by(areaid) %>%
  dplyr::mutate(leading_zeros = cumsum(cases)==0,
                num_nonzero = sum(cases > 0)) %>%
  ungroup()

# calculate slope and significance of log incidence (must have at least 3 non-zero incidence years)
# visualised at significance level of p<0.01
slopes = dd2 %>%
  dplyr::filter(num_nonzero >= 3 ) %>%
  dplyr::arrange(areaid, year) %>%
  #dplyr::filter(year != 2019) %>%
  dplyr::group_by(areaid) %>%
  dplyr::mutate(loginc = log(incidence+1),
                epoch = ifelse(year %in% 1998:2010, "One", "Two")) %>%
  dplyr::summarise(slope = lm(loginc ~ year)$coefficients[2],
                   pval = coef(summary(lm(loginc ~ year)))[ 2, 4],
                   sig = pval <= 0.05)

# offshores
shp_off = sf::st_read("./data/shapefiles/dengue_districts_shapefile.shp") %>%
  dplyr::filter(areaid %in% offshore_areas) 
ext2 = extent(shp_off)
ymin(ext2) = ymin(extent(shp))
shp_off = sf::st_crop(shp_off, ext2)

# big cities
munip = shp_prov %>%
  dplyr::filter(provincena %in% c("Ha Noi", "Hai Phong", "TP. Ho Chi Minh", "Can Tho", "Da Nang"))
munip_lb = shp %>%
  dplyr::filter(areanameen %in% c("Gia Lam", "Hai An", "Quan 9", "Son Tra", "Binh Thuy"))

lims = (exp(slopes$slope)-1)*100
lims = c(-max(abs(lims), na.rm=TRUE), max(abs(lims), na.rm=TRUE))

# visualise % change
p3 = shp %>%
  dplyr::left_join(
    slopes
  ) %>%
  dplyr::mutate(slope = replace(slope, sig==FALSE, NA)) %>%
  dplyr::mutate(facet = "facet") %>%
  ggplot() + 
  geom_sf(aes(fill= (exp(slope)-1)*100 ), col=NA) + 
  geom_sf(data=shp_off, fill=NA, col="grey20", size=0.2) + 
  geom_sf(data=shp_prov, fill=NA, col="grey70", alpha=0.2, size=0.2) + 
  geom_sf(data = shp_vt, fill=NA, col="grey20", size=0.2) +
  geom_sf(data = munip, fill=NA, col="grey10", size=0.2) +
  ggsflabel::geom_sf_text_repel(data = munip_lb, aes(label = areaprovin), force = 100, nudge_x = 4, seed = 10, col="black", segment.colour="grey20", segment.size=0.3, size=4.7) +
  maptheme + 
  facet_wrap(~facet) +
  #scale_fill_gradient2(high=scales::muted("red"), low=scales::muted("blue"), na.value="white", midpoint=1, name="Incidence\nslope\n(% change\nper year)") + 
  scale_fill_gradientn(colors=rev(pals::brewer.rdylbu(200)), na.value="white", limits=lims, name="Incidence\nslope\n(% change\nper year)", breaks=c(-40, -20, 0, 20, 40)) +
  theme(legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        strip.text = element_text(size=16, color="white"),
        legend.position=c(0.82, 0.5))




# ================== combine into map figure =====================

pc = gridExtra::grid.arrange(p2, p3, nrow=1, widths=c(2, 1))
pc = ggpubr::as_ggplot(pc)  +
  cowplot::draw_plot_label(label = c("a", "b"), 
                           fontface = "bold", size = 30, 
                           x = c(0, 0.66), y = c(0.89, 0.89))
ggsave(pc, file="./output/figures/Figure1_DengueTrends.jpg", device="jpeg", width=17, height=7.5, units="in", dpi=900)


