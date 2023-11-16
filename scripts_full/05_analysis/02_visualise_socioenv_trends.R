

# ======================= Visualise distributions, spatial and temporal trends in key covariates ====================

# Script produces MS Figure 2 and Supp. Figures visualising socio-env variables

# project root and dependencies
PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)
pacman::p_load("dplyr", "raster", "rgdal", "sf", "ecmwfr", "stringr", 
               "MetBrewer", "ggplot2", "lubridate", "magrittr", "vroom")

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
shp_vt = st_read("./data/shapefiles/vt_national.shp") %>%
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



# =================== Figure: Visualise socio-environmental covariates at regional scale ===================

# socioeco
sec = dd %>%
  dplyr::select(
    areaid, province, region, region3, year, population_census, popdens_census, pop_propurban_census, urbancover_km2,
    urban_pw, water_piped_year, flushtoilet_any_year, flux_grav1, traffic_milperskm,
    tmin_annualmean, tmin_coolestmonth, tmean_annualmean, precip_total
  ) %>%
  distinct() %>%
  dplyr::mutate(
    water_hh = water_piped_year * population_census,
    flushtoilet_hh = flushtoilet_any_year * population_census
  )

# area
shp$area = as.vector(st_area(shp))/10^6

# combine by region
reg = sec %>%
  left_join(shp %>% dplyr::select(areaid, area) %>% st_drop_geometry()) %>%
  dplyr::group_by(region, year) %>%
  dplyr::summarise(
    population = sum(population_census),
    #prop_urban = sum(population_census * pop_propurban_census) / population,
    area_km2 = sum(area),
    #popdens = population / area_km2,
    popdens = mean(popdens_census),
    urbancover = sum(urbancover_km2),
    water = sum(water_hh) / population,
    flushtoilet = sum(flushtoilet_hh) / population,
    gravity = mean(flux_grav1),
    #traffic_perinhab = sum(traffic_total) / population,
    tmin_coolestmonth = mean(tmin_coolestmonth),
    tmean_annual = mean(tmean_annualmean),
    precip_total = mean(precip_total),
  ) 

# calculate traffic separately bc province level source
traf = sec %>%
  dplyr::select(year, province, region, areaid, traffic_milperskm, population_census) %>%
  dplyr::group_by(province, year) %>%
  dplyr::summarise(traffic_km = head(traffic_milperskm, 1),
                   pop = sum(population_census),
                   region=head(region, 1)) %>%
  dplyr::group_by(region, year) %>%
  dplyr::summarise(
    traffic_km = sum(traffic_km),
    pop = sum(pop)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    traf_perinhab = (traffic_km * 10^6) / pop
  )

reg = left_join(reg, traf %>% dplyr::select(region, year, traf_perinhab))

# varnames
varnames = c("Mean population density\n(inhabitants per km2)", 
             "Flush toilet indoor/outdoor\n(% households)", 
             "Road traffic\n(km per inhabitant)", 
             "Built-up land area\n(km2)", 
             "Piped or drilled well water\n(% households)", 
             "Mean gravity flux",
             "Tmean annual mean\n(°C)",
             "Tmin coolest month\n(°C)", 
             "Precipitation annual\n(millimetres)")

region_facs = 
  c("Mekong River Delta",
    "Southeast",
    "Central Highlands", 
    "South Central Coast", 
    "North Central",
    "Red River Delta", 
    "Northwest",
    "Northeast")

# urb 2019
u19 = dd %>% dplyr::filter(year == 2019) %>% 
  dplyr::select(areaid, pop_propurban_census, urban_pw) %>%
  distinct()

# extra lines showing when data were available
line_dat = data.frame(
  var2 = c(rep("Mean population density\n(inhabitants per km2)", 3),
           rep("Flush toilet indoor/outdoor\n(% households)", 2),
           rep("Piped or drilled well water\n(% households)", 2),
           rep("Mean gravity flux", 3)),
  x_line = c(2000, 2009, 2019, 2009, 2019, 2009, 2019, 2000, 2009, 2019)
) %>%
  dplyr::mutate(var2 = factor(var2, levels=varnames, ordered=TRUE))

p1 = reg %>%
  reshape2::melt(id.vars = 1:2) %>%
  dplyr::filter(!variable %in% c("population", "area_km2")) %>%
  dplyr::left_join(
    data.frame(
      variable = c("popdens", 
                   "flushtoilet",
                   "traf_perinhab", 
                   "urbancover", 
                   "water", 
                   "gravity",
                   "tmean_annual",
                   "tmin_coolestmonth", 
                   "precip_total"),
      var2 = varnames)
  ) %>%
  dplyr::mutate(
    var2 = factor(var2, levels=varnames, ordered=TRUE),
    region = factor(region, levels=region_facs, ordered=TRUE)
  ) %>%
  ggplot() + 
  geom_line(aes(year, value, group=region, col=region), size=0.6, alpha=0.8) + 
  geom_vline(data = line_dat, aes(xintercept=x_line), lty=2, alpha=0.5, size=0.35, col="grey30") + 
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=11),
        axis.title.y = element_blank(),
        axis.text = element_text(size=10),
        axis.title.x = element_text(size=11)) +
  lemon::facet_rep_wrap(~var2, scales="free_y") +
  scale_color_discrete(type = met.brewer(name="Archambault", n=8)) +
  theme(legend.position="none", legend.title=element_blank(),
        legend.text = element_text(size=10.5)) +
  xlab("Year")

shp_prov2 = st_crop(shp_prov, shp)

p2 = shp %>%
  left_join(sec %>% dplyr::select(areaid, region) %>% distinct()) %>%
  dplyr::mutate(
    region = factor(region, levels=region_facs, ordered=TRUE)
  ) %>%
  ggplot() + 
  geom_sf(color=NA, aes(fill=region)) + 
  geom_sf(data=shp_prov2, fill=NA, col="grey50", alpha=0.6, size=0.3) +
  maptheme + 
  scale_fill_discrete(type = met.brewer(name="Archambault", n=8), guide=guide_legend(reverse=TRUE)) + 
  theme(legend.position=c(0.22, 0.45), legend.title=element_blank(), legend.text = element_text(size=11),
        legend.background = element_blank())

pc = gridExtra::grid.arrange(p1, p2, ncol=2, widths=c(1, 0.5))
ggsave(pc, file="./output/figures/Figure2_SocioEnvironmentTrends.jpg", device="jpg", units="in", dpi=600, width=11, height=7, scale=0.9)
ggsave(pc, file="./output/figures/Figure2_SocioEnvironmentTrends.pdf", device="pdf", units="in", width=11, height=7, scale=0.9)





# ===================== Supp: Climatic covariates over space and time =======================

plot_theme = theme(
  axis.text = element_text(size=10.5),
  axis.title = element_text(size=11), 
  strip.text = element_text(size=11),
  strip.background = element_blank(),
  plot.title=element_text(size=15, hjust=0.5)
)

# higher region grouping for visualisation
regionx = data.frame(
  region = 
    c("Mekong River Delta",
      "Southeast",
      "Central Highlands", 
      "South Central Coast", 
      "North Central",
      "Red River Delta", 
      "Northwest",
      "Northeast"),
  regionx = 
    c("South", "South", 
      "Central Highlands", 
      "South Central Coast", 
      "North Central",
      "Red River Delta",
      "Northeast/Northwest", 
      "Northeast/Northwest")
)

region_facs = c("South", 
                "Central Highlands", 
                "South Central Coast", 
                "North Central",
                "Red River Delta",
                "Northeast/Northwest")

# summarised by regeion
region_clim = dd %>%
  left_join(regionx) %>%
  dplyr::group_by(regionx, date) %>%
  dplyr::summarise(tmean = mean(tmean, na.rm=TRUE),
                   spei6 = mean(spei6_0m),
                   spei1 = mean(spei1_0m),
                   precip = mean(precip_0m)) %>%
  dplyr::mutate(regionx = factor(regionx, levels=rev(region_facs), ordered=TRUE))

# temperature by region
p1 = dd %>%
  left_join(regionx) %>%
  # dplyr::group_by(areaid, date) %>%
  # dplyr::summarise(tmean = mean(tmean), 
  #                  region = head(region, 1)) %>%
  dplyr::mutate(regionx = factor(regionx, levels=rev(region_facs), ordered=TRUE)) %>%
  ggplot() + 
  geom_line(aes(date, tmean, group=areaid), col="coral4", alpha=0.15, size=0.1) + 
  geom_line(data=region_clim, aes(date, tmean), col="grey10", size=0.4) + 
  theme_classic() + 
  plot_theme + facet_wrap(~regionx, ncol=1) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=13.5)) + 
  ylab("Tmean (°C)") + xlab("Month") +
  ggtitle("Tmean")

# precip by region
p2 = dd %>%
  left_join(regionx) %>%
  # dplyr::group_by(province, date) %>%
  # dplyr::summarise(tmean = mean(tmean), 
  #                  region = head(region, 1)) %>%
  dplyr::mutate(regionx = factor(regionx, levels=rev(region_facs), ordered=TRUE)) %>%
  ggplot() + 
  geom_line(aes(date, precip_0m, group=areaid), col="coral4", alpha=0.15, size=0.1) + 
  geom_line(data=region_clim, aes(date, precip), col="grey10", size=0.4) + 
  theme_classic() + 
  plot_theme + facet_wrap(~regionx, ncol=1) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=13.5)) + 
  ylab("Precipitation (mm/day)") + xlab("Month") +
  ggtitle("Precipitation")

# spei-6
p4 = dd %>%
  left_join(regionx) %>%
  # dplyr::group_by(province, date) %>%
  # dplyr::summarise(spei6 = mean(spei6_0m),
  #                  spei1 = mean(spei1_0m),
  #                  region = head(region, 1)) %>%
  dplyr::mutate(regionx = factor(regionx, levels=rev(region_facs), ordered=TRUE)) %>%
  ggplot() + 
  geom_line(aes(date, spei6_0m, group=areaid), col="skyblue4", alpha=0.15, size=0.1) + 
  geom_line(data=region_clim, aes(date, spei6), col="grey10", size=0.4) + 
  geom_hline(yintercept=0, lty=2) + 
  theme_classic() + 
  plot_theme + facet_wrap(~regionx, ncol=1) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=13.5)) + 
  ylab("SPEI-6") + xlab("Month") +
  ggtitle("SPEI-6")


# SPEI1 by region
p3 = dd %>%
  left_join(regionx) %>%
  # dplyr::group_by(province, date) %>%
  # dplyr::summarise(spei6 = mean(spei6_0m),
  #                  spei1 = mean(spei1_0m),
  #                  region = head(region, 1)) %>%
  dplyr::mutate(regionx = factor(regionx, levels=rev(region_facs), ordered=TRUE)) %>%
  ggplot() + 
  geom_line(aes(date, spei1_0m, group=areaid), col="cyan4", alpha=0.15, size=0.1) + 
  geom_line(data=region_clim, aes(date, spei1), col="grey10", size=0.4) + 
  geom_hline(yintercept=0, lty=2) + 
  theme_classic() + 
  plot_theme + facet_wrap(~regionx, ncol=1) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=13.5)) + 
  ylab("SPEI-1") + xlab("Month") +
  ggtitle("SPEI-1")

pc = gridExtra::grid.arrange(p1, p2, p3, p4, nrow=1)
ggsave(pc, file="./output/figures/SuppFigure_RegionalClimateRegimes.jpg", device="jpg", units="in", width=15, height=8, dpi=600, scale=0.85)







# ====================== Supp: Geography and correlation structure among socio-environmental factors ======================

cv = dd %>%
  dplyr::mutate(
    logpopdens = log(popdens_census)
  ) %>%
  dplyr::select(
    areaid, 
    year,
    logpopdens,
    urban_pw,
    flux_grav1,
    flux_rad,
    urbanexp_3yr_d,
    urbanexp_10yr_d,
    water_piped_year,
    flushtoilet_any_year, 
    traffic_thouskmperinhab,
    tmean_coolestmonth,
  ) %>%
  dplyr::mutate(
    gravityf_log = log(flux_grav1 + 1),
    radiationf_log = log(flux_rad + 1),
    traffic_kmperinhab_log = log(traffic_thouskmperinhab * 1000)
  ) %>% 
  dplyr::select(
    -flux_rad,
    -flux_grav1,
    -traffic_thouskmperinhab
  ) %>%
  dplyr::filter(!is.na(radiationf_log)) %>%
  distinct()



# ========== 1. Correlation plot across full time series ==============

# correlation plot
cx = as.data.frame(cor(cv %>% dplyr::select(-areaid, -year)))
cx$param = row.names(cx)

facs = unique(cx$param)
facs = data.frame(param=facs,
                     pname = c("Pop. density (log)",
                               "Built-up land",
                               "Urban exp. (3 yr)", 
                               "Urban exp. (10 yr)",
                               "Piped water",
                               "Hygienic toilet",
                               "Tmean coolest",
                               "Gravity flux (log)",
                               "Radiation flux (log)",
                               "Road travel (log)"))
fac_ord = c("Pop. density (log)",
            "Built-up land",
            "Urban exp. (3 yr)", 
            "Urban exp. (10 yr)",
            "Piped water",
            "Hygienic toilet",
            "Gravity flux (log)",
            "Radiation flux (log)",
            "Road travel (log)",
            "Tmean coolest")

cx = cx %>%
  dplyr::left_join(facs) %>%
  dplyr::select(-param) %>%
  reshape2::melt(id.vars = "pname") %>%
  dplyr::left_join(
    facs, by=c("variable"="param")
  ) %>%
  dplyr::rename("x"=pname.x, "y"=pname.y) %>%
  dplyr::mutate(x = factor(x, levels=fac_ord, ordered=TRUE),
                y = factor(y, levels=rev(fac_ord), ordered=TRUE))

cx$thresh = cx$value>0.7

# remove self-to-self
cx = cx[ -which(cx$x == cx$y), ]

p1 = cx %>%
  ggplot() + 
  geom_point(aes(x, y, fill=value, color=value, size=value), alpha=0.9, pch=21) +
  #geom_point(data = cx[ cx$thresh == TRUE, ], aes(x, y, size=value), fill=NA, color="black", pch=21, stroke=1) +
  theme_minimal() +
  scale_size(range=c(1,10), name="") +
  scale_x_discrete(position="top") +
  scale_fill_viridis_c(direction=-1, begin=0.05, name="Corr.") +
  scale_color_viridis_c(direction=-1, begin=0.05, name="Corr.") +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size=14, angle=90),
        axis.text.y = element_text(size=14, angle=0),
        legend.title=element_text(size=15),
        legend.text = element_text(size=12))




# ================ Maps in middle of time series ==================

comb = cv %>%
  dplyr::filter(year==2009) %>%
  dplyr::select(-year) %>%
  reshape2::melt(id.vars=1) %>%
  dplyr::left_join(facs, by=c("variable"="param"))

facs = facs[ c(1:6, 8:10, 7), ]

map_list = vector("list", length = 10)

for(i in 1:nrow(facs)){
  
  comb_i = comb[ comb$pname == facs$pname[i], ]
  
  p_i = shp %>%
    dplyr::left_join(comb_i) %>%
    ggplot() + 
    geom_sf(aes(fill=value), col=NA) + 
    geom_sf(data=shp_vt, fill=NA, col="grey20") +
    scale_fill_gradientn(colors=viridisLite::turbo(200)) +
    ggtitle(comb_i$pname[1]) +
    maptheme +
    theme(legend.title = element_blank(),
          legend.position = c(0.2, 0.45), 
          plot.title = element_text(size=15, hjust=0.5))
  
  map_list[[i]] = p_i
  
}

p2 = gridExtra::grid.arrange(
  grobs=map_list,
  nrow=2
)


# combine
pc = gridExtra::grid.arrange(
  grobs=list(p2, p1), ncol=2, nrow=3,
  layout_matrix = rbind(c(1, NA), c(1, 2), c(1, NA)),
  widths=c(0.825, 0.5), heights=c(0.15, 0.8, 0.15)
)
ggsave(pc, file="./output/figures/SuppFigure_SocioEnvironmentalMaps.jpg", device="jpg", units="in", width=20, height=9, dpi=600, scale=0.9)






