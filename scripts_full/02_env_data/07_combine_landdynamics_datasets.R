

# ==================== Combine and compare land cover change datasets ===================


# project root and dependencies
setwd("C:/Users/roryj/Documents/PhD/202007_LSHTM_Dengue/analysis/")
pacman::p_load("ggplot2", "rgdal", "raster", "dplyr", "sf", "exactextractr", "ggspatial")

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))

# viet districts shapefile
shp = st_read("./data/shapefiles/d_moss/vietnam_districts_final.shp") %>%
  filter(!areanameen %in% c("Truong Sa", "Hoang Sa")) 



# =============== Merged districts wtih lookup table to match ====================

# merged districts in shapefile
shpm = sf::st_read("./data/shapefiles/d_moss/vietnam_districts_merged.shp") %>%
  dplyr::filter(!areanameen %in% c("Truong Sa", "Hoang Sa"))

# lookup table cross-referencing areaid/areanames to updated ones
cross_ref = shpm %>%
  st_drop_geometry() %>%
  dplyr::select(areaid, areanameen) %>%
  dplyr::filter(!areaid %in% shp$areaid) 
cr2 = cross_ref %>%
  tidyr::separate_rows(areaid, sep="_") %>%
  dplyr::rename("areaid_orig" = areaid)
cross_ref = left_join(cr2, cross_ref)



# # ================ Commune-level data aggregated to district-level ====================
# 
# # population 
# pop = read.csv("./output/data_processed/population/VietnamAll_Population_Commune_WorldPop_Feb2021.csv", stringsAsFactors = FALSE) %>%
#   dplyr::select(year, GID_3, population)
# #table(pop$year[ is.na(pop$population) ])
# 
# # 1. forest loss: calculate cumulative losses in 3, 5, 10 years windows (short, med, long-ish time horizon, helps deal with uncertainty in year of conversion)
# # issue: Hansen data only present from 2000 so early years of time series require interpolation
# forest = read.csv("./output/data_processed/landcover/hansen/VietnamAll_ForestChange_Hansen_Communes_Feb2021.csv", stringsAsFactors = FALSE) %>%
#   dplyr::filter(!is.na(GID_3)) %>%
#   dplyr::arrange(GID_3, year) %>%
#   group_by(GID_3) %>%
#   dplyr::mutate(forestloss_3yr = data.table::frollsum(forestloss_area, 3, align="right"),
#                 forestloss_5yr = data.table::frollsum(forestloss_area, 5, align="right"),
#                 forestloss_10yr = data.table::frollsum(forestloss_area, 10, align="right"),
#                 forestloss_cumulative = cumsum(forestloss_area)) %>%
#   left_join(pop)
# 
# # replace areanames/areaids in forest with merged shp
# for(i in 1:nrow(cross_ref)){
#   a_i = cross_ref[ i, ]
#   forest$areanameen[ forest$areaid == a_i$areaid_orig ] = a_i$areanameen
#   forest$areaid[ forest$areaid == a_i$areaid_orig ] = a_i$areaid
# }
# 
# # aggregate to the district level (population-weighted mean of commune-level change)
# for_pop = forest %>%
#   dplyr::group_by(areaid, year) %>%
#   dplyr::mutate(pop_weight = population / sum(population, na.rm=TRUE)) %>%
#   dplyr::summarise(areanameen = head(areanameen, 1),
#                    areaprovin = head(areaprovin, 1),
#                    forestloss_3yr_cw = sum(forestloss_3yr * pop_weight, na.rm=TRUE),
#                    forestloss_5yr_cw = sum(forestloss_5yr * pop_weight, na.rm=TRUE),
#                    forestloss_10yr_cw = sum(forestloss_10yr * pop_weight, na.rm=TRUE),
#                    forestloss_cumulative_cw = sum(forestloss_cumulative * pop_weight, na.rm=TRUE))
# 
# # set years outside window at start of time series to NA
# for_pop$forestloss_10yr_cw[ for_pop$year < 2004 ] = NA
# for_pop$forestloss_5yr_cw[ for_pop$year < 1999 ] = NA
# for_pop$forestloss_3yr_cw[ for_pop$year < 1997 ] = NA
# 
# 
# 
# # 2. urban expansion rates 
# urban = read.csv("./output/data_processed/landcover/li_urbandynamics/VietnamAll_Urbanisation_Communes_Feb2021.csv", stringsAsFactors = FALSE) %>%
#   dplyr::filter(!is.na(GID_3)) %>%
#   dplyr::arrange(GID_3, year) %>%
#   group_by(GID_3) %>%
#   dplyr::mutate(u2 = replace(urbanexp_km2, is.na(urbanexp_km2), 0)) %>%
#   dplyr::rename("urbanexp" = urbanexp_km2) %>%
#   dplyr::mutate(urbanexp_3yr = data.table::frollsum(urbanexp, 3, align="right"),
#                 urbanexp_5yr = data.table::frollsum(urbanexp, 5, align="right"),
#                 urbanexp_10yr = data.table::frollsum(urbanexp, 10, align="right"),
#                 urbanexp_15yr = data.table::frollsum(urbanexp, 15, align="right"),
#                 urbanexp_cumulative = cumsum(u2)) %>%
#   left_join(pop) %>%
#   dplyr::select( -urbanexp_prop, -u2 )
# 
# # replace areanames/areaids with merged shp
# for(i in 1:nrow(cross_ref)){
#   a_i = cross_ref[ i, ]
#   urban$areanameen[ urban$areaid == a_i$areaid_orig ] = a_i$areanameen
#   urban$areaid[ urban$areaid == a_i$areaid_orig ] = a_i$areaid
# }
# 
# # calculate mean at the district level (population-weighted)
# urb_pop = urban %>%
#   dplyr::group_by(areaid, year) %>%
#   dplyr::mutate(pop_weight = population / sum(population, na.rm=TRUE)) %>%
#   dplyr::summarise(areanameen = head(areanameen, 1),
#                    urbanexp_3yr_cw = sum(urbanexp_3yr * pop_weight, na.rm=TRUE),
#                    urbanexp_5yr_cw = sum(urbanexp_5yr * pop_weight, na.rm=TRUE),
#                    urbanexp_10yr_cw = sum(urbanexp_10yr * pop_weight, na.rm=TRUE),
#                    urbanexp_15yr_cw = sum(urbanexp_15yr * pop_weight, na.rm=TRUE),
#                    urbanexp_cumulative_cw = sum(urbanexp_cumulative * pop_weight, na.rm=TRUE),
#                    areaprovin = head(areaprovin, 1))
# 
# # combine
# luc = left_join(urb_pop, for_pop) %>%
#   dplyr::select(-urbanexp_cumulative_cw, -forestloss_cumulative_cw)



# ====================== District-level data ========================

# forest change
forest = read.csv("./output/data_processed/landcover/hansen/VietnamAll_MergedSHP2_ForestChangeHansen_Districts_Feb2021.csv", stringsAsFactors = FALSE) %>%
  dplyr:: arrange(areaid, year) %>%
  group_by(areaid) %>%
  dplyr::mutate(forestloss_3yr_d = data.table::frollsum(forestloss_area, 3, align="right"),
                forestloss_5yr_d = data.table::frollsum(forestloss_area, 5, align="right"),
                forestloss_10yr_d = data.table::frollsum(forestloss_area, 10, align="right"),
                forestloss_cumulative_d = cumsum(forestloss_area))

# mean yearly change
forest2 = read.csv("./output/data_processed/landcover/hansen/VietnamAll_MergedSHP2_ForestChangeHansen_Districts_Feb2021.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(year >= 1998) %>%
  dplyr::arrange(areaid, year) %>%
  group_by(areaid) %>%
  dplyr::summarise(forestloss_meanannual = mean(forestloss_area, na.rm=TRUE))
forest = left_join(forest, forest2)

# urban change
urban = read.csv("./output/data_processed/landcover/li_urbandynamics/VietnamAll_MergedSHP2_UrbanisationRatesLiu_Districts_Feb2021.csv", stringsAsFactors = FALSE) %>%
  dplyr:: arrange(areaid, year) %>%
  group_by(areaid) %>%
  dplyr::mutate(urbanexp_3yr_d = data.table::frollsum(urbanexp_km2, 3, align="right"),
                urbanexp_5yr_d = data.table::frollsum(urbanexp_km2, 5, align="right"),
                urbanexp_10yr_d = data.table::frollsum(urbanexp_km2, 10, align="right"),
                urbanexp_15yr_d = data.table::frollsum(urbanexp_km2, 15, align="right"),
                urbanexp_cumulative_d = cumsum(urbanexp_km2)) %>%
  dplyr::filter(year >= 1990) %>%
  dplyr::mutate(urbanexp3_dev_d = as.vector(scale(urbanexp_3yr_d)),
                urbanexp10_dev_d = as.vector(scale(urbanexp_10yr_d)), 
                urbanexp3_dev_d = replace(urbanexp3_dev_d, is.nan(urbanexp3_dev_d), 0),
                urbanexp10_dev_d = replace(urbanexp10_dev_d, is.nan(urbanexp10_dev_d), 0))

# urban change mean annual
urban2 = read.csv("./output/data_processed/landcover/li_urbandynamics/VietnamAll_MergedSHP2_UrbanisationRatesLiu_Districts_Feb2021.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(year >= 1998) %>%
  dplyr:: arrange(areaid, year) %>%
  group_by(areaid) %>%
  dplyr::summarise(urbanexp_meanannual = mean(urbanexp_km2, na.rm=TRUE))
urban = left_join(urban, urban2)

# combine and save
lud = left_join(urban, forest) %>%
  dplyr::select(-liu_type, -urban_prop, -urbanexp_km2, -urbanexp_prop, -forestloss_area, -data_type, -urbanexp_cumulative_d, -forestloss_cumulative_d)
#comb = left_join(lud, luc)
comb = lud
write.csv(comb, "./code/viet_dengue_districts/output/covariates/LandUseDynamics_MergedSHP2_Feb2021.csv", row.names=FALSE)



# # ==================== Comparison of district and commune level data =======================
# 
# # comparisons 
# comb = left_join(lud, luc)
# #n_distinct(comb$areaid)
# 
# # generally correlated and closer in urban; 
# # in forest, notably population weighting substantially reduces the "exposure" to forest loss (as this tends to happen in more remote areas)
# ggplot(comb) + geom_point(aes(urbanexp_3yr_d, urbanexp_3yr_cw))
# ggplot(comb) + geom_point(aes(urbanexp_10yr_d, urbanexp_10yr_cw))
# ggplot(comb) + geom_point(aes(urbanexp_10yr_d, urbanexp_10yr_cw))
# ggplot(comb) + geom_point(aes(forestloss_3yr_d, forestloss_3yr_cw))
# ggplot(comb) + geom_point(aes(forestloss_10yr_d, forestloss_10yr_cw))
# 
# #
# xy = cbind(shpm[ , c("areaid", "areanameen")], st_coordinates(st_centroid(shpm))) %>%
#   st_drop_geometry() %>%
#   dplyr::arrange(desc(Y))
# comb$areaid = factor(comb$areaid, levels=xy$areaid, ordered=TRUE)
# 
# # plot urb
# colRamp = colorRampPalette(RColorBrewer::brewer.pal("YlGnBu", n=9))(400)
# comb$date = as.Date(paste(comb$year, "-01-01", sep=""))
# ggplot(comb[ comb$year >= 1998, ]) + 
#   geom_tile(aes(x=date, y=areaid, fill=urbanexp_10yr_d), width=366) +
#   #scale_fill_viridis_c(option="viridis", direction=-1, na.value = "grey60") +
#   scale_fill_gradientn(colors = rev(colRamp), na.value="grey90") +
#   theme_classic() +
#   #scale_x_date(date_breaks = "2 years", date_labels="%Y", limits=as.Date(c("1998-01-01", "2019-12-31")), position="bottom") +
#   scale_y_discrete(position="right") +
#   #scale_x_breaks=as.Date(paste(seq(1998, 2020, by=2), "-01-01", sep="")), labels=as.Date(paste(seq(1998, 2020, by=2), "-01-01", sep=""))) +
#   theme(legend.position="bottom",
#         legend.title=element_blank(),
#         axis.text.x = element_text(size=14.5),
#         axis.text.y = element_text(size=8.5),
#         axis.title.y = element_blank(),
#         axis.title.x = element_text(color="white", size=14),
#         axis.line.y = element_blank(),
#         axis.line.x = element_blank(),
#         axis.ticks.y = element_blank())



# ================== Viz =========================

# areas of districts
shpm$area = as.vector(st_area(shpm)/10^6)

# province level population density
pop = read.csv("./output/data_processed/population/VietPopulation_District_MergedSHP_GPW4_interpolated.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(year == 2019) %>%
  dplyr::left_join(shpm[ , c("areaid", "area") ] %>% st_drop_geometry()) %>%
  dplyr::mutate(popdens = population_gpw / area) %>%
  dplyr::group_by(areaprovin) %>%
  dplyr::summarise(popdens = mean(popdens, na.rm=TRUE)) %>%
  dplyr::arrange(desc(popdens))

# 
library(ggridges)
comb2 = comb[ , c("year", "areaid", "areaprovin", "urbanexp_10yr_d", "urbanexp_3yr_d", "urbanexp_10yr_cw", "forestloss_3yr_d")] %>%
  distinct() %>%
  dplyr::filter(year > 1997) 
comb2$areaprovin = factor(comb2$areaprovin, levels=rev(pop$areaprovin), ordered=TRUE)

# within 99th percentile of val
ggplot() + 
  geom_density_ridges(data = comb2[ comb2$urbanexp_3yr_d <= quantile(comb2$urbanexp_3yr_d, 0.99), ], scale=2,
                      aes(x=urbanexp_10yr_d, y=areaprovin), height=1, fill="skyblue3", col="skyblue4", size=0.7, alpha=0.6) +
  theme_classic()

# 
ggplot() + 
  geom_density_ridges(data = comb2[ comb2$forestloss_3yr_d <= quantile(comb2$forestloss_3yr_d, 0.95) & comb2$forestloss_3yr_d != 0, ], aes(x=forestloss_3yr_d, y=areaprovin), height=1, fill="skyblue4", alpha=0.6) +
  theme_classic()

# ggplot(comb[ comb$year >= 1998, ]) + 
#   geom_tile(aes(x=date, y=areaid, fill=urbanexp_10yr_cw), width=366) +
#   #scale_fill_viridis_c(option="viridis", direction=-1, na.value = "grey60") +
#   scale_fill_gradientn(colors = colRamp, na.value="grey50") +
#   theme_classic() +
#   #scale_x_date(date_breaks = "2 years", date_labels="%Y", limits=as.Date(c("1998-01-01", "2019-12-31")), position="bottom") +
#   scale_y_discrete(position="right") +
#   #scale_x_breaks=as.Date(paste(seq(1998, 2020, by=2), "-01-01", sep="")), labels=as.Date(paste(seq(1998, 2020, by=2), "-01-01", sep=""))) +
#   theme(legend.position="bottom",
#         legend.title=element_blank(),
#         axis.text.x = element_text(size=14.5),
#         axis.text.y = element_text(size=8.5),
#         axis.title.y = element_blank(),
#         axis.title.x = element_text(color="white", size=14),
#         axis.line.y = element_blank(),
#         axis.line.x = element_blank(),
#         axis.ticks.y = element_blank())



# =============== simple pop dens vs land use ================

# # 
# pop = read.csv("./output/data_processed/population/VietPopulation_District_MergedSHP_GPW4_interpolated.csv", stringsAsFactors = FALSE) 
# comb = left_join(comb, pop)
# comb = comb[ comb$year > 1997, ]
# 
# ggplot(comb) + geom_point(aes(population_gpw, urbanexp_10yr_d))
# ggplot(comb) + geom_point(aes(population_gpw, urbanexp_3yr_d))
# ggplot(comb) + geom_point(aes(population_gpw, urbanexp_3yr_cw))
# ggplot(comb) + geom_point(aes(population_gpw, urbanexp_10yr_cw))
# ggplot(comb) + geom_point(aes(population_gpw, forestloss_3yr_d))
# ggplot(comb) + geom_point(aes(population_gpw, forestloss_3yr_cw))
# 
# 
# cor(c2$population_gpw, c2$urbanexp_10yr_d)
