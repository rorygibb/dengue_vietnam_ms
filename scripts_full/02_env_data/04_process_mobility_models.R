


# ======================= Connectivity models using travel-time based distance estimates ==============================

# uses travel time estimates between districts estimated using MAP friction surfaces (see process_traveltime_friction.R)
# and annual population estimates as weighting from HRW (focal provinces) and WorldPop (rest of Vietnam)
# estimate population flux between districts using naive gravity and radiation models, and comparison with overall traffic data from Vietnam GSO

setwd("C:/Users/roryj/Documents/PhD/202007_LSHTM_Dengue/analysis/")
pacman::p_load("gdistance", "abind", "rje", "ggplot2", "malariaAtlas", "rgdal", "raster", "dplyr", "data.table", 
               "sf", "exactextractr", "doParallel")

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))



# =============================== data ====================================

# polygons for entire of Vietnam (same areaid matches)
# specify unique field identifier (called "areaid") for processing
shp = st_read("./data/shapefiles/d_moss/vietnam_districts_merged.shp")
shp = shp[ ! shp$areanameen %in% c("Truong Sa", "Hoang Sa"), ] # remove offshore districts
shp$areaid = shp$areaid

# travel time matrix 
load("./output/data_processed/connectivity/VietnamAll_MergedSHP2_PairwiseTravelTime_Centroid_MAP_offshoreadjusted.R")
tmatrix = tt_viet$travelmatrix

# annual population for all districts (GPW)
# ppd = read.csv("./output/data_processed/population/VietPopulation_District_MergedSHP_GPW4_interpolated.csv", stringsAsFactors = FALSE) %>%
#   dplyr::rename("population"=population_gpw) %>%
#   dplyr::filter(areaid %in% shp$areaid)

# annual population (census)
ppd = read.csv("./code/viet_dengue_districts/output/covariates/Population_Census_2009_2019.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(areaid %in% shp$areaid) %>%
  dplyr::select(-urban, -rural) %>%
  dplyr::rename("population"=total)

# create table of pairwise distances, populations per year (711 * 711 * 22 years), and remove NAs
ppd_i = ppd %>% rename("area_i"=areaid, "population_i"=population) %>% dplyr::select(area_i, year, population_i)
ppd_j = ppd %>% rename("area_j"=areaid, "population_j"=population) %>% dplyr::select(area_j, year, population_j)
dd = reshape2::melt(tmatrix) %>%
  rename("area_i"=Var1, "area_j"=Var2, "dist_ij"=value) %>%
  left_join(ppd_i) %>%
  left_join(ppd_j) %>%
  filter(!is.na(dist_ij))

# create data.table 
dd = data.table(dd)




# ===================== connectivity estimates: naive gravity and radiation model ====================

# naive gravity = unparameterised gravity/diffusion model without weights on source/sink and distance decay https://www.pnas.org/content/112/38/11887
# T_ij = (pop_i^alpha * pop_j^beta) / dist_ij^lambda
# where alpha, beta and lambda are all == 1

# radiation = parameter-free general model for commuting accounting for attractiveness of other population centres in same distance radius
# doi.org/10.1038/nature10856
# T_ij =  T_i * ( (pop_i * pop_j) / ((pop_i + s_ij) * (pop_i + pop_j + s_ij)) )
# where T_i is the proportion of the population of district i that commute (can be fixed if unknown)
# and s_ij is the total population of the circle of radius dist_ij centred at district i (can be Euclidean distance or travel distance)

# n.b. Wesolowski 2015, PLOS Comp Biol for Sub Saharan Africa (i.e. low income setting)
# Radiation model does much better at predicting low volumes of travel, and gravity at high volumes of travel (e.g. urban-to-urban)
# "Scenarios to use a gravity model over a radiation model include: travel to and from a major population centre, over short distances, and when predicting large volumes of travel"
# "A radiation model should be used over a gravity model when describing travel between rural areas and low volumes of travel."
# "Caution should be taken using either model if the travel is between locations of intermediate rural population and over short distances"
# schematic suggests: 
# gravity = localised urban spreading; urban commuting
# radation = migration into cities; national spread of highly infectious disease; family visits
# in general areas with higher gravity factors (i.e. naive gravity) are better suited to gravity model

# n.b. Rabaa et al 2010; 2013: evidence that DENV diffusion within Vietnam follows gravity-like pathways around urban centres (incl. within HCMC)
# that HCMC acts as source population for majority of DENV diversity in country, including to SE (Dong Nai), north (incl. Ha Noi) and central
# evidence for year round persistence in south, central highlands and south central coast, and annual long-range introductions from south to north
# little evidence for persistence in north (but not specific to Hanoi and some evidence that some years show interannual persistence in more urbanised areas)
# Economic migration to HCMC and SouthEast (inc. Dong Nai - increasingly industrail) is especially common among young adults from Mekong Delta, Red River Delta, Central Coast
# So these links could act as diffusion pathways of infection across country
# Transmission pathways across regions (and more locally) explained by relative endemicity difference between donor and recipient area; and gravity-like connections


# -------------- gravity model --------------

# 1. naive gravity model (scale pops by 1000 just to reduce size)
k = 1
g_alpha = 1
g_beta = 1
g_lambda = 1

# gravity model 
gravityMod = function(k, pop_i, pop_j, dist_ij, alpha, beta, lambda){
  t_ij = k * ((pop_i^alpha * pop_j^beta) / (dist_ij^lambda))
  return(t_ij)
}

dd$Tij_grav1 = gravityMod(k = 1, pop_i=dd$population_i, pop_j=dd$population_j, dist_ij=dd$dist_ij, alpha=g_alpha, beta=g_beta, lambda=g_lambda)
dd$Tij_grav1 = dd$Tij_grav1 / 10^6

# 2. parameterised using estimates from mobile phone trip data and travel time estimates in Thailand 
# https://www.nature.com/articles/s41598-020-79438-0#Sec191
# provided in SI of article
k = 1
g_alpha = 0.751
g_beta = 0.683
g_lambda = 1.947
dd$Tij_grav2 = gravityMod(k = k, pop_i=dd$population_i, pop_j=dd$population_j, dist_ij=dd$dist_ij, alpha=g_alpha, beta=g_beta, lambda=g_lambda)
dd$Tij_grav2 = dd$Tij_grav2 / 1000






# ------------- radiation ---------------

# workflow:
# calculate S_ij (total population within travel radius between 2 focal districts, excluding source and destination)
# combine into df and calculate radiation estimate

# set lookup keys in data.table for fast lookup
# distances between all pairs of districts
distances = dd[ dd$year == 1998, c("area_i", "area_j", "dist_ij")]
distances$lookup = distances$area_i
distances = data.table::setkey(distances, lookup)

# population per year, per district
dpop = as.data.table(ppd)
dpop$lookup = dpop$areaid
dpop = data.table::setkey(dpop, lookup)

# for each pair of districts, get other districts within travel time radius, excluding source and destination
# applied to the xth row of distances df
getDistrictsInRadius = function(x){

  focal = distances[ x ,]
  radius = distances[ .(focal$lookup) ] %>%
    dplyr::filter(dist_ij <= focal$dist_ij) %>%
    dplyr::filter(!area_j %in% c(focal$area_i, focal$area_j))

  return(list(
    area_i = focal$area_i,
    area_j = focal$area_j,
    radius_districts = radius$area_j
  ))
}

# calculate total population within radius districts per year
# for each xth element of list produced by getDistrictsInRadius function
calculate_Sij = function(x){

  if(length(x$radius_districts) == 0){

    return(data.frame(area_i = x$area_i,
                      area_j = x$area_j,
                      year = 1998:2020,
                      S_ij = 0))
  } else{

    px = dpop[ .(x$radius_districts) ] %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(S_ij = sum(population)) %>%
      dplyr::mutate(area_i = x$area_i,
                    area_j = x$area_j)
    return(px)
  }
}

# run S_ij calculation then save
radius = lapply(1:nrow(distances), getDistrictsInRadius)
sij = lapply(radius, calculate_Sij)
sij_df = do.call(rbind.data.frame, sij)
write.csv(sij_df, "./output/data_processed/connectivity/VietnamAll_Radiation_Sij_Annual_MergedSHP2.csv", row.names=FALSE)

# combine with dd
#sij_df = read.csv("./output/data_processed/connectivity/VietnamAll_Radiation_Sij_Annual_MergedSHP2.csv", stringsAsFactors = FALSE)
dd = left_join(dd, sij_df)

# calculate radiation
# r_Ti = proportion persons commuting, assume all (constant)
radiationMod = function(r_Ti, pop_i, pop_j, s_ij){
  T_i = r_Ti * pop_i
  T_ij = T_i * ( (pop_i * pop_j) / ((pop_i + s_ij) * (pop_i + pop_j + s_ij)) )
  return(T_ij)
}
dd$Tij_rad = radiationMod(r_Ti = 1, pop_i=dd$population_i, pop_j = dd$population_j, s_ij = dd$S_ij)
save(dd, file="./output/data_processed/connectivity/Vietnam_MergedSHP2_ConnectivityModelsMatrix_TravelTime.R")





# ================= calculate summary mobility network metrics for model covariates ===============

# connectedness matrix
load("./output/data_processed/connectivity/Vietnam_MergedSHP2_ConnectivityModelsMatrix_TravelTime.R")

# remove same-to-same
dd = dd[ -which(dd$area_i == dd$area_j), ]

# # set same-to-same as NA
# dd$Tij_grav[ dd$area_i == dd$area_j ] <- NA
# dd$Tij_rad[ dd$area_i == dd$area_j ] <- NA


# ------------- 1. Mean flux per district per year (equivalent to sum of weights on all network edges, i.e. a "weighted degree" metric) -----------

# gravity model
flux_g = dd %>%
  #dplyr::filter(year != 2020) %>%
  group_by(year, area_i) %>%
  dplyr::summarise(flux_grav1 = mean(Tij_grav1, na.rm=TRUE),
                   flux_grav2 = mean(Tij_grav2, na.rm=TRUE),
                   population = head(population_i, 1)) %>%
  dplyr::rename("areaid" = area_i)

# radiation model: n.b. with area j as directed, so measuring influx
flux_r = dd %>%
  #dplyr::filter(year != 2020) %>%
  group_by(year, area_j) %>% 
  dplyr::summarise(flux_rad = mean(Tij_rad, na.rm=TRUE)) %>%
  dplyr::rename("areaid" = area_j)

flux = left_join(flux_g, flux_r) %>%
  left_join(st_drop_geometry(shp[ , c("areaid", "areanameen", "areaprovin") ]))




# --------------- 2. Degree per district per year -----------------------

# calculate modified degree (see Dhaka paper that takes this approach)
# traditional degree (number of edges connected to a node) is impossible for a fully-connected network
# instead define an edge-weight threshold beyond which that edge is considered to be "connected"
# then calculate degree based on those edges that satisfy that criterion
# in Dhaka study this threshold was defined as the mean weight across all edges (can also sensitivity analyse to median)

# run calculation
degree = dd %>%
  #dplyr::filter(year != 2020) %>%
  dplyr::mutate(grav1_exceeds_mean = Tij_grav1 >= mean(Tij_grav1, na.rm=TRUE),
                grav2_exceeds_mean = Tij_grav2 >= mean(Tij_grav2, na.rm=TRUE),
                rad_exceeds_mean = Tij_rad >= mean(Tij_rad, na.rm=TRUE)) %>%
  dplyr::group_by(area_i, year) %>%
  dplyr::summarise(population_i = head(population_i, 1),
                   degree_grav1 = sum(grav1_exceeds_mean),
                   degree_grav2 = sum(grav2_exceeds_mean),
                   degree_rad = sum(rad_exceeds_mean)) %>%
  dplyr::rename("areaid" = area_i)

# compare to setting threshold as median
# relationship is broadly logarithmic; relaxing threshold means that many more districts become much more connected
# but with some variability; may choose to explore this
# degreem = dd %>%
#   dplyr::filter(year != 2020) %>%
#   dplyr::mutate(grav1_exceeds_mean = Tij_grav1 >= median(Tij_grav1, na.rm=TRUE),
#                 grav2_exceeds_mean = Tij_grav2 >= median(Tij_grav2, na.rm=TRUE),
#                 rad_exceeds_mean = Tij_rad >= median(Tij_rad, na.rm=TRUE)) %>%
#   dplyr::group_by(area_i, year) %>%
#   dplyr::summarise(population_i = head(population_i, 1),
#                    degree_grav1 = sum(grav1_exceeds_mean),
#                    degree_grav2 = sum(grav2_exceeds_mean),
#                    degree_rad = sum(rad_exceeds_mean))
# plot(degree$degree_grav1, degreem$degree_grav1)
# #plot(degree$degree_grav2, degreem$degree_grav2)
# plot(degree$degree_rad, degreem$degree_rad)


# ------------ Combine metrics and save ------------

# combine and save
mobility = left_join(degree, flux_g)
mobility = left_join(mobility, flux_r)
mobility = mobility %>% dplyr::select(-population)
write.csv(mobility, "./code/viet_dengue_districts/output/covariates/MobilityMetrics_MergedSHP2_Feb2021.csv", row.names=FALSE)

# shp %>%
#   dplyr::full_join(
#     mobility %>% dplyr::filter(year %in% c(1999, 2009, 2019))
#     ) %>%
#   sf::st_crop(extent(shp %>% filter(areaprovin == "Ha Noi")) + 1) %>%
#   ggplot() +
#   geom_sf(aes(fill=flux_rad), col=NA) +
#   facet_wrap(~year) +
#   scale_fill_gradientn(colors=viridisLite::turbo(200))

