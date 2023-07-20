
# ==========================================================================================#
#                                                                                           #
#         EXAMPLE PIPELINE 2: COMBINES DENGUE, SOCIO, CLIMATE, MOBILITY COVARIATES          #
#                         OUTPUTS SAVED TO "OUTPUT/MODEL_DATA"                              #
#                                                                                           #
# ==========================================================================================#



# ================ Combine dengue data and covariates from socioeconomic / earth observation into dataset for modelling ===============

# script combines all covariates produced using "process" scripts with dengue incidence data
# produces final dataset for use in statistical models
# scripts that produced each constituent dataset/covariate are named below

# project root and dependencies
PATH = dirname(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(PATH)
pacman::p_load("dplyr", "raster", "rgdal", "sf", "ecmwfr", "stringr", "ggplot2", "lubridate", "magrittr", "vroom")

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))



# ================= key objects =================

# districts to be excluded (offshore)
offshore_areas = c(70154, 70339, 70273, 70355, 70698)

# districts shapefile for Vietnam
shp = st_read("./data/shapefiles/dengue_districts_shapefile.shp") %>%
  dplyr::filter(!areaid %in% offshore_areas)
shp$area = st_area(shp)/10^6




# ================ dengue data and derive metrics ===================

# dengue and set year/month to dengue year/month
dd = read.csv("./data/dengue/dengue_districts_19982020.csv") %>%
  dplyr::filter(!areaid %in% offshore_areas) %>%
  dplyr::mutate(year = yeardengue, month = monthdengue) %>%
  dplyr::select(-yeardengue, -monthdengue) %>%
  dplyr::mutate(cases = as.numeric(cases),
                date = as.Date(date),
                year_month = paste(year, month, sep="_")) %>%
  dplyr::arrange(province, areaid, date) %>%
  dplyr::select(-population, -population_urban, -incidence)

# test
if(all(dd$areaid %in% shp$areaid)){ print("All districts matched in shapefile = TRUE") }

# human population
pop = read.csv("./output/covariates/Population_Census_2009_2019.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(!areaid %in% offshore_areas) %>%
  dplyr::rename("population_census"=4, "urbanpop_census"=5, "ruralpop_census"=6) %>%
  dplyr::left_join(shp %>% dplyr::select(areaid, area) %>% st_drop_geometry()) %>%
  dplyr::mutate(popdens_census = as.vector(population_census / area),
                pop_propurban_census = urbanpop_census / population_census) %>%
  dplyr::select(-area, -areanameen, -province)

# test and combine
if(all(dd$areaid %in% pop$areaid)){ print("Pop data for all districts = TRUE") }
dd = dplyr::left_join(dd, pop)

# calculate crude incidence (cases/100,000)
dd$incidence = (dd$cases / dd$population_census) * 100000

# fix for name
dd$district[ dd$areaid == 70266 ] = "Ninh Binh"



# =================== socio-environmental metrics =====================

# land change dynamics
fx1 = read.csv("./output/covariates/LandUseDynamics_MergedSHP2_Feb2021.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(!areaid %in% offshore_areas) %>%
  dplyr::select(-areanameen, -areaprovin, -urbanexp3_dev_d, -urbanexp10_dev_d) 

# agriculture and other land cover types
fx2 = read.csv("./output/covariates/LandUseCover_ESA_MergedSHP2.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(!areaid %in% offshore_areas) %>%
  dplyr::select(-id, -areanameen)

# socioeconomic variables 
fx3 = read.csv("./output/covariates/Socioeconomic_PCAs_TimeSeriesProj.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(!areaid %in% offshore_areas) %>%
  dplyr::select(areaid, year, water_piped_well, sanitation_flushtoilet_indoor, sanitation_flushtoilet_any, sanitation_flushtoilet_outdoor) %>%
  dplyr::rename("water_piped_year"=3, "flushtoilet_indoor_year"=4, "flushtoilet_any_year"=5, "flushtoilet_outdoor_year"=6)

# human mobility/connectivity estimates using mobility models and province-level traffic estimates
fx4 = read.csv("./output/covariates/MobilityMetrics_MergedSHP2_Feb2021.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(!areaid %in% offshore_areas) %>%
  dplyr::select(-population_i, -degree_grav2, -flux_grav2) %>%
  dplyr::left_join(shp %>% st_drop_geometry() %>% dplyr::select(areaid, areaprovin)) %>% 
  dplyr::rename("province"=areaprovin) %>%
  dplyr::left_join(read.csv("./output/covariates/Mobility_TrafficMetrics_VietGSO.csv")) %>%
  dplyr::select(-province, -provinceid)

# travel time to city/airport
fx5 = read.csv("./output/covariates/VietnamAll_MergedSHP2_CityAirportTravelTimes_MAP.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(!areaid %in% offshore_areas) %>%
  dplyr::select(-areanameen, -areaprovin)
  
# combine and include only dengue years
fx = left_join(fx1, left_join(fx2, left_join(fx3, left_join(fx4, fx5))))
fx = fx[ fx$year %in% dd$year, ]

# test
if(all(dd$areaid %in% fx$areaid)){ print("Annual covar data for all districts = TRUE") }



# ===================== climate data (ERA5-LAND temperature; WFDE5/ERA5 precip): calculate lags ========================

# other variables from ERA5-LAND (don't use ERA5 precip for full time series; only for final years)
# combine with precip from "prec" above
# data produced using process_clim_era5land.R
clim = read.csv("./output/covariates/VietDistricts_ERA5Land_monthly_MergedSHP2.csv", stringsAsFactors = FALSE) %>%
  dplyr::select(areaid, Date, variable, value) %>%
  dplyr::filter(!areaid %in% offshore_areas) %>%
  reshape2::dcast(areaid + Date ~ variable) %>%
  dplyr::rename_all(tolower) %>%
  mutate(date = as.Date(date)) %>%
  dplyr::select(-precipmean) %>%
  distinct()
 
# precipitation indicators
precip = read.csv("./output/covariates/PrecipitationIndicators_VietAll_SPEI_MergedSHP2.csv", stringsAsFactors  = FALSE) %>%
  dplyr::select(-areanameen, -areaprovin) %>%
  dplyr::mutate(date = as.Date(date)) %>%
  dplyr::filter(lubridate::year(date) %in% lubridate::year(clim$date)) %>%
  dplyr::filter(!areaid %in% offshore_areas)

# fix single issue value in spei
precip$spei_1[ precip$areaid == "70495" & precip$date == as.Date("2004-10-01") ] = 
  mean(
    precip$spei_1[ precip$areaid == "70495" & precip$date == as.Date("2004-09-01") ],
    precip$spei_1[ precip$areaid == "70495" & precip$date == as.Date("2004-11-01") ]
  )

# combine
clim = clim %>% dplyr::left_join(precip)


# ------------------ create annual mean variables -------------------

# annual means per dengue year
amc = clim %>%
  dplyr::left_join(dd[ , c("areaid", "date", "year")] %>% distinct()) %>%
  dplyr::mutate(precip_tot = precip * lubridate::days_in_month(date)) %>%
  dplyr::filter(!is.na(year)) %>%
  dplyr::filter(!is.na(tmean)) %>%
  dplyr::group_by(areaid, year) %>%
  dplyr::summarise(tmean_annualmean = mean(tmean, na.rm=TRUE),
                   tmin_annualmean = mean(tmin, na.rm=TRUE),
                   tmean_coolestmonth = min(tmean, na.rm=TRUE),
                   tmin_coolestmonth = min(tmin, na.rm=TRUE),
                   tmax_warmestmonth = max(tmax, na.rm=TRUE),
                   precip_total = sum(precip_tot, na.rm=TRUE))


# ----------------- create lagged climate variables ----------------

# calcuate climate in lagged window of n months up to and including focal month (set to 3)
window_size = 2
clim = clim %>%  
  dplyr::arrange(areaid, date) %>%
  dplyr::group_by(areaid) %>%
  dplyr::mutate(tdrange_0m = data.table::frollmean(tdrange, window_size, align="right"),
                tmax_0m = data.table::frollmean(tmax, window_size, align="right"),
                tmean_0m = data.table::frollmean(tmean, window_size, align="right"),
                tmin_0m = data.table::frollmean(tmin, window_size, align="right"),
                windspeed_0m = data.table::frollmean(windspeed, window_size, align="right"),
                precip_0m = data.table::frollmean(precip, window_size, align="right"),
                spei1_0m = data.table::frollmean(spei_1, window_size, align="right"),
                spei3_0m = data.table::frollmean(spei_3, window_size, align="right"),
                spei6_0m = data.table::frollmean(spei_6, window_size, align="right"),
                spei12_0m = data.table::frollmean(spei_12, window_size, align="right")) %>%
  dplyr::select(areaid, date, tmean, tmin, tmax, tdrange,
                tmean_0m, tmin_0m, tmax_0m, tdrange_0m, precip_0m, spei1_0m, spei3_0m, spei6_0m, spei12_0m, windspeed_0m)
#if(all(as.vector(climl_3[ 3, 9:14] == apply(climl_3[ 1:3, 3:8], 2, mean)))){ print("Rollmean test passed")}

# create lags
# subset to only env data, rename, add X months to date (to bring up to lag time) and then left_join to clim
c2 = clim
for(lag in 1:6){
  lag1 = c2 %>% dplyr::select(-tmean, -tmin, -tmax, -tdrange)
  new_names = unlist(lapply(strsplit(names(lag1)[3:ncol(lag1)], "_"), "[", 1))
  names(lag1)[3:ncol(lag1)] = paste0( new_names, "_", lag, "m", sep="")
  lag1$date = ymd(lag1$date) %m+% months(lag)
  clim = left_join(clim, lag1)
}

# # create future temperature lags for causality testing (1-6 months in future)
# c2 = clim
# for(lag in c(1, 2, 3, 6)){
#   lag1 = c2 %>% dplyr::select(areaid, date, tmean_0m, tmin_0m, tmax_0m)
#   new_names = unlist(lapply(strsplit(names(lag1)[3:ncol(lag1)], "_"), "[", 1))
#   names(lag1)[3:ncol(lag1)] = paste0( new_names, "_future_", lag, "m", sep="")
#   lag1$date = ymd(lag1$date) %m-% months(lag)
#   clim = left_join(clim, lag1)
# }

# add annual mean estimates
clim = clim %>%
  dplyr::left_join(dd[ , c("areaid", "date", "year")] %>% distinct()) %>%
  left_join(amc)

# climate for offshore islands and coastal peninsulas without climate data: replace with nearest mainland location
# Son Tra (70361) replace with Hai Chau (70363)
# Ly Son (areaid = 70382) replace with Binh Son (very close geographically) (areaid 70384)
# Phu Quy (areaid = 70596) replae with Bac Binh district areaid = 70517
# Kien Hai (areaid = 70671) replace with An Minh areaid = 70680
# Con Dao (areaid = 70711) replace with Vinh Chau (70699)

clim = rbind(clim[ clim$areaid != 70361, ], clim[ clim$areaid == 70363, ] %>% dplyr::mutate(areaid = "70361", areanameen = "Son Tra"))
clim = rbind(clim[ clim$areaid != 70382, ], clim[ clim$areaid == 70384, ] %>% dplyr::mutate(areaid = "70382", areanameen = "Ly Son"))
clim = rbind(clim[ clim$areaid != 70596, ], clim[ clim$areaid == 70517, ] %>% dplyr::mutate(areaid = "70596", areanameen = "Phu Quy"))
clim = rbind(clim[ clim$areaid != 70671, ], clim[ clim$areaid == 70680, ] %>% dplyr::mutate(areaid = "70671", areanameen = "Kien Hai"))
clim = rbind(clim[ clim$areaid != 70711, ], clim[ clim$areaid == 70699, ] %>% dplyr::mutate(areaid = "70711", areanameen = "Con Dao"))




# ===================== table of climatic subregions of Vietnam =========================

# three classification regimes:
# 1: Vietnam subregions classification as defined nationally, but with Southeast and Mekong Delta combined because of climatic similarity
# 2. "Climatic regimes": combines several clusters in north and north central together to more closely reflect temperature regimes 
# 3. North-South-Central climatic regimes

# table of regions per province
rr = read.csv("./data/shapefiles/regions_lookup.csv") %>%
  dplyr::rename("province"=provincena)

# Region 1: combine southeast and Mekong delta (n=7 zones)
rr$region1 = rr$region
rr$region1[ rr$region1 %in% c("Southeast", "Mekong River Delta")] = "South"

# Region 2: "Climatic regimes" 
rr$region2 = rr$region1
rr$region2[ rr$region2 %in% c("Southeast", "Mekong River Delta")] = "South"
rr$region2[ rr$province == "Binh Thuan" ] = "South"
rr$region2[ rr$region2 %in% c("Northeast", "Northwest", "Red River Delta") ] = "Red River Delta"
rr$region2[ rr$province %in% c("Lang Son", "Yen Bai", "Bac Kan", "Cao Bang", "Lao Cai", "Ha Giang") ] = "North/Northeast"
rr$region2[ rr$province %in% c("Dien Bien", "Son La", "Lai Chau") ] = "Northwest"
rr$region2[ rr$province %in% c("Thanh Hoa") ] = "Red River Delta"

# Region 3: Broad "transmission settings": north, south, central
rr$region3 = rr$region2
rr$region3[ rr$region3 %in% c("North/Northeast", "Northwest", "Red River Delta") ] = "North"
rr$region3[ rr$region3 %in% c("North Central", "Central Highlands", "South Central Coast") ] = "Central"



# ===================== save everything in a modelling materials folder ======================

# combine these in modelling script (reduces space issues for git)
# set date to the middle of the month for dd, clim, anom and save

dd$date = dd$date + 15
write.csv(dd, "./output/model_data/ModelData_Dengue_VietAll.csv", row.names=FALSE)

write.csv(fx, "./output/model_data/ModelData_SocioEcologicalCovar_VietAll.csv", row.names=FALSE)

clim$date = clim$date + 15
vroom::vroom_write(clim, "./output/model_data/ModelData_ClimateLags_VietAll.csv.gz", delim=",")
write.csv(clim, "./output/model_data/ModelData_ClimateLags_VietAll.csv", row.names=FALSE)

write.csv(rr, "./output/model_data/ModelData_ClimaticRegions.csv")


