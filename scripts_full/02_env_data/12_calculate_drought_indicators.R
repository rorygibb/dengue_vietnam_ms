

# ================ Calculate precipitation indices (SPEI; mean monthly; number of dry days) ===============

# project dir
PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

# dependencies
library(raster)
library(rgdal)
library(sf)
library(ecmwfr)
library(stringr)
library(lubridate)
library(SPEI)
library(dplyr)

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))

# districts shapefile for Vietnam
shp = st_read("./data/shapefiles/dengue_districts_shapefile.shp")
shp = shp[ ! shp$areanameen %in% c("Truong Sa", "Hoang Sa"), ] # remove offshore districts




# ============== calculate SPEI using long-term time series ================


# ------------------ create dataset: precip, tmean, longitude, latitude, pev, since 1981 -------------------

# monthly precipitation, mean daily metres (WFDE5 to 2016, ERA5-LAND 2017-2010) 
prec1 = raster::brick("./data/climate/wfde5/PrecipMean_1979-01-01_2019-12-01_monthly_wfde5gpcc_v2.1_Vietnam.grd")
names(prec1) = paste("WFDE5GPCC", unlist(lapply(strsplit(names(prec1), "_"), "[", 2)), sep="_")

prec2a = terra::rast("./data/climate/era5-land/total_precipitation/meanprecip_month_viet_2020.nc")
prec2b = terra::rast("./data/climate/era5-land/total_precipitation/meanprecip_month_viet_2021.nc")
prec2 = c(prec2a, prec2b)
names(prec2) = paste("ERA5LAND", as.Date(terra::time(prec2)), sep="_")

# WFDE5 CRUGPCC: units are m/day so multiply by 1000 to daily mm
# convert date to 16th of the month (mid-month)
e1 = exactextractr::exact_extract(prec1, shp, fun="mean", progress=FALSE) %>%
  cbind(st_drop_geometry(shp[ , c("areaid", "areanameen", "areaprovin")])) %>%
  reshape2::melt(id.vars=c("areaid", "areanameen", "areaprovin")) %>%
  dplyr::mutate(date = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 2)),
                dataset = "WFDE5-CRU+GPCC",
                value = value * 1000) %>%
  dplyr::select(-variable) %>%
  dplyr::mutate(date = paste(substr(date, 1, 7), "01", sep="."),
                date = as.Date(date, "%Y.%m.%d"))

# ERA5-LAND: units are m/day so multiply by 1000 to daily mm
# convert date to 16th of the month (mid-month)
e2 = exactextractr::exact_extract(prec2, shp, fun="mean", progress=FALSE) %>%
  cbind(st_drop_geometry(shp[ , c("areaid", "areanameen", "areaprovin")])) %>%
  reshape2::melt(id.vars=c("areaid", "areanameen", "areaprovin")) %>%
  dplyr::mutate(date = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 2)),
                dataset = "ERA5-LAND",
                value = value * 1000) %>%
  dplyr::select(-variable) %>%
  dplyr::mutate(date = paste(substr(date, 1, 7), "01", sep="-"), 
                date = as.Date(date,  "%Y-%m%d"))

# combine and keep only ERA5 after WDFE5 ends
prec = do.call(rbind.data.frame, list(e1, e2)) 
table(prec$dataset, lubridate::year(prec$date))
prec = prec %>% dplyr::rename("precip"=value)
prec$precip = prec$precip * lubridate::days_in_month(prec$date) # convert to total accumulated over month (* num days)

# potential evapotranspiration (ERA5-LAND) in metres
pev = terra::rast("./data/climate/era5-land/potential_evaporation/meanpotentialevap_month_viet_19812020.nc")
pevb = terra::rast("./data/climate/era5-land/potential_evaporation/meanpotentialevap_month_viet_2020.nc")
pevc = terra::rast("./data/climate/era5-land/potential_evaporation/meanpotentialevap_month_viet_2021.nc")
pev = c(pev, pevb, pevc)
names(pev) = paste("PEV", as.Date(terra::time(pev)), sep="_")

pevs = exactextractr::exact_extract(pev, shp, fun="mean", progress=FALSE) %>%
  cbind(st_drop_geometry(shp[ , c("areaid", "areanameen", "areaprovin")])) %>%
  reshape2::melt(id.vars=c("areaid", "areanameen", "areaprovin")) %>%
  dplyr::mutate(date = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 2)),
                dataset = "ERA5-LAND") %>%
  dplyr::select(-variable) %>%
  dplyr::rename("pev_era5"=value)
pevs$date = as.Date(pevs$date, "%Y-%m-%d")
pevs$pev_era5 = pevs$pev_era5 * 1000 # convert to mean mm/day
pevs$pev_era5 = pevs$pev_era5 * lubridate::days_in_month(pevs$date) # convert to total accumulated over month (* num days)
pevs$pev_era5 = -pevs$pev_era5 # transform because convention in ERA5 is to express evaporation as negative (i.e. invert)

# mean monthly temperature (from ERA5-LAND to approximate PET; convert from K to Celcius)
# tmean = brick("./data/climate/era5-land/Vietnam/tmean_monthly/meantemp_month_viet_19812020.nc")
# tmeanb = brick("./data/climate/era5-land/Vietnam/tmean_monthly/meantemp_month_viet_20202021.nc")
# tmean = stack(tmean, tmeanb)
# names(tmean) = paste("Tmean_", substr(names(tmean), 2, 40), sep="")
# tmean = tmean - 273.15
# tmeans = exactextractr::exact_extract(tmean, shp, fun="mean", progress=FALSE) %>%
#   cbind(st_drop_geometry(shp[ , c("areaid", "areanameen", "areaprovin")])) %>%
#   reshape2::melt(id.vars=c("areaid", "areanameen", "areaprovin")) %>%
#   dplyr::mutate(date = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 2)),
#                 dataset = "ERA5-LAND") %>%
#   dplyr::select(-variable) %>%
#   dplyr::rename("tmean_era5"=value)
# tmeans$date = as.Date(tmeans$date, "%Y.%m.%d")

# spei_dat
dd = left_join(prec[ , 1:5], pevs[ , 1:5])
#dd = left_join(dd, tmeans[ , 1:5])

# access latitude
coords = cbind(shp$areaid, as.data.frame(st_coordinates(st_centroid(shp)))) %>%
  dplyr::rename("areaid"=1, "long"=2, "lat"=3)
dd = left_join(dd, coords)
dd = dd %>%
  dplyr::arrange(areaid, date)



# ------------------- coastal areas with NA values because do not cross land -------------------

# na_areas are almost all offshore
na_areas = unique(dd$areaid[ is.nan(dd$pev_era5) ])
# plot(shp$geometry, border="grey90")
# plot(shp$geometry[ shp$areaid %in% na_areas  ], col="red", border="black", add=TRUE)

# main issues are the coastal cities; da nang; for now leave others out as offshore are not an issue
# 2nd is Son Tra in Da Nang, areaid 70361l replace with Thanh Khe district
dd$precip[ dd$areaid == 70361 ] = dd$precip[ dd$areaid == 70364  ]
dd$pev_era5[ dd$areaid == 70361 ] = dd$pev_era5[ dd$areaid == 70364  ]
#dd$tmean_era5 [ dd$areaid == 70361 ] = dd$tmean_era5 [ dd$areaid == 70364  ]

# remove NA districts
dd = dd[ !is.nan(dd$pev_era5), ]

# all years with complete climate data
dd = dd[ lubridate::year(dd$date)>1980, ]


# --------------- calculate components of spei using 'SPEI' package -----------------

# per district
calcSPEI = function(x){
  
  # reporting
  print(x)
  
  # data for specific district
  dx = dd[ dd$areaid == unique(dd$areaid)[x], ] %>%
    dplyr::arrange(date)
  
  # calculate PET using thornwaite method to compare to ERA5 estimates
  # th = thornthwaite(Tave=dx$tmean_era5, lat=dx$lat[1])
  # dx$pet_thorn = as.vector(th)
  # #comparison w/ era5 estimate; much more variability with ERA5
  # ggplot(dx) + geom_line(aes(date, pet_thorn), col="blue") + geom_line(aes(date, pev_era5), col="red")

  # SPEI is calculated on Pi - PETi (i.e. monthly precip - monthly PET, which is a simple measure of water surplus/deficit in a given month)
  # i.e. the amount of water lost to evaporation that is not replenished by precip; or the amount of water falling that is not lost by evaporation
  # https://journals.ametsoc.org/view/journals/clim/23/7/2009jcli2909.1.xml
  # values are then calculated at different time scales (i.e. number of months ), as with the SPI
  
  # calculate water imbalance (D) and then calculate SPEI in 1, 3, 6, and 12 month windows
  dx$Di = dx$precip - dx$pev_era5
  dx$spei_1 = as.vector(SPEI::spei(dx$Di, scale = 1, fit="max-lik")$fitted)
  dx$spei_3 = as.vector(SPEI::spei(dx$Di, scale = 3, fit="max-lik")$fitted)
  dx$spei_6 = as.vector(SPEI::spei(dx$Di, scale = 6, fit="max-lik")$fitted)
  dx$spei_12 = as.vector(SPEI::spei(dx$Di, scale = 12, fit="max-lik")$fitted)
  
  # p1 = dx %>%
  #   dplyr::select(date, precip, pev_era5, Di) %>%
  #   dplyr::rename("Precipitation_WFDE5"=2, "PotentialEvapotranspiration_ERA5"=3, "Precip_Minus_PET"=4) %>%
  #   reshape2::melt(id.vars=1) %>%
  #   dplyr::mutate(variable = factor(variable, levels=c("Precipitation_WFDE5", "PotentialEvapotranspiration_ERA5", "Precip_Minus_PET"), ordered=TRUE)) %>%
  #   ggplot() + 
  #   geom_line(aes(date, value, col=variable), size=0.5) + 
  #   facet_wrap(~variable, ncol=1, scales="free_y") + 
  #   theme_minimal() +
  #   theme(legend.position="none", strip.text=element_text(size=14)) + 
  #   scale_color_viridis_d(end=0.75)
  # ggsave(p1, file="./spei_hoankiem_precipindicators.jpeg", dpi=600, units="in", width=9, height=7.5, scale=0.95)
  # 
  # p2 = dx %>%
  #   dplyr::select(date, spei_2, spei_3, spei_6, spei_12) %>%
  #   dplyr::rename_all(toupper) %>%
  #   dplyr::rename("date"=DATE) %>%
  #   reshape2::melt(id.vars=1) %>%
  #   ggplot() + 
  #   geom_line(aes(date, value, col=variable), size=0.5) + 
  #   facet_wrap(~variable, ncol=1, scales="free_y") + 
  #   theme_minimal() +
  #   theme(legend.position="none", strip.text=element_text(size=14)) + 
  #   scale_color_viridis_d(end=0.7, option="magma") +
  #   geom_hline(yintercept=0, lty=2)
  # ggsave(p2, file="./spei_hoankiem_spei.jpeg", dpi=600, units="in", width=9, height=9, scale=0.95)
  
  # return
  return(dx[ , c("areaid", "areanameen", "areaprovin", "date", "spei_1", "spei_3", "spei_6", "spei_12")])
}

# run and return
dds = do.call(rbind.data.frame, lapply(1:n_distinct(dd$areaid), calcSPEI))

distx = "Chau Thanh A"
pa = ggplot(dds[ dds$areanameen == distx, ]) +
  geom_line(aes(date, spei_1), col="grey20") +
  geom_line(aes(date, spei_6), col="red", size=1) +
  geom_hline(lty=2, yintercept=0) + theme_minimal()
pb = ggplot(dd[ dd$areanameen == distx, ]) +
  geom_line(aes(date, precip)) +
  geom_hline(lty=2, yintercept=0) + theme_minimal()
gridExtra::grid.arrange(pa, pb, ncol=1)

pa = ggplot(dds[ dds$areanameen == distx, ]) +
  geom_line(aes(date, spei_1), col="grey20") +
  geom_line(aes(date, spei_12), col="red", size=1) +
  geom_hline(lty=2, yintercept=0) + theme_minimal()
pb = ggplot(dd[ dd$areanameen == distx, ]) +
  geom_line(aes(date, precip)) +
  geom_hline(lty=2, yintercept=0) + theme_minimal()
gridExtra::grid.arrange(pa, pb, ncol=1)




# =================== crude precipitation: within focal month and in 2-3 month windows prior ====================

# rescale precip to mean monthly
prec$precip = prec$precip / lubridate::days_in_month(prec$date) # convert to total accumulated over month (* num days)

# group by areaid and calculate across lagged window including focal month
# prec = prec %>%  
#   dplyr::arrange(areaid, date) %>%
#   dplyr::group_by(areaid) %>%
#   dplyr::mutate(precip_2mwindow = data.table::frollmean(precip, 2, align="right"),
#                 precip_3mwindow = data.table::frollmean(precip, 3, align="right"))

# combine with dds into preciptiation indicators df
#indicators = left_join(dds, prec[ , c("areaid", "date", "precip", "precip_2mwindow", "precip_3mwindow")])
indicators = left_join(dds, prec[ , c("areaid", "date", "precip")])
write.csv(indicators, "./output/covariates/PrecipitationIndicators_VietAll_SPEI_MergedSHP2.csv", row.names=FALSE)




