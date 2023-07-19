


# ================= R script to access and process era5 reanalysis data =================

# project dir
PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

# # dependencies
library(raster)
library(rgdal)
library(sf)
library(ecmwfr)
library(stringr)
library(lubridate)


# ================= key objects =================

# districts shapefile for Vietnam
shp = sf::st_read("./data/shapefiles/dengue_districts_shapefile.shp")




# =============== process temperature data to derive daily variables ===================

# location for accessing/saving hourly temperature data (external)
temp_loc = "E:/ResearchProjects/202007_lshtm_dengue/analysis/data/climate/era5-land/Vietnam/2m_temperature/"
files = list.files(paste(temp_loc, "hourly/", sep=""), pattern=".nc", full.names=TRUE)
files = files[ -grep("19891996|19801988", files)] 
files =  files[25]

# function to convert hourly data into daily Tmean, Tmin, Tmax and diurnal range
deriveDailyTemperatureMetrics = function(x){
  
  # read-out
  print(paste("Processing...", x, sep=""))
  
  # read annual raster, convert to celcius, and create metadata table
  rx = stack(x)
  meta = do.call(rbind.data.frame, strsplit(substr(names(rx), 2, 30), "[.]"))
  names(meta) = c("year", "month", "day", "hour", "minute", "second")
  meta$date = paste(meta$year, meta$month, meta$day, sep="-")
  meta$day_id = as.numeric(as.factor(meta$date))
  
  # calculate daily mean, minimum, maximum
  dailyTempCalc = function(day_id, func=mean){ 
    if(day_id %in% seq(10, 10^6, by=10)){ cat(paste(day_id, "...", sep=""))}
    dd = rx[[ which(meta$day_id == day_id) ]]
    dd = raster::calc(dd, func)
    names(dd) = paste("X", meta$date[ meta$day_id == day_id][1], sep="")
    return(dd)
  }
  dailyMean = do.call(raster::stack, lapply(unique(meta$day_id), dailyTempCalc, func=mean)); names(dailyMean) = paste("Tmean", names(dailyMean), sep="_")
  dailyMin = do.call(raster::stack, lapply(unique(meta$day_id), dailyTempCalc, func=min)); names(dailyMin) = paste("Tmin", names(dailyMin), sep="_")
  dailyMax = do.call(raster::stack, lapply(unique(meta$day_id), dailyTempCalc, func=max)); names(dailyMax) = paste("Tmax", names(dailyMax), sep="_")
  dailyDiurnalRange = dailyMax - dailyMin; names(dailyDiurnalRange) = paste("Tdrange", substr(names(dailyMean), 7, 30), sep="_")
  
  # result
  resx = list(
    dailyMean = dailyMean,
    dailyMin = dailyMin,
    dailyMax = dailyMax,
    dailyDiurnalRange = dailyDiurnalRange
  )
  resx
}

# derive temperature metrics and create separate stacks
erat = lapply(files, deriveDailyTemperatureMetrics)
tmean = do.call(raster::stack, lapply(erat, function(x) x$dailyMean))
tmin = do.call(raster::stack, lapply(erat, function(x) x$dailyMin))
tmax = do.call(raster::stack, lapply(erat, function(x) x$dailyMax))
tdrange = do.call(raster::stack, lapply(erat, function(x) x$dailyDiurnalRange))

# save 
# writeRaster(tmean, file=paste(temp_loc, "daily/Tmean_1997-01-01_2020-02-01_daily_era5land_Vietnam.grd", sep=""), format="raster")
# writeRaster(tmin, file=paste(temp_loc, "daily/Tmin_1997-01-01_2020-02-01_daily_era5land_Vietnam.grd", sep=""), format="raster")
# writeRaster(tmax, file=paste(temp_loc, "daily/Tmax_1997-01-01_2020-02-01_daily_era5land_Vietnam.grd", sep=""), format="raster")
# writeRaster(tdrange, file=paste(temp_loc, "daily/Tdrange_1997-01-01_2020-02-01_daily_era5land_Vietnam.grd", sep=""), format="raster")

# writeRaster(tmean, file=paste(temp_loc, "daily/Tmean_2020_daily_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
# writeRaster(tmin, file=paste(temp_loc, "daily/Tmin_2020_daily_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
# writeRaster(tmax, file=paste(temp_loc, "daily/Tmax_2020_daily_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
# writeRaster(tdrange, file=paste(temp_loc, "daily/Tdrange_2020_daily_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)

# writeRaster(tmean, file=paste(temp_loc, "daily/Tmean_19801996_daily_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
# writeRaster(tmin, file=paste(temp_loc, "daily/Tmin_19801996_daily_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
# writeRaster(tmax, file=paste(temp_loc, "daily/Tmax_19801996_daily_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
# writeRaster(tdrange, file=paste(temp_loc, "daily/Tdrange_19801996_daily_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)

writeRaster(tmean, file=paste(temp_loc, "daily/Tmean_2021_daily_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
writeRaster(tmin, file=paste(temp_loc, "daily/Tmin_2021_daily_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
writeRaster(tmax, file=paste(temp_loc, "daily/Tmax_2021_daily_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
writeRaster(tdrange, file=paste(temp_loc, "daily/Tdrange_2021_daily_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)



# ================== create monthly temperature variables from daily stacks ==================

# function creates monthly mean rasters from daily stacks produced above
#' @param rr raster stack

createMonthlyMeans = function(rr){
  
  # metric
  metric = lapply(strsplit(names(rr), "_"), "[", 1)[[1]]
  
  # month id
  meta = do.call(rbind.data.frame, strsplit(substr(unlist(lapply(strsplit(names(rr), "_"), "[", 2)), 2, 30), "[.]"))
  names(meta) = c("year", "month", "day")
  meta$month_id = as.numeric(as.factor(paste(meta$year, meta$month)))
  
  # calculate monthly mean, convert to celcius
  calcMonthRas = function(month_id){
    rx = mean(rr[[ which(meta$month_id == month_id) ]])
    if(metric != "Tdrange"){ rx = rx - 273.15 }
    names(rx) = paste0(metric, "_", meta$year[ meta$month_id == month_id ][1], "-", meta$month[ meta$month_id == month_id ][1], "-01", sep="")
    rx
  }
  res = do.call(raster::stack, lapply(unique(meta$month_id), calcMonthRas))
  return(res)
}

# create monthly means and save as .grd files
tmeanm = createMonthlyMeans(tmean)
tminm = createMonthlyMeans(tmin)
tmaxm = createMonthlyMeans(tmax)
tdrangem = createMonthlyMeans(tdrange)

# save to repo
# writeRaster(tmeanm, file=paste(temp_loc, "monthly/Tmean_1997-01-01_2020-02-01_month_era5land_Vietnam.grd", sep=""), format="raster")
# writeRaster(tminm, file=paste(temp_loc, "monthly/Tmin_1997-01-01_2020-02-01_month_era5land_Vietnam.grd", sep=""), format="raster")
# writeRaster(tmaxm, file=paste(temp_loc, "monthly/Tmax_1997-01-01_2020-02-01_month_era5land_Vietnam.grd", sep=""), format="raster")
# writeRaster(tdrangem, file=paste(temp_loc, "monthly/Tdrange_1997-01-01_2020-02-01_month_era5land_Vietnam.grd", sep=""), format="raster")

# writeRaster(tmeanm, file=paste(temp_loc, "monthly/Tmean_2020_month_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
# writeRaster(tminm, file=paste(temp_loc, "monthly/Tmin_2020_month_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
# writeRaster(tmaxm, file=paste(temp_loc, "monthly/Tmax_2020_month_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
# writeRaster(tdrangem, file=paste(temp_loc, "monthly/Tdrange_2020_month_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)

# writeRaster(tmeanm, file=paste(temp_loc, "monthly/Tmean_19801996_month_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
# writeRaster(tminm, file=paste(temp_loc, "monthly/Tmin_19801996_month_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
# writeRaster(tmaxm, file=paste(temp_loc, "monthly/Tmax_19801996_month_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)
# writeRaster(tdrangem, file=paste(temp_loc, "monthly/Tdrange_19801996_month_era5land_Vietnam.grd", sep=""), format="raster", overwrite=TRUE)

writeRaster(tmeanm, file="./data/climate/era5_land/2m_temperature/Tmean_2021_month_era5land_Vietnam.grd", format="raster", overwrite=TRUE)
writeRaster(tminm, file="./data/climate/era5_land/2m_temperature/Tmin_2021_month_era5land_Vietnam.grd", format="raster", overwrite=TRUE)
writeRaster(tmaxm, file="./data/climate/era5_land/2m_temperature/Tmax_2021_month_era5land_Vietnam.grd", format="raster", overwrite=TRUE)
writeRaster(tdrangem, file="./data/climate/era5_land/2m_temperature/Tdrange_2021_month_era5land_Vietnam.grd", format="raster", overwrite=TRUE)


# as.data.frame(tmaxm, xy=TRUE) %>%
#   reshape2::melt(id.vars = 1:2) %>%
#   dplyr::mutate(variable = substr(variable, 7, 30)) %>%
#   ggplot() + 
#   geom_raster(aes(x, y, fill=value)) + 
#   theme_bw() + 
#   facet_wrap(~variable) + 
#   scale_fill_gradientn(colors=viridisLite::turbo(200), na.value="white")
  
# some initial viz of temperature
# spplot(tmeanm[[grep("1998", names(tmeanm))]])
# spplot(tmeanm[[grep("2018", names(tmeanm))]])
# spplot(tminm[[grep("1998", names(tminm))]])
# spplot(tminm[[grep("2018", names(tminm))]])
# spplot(tdrangem[[grep("2018", names(tdrangem))]])

# r1 = stack(paste(temp_loc, "monthly/Tmin_1997-01-01_2020-02-01_month_era5land_Vietnam.grd", sep=""))
# r2 = stack(paste(temp_loc, "monthly/Tdrange_1997-01-01_2020-02-01_month_era5land_Vietnam.grd", sep=""))
# 
# aa = as.data.frame(r2[[grep("1998|1999", names(r2))]], xy=TRUE) %>%
#   reshape2::melt(id.vars=1:2)
# ggplot(aa) + 
#   geom_raster(aes(x, y, fill=value)) + 
#   facet_wrap(~variable) + 
#   scale_fill_viridis_c(option="magma", direction=-1) + 
#   coord_fixed()
# 
# aa = as.data.frame(r1[[grep("1998", names(r1))]], xy=TRUE) %>%
#   reshape2::melt(id.vars=1:2)
# ggplot(aa) + 
#   geom_raster(aes(x, y, fill=value)) + 
#   facet_wrap(~variable) + 
#   scale_fill_viridis_c(option="magma", direction=-1) + 
#   coord_fixed()



# =================== extract monthly climate variables for Vietnamese districts ====================

# read and initial formatting of climate (precip, wind, temperature)
clim_loc = getwd()

# prec
pp = list.files(paste(clim_loc, "/data/climate/era5-land/total_precipitation/", sep=""), full.names=TRUE)
prec = do.call(c, lapply(pp, terra::rast))
names(prec) = paste("PrecipMean", as.character(as.Date(terra::time(prec))), sep="_")
shp2 = sf::st_transform(shp, st_crs(prec))
px = exactextractr::exact_extract(prec, shp2, fun="mean")
px = cbind(shp %>% st_drop_geometry(), px)
px = reshape2::melt(px, id.vars=1:5)
px$Date = as.Date(unlist(lapply(strsplit(substr(px$variable, 6, 30), "_"), "[", 2)), "%Y-%m-%d")
px$variable = unlist(lapply(strsplit(substr(px$variable, 6, 30), "_"), "[", 1))
px = px %>% dplyr::select(-areanamelo, -areatypeid, -areaprovin)

# wind speed
wind_u = list.files(paste(clim_loc, "/data/climate/era5-land/10m_u_component_of_wind/", sep=""), full.names=TRUE)
wind_u = do.call(c, lapply(wind_u, terra::rast))
wind_v = list.files(paste(clim_loc, "/data/climate/era5-land/10m_v_component_of_wind/", sep=""), full.names=TRUE)
wind_v = do.call(c, lapply(wind_v, terra::rast))
wspeed = sqrt(wind_u^2 + wind_v^2)
names(wspeed) = paste("WindSpeed", as.character(as.Date(terra::time(wind_u))), sep="_")
wx = exactextractr::exact_extract(wspeed, shp2, fun="mean")
wx = cbind(shp %>% st_drop_geometry(), wx)
wx = reshape2::melt(wx, id.vars=1:5)
wx$Date = as.Date(unlist(lapply(strsplit(substr(wx$variable, 6, 30), "_"), "[", 2)), "%Y-%m-%d")
wx$variable = unlist(lapply(strsplit(substr(wx$variable, 6, 30), "_"), "[", 1))
wx = wx %>% dplyr::select(-areanamelo, -areatypeid, -areaprovin)

# temperature
tmean = do.call(raster::stack,
                lapply(list(paste(clim_loc, "data/climate/era5-land/2m_temperature/Tmean_1997-01-01_2020-02-01_month_era5land_Vietnam.grd", sep=""), 
                            paste(clim_loc, "data/climate/era5-land/2m_temperature/Tmean_2020_month_era5land_Vietnam.grd", sep=""), 
                            paste(clim_loc, "data/climate/era5-land/2m_temperature/Tmean_2021_month_era5land_Vietnam.grd", sep="")),
                       raster::brick))

tmin = do.call(raster::stack,
               lapply(list(paste(clim_loc, "data/climate/era5-land/2m_temperature/Tmin_1997-01-01_2020-02-01_month_era5land_Vietnam.grd", sep=""), 
                           paste(clim_loc, "data/climate/era5-land/2m_temperature/Tmin_2020_month_era5land_Vietnam.grd", sep=""), 
                           paste(clim_loc, "data/climate/era5-land/2m_temperature/Tmin_2021_month_era5land_Vietnam.grd", sep="")),
                      raster::brick))


tmax = do.call(raster::stack,
               lapply(list(paste(clim_loc, "data/climate/era5-land/2m_temperature/Tmax_1997-01-01_2020-02-01_month_era5land_Vietnam.grd", sep=""), 
                           paste(clim_loc, "data/climate/era5-land/2m_temperature/Tmax_2020_month_era5land_Vietnam.grd", sep=""), 
                           paste(clim_loc, "data/climate/era5-land/2m_temperature/Tmax_2021_month_era5land_Vietnam.grd", sep="")),
                      raster::brick))

tdrange = do.call(raster::stack,
                  lapply(list(paste(clim_loc, "data/climate/era5-land/2m_temperature/Tdrange_1997-01-01_2020-02-01_month_era5land_Vietnam.grd", sep=""), 
                              paste(clim_loc, "data/climate/era5-land/2m_temperature/Tdrange_2020_month_era5land_Vietnam.grd", sep=""), 
                              paste(clim_loc, "data/climate/era5-land/2m_temperature/Tdrange_2021_month_era5land_Vietnam.grd", sep="")),
                         raster::brick))

# extract for temperature
climr = raster::stack(tmean, tmin, tmax, tdrange)
shp2 = sf::st_transform(shp, st_crs(climr))
ee = exactextractr::exact_extract(climr, shp2, fun="mean")
ee = cbind(shp %>% st_drop_geometry(), ee)
ee2 = reshape2::melt(ee, id.vars=1:5)
ee2$Date = as.Date(unlist(lapply(strsplit(substr(ee2$variable, 6, 30), "_"), "[", 2)), "%Y.%m.%d")
ee2$variable = unlist(lapply(strsplit(substr(ee2$variable, 6, 30), "_"), "[", 1)) 
ee2 = ee2 %>% dplyr::select(-areanamelo, -areatypeid, -areaprovin)

# combine
clim_data = do.call(rbind.data.frame, list(px, wx, ee2))

write.csv(clim_data, "./output/covariates/VietDistricts_ERA5Land_monthly_MergedSHP2.csv", row.names=FALSE)

# ee2 %>%
#   dplyr::filter(areanameen == "Hoan Kiem") %>%
#   ggplot() +
#   geom_line(aes(Date, value)) +
#   facet_wrap(~variable, scales="free_y")

