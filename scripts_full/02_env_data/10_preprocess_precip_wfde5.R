

# ================ Process WFDE5 hourly precipitation data into monthly and daily summaries ===============

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

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))

# districts shapefile for Vietnam
shp = st_read("./data/shapefiles/dengue_districts_shapefile.shp")
shp = shp[ ! shp$areanameen %in% c("Truong Sa", "Hoang Sa"), ] # remove offshore districts


# =============== process WFDE5 precipitation data (both CRU and CRU+GPCC-corrected) ===================

# function to convert hourly data into monthly precipitation values
#' @param x location of .nc file of hourly data
 
deriveMonthlyPrecip = function(x){
  
  # read-out
  print(paste("Processing...", x, sep=""))

  # read monthly stack of hourly precip
  # native units kg/m2/sec, so multiply by num sec to convert to total mm in hour-long window
  rx = terra::rast(x)
  #ry = raster::stack(x)
  #rx = crop(rx, extent(shp)+0.25)
  rx = (rx * 3600)/1000 # multiply by num seconds in an hour, and divide by 1000 (to give total precip in metres per hour)
  
  # metadata
  # meta = do.call(rbind.data.frame, strsplit(substr(names(rx), 2, 30), "[.]"))
  # names(meta) = c("year", "month", "day", "hour", "minute", "second")
  # meta$date = paste(meta$year, meta$month, meta$day, sep="-")
  # meta$day_id = as.numeric(as.factor(meta$date))
  
  # # # calculate daily total precip (mm)
  # dailyPrecipCalc = function(day_id, func=sum){
  #   if(day_id %in% seq(10, 10^6, by=10)){ cat(paste(day_id, "...", sep=""))}
  #   dd = rx[[ which(meta$day_id == day_id) ]]
  #   dd = raster::calc(dd, func)
  #   names(dd) = paste("X", meta$date[ meta$day_id == day_id][1], sep="")
  #   return(dd)
  # }
  # dailyPrecip = do.call(raster::stack, lapply(unique(meta$day_id), dailyPrecipCalc, func=sum))
  
  # calculate monthly total then average
  monthlyPrecip = sum(rx)
  monthlyPrecip = monthlyPrecip / dplyr::n_distinct(as.Date(terra::time(rx)))
  names(monthlyPrecip) = paste("PrecipMean_", as.Date(terra::time(rx[[1]])), sep="")
  
  monthlyPrecip
}


# 1. derive for CRU-GPCC bias-corrected data (up to 2019; stored externally)
prec_loc = "C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/data/climate/wfde5/precipitation_biascorrected/cru_gpcc/"
files = list.files(paste(prec_loc, "hourly/", sep=""), pattern=".nc", full.names=TRUE)
prec = lapply(files, deriveMonthlyPrecip)
precr = do.call(c, prec)
precr_ras = raster::stack(precr)
precr_ras = precr_ras[[ order(names(precr_ras)) ]]
writeRaster(precr_ras, file="./data/climate/wfde5/PrecipMean_1979-01-01_2019-12-01_monthly_wfde5gpcc_v2.1_Vietnam.grd", format="raster", overwrite=TRUE)



# ============ OLD CODE ==============

# # ==================== run extraction for districts ========================
# 
# # extract for districts
# ee = exactextractr::exact_extract(precr, shp, fun="mean")
# ee = cbind(shp %>% st_drop_geometry(), ee)
# ee = reshape2::melt(ee, id.vars=1:5)
# ee$Date = as.Date(unlist(lapply(strsplit(substr(ee$variable, 6, 30), "_"), "[", 2)), "%Y.%m.%d")
# ee$variable = unlist(lapply(strsplit(substr(ee$variable, 6, 30), "_"), "[", 1))                
# ee$Dataset = "WFDE5"
# ee = ee %>% dplyr::select(-areanamelo, -areatypeid, -areaprovin)
# ee$value[ is.nan(ee$value) ] = NA
# write.csv(ee, "./code/viet_dengue_districts/output/covariates/VietDistricts_WFDE5precip_monthly.csv", row.names=FALSE)
# 


# ==================== derive daily precipitation estimates (used for num/num consecutive dry days) ======================

# # calculate daily precip for num dry days
# deriveDailyPrecip = function(x){
#   
#   # read-out
#   print(paste("Processing...", x, sep=""))
#   
#   # read monthly stack of hourly precip
#   # native units kg/m2/sec, so multiply by num sec to convert to total mm in hour-long window
#   rx = stack(x)
#   rx = crop(rx, extent(shp)+0.25)
#   rx = (rx * 3600) # multiply by num seconds in an hour (total precip in mm)
#   
#   # metadata
#   meta = do.call(rbind.data.frame, strsplit(substr(names(rx), 2, 30), "[.]"))
#   names(meta) = c("year", "month", "day", "hour", "minute", "second")
#   meta$date = paste(meta$year, meta$month, meta$day, sep="-")
#   meta$day_id = as.numeric(as.factor(meta$date))
#   
#   # calculate daily total precip (m)
#   dailyPrecipCalc = function(day_id, func=sum){ 
#     if(day_id %in% seq(10, 10^6, by=10)){ cat(paste(day_id, "...", sep=""))}
#     dd = rx[[ which(meta$day_id == day_id) ]]
#     dd = raster::calc(dd, func)
#     names(dd) = paste("X", meta$date[ meta$day_id == day_id][1], sep="")
#     return(dd)
#   }
#   dailyPrecip = do.call(raster::stack, lapply(unique(meta$day_id), dailyPrecipCalc, func=sum))
#   names(dailyPrecip) = paste("Precip_mm_", substr(names(dailyPrecip), 2, 30), sep="")
#   return(dailyPrecip)
#   
#   # calculate mean monthly and return
#   # monthlyPrecip = mean(dailyPrecip); names(monthlyPrecip) = paste("PrecipMean_", meta$date[1], sep="")
#   # monthlyPrecip
# }
# 
# prec_loc = "./data/climate/wfde5/precipitation_biascorrected/cru_gpcc/"
# files = list.files(paste(prec_loc, "hourly_1981/", sep=""), pattern=".nc", full.names=TRUE)
# precd = lapply(files, deriveDailyPrecip)
# precdr = do.call(raster::stack, precd)
# writeRaster(precdr, file=paste(prec_loc, "daily/PrecipDaily_1981-01-01_2018-12-01_daily_wfde5gpcc_Vietnam_20201020.grd", sep=""), format="raster")


# # ================== visualise =================
# 
# library(ggplot2)
# plotClimTS = function(dist_name){
#   ggplot(ee[ ee$areanameen == dist_name, ]) +
#     geom_line(aes(Date, value, col=variable, group=variable)) +
#     theme_minimal() + 
#     facet_wrap(~variable, scales="free_y") +
#     ggtitle(dist_name) + 
#     theme(plot.title=element_text(size=13, hjust=0.5))
# }
# plotClimTS("Dong Da")
# plotClimTS("Ha Dong")
# plotClimTS("Chuong My")
# plotClimTS("Buon Ma Thuot")
# plotClimTS("Cu Mgar")
# plotClimTS("Nha Trang")
# plotClimTS("Ninh Hoa")
# plotClimTS("Bien Hoa")
# plotClimTS("Long Thanh")

