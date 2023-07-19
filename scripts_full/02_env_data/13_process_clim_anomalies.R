
# ======================= Process climate anomaly indices (ONI and DMI) =========================

#
setwd("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/")
pacman::p_load("rgdal", "sp", "spdep", "INLA", "corrplot", "usdm", "dplyr", "reshape2", "ggplot2", "ncdf4")
source("./code/viet_dengue_districts/scripts/00_functions.R")

# Oceanic Nino Index (1950-present), which is calculated as rolling mean of 3 month windows
# Set "date" as final month of time window (because this means that dengue has occurred across all 3 months of this window, i.e. no forward lag)
oni = read.table("./data/climate/anomalies/oceanicnino_cpc_noaa/oni.ascii_15102020.txt", header=TRUE, sep="") %>%
  rename("season" = SEAS, "year" = YR, "total" = TOTAL, "oni_3.4" = ANOM)
foo = oni[ 1:12, "season", drop=FALSE]
#foo$month = c(2:12, 1)
#foo$startmonth = c(12, 1:11)
foo$month = 1:12
oni = left_join(oni, foo, by="season")
oni$date = as.Date(paste(oni$year, oni$month, "01", sep="-"), "%Y-%m-%d")
oni$date = oni$date %m+% months(1)

# indian ocean dipole; dipole mode index
# https://stateoftheocean.osmc.noaa.gov/sur/ind/dmi.php
# read nc and convert days since 1900-01-01 to date
nc = nc_open("./data/climate/anomalies/indianoceandipole_noaa/dmi.nc")
dmi = data.frame(day = ncvar_get(nc, "WEDCEN2"),
                 dmi = ncvar_get(nc, "DMI")) %>%
  dplyr::mutate(date = as.Date(as.POSIXct(day*24*60*60, origin = "1900-01-01 00:00:00", tz="UTC"))) %>%
  dplyr::mutate(year = lubridate::year(date),
                month = lubridate::month(date),
                yearmonth = paste(year, month, sep="_"))
nc_close(nc)

# calculate monthly means and then 3 month rolling mean covering preceding months (i.e. to harmonise with ONI)
dmi_m = dmi %>%
  group_by(yearmonth) %>%
  dplyr::summarise(year = year[1],
                   month = month[1],
                   date = as.Date(paste(year, month, "01", sep="-")),
                   dmi = mean(dmi)) %>%
  dplyr::arrange(date) %>%
  dplyr::mutate(dmi_3m = data.table::frollmean(dmi, 3, align="right"))

# combine into anomalies dataframe with key variables (startdate, season)
anom = left_join(oni[ , c("date", "season", "oni_3.4")],
                 dmi_m[ , c("date", "dmi_3m")],
                 by="date")
write.csv(anom, "./code/viet_dengue_districts/output/covariates/ClimateAnomalies_ONI_DMI.csv", row.names=FALSE)




# some visualisation
# ggplot(anom[ anom$date > as.Date("1980-01-01"), ]) + 
#   geom_line(aes(date, oni_3.4), col="blue") +
#   geom_line(aes(date, dmi_3m), col="red")
  


