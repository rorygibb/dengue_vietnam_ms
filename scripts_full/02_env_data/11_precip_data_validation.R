

# ======================== Some rainfall data validation: compare ERA5 and WFDE5 GPCC+CRU datasets with station data =====================

# take home from this subanalysis is that precipitation data from WFDE5 are very accurate with respect to local station data
# use these to calculate rainfall indicators for models; in next script calculate both precipitation crude and SPEI (accounts for PEV)

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
library(ggplot2)
library(dplyr)
source("./scripts/00_plot_themes.R")
offshore_areas = c(70154, 70339, 70273, 70355, 70698)

# districts shapefile for Vietnam
shp = sf::st_read("./data/shapefiles/dengue_districts_shapefile.shp") %>%
  dplyr::filter(!areaid %in% offshore_areas)



# =============== Monthly rainfall at stations in cities/provinces from the Vietnam GSO ================

rainfall_stations = read.csv("./data/climate/vietnam_stations/rainfall_mm_monthly_stations_gso.csv")
names(rainfall_stations) = c("Year", "Locale", paste("X", 1:12, sep=""))
rainfall_stations = rainfall_stations %>%
  reshape2::melt(id.vars = 1:2) %>%
  dplyr::mutate(month = substr(variable, 2, 3))
rainfall_stations$date = as.Date(paste(rainfall_stations$Year, rainfall_stations$month, "01", sep="-"))
rain = rainfall_stations %>% 
  dplyr::select(-variable) %>%
  dplyr::rename("Precip_mm" = value)
rain$Precip_mm[ rain$Precip_mm == ".." ] = NA
rain$Precip_mm = as.numeric(rain$Precip_mm)




# ============== Name matching to polygons and access centroids ==================

# names
nm = data.frame(rain_locale = unique(rain$Locale))

# any that match exactly and combine with centroids
nm_match = nm[ nm$rain_locale %in% shp$areanameen, , drop=FALSE ]
nm_match$areanameen = nm_match$rain_locale

# non matched (set to central)
nm_nomatch = nm[ !nm$rain_locale %in% shp$areanameen, , drop=FALSE ]
nm_nomatch$areanameen = NA
nm_nomatch$areanameen[ nm_nomatch$rain_locale == "Ha Noi" ] = "Hoan Kiem" # central Ha Noi
nm_nomatch$areanameen[ nm_nomatch$rain_locale == "Bai Chay" ] = "Ha Long" # district containing Bai Chay
nm_nomatch$areanameen[ nm_nomatch$rain_locale == "Da Nang" ] = "Hai Chau" # central Da Nang
nm_nomatch$areanameen[ nm_nomatch$rain_locale == "Qui Nhon" ] = "Quy Nhon" # different spelling
nm_nomatch$areanameen[ nm_nomatch$rain_locale == "Playku" ] = "Pleiku" # different spelling
nm_match = rbind(nm_match, nm_nomatch)

# access centroids
shpm = shp[ shp$areanameen %in% nm_match$areanameen, c("areaid", "areanameen")]
coords = cbind(shpm[ , c("areaid", "areanameen")] %>% st_drop_geometry(), st_coordinates(st_centroid(shpm)))
nm_match = left_join(nm_match, coords, by=c("areanameen"))

# visualise distribution of points: pan-Vietnam
p1 = ggplot() + 
  geom_sf(data=shp, fill="grey80", col=NA) +
  geom_point(data=nm_match, aes(X, Y), col="darkred", pch=21, fill="coral2", size=3) + 
  maptheme
  



# ============= Extract ERA5 and WFDE5 precipitation data =================

# era5-land and wfde5
prec_era = terra::rast("./data/climate/era5-land/total_precipitation/meanprecip_month_viet_19972019.nc")
names(prec_era) = paste("X", as.character(as.Date(terra::time(prec_era))), sep="")
prec_wfde = raster::brick("./data/climate/wfde5/PrecipMean_1979-01-01_2019-12-01_monthly_wfde5gpcc_v2.1_Vietnam.grd")

# extract for era (sum of daily rainfall)
# units are m/day so multiply by 1000 for mm, and then multiply by number of days in the month
shpm = shp[ shp$areanameen %in% nm_match$areanameen, c("areaid", "areanameen")] %>%
  st_transform(st_crs(prec_era))
era = exactextractr::exact_extract(prec_era, shpm, fun="mean")
era = cbind(shpm %>% st_drop_geometry(), era)
era = era %>%
  reshape2::melt(id.vars = 1:2) %>%
  dplyr::mutate(date = as.Date(substr(variable, 7, 30), "%Y-%m-%d"),
                value = value * 1000) 
era$precip_mm = era$value * as.vector(lubridate::days_in_month(era$date))
era = era %>%
  dplyr::select(areaid, areanameen, date, precip_mm) %>%
  dplyr::mutate(dataset="ERA5-LAND")

# wfde5 cru+gpcc; same units so run same conversion
wf = exactextractr::exact_extract(prec_wfde, shpm, fun="mean")
wf = cbind(shpm %>% st_drop_geometry(), wf)
wf = wf %>%
  reshape2::melt(id.vars = 1:2) %>%
  dplyr::mutate(date = as.Date(substr(variable, 17, 30), "%Y.%m.%d"),
                value = value * 1000) 
wf$precip_mm = wf$value * as.vector(lubridate::days_in_month(wf$date))
wf = wf %>%
  dplyr::select(areaid, areanameen, date, precip_mm) %>%
  dplyr::mutate(dataset="WFDE5")



# =============== Combine with station data ================

rain = left_join(rain, nm_match, by=c("Locale" = "rain_locale"))
rainx = rain %>% 
  dplyr::select(areaid, areanameen, date, Precip_mm) %>%
  dplyr::rename("precip_mm" = Precip_mm) %>%
  dplyr::mutate(dataset = "Stations")

era = era %>% dplyr::filter(date %in% rainx$date)
wf = wf %>% dplyr::filter(date %in% rainx$date)

dd = rbind(rainx, era, wf)




# ============= Calculate error estimates (MAE) ===============

era2 = era %>% dplyr::rename("era5" = precip_mm) %>% dplyr::select(-dataset)
wf2 = wf %>% dplyr::rename("wfde5" = precip_mm) %>% dplyr::select(-dataset)
rainx2 = rainx %>% dplyr::rename("stations" = precip_mm) %>% dplyr::select(-dataset)
err_dat = left_join(era2, wf2)
err_dat = left_join(err_dat, rainx2)
err_dat = err_dat[ !is.na(err_dat$stations), ]

# error
err_dat$err_wfde5 = err_dat$wfde5  - err_dat$stations
err_dat$err_era5 = err_dat$era5  - err_dat$stations

# calculate MAE on two different datasets
mae_wf = mean(abs(err_dat$err_wfde5), na.rm=TRUE)
mae_era5 = mean(abs(err_dat$err_era5), na.rm=TRUE)

# and at different locales
mae_loc = err_dat %>%
  dplyr::group_by(areanameen) %>%
  dplyr::summarise(WFDE5 = mean(abs(err_wfde5), na.rm=TRUE),
                   ERA5 = mean(abs(err_era5), na.rm=TRUE))
mae_loc = rbind(data.frame(areanameen = "Vietnam", WFDE5 = mae_wf, ERA5 = mae_era5), mae_loc)
mae_loc = reshape2::melt(mae_loc, id.vars=1) %>%
  dplyr::mutate(variable = as.vector(variable),
                variable = replace(variable, variable=="ERA5", "ERA5-Land"))

fac_order = c(nm_match$areanameen[ order(nm_match$Y, decreasing = FALSE) ], "Vietnam" )
mae_loc$areanameen = factor(mae_loc$areanameen, levels=fac_order, ordered=TRUE)
p2 = ggplot(mae_loc) + 
  geom_bar(aes(areanameen, value, group=variable, fill=variable), stat="identity", position=position_dodge(), col="grey50", width=0.5) +
  theme_classic() + 
  xlab("Locale") + 
  ylab("Mean absolute difference\nfrom station observations (mm)") +
  scale_fill_viridis_d(option="viridis", begin=0.3, end=0.85, name="Dataset", guide=guide_legend(reverse=TRUE)) + 
  coord_flip() +
  theme(legend.position=c(0.9, 0.9))
# ggsave(px, file="./precipvalidation_denv_vietnam.png", dpi=600, units="in", width=6, height=3, device="png")


gridExtra::grid.arrange(p1, p2, ncol=2, widths=c(0.65, 1))


# ================ Visualise comparisons for some locations =================

dd$dataset = as.vector(dd$dataset)
dd$dataset[ dd$dataset == "ERA5-LAND" ] = "ERA5-Land"
dd$dataset = factor(dd$dataset, level=c("Stations", "ERA5-Land", "WFDE5"), ordered=TRUE)

p3 = dd %>%
  dplyr::filter(dataset == "WFDE5") %>%
  left_join(
    dd %>% dplyr::filter(dataset == "Stations") %>% dplyr::rename("p2"=precip_mm) %>% dplyr::select(-dataset)
  ) %>% 
  dplyr::filter(
    areanameen %in% c("Lai Chau", "Hoan Kiem", "Vinh", "Hue", "Nha Trang", "Da Lat")
  ) %>%
  dplyr::mutate(areanameen = factor(areanameen, levels=c("Lai Chau", "Hoan Kiem", "Vinh", "Hue", "Nha Trang", "Da Lat"), ordered=TRUE)) %>%
  dplyr::filter(date >= as.Date("2002-01-01")) %>%
  ggplot() + 
  geom_bar(aes(date, p2), stat="identity", fill="grey75", col="grey75", width=4) +
  geom_line(aes(date, precip_mm), alpha=1, size=0.5, col=viridis::viridis(200)[130]) +
  facet_wrap(~areanameen, ncol=1, scales="free_y") +
  theme_minimal() +
  xlab("Month") + ylab("Total precipitation (mm)") + 
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=8),
        strip.text = element_text(size=10))

# combine
pc = gridExtra::grid.arrange(p1, p2, p3, ncol=3, widths=c(0.65, 1, 1.2))
pc = ggpubr::as_ggplot(pc)  +
  cowplot::draw_plot_label(label = c("a", "b", "c"),
                           fontface = "bold", size = 22.5,
                           x = c(0.01, 0.24, 0.62), y = c(0.99, 0.99, 0.99))
ggsave(pc, file="./output/figures/SuppFigure_PrecipGroundTruthing.jpg", device="jpg", units="in", width=12, height=6, dpi=600)





# # ============== Some comparisons: initial viz shows really good performance of WFDE5 ==================
# 
# name = "Nha Trang"
# dd$dataset = factor(dd$dataset, level=c("Stations", "ERA5-LAND", "WFDE5"), ordered=TRUE)
# p1 = ggplot(dd[ dd$areanameen == name & dd$dataset != "WFDE5", ]) + 
#   geom_line(aes(date, precip_mm, group=dataset, col=dataset), alpha=0.7, size=1) +
#   facet_wrap(~areanameen) +
#   theme_minimal() +
#   scale_color_viridis_d(option="viridis", begin=0, end=0.5)
# p2 = ggplot(dd[ dd$areanameen == name & dd$dataset != "ERA5-LAND", ]) + 
#   geom_line(aes(date, precip_mm, group=dataset, col=dataset), alpha=0.7, size=1) +
#   facet_wrap(~areanameen) +
#   theme_minimal() +
#   scale_color_viridis_d(option="viridis", begin=0, end=0.7)
# gridExtra::grid.arrange(p1, p2, nrow=2)
# 
# 
# 
# 
# 
# 
# # ========================== TEMPERATURE: comparison of ERA5-LAND versus WFDE5 temperature data ======================
# 
# temp_stations = read.csv("./data/climate/vietnam_stations/temperature_annual_stations_gso.csv")
# names(temp_stations)[1] = "Locale"
# temp_stations = temp_stations %>%
#   reshape2::melt(id.vars = 1) %>%
#   dplyr::mutate(year = substr(variable, 2, 5))
# temp = temp_stations %>% 
#   dplyr::select(-variable) %>%
#   dplyr::rename("Tmean" = value)
# 
# # names
# nm = data.frame(temp_locale = unique(temp$Locale))
# nm_match = nm[ nm$temp_locale %in% shp$areanameen, , drop=FALSE ]
# nm_match$areanameen = nm_match$temp_locale
# 
# # non matched (set to central)
# nm_nomatch = nm[ !nm$temp_locale %in% shp$areanameen, , drop=FALSE ]
# nm_nomatch$areanameen = NA
# nm_nomatch$areanameen[ nm_nomatch$temp_locale == "Ha Noi" ] = "Hoan Kiem" # central Ha Noi
# nm_nomatch$areanameen[ nm_nomatch$temp_locale == "Bai Chay" ] = "Ha Long" # district containing Bai Chay
# nm_nomatch$areanameen[ nm_nomatch$temp_locale == "Da Nang" ] = "Hai Chau" # central Da Nang
# nm_nomatch$areanameen[ nm_nomatch$temp_locale == "Qui Nhon" ] = "Quy Nhon" # different spelling
# nm_match = rbind(nm_match, nm_nomatch)
# 
# # access centroids
# shpm = shp[ shp$areanameen %in% nm_match$areanameen, c("areaid", "areanameen")]
# coords = cbind(shpm[ , c("areaid", "areanameen")] %>% st_drop_geometry(), st_coordinates(st_centroid(shpm)))
# nm_match = left_join(nm_match, coords, by=c("areanameen"))
# 
# # visualise distribution of points: pan-Vietnam
# # ggplot() + 
# #   geom_sf(data=shp, fill="grey70", col=NA) +
# #   geom_point(data=nm_match, aes(X, Y), col="darkred", pch=21, fill="coral2", size=2.5)
# 
# # extract annual Tmean for ERA5 and WDFE5
# er = raster::brick("./data/climate/era5-land/Vietnam/2m_temperature/monthly/Tmean_1997-01-01_2020-02-01_month_era5land_Vietnam.grd")
# wf = raster::brick("./data/climate/wfde5/temperature_biascorrected/monthly/Tmean_1997-01-01_2018-12-01_month_wfde5_Vietnam.grd")
# 
# # era5
# shpm = shp[ shp$areanameen %in% nm_match$areanameen, c("areaid", "areanameen")]
# era = exactextractr::exact_extract(er, shpm, fun="mean")
# era = cbind(shpm %>% st_drop_geometry(), era)
# era = era %>%
#   reshape2::melt(id.vars = 1:2) %>%
#   dplyr::mutate(date = as.Date(substr(variable, 12, 30), "%Y.%m.%d")) %>%
#   dplyr::mutate(year = lubridate::year(date)) %>%
#   dplyr::group_by(areanameen, year) %>%
#   dplyr::summarise(Tmean = mean(value)) %>%
#   dplyr::mutate(dataset = "ERA5-LAND")
# 
# wfd = exactextractr::exact_extract(wf, shpm, fun="mean")
# wfd = cbind(shpm %>% st_drop_geometry(), wfd)
# wfd = wfd %>%
#   reshape2::melt(id.vars = 1:2) %>%
#   dplyr::mutate(date = as.Date(substr(variable, 12, 30), "%Y.%m.%d")) %>%
#   dplyr::mutate(year = lubridate::year(date)) %>%
#   dplyr::group_by(areanameen, year) %>%
#   dplyr::summarise(Tmean = mean(value)) %>%
#   dplyr::mutate(dataset = "WFDE5")
# 
# 
# # ============ combine datasets together ===============
# 
# temp2 = left_join(temp, nm_match[ , c("temp_locale", "areanameen")], by=c("Locale"="temp_locale")) %>%
#   dplyr::select(-Locale) %>%
#   dplyr::mutate(dataset="Stations")
# 
# dat = do.call(rbind.data.frame, list(temp2, wfd, era))
# 
# # plot for a location
# name = "Nha Trang"
# dd = dat
# dd$dataset = factor(dd$dataset, level=c("Stations", "ERA5-LAND", "WFDE5"), ordered=TRUE)
# p1 = ggplot(dd[ dd$areanameen == name & dd$dataset != "WFDE5", ]) + 
#   geom_line(aes(year, Tmean, group=dataset, col=dataset), alpha=0.7, size=1) +
#   facet_wrap(~areanameen) +
#   theme_minimal() +
#   scale_color_viridis_d(option="viridis", begin=0, end=0.5)
# p2 = ggplot(dd[ dd$areanameen == name & dd$dataset != "ERA5-LAND", ]) + 
#   geom_line(aes(year, Tmean, group=dataset, col=dataset), alpha=0.7, size=1) +
#   facet_wrap(~areanameen) +
#   theme_minimal() +
#   scale_color_viridis_d(option="viridis", begin=0, end=0.7)
# gridExtra::grid.arrange(p1, p2, nrow=2)
# 
# ggplot(dd) + 
#   geom_line(aes(as.numeric(year), Tmean, group=dataset, col=dataset), alpha=0.7, size=1) +
#   facet_wrap(~areanameen) +
#   theme_minimal() 
# 
# 
# # ============ error estimates ===============
# 
# era2 = era %>% dplyr::rename("era5" = Tmean) %>% dplyr::select(-dataset)
# wf2 = wfd %>% dplyr::rename("wfde5" = Tmean) %>% dplyr::select(-dataset)
# tempx = temp2 %>% dplyr::rename("stations" = Tmean) %>% dplyr::select(-dataset) %>% dplyr::mutate(year = as.numeric(year))
# err_dat = left_join(era2, wf2)
# err_dat = left_join(err_dat, tempx)
# err_dat = err_dat[ !is.na(err_dat$stations), ]
# 
# # error
# err_dat$err_wfde5 = err_dat$wfde5  - err_dat$stations
# err_dat$err_era5 = err_dat$era5  - err_dat$stations
# 
# # calculate MAE on two different datasets
# mae_wf = mean(abs(err_dat$err_wfde5), na.rm=TRUE)
# mae_era5 = mean(abs(err_dat$err_era5), na.rm=TRUE)
# 
# # and at different locales
# mae_loc = err_dat %>%
#   dplyr::group_by(areanameen) %>%
#   dplyr::summarise(WFDE5 = mean(abs(err_wfde5), na.rm=TRUE),
#                    ERA5 = mean(abs(err_era5), na.rm=TRUE))
# mae_loc = rbind(data.frame(areanameen = "Vietnam", WFDE5 = mae_wf, ERA5 = mae_era5), mae_loc)
# mae_loc = reshape2::melt(mae_loc, id.vars=1)
# 
# # for temperature (compared to yearly data) error strongly varies between locales
# # in some areas there is a very strong discrepancy; this is likely to do with altitude and local weather conditions
# # strongly suggests that ERA5-Land is more appropriate due to higher spatial res
# fac_order = c("Vietnam", nm_match$areanameen[ order(nm_match$Y, decreasing = TRUE) ])
# mae_loc$areanameen = factor(mae_loc$areanameen, levels=fac_order, ordered=TRUE)
# px = ggplot(mae_loc) + 
#   geom_bar(aes(areanameen, value, group=variable, fill=variable), stat="identity", position=position_dodge(), col="grey50", width=0.5) +
#   theme_classic() + 
#   theme(axis.text.x = element_text(angle=90)) + 
#   xlab("Locale") + 
#   ylab("Mean absolute error") +
#   scale_fill_viridis_d(begin=0.3, end=0.7, name="Dataset")
# 
