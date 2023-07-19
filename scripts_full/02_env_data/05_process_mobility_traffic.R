

# ==================== Processing of travel data from Vietnam Office of Statistics ===================

# ------------- dependencies -----------

setwd("C:/Users/roryj/Documents/PhD/202007_LSHTM_Dengue/analysis/")

library(ggplot2)
library(rgdal)
library(raster)
library(dplyr)
library(sf)
library(animation)

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))



# ================= 1. Millions of person kilometres by road ======================

# persons traffic by road (millions of person kilometres)
# measure of changes in human movement at province-level over time (does not report on within/between province)
stats1 = read.csv("./data/viet_statistics/transport/numpassengers_traffic_byprovince_millionpersonkm.csv", stringsAsFactors = FALSE)
stats1 = stats1[ 1:71, ]

# initial processing:
sx = reshape2::melt(stats1, id.vars=1:7)
names(sx)[1] = "Metric"
names(sx)[8:9] = c("Year", "value")
sx$Year = as.numeric(substr(sx$Year, 2, 5))
sx = sx[ order(sx$GeoUnit, sx$Year), ]
sx$value = as.numeric(sx$value)

# adjustment of statistics for specific years due to reporting
# combine provinces that were later combined
# conduct simple interpolation and adjustment of values for merged provinces in early part of time series
# extend time series in either direction (to 1998 and 2019) using linear interpolation (assume change rate is same as prev 5 years)

# (1) Ha Noi = Ha Noi + Ha Tay for 2000 to 2007 (states merge in 08)
sx$value[ sx$GeoUnit == "Ha Noi" & sx$Year %in% 2000:2007 ] = sx$value[ sx$GeoUnit == "Ha Noi" & sx$Year %in% 2000:2007 ] + sx$value[ sx$GeoUnit == "Ha Tay" & sx$Year %in% 2000:2007 ]
sx = sx[ sx$GeoUnit != "Ha Tay", ]

# (2) Dak Lak statistics contain Dak Nong for years 2000 to 2002
# correction: assume annual rate of change in Dak Nong is same from 2000 to 2002 as in 2003:2004, then subtract from Dak Lak
# reasonable since traffic starts to significantly increases country wide around 2005
dnong1 = sx$value[ sx$GeoUnit == "Dak Nong" ]
rates = dnong1[ 5:length(dnong1)] / dnong1[ 4:(length(dnong1)-1)]
rate_change = mean(rates[1:3])
imputed_vals = rep(NA, 3)
imputed_vals[3] = dnong1[4] / rate_change
imputed_vals[2] = imputed_vals[3]  / rate_change
imputed_vals[1] = imputed_vals[2]  / rate_change

# impute to Dak Nong and subtract from Dak Lak
sx$value[ sx$GeoUnit == "Dak Nong" ][1:3] = imputed_vals
sx$value[ sx$GeoUnit == "Dak Lak" ][1:3] = sx$value[ sx$GeoUnit == "Dak Lak" ][1:3] - imputed_vals

# (3) Dien Bien statistics contain Lai Chau for years 2000 to 2002
# same correction as above
laichau1 = sx$value[ sx$GeoUnit == "Lai Chau" ]
rates = laichau1[ 5:length(laichau1)] / laichau1[ 4:(length(laichau1)-1)]
rate_change = mean(rates[1:3])
imputed_vals = rep(NA, 3)
imputed_vals[3] = laichau1[4] / rate_change
imputed_vals[2] = imputed_vals[3]  / rate_change
imputed_vals[1] = imputed_vals[2]  / rate_change

# impute to Lai Chau and subtract from Dien Bien
sx$value[ sx$GeoUnit == "Lai Chau" ][1:3] = imputed_vals
sx$value[ sx$GeoUnit == "Dien Bien" ][1:3] = sx$value[ sx$GeoUnit == "Dien Bien" ][1:3] - imputed_vals

# (4) Can Tho statistics contain Hau Giang for years 2000 to 2002
# same correction as above
hg1 = sx$value[ sx$GeoUnit == "Hau Giang" ]
rates = hg1[ 5:length(hg1)] / hg1[ 4:(length(hg1)-1)]
rate_change = mean(rates[1:3])
imputed_vals = rep(NA, 3)
imputed_vals[3] = hg1[4] / rate_change
imputed_vals[2] = imputed_vals[3]  / rate_change
imputed_vals[1] = imputed_vals[2]  / rate_change

# impute to Lai Chau and subtract from Dien Bien
sx$value[ sx$GeoUnit == "Hau Giang" ][1:3] = imputed_vals
sx$value[ sx$GeoUnit == "Can Tho" ][1:3] = sx$value[ sx$GeoUnit == "Can Tho" ][1:3] - imputed_vals

# conduct imputation of start/end for each time series
imputeTimeSeriesEnds = function(geounit){
  
  dfx = sx[ sx$GeoUnit == geounit, ]
  dfx = do.call(rbind.data.frame, list(dfx[ 1, ], dfx[ 1, ], dfx))
  #dfx = do.call(rbind.data.frame, list(dfx, dfx[ nrow(dfx), ], dfx[ nrow(dfx), ], dfx[ nrow(dfx), ]))
  dfx$Year = 1998:2020
  
  # change rates start and end
  rate_start = mean(dfx$value[ dfx$Year %in% 2001:2006] / dfx$value[ dfx$Year %in% 2000:2005])
  #rate_end = mean(dfx$value[ dfx$Year %in% 2013:2017] / dfx$value[ dfx$Year %in% 2012:2016])
  
  # calculate imputed
  dfx$value[ dfx$Year == 1999 ] = dfx$value[ dfx$Year == 2000 ] / rate_start
  dfx$value[ dfx$Year == 1998 ] = dfx$value[ dfx$Year == 1999 ] / rate_start

  #ggplot(dfx, aes(Year, value)) + geom_line() + theme_classic() + facet_wrap(~GeoUnit, scales="free_y")
  
  # return
  return(dfx)
}

# run imputation
sx_imp = do.call(rbind.data.frame, lapply(unique(sx$GeoUnit), imputeTimeSeriesEnds))

# combine with lookup table for shapefile 
shp_lookup = read.csv("./output/data_processed/viet_statistics/lookup_tables/prov_shp_names.csv", stringsAsFactors = FALSE)
names(sx_imp)[2] = "Name_VSO"
sx_imp = left_join(sx_imp, shp_lookup)

# save
#write.csv(sx_imp, "./output/data_processed/viet_statistics/data_tables/NumPassengers_Traffic_byProvince_19982019.csv", row.names=FALSE)

# processed traffic data; millions of person kilometres
trav1 = sx_imp %>%
  dplyr::select(provinceid, provincena, Year, value) %>%
  dplyr::rename("province"=2, "year"=3, "traffic_milperskm"=4) %>%
  dplyr::filter(!is.na(provinceid))



# 
# # =========================== 2. Number of passengers travelling by road =======================
# 
# # persons traffic by road (millions of person kilometres)
# # measure of changes in human movement at province-level over time (does not report on within/between province)
# stats1 = read.csv("./data/viet_statistics/transport/numpassengers_road_byprovince.csv", stringsAsFactors = FALSE) %>%
#   dplyr::rename("Metric"=1, "GeoUnit"=2) %>%
#   dplyr::select(-Notes_Vstats, -Notes_RG)
# 
# # initial processing:
# sx = reshape2::melt(stats1, id.vars=1:2)
# names(sx)[3:4] = c("Year", "value")
# sx$Year = as.numeric(substr(sx$Year, 2, 5))
# sx = sx[ order(sx$GeoUnit, sx$Year), ]
# sx$value = as.numeric(sx$value)
# 
# # adjustment of statistics for specific years due to reporting
# # combine provinces that were later combined
# # conduct simple interpolation and adjustment of values for merged provinces in early part of time series
# # extend time series in either direction (to 1998 and 2019) using linear interpolation (assume change rate is same as prev 5 years)
# 
# # (1) Ha Noi = Ha Noi + Ha Tay for 2000 to 2007 (states merge in 08)
# sx$value[ sx$GeoUnit == "Ha Noi" & sx$Year %in% 2000:2007 ] = sx$value[ sx$GeoUnit == "Ha Noi" & sx$Year %in% 2000:2007 ] + sx$value[ sx$GeoUnit == "Ha Tay" & sx$Year %in% 2000:2007 ]
# sx = sx[ sx$GeoUnit != "Ha Tay", ]
# 
# # (2) Dak Lak statistics contain Dak Nong for years 2000 to 2002
# # correction: assume annual rate of change in Dak Nong is same from 2000 to 2002 as in 2003:2004, then subtract from Dak Lak
# # reasonable since traffic starts to significantly increases country wide around 2005
# dnong1 = sx$value[ sx$GeoUnit == "Dak Nong" ]
# rates = dnong1[ 5:length(dnong1)] / dnong1[ 4:(length(dnong1)-1)]
# rate_change = mean(rates[1:3])
# imputed_vals = rep(NA, 3)
# imputed_vals[3] = dnong1[4] / rate_change
# imputed_vals[2] = imputed_vals[3]  / rate_change
# imputed_vals[1] = imputed_vals[2]  / rate_change
# 
# # impute to Dak Nong and subtract from Dak Lak
# sx$value[ sx$GeoUnit == "Dak Nong" ][1:3] = imputed_vals
# sx$value[ sx$GeoUnit == "Dak Lak" ][1:3] = sx$value[ sx$GeoUnit == "Dak Lak" ][1:3] - imputed_vals
# 
# # (3) Dien Bien statistics contain Lai Chau for years 2000 to 2002
# # same correction as above
# laichau1 = sx$value[ sx$GeoUnit == "Lai Chau" ]
# rates = laichau1[ 5:length(laichau1)] / laichau1[ 4:(length(laichau1)-1)]
# rate_change = mean(rates[1:3])
# imputed_vals = rep(NA, 3)
# imputed_vals[3] = laichau1[4] / rate_change
# imputed_vals[2] = imputed_vals[3]  / rate_change
# imputed_vals[1] = imputed_vals[2]  / rate_change
# 
# # impute to Lai Chau and subtract from Dien Bien
# sx$value[ sx$GeoUnit == "Lai Chau" ][1:3] = imputed_vals
# sx$value[ sx$GeoUnit == "Dien Bien" ][1:3] = sx$value[ sx$GeoUnit == "Dien Bien" ][1:3] - imputed_vals
# 
# # (4) Can Tho statistics contain Hau Giang for years 2000 to 2002
# # same correction as above
# hg1 = sx$value[ sx$GeoUnit == "Hau Giang" ]
# rates = hg1[ 5:length(hg1)] / hg1[ 4:(length(hg1)-1)]
# rate_change = mean(rates[1:3])
# imputed_vals = rep(NA, 3)
# imputed_vals[3] = hg1[4] / rate_change
# imputed_vals[2] = imputed_vals[3]  / rate_change
# imputed_vals[1] = imputed_vals[2]  / rate_change
# 
# # impute to Hau Giang and subtract from Can Tho
# sx$value[ sx$GeoUnit == "Hau Giang" ][1:3] = imputed_vals
# sx$value[ sx$GeoUnit == "Can Tho" ][1:3] = sx$value[ sx$GeoUnit == "Can Tho" ][1:3] - imputed_vals
# 
# # conduct imputation of start/end for each time series
# imputeTimeSeriesEnds = function(geounit){
#   
#   dfx = sx[ sx$GeoUnit == geounit, ]
#   dfx = do.call(rbind.data.frame, list(dfx[ 1, ], dfx[ 1, ], dfx))
#   dfx = do.call(rbind.data.frame, list(dfx, dfx[ nrow(dfx), ], dfx[ nrow(dfx), ], dfx[ nrow(dfx), ]))
#   dfx$Year = 1998:2020
#   
#   # change rates start and end
#   rate_start = mean(dfx$value[ dfx$Year %in% 2001:2006] / dfx$value[ dfx$Year %in% 2000:2005])
#   rate_end = mean(dfx$value[ dfx$Year %in% 2013:2017] / dfx$value[ dfx$Year %in% 2012:2016])
#   
#   # calculate imputed
#   dfx$value[ dfx$Year == 1999 ] = dfx$value[ dfx$Year == 2000 ] / rate_start
#   dfx$value[ dfx$Year == 1998 ] = dfx$value[ dfx$Year == 1999 ] / rate_start
#   dfx$value[ dfx$Year == 2018 ] = dfx$value[ dfx$Year == 2017 ] * rate_end
#   dfx$value[ dfx$Year == 2019 ] = dfx$value[ dfx$Year == 2018 ] * rate_end
#   dfx$value[ dfx$Year == 2020 ] = dfx$value[ dfx$Year == 2019 ] * rate_end
#   
#   #ggplot(dfx, aes(Year, value)) + geom_line() + theme_classic() + facet_wrap(~GeoUnit, scales="free_y")
#   
#   # return
#   return(dfx)
# }
# 
# # run imputation
# sx_imp = do.call(rbind.data.frame, lapply(unique(sx$GeoUnit), imputeTimeSeriesEnds))
# 
# # combine with lookup table for shapefile 
# shp_lookup = read.csv("./output/data_processed/viet_statistics/lookup_tables/prov_shp_names.csv", stringsAsFactors = FALSE)
# names(sx_imp)[2] = "Name_VSO"
# sx_imp = left_join(sx_imp, shp_lookup)
# 
# # save
# #write.csv(sx_imp, "./output/data_processed/viet_statistics/data_tables/NumPassengers_Traffic_byProvince_19982019.csv", row.names=FALSE)
# 
# # processed traffic data; millions of person kilometres
# trav2 = sx_imp %>%
#   dplyr::select(provinceid, provincena, Year, value) %>%
#   dplyr::rename("province"=2, "year"=3, "traffic_milpassengers"=4) %>%
#   dplyr::filter(!is.na(provinceid))



# ========================== combine both data sources ========================

#trav = left_join(trav1, trav2) #%>%
  #dplyr::mutate(traffic_jlengthkm = traffic_milperskm / traffic_milpassengers)

trav = trav1

# populations
pop = read.csv("./code/viet_dengue_districts/output/covariates/Population_Census_2009_2019.csv") %>%
  dplyr::group_by(province, year) %>%
  dplyr::summarise(population = sum(total))

trav = left_join(trav, pop) %>%
  dplyr::mutate(
    traffic_thouskmperinhab = (traffic_milperskm / (population / 10^6)) / 1000
  )

# save covariates
write.csv(trav, "./code/viet_dengue_districts/output/covariates/Mobility_TrafficMetrics_VietGSO.csv", row.names=FALSE)

ggplot(trav) +
  geom_line(aes(year, traffic_thouskmperinhab)) +
  facet_wrap(~province, nrow=9)
ggplot(trav) +
  geom_line(aes(year, traffic_milperskm)) +
  facet_wrap(~province, nrow=9)



# # ============== subset to only provinces and map ===============
# 
# shp = st_read("./code/viet_dengue_districts/data/shapefiles/provinces.shp")
# shp %>%
#   full_join(trav) %>%
#   dplyr::filter(year %in% c(1998, 2008, 2018)) %>%
#   ggplot() + 
#   geom_sf(aes(fill=log(traffic_thouskmperinhab+1)), col=NA) + 
#   facet_wrap(~year) +
#   scale_fill_gradientn(colors=viridisLite::turbo(200))
# 




# # -------------- initial visualising -----------
# 
# # stats
# sx_imp = read.csv("./output/data_processed/viet_statistics/data_tables/NumPassengers_Traffic_byProvince_19982019.csv", stringsAsFactors = FALSE)
# sx_imp = sx_imp[ !is.na(sx_imp$provinceid), ]
# 
# # plot for focal provinces
# prov_plot = ggplot(sx_imp, aes(Year, value)) +
#   geom_line(col="darkred") + 
#   theme_minimal() + 
#   facet_wrap(~provincena) + 
#   theme(panel.grid = element_blank(), 
#         strip.text = element_text(size=13),
#         plot.title = element_text(size=15, hjust=0.5)) + 
#   ylab("Millions passenger km (road)")
# 
# # comparison with population data; does travel per capita increase/decrease/same?
# # travel per capita increases substantially, with similar rapidity in most provinces except Dak Lak which is slower
# pop = read.csv("./output/data_processed/population/vietpop_perdistrict_HRWimputed_19982019.csv", stringsAsFactors = FALSE) %>%
#   group_by(Province.ID, Province.name, Year) %>%
#   dplyr::summarise(TotalPop = sum(tsvalue, na.rm=TRUE)) %>%
#   rename("provincena" = Province.name, "provinceid" = Province.ID)
# pt = left_join(pop, sx_imp)
# pt$km_per_person = (pt$value * 10^6) / pt$TotalPop # convert millions to passenge km and caculate km per person
# ggplot(pt) + 
#   geom_point(aes(TotalPop/10^6, value, col=Year, size=Year)) + 
#   facet_wrap(~provincena, scales="free") + 
#   theme_minimal() + 
#   scale_colour_viridis_c(begin=0, end=0.8)
# ggplot(pt) + 
#   geom_point(aes(Year, km_per_person)) + 
#   facet_wrap(~provincena) + 
#   theme_minimal() 
# 
# # read provinces shapefile and combine with 
# shp1 = st_read("./data/shapefiles/d_moss/provinces.shp")
# shp1 = st_crop(shp1, ext)
# shp1 = full_join(shp1, sx_imp[ , c("provinceid", "Year", "value")])
# 
# # viz over time (natural scale)
# colScale = colorRampPalette(RColorBrewer::brewer.pal(9, "PuBuGn"))(200)
# plotx = ggplot() +
#   geom_sf(data=shp1, aes(fill=value), colour="grey30", size=0.1) +
#   facet_wrap(~Year, nrow=2) + 
#   theme_minimal() + 
#   ggtitle("Annual road travel (millions passenger km)") +
#   #scale_fill_gradientn(colors=colScale, na.value = "white") +
#   scale_fill_viridis_c(option="magma", direction = -1) +
#   theme(axis.text = element_blank(),
#         panel.grid = element_blank(), 
#         strip.text = element_text(size=13),
#         plot.title = element_text(size=15, hjust=0.5))
# ggsave(plotx, file="./output/figures/prelim/RoadJourneys_Viet_natural.png", device="png", dpi=300, units="in", width=12, height=8, scale=0.9)

# # video
# index = unique(shp1$Year)
# vidname = "Vietnam_RoadTravel_timeseries"  
# 
# saveVideo({
#   for(i in index){
#     
#     # plot map
#     map = ggplot() +
#       geom_sf(data=shp1[ shp1$Year == i, ], aes(fill=value), colour="grey30", size=0.1) +
#       facet_wrap(~Year) +
#       theme_minimal() + 
#       ggtitle("Annual road travel (millions passenger km)") +
#       scale_fill_viridis_c(option="magma", direction = -1, limits=range(shp1$value)) +
#       theme(axis.text = element_blank(),
#             panel.grid = element_blank(), 
#             strip.text = element_text(size=13),
#             strip.text.y.left = element_text(angle=0),
#             legend.title = element_blank(),
#             plot.title = element_text(size=15, hjust=0.5))
#     
#     print(map)
#   }
# }, ani.width = 500, ani.height= 600, video.name = paste(vidname, ".mp4", sep=""), interval=1, other.opts = "-b 500k")  # higher bitrate, better quality



# 
# # ================ distribute across districts within provinces ================= 
# 
# # possible scenarios
# # 1: per capita (i.e. kilometres distributed equally among all people within province per year)
# # 2: urban/rural (urban residents take more and shorter trips; rural residents take fewer and longer, following Hanoi study) https://www.researchgate.net/publication/284323148_Tracking_Sustainable_Transport_in_Vietnam_Data_and_Policy_Review_for_Energy_Efficiency_and_Climate_Change_2015?enrichId=rgreq-ecc29d8703c7c7bba8462253682a0259-XXX&enrichSource=Y292ZXJQYWdlOzI4NDMyMzE0ODtBUzoyOTg1MDU0NzEwNTM4MjhAMTQ0ODE4MDY1ODQ4Mg%3D%3D&el=1_x_2&_esc=publicationCoverPdf
# # 3: road density (distributed to districts depending on roads density - i.e. all roads are travelled in equal proportion to their density)
# # 4: connectivity (distributed to districts depending on their travel time to all other districts (i.e. most connected districts have more travel))
# 
# pop = read.csv("./output/data_processed/population/vietpop_perdistrict_HRWimputed_19982019.csv", stringsAsFactors = FALSE)
# traffic = read.csv("./output/data_processed/viet_statistics/data_tables/NumPassengers_Traffic_byProvince_19982019.csv", stringsAsFactors = FALSE)
# traffic = traffic[ !is.na(traffic$provinceid), ]
# 
# # scenario 1: per capita
# # distribute kilometres equally among all people per year (i.e. everyone has same per-capita movement annually)
# # however, means that travel is perfectly collinear with population at the province level; may not be appropriate
# pop$yp = paste(pop$Province.name, pop$Year, sep="_")
# scen1 = data.frame()
# for(i in unique(pop$yp)){
#   
#   yx = pop[ pop$yp == i, ]
#   tx = traffic[ traffic$Year == yx$Year[1] & traffic$provinceid == yx$Province.ID[1], ]
#   yx$weight = yx$tsvalue / sum(yx$tsvalue)
#   yx$traffic_millkm_ds1 = (tx$value) * yx$weight
#   yx = yx[ , c("areaid", "tsvalue", "Province.name", "Year", "traffic_millkm_ds1")] %>%
#     rename("population" = tsvalue, "areaprovin"=Province.name)
#   scen1 = rbind(scen1, yx)
# }
# 
# ggplot(scen1[ scen1$Year %in% 2012:2018,]) + 
#   geom_point(aes(population, traffic_millkm_ds1)) + 
#   facet_wrap(Year~areaprovin)
# 
# # ggplot(scen1[ scen1$areaprovin == "Ha Noi", ]) + 
# #   geom_point(aes(population, traffic_millkm_ds1, col=Year)) + 
# #   facet_wrap(~areaid)
# 
# write.csv(scen1, "./output/data_processed/connectivity/VietnamGSO_Traffic_perdistricts.csv", row.names=FALSE)
