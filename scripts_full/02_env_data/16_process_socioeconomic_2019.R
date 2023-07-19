


# =================== Time-series projection of socioeconomic indicators =====================

# Use district-level indicators from 2009 and 2019 census for major socioeconomic axes
# Project backwards assuming consistent rates

setwd("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/")

library(gdistance)
library(abind)
library(rje)
library(ggplot2)
library(malariaAtlas)
library(ggplot2)
library(rgdal)
library(raster)
library(dplyr)
library(sf)
library(exactextractr)
library(ggfortify)

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))

# theme for mapping
maptheme = theme_classic() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.ticks = element_blank(),
        plot.title = element_text(hjust=0.5, size=12),
        legend.title = element_text(size=10))

# merged shapefile
shpm = st_read("./code/viet_dengue_districts/data/shapefiles/vietnam_districts_merged.shp") %>%
  dplyr::filter(!areanameen %in% c("Truong Sa", "Hoang Sa"))




# ======================== district-level 2009 census socioeconomic data (already merged) ==================

wb09 = read.csv("./code/viet_dengue_districts/output/covariates/Socioeconomic_WorldBank_MAPVietnam_2009.csv") %>%
  dplyr::left_join(st_read("./code/viet_dengue_districts/data/shapefiles/vietnam_districts_merged.shp") %>%
                     st_drop_geometry() %>%
                     dplyr::select(areaid, areaprovin))  %>%
  dplyr::rename("province"=areaprovin) %>%
  dplyr::filter(areaid != 70273) # remove Bach Long Vi as has no data

# sanitation any
wb09$sanitation_flushtoilet_any = wb09$sanitation_flushtoilet_indoor + wb09$sanitation_flushtoilet_outdoor

# select vars
wb09 = wb09 %>%
  dplyr::select(areaid, province, population, water_piped_well, sanitation_flushtoilet_indoor, sanitation_flushtoilet_outdoor, sanitation_flushtoilet_any)




# ===================== district-level 2019 census data =====================

# 2019 census data
wb19 = read.csv("./code/viet_dengue_districts/data/census_2019/census2019_sociovariables_dec2021.csv") %>%
  dplyr::filter(!district_dmoss %in% c("Truong Sa", "Hoang Sa")) %>%
  dplyr::mutate(
    Water_Tap = replace(Water_Tap, is.na(Water_Tap), 0),
    Water_DrilledWell = replace(Water_DrilledWell, is.na(Water_DrilledWell), 0),
    Sanitation_FlushToilet_Indoor = replace(Toilet_IndoorHygeinic, is.na(Toilet_IndoorHygeinic), 0),
    Sanitation_FlushToilet_Outdoor = replace(Toilet_OutdoorHygeinic, is.na(Toilet_OutdoorHygeinic), 0),
    NetMigrationRate_2019 = replace(NetMigrationRate_ImmigMinusEmig, is.na(NetMigrationRate_ImmigMinusEmig), 0),
    Water_Piped_Well = Water_Tap + Water_DrilledWell,
    Sanitation_FlushToilet_Any = Sanitation_FlushToilet_Indoor + Sanitation_FlushToilet_Outdoor
  ) %>%
  dplyr::select(areaid, Water_Piped_Well, Sanitation_FlushToilet_Indoor, Sanitation_FlushToilet_Outdoor, Sanitation_FlushToilet_Any, NetMigrationRate_2019) %>%
  dplyr::rename_all(tolower) %>%
  dplyr::mutate(year = 2019)

# 2019 population data
pop = read.csv("./code/viet_dengue_districts/data/census_2019/census2019_population.csv") %>%
  dplyr::filter(areaid %in% wb19$areaid) %>%
  dplyr::select(areaid, Total) %>%
  dplyr::mutate(
    Total = as.numeric(stringr::str_replace_all(Total, stringr::fixed(" "), ""))
  ) %>%
  dplyr::rename("population" = Total) 
wb19 = left_join(wb19, pop)

# merge to combined districts
tomerge = wb19 %>% dplyr::filter(!areaid %in% shpm$areaid)

# new areaids
merge_lk = shpm[ !shpm$areaid %in% wb19$areaid, ] %>%
  st_drop_geometry() %>%
  dplyr::mutate(idx = areaid) %>%
  tidyr::separate_rows(areaid, sep="_") %>%
  dplyr::mutate(areaid = as.numeric(areaid))

# combine
tomerge = left_join(tomerge, merge_lk[ , c("areaid", "idx")], by=c("areaid"))

# calculate pop weighted % for %
merged = tomerge %>%
  dplyr::group_by(idx) %>%
  dplyr::summarise(water_piped_well  = sum(water_piped_well *population)/sum(population),
                   sanitation_flushtoilet_indoor  = sum(sanitation_flushtoilet_indoor *population)/sum(population), 
                   sanitation_flushtoilet_outdoor  = sum(sanitation_flushtoilet_outdoor *population)/sum(population), 
                   sanitation_flushtoilet_any  = sum(sanitation_flushtoilet_any *population)/sum(population), 
                   netmigrationrate_2019  = mean(as.numeric(netmigrationrate_2019), na.rm=TRUE),
                   population = sum(population)) %>%
  dplyr::rename("areaid" = idx) %>%
  dplyr::mutate(year = 2019)

# subset to required vars and remove Bach Long Vi as missing from 2009 data and offshore
wb19 = wb19 %>% 
  filter(areaid %in% shpm$areaid) %>%
  rbind(merged) %>%
  dplyr::filter(areaid != 70273) %>%
  dplyr::mutate(
    water_piped_well = water_piped_well/100,
    sanitation_flushtoilet_indoor = sanitation_flushtoilet_indoor/100,
    sanitation_flushtoilet_outdoor = sanitation_flushtoilet_outdoor/100,
    sanitation_flushtoilet_any = sanitation_flushtoilet_any/100
  )

# visualse changes
# wb19 %>%
#   dplyr::left_join(wb09 %>% dplyr::select(areaid, province)) %>%
#   dplyr::select(-netmigrationrate_2019) %>%
#   rbind(
#     wb09 %>% dplyr::mutate(year = 2009, water_piped_well = water_piped_well * 100,
#                            sanitation_flushtoilet_any= sanitation_flushtoilet_any * 100) 
#   ) %>%
#   ggplot() + 
#   geom_point(aes(factor(year), water_piped_well, group=areaid)) + 
#   geom_line(aes(factor(year), water_piped_well, group=areaid)) + 
#   facet_wrap(~province) + 
#   theme_bw()
#   
# # visualse changes
# wb19 %>%
#   dplyr::left_join(wb09 %>% dplyr::select(areaid, province)) %>%
#   dplyr::select(-netmigrationrate_2019) %>%
#   rbind(
#     wb09 %>% dplyr::mutate(year = 2009, water_piped_well = water_piped_well * 100,
#                            sanitation_flushtoilet_any= sanitation_flushtoilet_any * 100) 
#   ) %>%
#   ggplot() + 
#   geom_point(aes(factor(year), sanitation_flushtoilet_any, group=areaid)) + 
#   geom_line(aes(factor(year), sanitation_flushtoilet_any, group=areaid)) + 
#   facet_wrap(~province) + 
#   theme_bw()





# ==================== project for interim and future years ================

projectIndicators = function(x){
  
  # linear between 
  w09x = wb09[ wb09$areaid == x, ] %>% dplyr::mutate(year = 2009) %>% dplyr::select(-province, -population, -areaid)
  w19x = wb19[ wb19$areaid == x, ] %>% dplyr::mutate(year = 2019) %>% dplyr::select(-population, -netmigrationrate_2019, -areaid)
  wx = rbind(w09x, w19x)
  rates = (wx[ 2, 1:4] - wx[ 1, 1:4]) / 10 
  result = wx[ 1, ]
  for(y in 2010:2018){
    yy = result[ nrow(result), ] %>%
      dplyr::mutate(year = y)
    yy[ , 1:4 ] = yy[ , 1:4] + rates
    result = rbind(result, yy)
  }
  result = rbind(result, wx[ 2, ])
  row.names(result) = c()
  
  # change rates (mean over 5 years) and project backward
  calcChangeRates = function(x){  mean( x[ 1:(length(x)-1)] / x[ 2:length(x)] ) }
  chr = apply(result[ 1:6, 1:4], 2, calcChangeRates)
  result2 = wx[ 1, ]
  for(y in 2008:1998){
    yy = result2[ 1, ] %>%
      dplyr::mutate(year = y)
    yy[ , 1:4 ] = yy[ , 1:4] * chr
    result2 = rbind(yy, result2)
  }
  
  # combine
  result = rbind(
    result2 %>% dplyr::filter(year != 2009),
    result
  )
  row.names(result) = c()
 
  # set any above 1 to 1, and any below 0 to zero
  constrain = function(x){
    x[ x>1 ] = 1
    x[ x<0 ] = 0
    return(x)
  }
  result[ , 1:4] = apply(result[ , 1:4], 2, constrain)
  result$areaid = x
  
  # set 2020 to 2019
  result = rbind(result, result %>% dplyr::filter(year == 2019) %>% dplyr::mutate(year = 2020))
  
  # test viz
  # result %>%
  #   reshape2::melt(id.vars = 5:6) %>%
  #   ggplot() + 
  #   geom_line(aes(year, value)) + 
  #   facet_wrap(~variable)
  
  return(result)
}

# project
projected = do.call(
  rbind.data.frame, 
  lapply(
    wb19$areaid, projectIndicators
  )
)

projected$sanitation_flushtoilet_outdoor = projected$sanitation_flushtoilet_any - projected$sanitation_flushtoilet_indoor
projected$water_piped_well[ is.na(projected$water_piped_well) ] = 0


# =============== save output ===============

write.csv(projected, "./code/viet_dengue_districts/output/covariates/Socioeconomic_TimeSeriesProj.csv", row.names=FALSE)


