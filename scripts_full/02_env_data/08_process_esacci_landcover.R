

# ====================== processing and extraction of ESA-CCI land cover data ==========================

# extracts proportional and population-weighted land cover from ESA-CCI land cover rasters across study time series

setwd("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/")

# dependencies
library(ggplot2)
library(rgdal)
library(raster)
library(dplyr)
library(sf)
library(exactextractr)

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))

# viet shapefile
shp = st_read("./data/shapefiles/d_moss/vietnam_districts_merged.shp")
shp = shp[ ! shp$areanameen %in% c("Truong Sa", "Hoang Sa"), ] # remove offshore districts



# ============================== raster data ===========================

# ESA-CCI rasters (replicate 2018 for 2019)
loc = "./data/landcover/esacci_viet/"
rr = stack(list.files(loc, pattern=".tif", full.names=TRUE))
rr = rr[[ 7:nlayers(rr) ]]
rr = stack(rr, rr[[ nlayers(rr) ]], rr[[ nlayers(rr) ]])
names(rr) = paste("esacci", 1998:2020, sep="_")

# world pop
wp = raster::stack(list.files("./data/population/viet_population_worldpop/", full.names=TRUE))
wp = stack(wp[[1]], wp[[1]], wp, wp[[nlayers(wp)]])
names(wp) = paste("population_", 1998:2020, sep="")

# area raster
areax = raster::area(rr[[1]]); names(areax) = "area"



# ======================= for each year create stack and extract metrics of interest ===========================

# for storing results
result = data.frame()

# run for each year
for(year in 1998:2020){
  
  print(year)
  
  # rasters
  esax = rr[[ grep(year, names(rr)) ]]
  wpx = wp[[ grep(year, names(wp)) ]]
  
  # aggregate and resample pop to same res as esa
  wpx = aggregate(wpx, fact=3, fun=sum)
  wpx = resample(wpx, esax, fun="ngb")
  
  # create stack
  stackx = raster::stack(esax, wpx, areax)
  names(stackx) = c("land", "population", "area")
  
  # extract
  ext = exactextractr::exact_extract(stackx, shp)
  
  # reclassify land into areas of interest
  getLCMetrics = function(x){
    
    poly = ext[[x]]
    total_area = sum(poly$area * poly$coverage_fraction)
    total_population = sum(poly$population * poly$coverage_fraction, na.rm=TRUE)
    
    # remove NAs and classify LC classes of interest
    poly = poly[ !is.na(poly$land), ]
    poly$forest = ifelse(poly$land %in% c(50:100, 160, 170), 1, 0)
    poly$cropland = ifelse(poly$land %in% c(10:40), 1, 0)
    poly$cropland_irr = ifelse(poly$land == 20, 1, 0)
    poly$urban = ifelse(poly$land == 190, 1, 0)
    
    # calculate
    resx = data.frame(id = x,
                      year = year,
                      
                      # proportion cover
                      forest_prop = sum(poly$area[ poly$forest == 1 ] * poly$coverage_fraction[ poly$forest == 1 ]) / total_area,
                      cropland_prop = sum(poly$area[ poly$cropland == 1 ] * poly$coverage_fraction[ poly$cropland == 1 ]) / total_area,
                      croplandirrigated_prop = sum(poly$area[ poly$cropland_irr == 1 ] * poly$coverage_fraction[ poly$cropland_irr == 1 ]) / total_area,
                      urban_prop = sum(poly$area[ poly$urban == 1 ] * poly$coverage_fraction[ poly$urban == 1 ]) / total_area,
                      
                      # area
                      urbanarea_esa = sum(poly$area[ poly$urban == 1 ] * poly$coverage_fraction[ poly$urban == 1 ]),
                      
                      # population weighted cover
                      forest_pw = sum(poly$population[ poly$forest == 1 ] * poly$coverage_fraction[ poly$forest == 1], na.rm=TRUE) / total_population,
                      cropland_pw = sum(poly$population[ poly$cropland == 1 ] * poly$coverage_fraction[ poly$cropland == 1], na.rm=TRUE) / total_population,
                      croplandirrigated_pw = sum(poly$population[ poly$cropland_irr == 1 ] * poly$coverage_fraction[ poly$cropland_irr == 1], na.rm=TRUE) / total_population,
                      urban_pw = sum(poly$population[ poly$urban == 1 ] * poly$coverage_fraction[ poly$urban == 1], na.rm=TRUE) / total_population
                      )
    return(resx)
  }
  
  # run for all districts and append
  landcover = do.call(rbind.data.frame, lapply(1:length(ext), getLCMetrics))
  landcover = left_join(data.frame(id = 1:nrow(shp), areaid = shp$areaid, areanameen = shp$areanameen), landcover)
  result = rbind(result, landcover)
}

# save
write.csv(result, "./code/viet_dengue_districts/output/covariates/LandUseCover_ESA_MergedSHP2.csv", row.names=FALSE)




