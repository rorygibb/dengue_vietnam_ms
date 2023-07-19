

# ============== Population at Vietnam district-level from GPW ===================

# UN-adjusted GPWv4 in 2000, 2005, 2010, 2015, 2020
# Linear interpolate between values for each district

# project root and dependencies
setwd("C:/Users/roryj/Documents/PhD/202007_LSHTM_Dengue/analysis/")
pacman::p_load("ggplot2", "rgdal", "raster", "dplyr", "sf", "exactextractr", "ggspatial")

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))

# districts shape
shp = st_read("./code/dengue_vietnam/data/shapefiles/dengue_districts_shapefile.shp")




# ============= population raster data: extract for all years/districts ================

# gridded pw 2000-2020 in 5 year increments
gpw = raster::stack(list.files("./data/population/gpw_v4/gpw_v4_rasters/", pattern=".tif", full.names=TRUE))
names(gpw) = paste("pop", seq(2000, 2020, by=5), sep="_")

# plot example
# pp = crop(gpw[[4]], shp)
# plot(log(pp+1), col=viridis::inferno(200))

# extract for all districts
px = exactextractr::exact_extract(gpw, shp, fun="sum")
px = cbind(shp$areaid, px) 
names(px) = c("areaid", paste("pop", seq(2000, 2020, by=5), sep="_"))

# linear interpolate between years for all districts
popInterp = function(x){
  
  # linearly interpolate between years and back to 1998
  popx = px[ x, ]
  a = popx$pop_2000 + ((popx$pop_2005 - popx$pop_2000)/5) * -2:5
  b = popx$pop_2005 + ((popx$pop_2010 - popx$pop_2005)/5) * 1:5
  c = popx$pop_2010 + ((popx$pop_2015 - popx$pop_2010)/5) * 1:5
  d = popx$pop_2015 + ((popx$pop_2020 - popx$pop_2015)/5) * 1:5
  resx = data.frame(areaid = popx$areaid, 
                    year = 1998:2020,
                    population_gpw = c(a, b, c, d))
  
  # test
  if(all(popx$pop_2000 == resx$population_gpw[ resx$year == 2000 ], 
         popx$pop_2010 == resx$population_gpw[ resx$year == 2010 ],
         popx$pop_2020 == resx$population_gpw[ resx$year == 2020 ])){
    print("Interpolation correct")
  } else{
    print("Interpolation error")
  }
  
  # return for district
  return(resx)
}

# result
result = do.call(rbind.data.frame, lapply(1:nrow(px), popInterp))
result = left_join(result, shp[ , c("areaid", "areanameen", "areaprovin")] %>% st_drop_geometry())

# # viz
# ggplot(result[ result$areaprovin == "Ha Noi",]) +
#   geom_line(aes(year, population_gpw)) +
#   facet_wrap(~areanameen)

# save
write.csv(result, "./output/data_processed/population/VietPopulation_District_MergedSHP_GPW4_interpolated.csv", row.names=FALSE)
