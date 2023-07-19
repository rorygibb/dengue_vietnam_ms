

# =========================== Extract and process human population data =============================

# 1. Extract annual district-level estimates from WorldPop, and use linear interpolation to impute to 1998
# These include total population, and different forms of population density (pop dens by area, pop weighted density)
# 2. Examine annual district-level population from HRW and use linear interpolation to project to 2019
# 3. Some basic cross-comparison of each source


# ------------- dependencies -----------

library(ggplot2)
library(rgdal)
library(raster)
library(dplyr)
library(sf)
library(exactextractr)
setwd("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/")

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))


# ------------ load data ------------

# vietnam polygons and harmonise crs
# combine with metadata on district names, and remove truong sa (offshore)
shp = st_read("./data/shapefiles/d_moss/vietnam_districts_merged.shp")

# population rasters and calculate density
pop = stack(list.files("./data/population/viet_population_worldpop/", pattern=".tif", full.names=TRUE))
areax = raster::area(pop[[1]])



# ================== extract population information for each district over time ===============

# extract area and populations for each district
ex.area = exactextractr::exact_extract(areax, shp)
ex.pop = exactextractr::exact_extract(pop, shp)

# # total population (sum)
# popTotal = function(x){
#   cat(paste(x, "...", sep=""))
#   px = ex.pop[[ x ]]
#   px[ , 1:(ncol(px)-1 ) ] = px[ , 1:(ncol(px)-1 ) ] * px$coverage_fraction
#   result = colSums(px[ , 1:(ncol(px)-1 ) ], na.rm=TRUE)
#   return(as.data.frame(t(as.data.frame(result))) %>% mutate())
# }

# mean or sd population density by area (everywhere weighted equally, i.e. mean pop dens per pixel)
popDensArea = function(x, func="mean"){
  cat(paste(x, "...", sep=""))
  px = ex.pop[[ x ]]
  px$area = ex.area[[ x ]]$value
  px[ , grep("vnm_ppp", names(px)) ] = px[ , grep("vnm_ppp", names(px)) ] / px$area # calculate density
  px[ , grep("vnm_ppp", names(px)) ] = px[ , grep("vnm_ppp", names(px)) ] * px$coverage_fraction # by coverage fraction
  if(func == "mean"){ result = colMeans(px[ , grep("vnm_ppp", names(px)) ], na.rm=TRUE) }
  if(func == "sd"){ result = apply(px[ , grep("vnm_ppp", names(px)) ], 2, sd, na.rm=TRUE) }
  return(as.data.frame(t(as.data.frame(result))) %>% mutate())
}

# population weighted mean population density (i.e. average density per inhabitant rather than per area)
popDensPW = function(x){
  cat(paste(x, "...", sep=""))
  px = ex.pop[[ x ]]
  px$area = ex.area[[ x ]]$value
  popx = px[ , grep("vnm_ppp", names(px)) ]
  popweights = apply(popx, 2, function(vec){ vec / sum(vec, na.rm=TRUE) }) # calculate population weights
  popx = popx / px$area # calculate density
  popx = popx * popweights # weight by population
  result = colSums(popx * px$coverage_fraction, na.rm=TRUE) # by coverage fraction
  return(as.data.frame(t(as.data.frame(result))) %>% mutate())
}

# run extractions for different density metrics
#r1 = do.call(rbind.data.frame, lapply(1:length(ex.pop), popTotal))
r2 = do.call(rbind.data.frame, lapply(1:length(ex.pop), popDensArea, func="mean"))
r3 = do.call(rbind.data.frame, lapply(1:length(ex.pop), popDensArea, func="sd"))
r4 = do.call(rbind.data.frame, lapply(1:length(ex.pop), popDensPW))

# wrapper function for combining with districts info
combx = function(resx, variable_name){
  shp %>% 
    sf::st_drop_geometry() %>%
    cbind(resx) %>%
    dplyr::select(-areanamelo) %>%
    reshape2::melt(id.vars=1:4) %>%
    dplyr::mutate(Year = as.numeric(unlist(lapply(strsplit(as.vector(variable), "_"), "[", 3)))) %>%
    dplyr::rename(!!variable_name := value) %>%
    dplyr::select(-variable)
}

# save
#res1 = combx(r1, "Population_WorldPop")
res2 = combx(r2, "PopDens_Mean_WorldPop")
res3 = combx(r3, "PopDens_SD_WorldPop")
res4 = combx(r4, "PopDens_PW_WorldPop")

# combine
result = left_join(res2, left_join(res3, res4))
write.csv(result, "./output/data_processed/population/PopulationMetrics_PerDistrict_WorldPop_Nov2020.csv", row.names=FALSE)



# impute from 1998-1999 for each district
# using same rate as subseqeuent 3 years

# conduct imputation for each time series
imputeTimeSeriesPop = function(district){
  
  dfx = result[ result$areaid == district, ]
  dfx = do.call(rbind.data.frame, list(dfx[ 1, ], dfx[ 1, ], dfx))
  dfx$Year[1:2] = 1998:1999
  
  # change rates start and end
  #rate_start1 = mean(dfx$Population_WorldPop[ dfx$Year %in% 2001:2005] / dfx$Population_WorldPop[ dfx$Year %in% 2000:2004])
  rate_start2 = mean(dfx$PopDens_Mean_WorldPop[ dfx$Year %in% 2001:2005] / dfx$PopDens_Mean_WorldPop[ dfx$Year %in% 2000:2004])
  rate_start3 = mean(dfx$PopDens_PW_WorldPop[ dfx$Year %in% 2001:2005] / dfx$PopDens_PW_WorldPop[ dfx$Year %in% 2000:2004])

  # calculate imputed
  # dfx$Population_WorldPop[ dfx$Year == 1999 ] = dfx$Population_WorldPop[ dfx$Year == 2000 ] / rate_start1
  # dfx$Population_WorldPop[ dfx$Year == 1998 ] = dfx$Population_WorldPop[ dfx$Year == 1999 ] / rate_start1
  dfx$PopDens_Mean_WorldPop[ dfx$Year == 1999 ] = dfx$PopDens_Mean_WorldPop[ dfx$Year == 2000 ] / rate_start2
  dfx$PopDens_Mean_WorldPop[ dfx$Year == 1998 ] = dfx$PopDens_Mean_WorldPop[ dfx$Year == 1999 ] / rate_start2
  dfx$PopDens_PW_WorldPop[ dfx$Year == 1999 ] = dfx$PopDens_PW_WorldPop[ dfx$Year == 2000 ] / rate_start3
  dfx$PopDens_PW_WorldPop[ dfx$Year == 1998 ] = dfx$PopDens_PW_WorldPop[ dfx$Year == 1999 ] / rate_start3
  dfx$PopDens_SD_WorldPop[ dfx$Year == 1998 ] = NA
  dfx$PopDens_SD_WorldPop[ dfx$Year == 1999 ] = NA

  # return
  return(dfx)
}

# run imputation
result_imp = do.call(rbind.data.frame, lapply(unique(result$areaid), imputeTimeSeriesPop))

# view
ggplot(result_imp[ result_imp$areanameen == "Dong Da", ], aes(Year, PopDens_PW_WorldPop, group=areaid)) +
  geom_line() 

# save
write.csv(result_imp, "./output/data_processed/population/PopulationMetrics_PerDistrict_WorldPopimputed_Nov2020.csv", row.names=FALSE)




# ========================= GPW population estimates used to calculate population density over area =======================

# gpw estimates
gpw = read.csv("./output/data_processed/population/VietPopulation_District_MergedSHP_GPW4_interpolated.csv")

# get area for each and calculate crude population denstiy (pop/area)
shp$area = as.vector(st_area(shp))/10^6
gpw = left_join(gpw, shp[ , c("areaid", "area")] %>% st_drop_geometry())
gpw$popdens_gpw = gpw$population_gpw / gpw$area

# combine with weighted density estimates from WP and save
gpw = left_join(gpw, 
                result_imp[ , c("areaid", "Year", "PopDens_Mean_WorldPop", "PopDens_PW_WorldPop")] %>%
                  dplyr::rename("year"=Year, "popdens_area"=PopDens_Mean_WorldPop, "popdens_pw"=PopDens_PW_WorldPop))
write.csv(gpw, "./code/viet_dengue_districts/output/covariates/Population_GPW_WP.csv", row.names=FALSE)


