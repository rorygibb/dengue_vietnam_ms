

# ============================= Land cover change: extract data at district-level ==============================

# project root and dependencies
setwd("C:/Users/roryj/Documents/PhD/202007_LSHTM_Dengue/analysis/")
pacman::p_load("ggplot2", "rgdal", "raster", "dplyr", "sf", "exactextractr", "ggspatial")

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))

# viet districts shapefile
shp = st_read("./data/shapefiles/d_moss/vietnam_districts_merged.shp") %>%
  filter(!areanameen %in% c("Truong Sa", "Hoang Sa")) 






# ====================== LIU URBAN DYNAMICS: urban land cover change per commune ==========================

# urban dynamics raster: grid cell entry == year of urban conversion (base year == 1985)
urb = raster("./output/data_processed/landcover/li_urbandynamics/urban_grid_combined/urban_grid_combined.tif")
areax = raster::area(urb); names(areax) = "area"
urban = stack(urb, areax)
names(urban) = c("urban", "area")
st_crs(shp) = crs(urb)

# function to extract urban cover across all years for the xth polygon
extractUrbanChange = function(x){
  
  # setup
  poly = shp[ x, ]
  poly_area = st_area(poly)/10^6
  print(poly$areaid)
  
  # extraction (add dummy variable so areas with no urban are still counted but do not add to total)
  ex.urb = exactextractr::exact_extract(urban, poly)[[1]] 
  ex.urb$area = ex.urb$area * ex.urb$coverage_fraction
  dummy = data.frame(urban = 1985:2015, area=0, coverage_fraction=0)
  ex.urb = rbind(ex.urb, dummy)
  ex.urb = ex.urb[ ex.urb$urban != 0, ]

  # extract total urban cover per year from 1985 to 2015
  extractUrbTotal = function(year){ sum(ex.urb$area[ ex.urb$urban <= year]) }
  resx = data.frame(areaid = poly$areaid,
                    areanameen = poly$areanameen,
                    areaprovin = poly$areaprovin,
                    year = 1985:2015,
                    urbancover_km2 = sapply(1985:2015, extractUrbTotal))
  
  # impute forward 2016-2019 assuming same rate of change as prev 3 years (consistent with ESACCI)
  years_to_calc1 = 2012:2014
  change_rate = mean(resx$urbancover_km2[ resx$year %in% (years_to_calc1+1) ] - resx$urbancover_km2[ resx$year %in% (years_to_calc1) ])
  foo = resx[ rep(nrow(resx), 5), ] %>%
    dplyr::mutate(year = 2016:2020)
  foo$urbancover_km2 = foo$urbancover_km2 + (change_rate * 1:5)
  resx = rbind(resx, foo)
  
  # impute backward 1985-1982 assuming same rate of change 
  years_to_calc2 = 1985:1987
  change_rate = mean(resx$urbancover_km2[ resx$year %in% (years_to_calc2+1) ] - resx$urbancover_km2[ resx$year %in% (years_to_calc2) ])
  foo = resx[ rep(1, 3), ] %>%
    dplyr::mutate(year = 1982:1984)
  foo$urbancover_km2 = foo$urbancover_km2 - (change_rate * 3:1)
  resx = rbind(foo, resx)
  
  # label interpolated
  resx$liu_type = ifelse(resx$year %in% c(1982:1984, 2016:2020), "lin_interp", "landsat")
  
  # calculate urban change in km2 between years
  # calculate proportion urban per year
  resx$urban_prop = as.vector(resx$urbancover_km2/poly_area)
  resx$urbanexp_km2 = c(NA, resx$urbancover_km2[ 2:nrow(resx)] - resx$urbancover_km2[ 1:(nrow(resx)-1)])
  
  # calculate urban change as a proportion of existing urban cover
  resx$urbanexp_prop = log(c(NA, resx$urbancover_km2[ 2:nrow(resx)] / resx$urbancover_km2[ 1:(nrow(resx)-1)]))
  
  # replace NaN values
  resx$urbancover_km2[ is.nan(resx$urbancover_km2) ] = 0
  resx$urban_prop[ is.nan(resx$urban_prop) ] = 0
  resx$urbanexp_km2[ is.nan(resx$urbanexp_km2) ] = 0
  resx$urbanexp_prop[ is.nan(resx$urbanexp_prop) ] = 0

  # return
  return(resx)
}

# run extraction
result = lapply(1:nrow(shp), extractUrbanChange)
result = do.call(rbind.data.frame, result)
write.csv(result, "./output/data_processed/landcover/li_urbandynamics/VietnamAll_MergedSHP2_UrbanisationRatesLiu_Districts_Feb2021.csv", row.names=FALSE)

# sss = full_join(shpc, result)
# ggplot(sss[ sss$year > 1999, ]) +
#   geom_sf(aes(fill=urbanexp_km2), col=NA) +
#   facet_wrap(~year) +
#   scale_fill_viridis_c(option="magma", direction=-1)
# ggplot(result) + 
#   geom_line(aes(year, urbanexp_km2)) +
#   facet_wrap(~areanameen)





# =================== HANSEN FOREST LOSS  =======================


# ----------- combine separate rasters into full Vietnam layer -------------

# baseline (2000)
# b1 = raster("./data/landcover/hansen_forestloss/v1.7/Hansen_GFC-2019-v1.7_treecover2000_10N_100E.tif")
# b2 = raster("./data/landcover/hansen_forestloss/v1.7/Hansen_GFC-2019-v1.7_treecover2000_20N_100E.tif")
# b3 = raster("./data/landcover/hansen_forestloss/v1.7/Hansen_GFC-2019-v1.7_treecover2000_30N_100E.tif")
# baseline = raster::merge(b1, b2)
# baseline = raster::merge(baseline, b3)
# writeRaster(baseline, file="./output/data_processed/landcover/hansen/vietnam_forestcover_baseline_2000.tif", format="GTiff", overwrite=TRUE)

# loss by year
# l1 = raster("./data/landcover/hansen_forestloss/v1.7/Hansen_GFC-2019-v1.7_lossyear_10N_100E.tif")
# l2 = raster("./data/landcover/hansen_forestloss/v1.7/Hansen_GFC-2019-v1.7_lossyear_20N_100E.tif")
# l3 = raster("./data/landcover/hansen_forestloss/v1.7/Hansen_GFC-2019-v1.7_lossyear_30N_100E.tif")
# loss = raster::merge(l1, l2)
# loss = raster::merge(loss, l3)
# writeRaster(loss, file="./output/data_processed/landcover/hansen/vietnam_forestcover_lossyear_20012019.tif", format="GTiff", overwrite=TRUE)
# 
# # area
# areax = raster::area(baseline)
# names(areax) = "area"
# writeRaster(areax, file="./output/data_processed/landcover/hansen/vietnam_arearaster_derived.tif", format="GTiff")



# ---------------- extract annual forest change per commune polygon --------------------

# forest rasters: baseline (2000 prop cover); loss year; area
baseline = raster("./output/data_processed/landcover/hansen/vietnam_forestcover_baseline_2000.tif")
loss = raster("./output/data_processed/landcover/hansen/vietnam_forestcover_lossyear_20012019.tif")
areax = raster("./output/data_processed/landcover/hansen/vietnam_arearaster_derived.tif")
forest = stack(baseline, loss, areax)
names(forest) = c("baseline", "loss", "area")

# harmonise crs
st_crs(shp) = crs(loss)

# function to extract forest change across all years for the xth polygon
extractForestChange = function(x){
  
  # calculate column of baseline forested area within grid cell (proportion forest * area * coverage)
  poly = shp[ x, ]
  print(x)
  ex.for = exactextractr::exact_extract(forest, poly)[[1]] 
  ex.for$baseline_area = (ex.for$baseline/100) * ex.for$area * ex.for$coverage_fraction
  
  # test: is area < total area
  #sum(ex.for$baseline_area) <= sum(ex.for$area)
  
  # for areas with no loss
  if(sum(ex.for$loss) == 0){
    resx = data.frame(year = 1995:2019, 
                      forestloss_area = 0,
                      areaid = poly$areaid,
                      areanameen = poly$areanameen,
                      areaprovin = poly$areaprovin,
                      data_type = "hansen_landsat",
                      forestcover_total = sum(ex.for$baseline_area, na.rm=TRUE), 
                      forestcover_prop =  sum(ex.for$baseline_area, na.rm=TRUE) / sum(ex.for$area, na.rm=TRUE))
  } else{
    
    # add dummy rows for all years
    dummy = data.frame(baseline = 0, 
                       loss = 1:19,
                       area = 0,
                       coverage_fraction = 0, 
                       baseline_area = 0)
    ex.for = rbind(ex.for, dummy)
    
    # summarise annual loss
    resx = ex.for[ ex.for$loss != 0 & !is.na(ex.for), ] %>%
      dplyr::group_by(loss) %>%
      dplyr::summarise(forestloss_area = sum(baseline_area, na.rm=TRUE)) %>%
      dplyr::rename("year" = loss) %>%
      dplyr::mutate(year = year + 2000,
                    areaid = poly$areaid,
                    areanameen = poly$areanameen,
                    areaprovin = poly$areaprovin,
                    GID_3 = poly$GID_3,
                    data_type = "hansen_landsat")
    
    # interpolate to 1996 assuming rate of forest cover loss is average of past 3 years (this is problematic because rates can change v quickly)
    res.int = resx[ 1:6, ] %>%
      dplyr::mutate(year = 1995:2000,
                    forestloss_area = mean(forestloss_area[1:3]),
                    data_type = "lin_interp")
    resx = rbind(res.int, resx)
    resx$forestcover_total = sum(c(ex.for$baseline_area, res.int$forestloss_area)) - cumsum(resx$forestloss_area)
    resx$forestcover_prop = resx$forestcover_total / sum(ex.for$area, na.rm=TRUE)
    
    # replace NaN values
    resx$forestcover_total[ is.nan(resx$forestcover_total) ] = 0
    resx$forestcover_prop[ is.nan(resx$forestcover_prop) ] = 0
    resx$forestloss_area[ is.nan(resx$forestloss_area) ] = 0
  }
  
  # ggplot(resx) + geom_line(aes(year, forestloss_area))
  # ggplot(resx) + geom_line(aes(year, forestcover_total))
  
  # return
  return(resx)
}

# run
res = lapply(1:nrow(shp), extractForestChange)
result = do.call(rbind.data.frame, res)

# remove NAs and any extrapolation where forest cover < 0 
result = result[ !is.na(result$year), ]
result$forestcover_prop[ result$forestcover_prop < 0 ] = 0

# save
write.csv(result, "./output/data_processed/landcover/hansen/VietnamAll_MergedSHP2_ForestChangeHansen_Districts_Feb2021.csv", row.names=FALSE)




# # ======================= ESA-CCI LAND COVER CHANGE METRICS (lower resolution) =========================
# 
# # rasters
# loc = "./data/landcover/esacci_viet/"
# esa = stack(list.files(loc, pattern=".tif", full.names=TRUE))
# names(esa) = paste("esa", 1992:2018, sep="_")
# areax = area(esa[[1]]); names(areax) = "area"
# esa = stack(esa, areax)
# st_crs(shpc) = crs(esa)
# 
# # function to extract urban cover across all years for the xth polygon
# # vals is a vector of land cover class ids from ESA-CCI data (e.g. 190 for urban)
# # var_name specifies name of variable for column
# extractESAChange = function(x, class_ids, var_name){
# 
#   # setup
#   poly = shpc[ x, ]
#   poly_area = st_area(poly)/10^6
#   print(x)
# 
#   # extraction (add dummy variable so areas with no urban are still counted but do not add to total)
#   ext = exactextractr::exact_extract(esa, poly)[[1]]
#   ext$area = ext$area * ext$coverage_fraction
# 
#   # function to calculate cover of specified 
#   funcx = function(x){  ifelse(x %in% class_ids, 1, 0) * ext$area }
#   resx = data.frame(areaid = poly$areaid,
#                     areanameen = poly$areanameen,
#                     areaprovin = poly$areaprovin,
#                     GID_3 = poly$GID_3,
#                     year = 1992:2018,
#                     cover_km2 = colSums(apply(ext[ , grep("esa", names(ext)) ], 2, funcx), na.rm=TRUE))
#   
#   # impute forward to 2019 assuming same rate of change as prev 3 years
#   years_to_calc1 = 2015:2017
#   change_rate = mean(resx$cover_km2[ resx$year %in% (years_to_calc1+1) ] / resx$cover_km2[ resx$year %in% (years_to_calc1) ])
#   foo = resx[ nrow(resx), ] %>%
#     dplyr::mutate(year = 2019)
#   foo$cover_km2 = foo$cover_km2 * change_rate
#   resx = rbind(resx, foo)
#   
#   # impute backward to 1990 assuming same rate of change as earliest 3 years
#   years_to_calc2 = 1992:1994
#   change_rate = mean(resx$cover_km2[ resx$year %in% (years_to_calc2+1) ] / resx$cover_km2[ resx$year %in% (years_to_calc2) ])
#   foo = resx[ rep(1, 2), ] %>%
#     dplyr::mutate(year = 1990:1991)
#   foo$cover_km2[2] = foo$cover_km2[2] / change_rate
#   foo$cover_km2[1] = foo$cover_km2[2] / change_rate
#   resx = rbind(foo, resx)
#   row.names(resx) = c()
#   
#   # replace NaNs
#   resx$cover_km2 = replace(resx$cover_km2, is.nan(resx$cover_km2), 0)
# 
#   # calculate  change in km2 between years
#   # calculate proportion urban per year
#   resx$prop = as.vector(resx$cover_km2/sum(ext$area, na.rm=TRUE))
#   resx$prop[ resx$prop > 1] = 1
#   resx$exp_km2 = c(NA, resx$cover_km2[ 2:nrow(resx)] - resx$cover_km2[ 1:(nrow(resx)-1)])
# 
#   # change as a proportion of existing urban cover
#   resx$exp_prop = log(c(NA, resx$cover_km2[ 2:nrow(resx)] / resx$cover_km2[ 1:(nrow(resx)-1)]))
#   resx$exp_prop = replace(resx$exp_prop, is.nan(resx$exp_prop), 0)
#   
#   # rename columns to land cover type
#   names(resx)[ (ncol(resx)-3):ncol(resx) ] = paste(var_name, names(resx)[ (ncol(resx)-3):ncol(resx) ], sep="")
#   
#   # sense check
#   #ggplot(resx) + geom_line(aes(year, exp_km2))
#   return(resx)
# }
# 
# # run and save for forest cover
# result = lapply(1:nrow(shpc), extractESAChange, class_ids = c(50:100, 160, 170), var_name = "forest")
# res2 = do.call(rbind.data.frame, result)
# #write.csv(res2, "./output/data_processed/landcover/esa_cci/vietnam_esacci_forest_percommune.csv", row.names=FALSE)
# 
# # run and save for urban cover
# result = lapply(1:nrow(shpc), extractESAChange, class_ids = c(190), var_name = "urban")
# res2 = do.call(rbind.data.frame, result)
# #write.csv(res2, "./output/data_processed/landcover/esa_cci/vietnam_esacci_urban_percommune.csv", row.names=FALSE)


# 
# 
# 
# # ================== calculate population-weighted land cover change metrics at the district level (across all communes) ===============
# 
# # population 
# pop = read.csv("./output/data_processed/population/VietnamAll_Population_Commune_WorldPop_Feb2021.csv", stringsAsFactors = FALSE) %>%
#   dplyr::select(year, GID_3, population)
# table(pop$year[ is.na(pop$population) ])
# 
# # 1. forest loss: calculate cumulative losses in 3, 5, 10 years windows (short, med, long-ish time horizon, helps deal with uncertainty in year of conversion)
# # issue: Hansen data only present from 2000 so early years of time series require interpolation
# forest = read.csv("./output/data_processed/landcover/hansen/VietnamAll_ForestChange_Hansen_Communes_Feb2021.csv", stringsAsFactors = FALSE) %>%
#   dplyr::filter(!is.na(GID_3)) %>%
#   dplyr::arrange(GID_3, year) %>%
#   group_by(GID_3) %>%
#   dplyr::mutate(forestloss_3yr = data.table::frollsum(forestloss_area, 3, align="right"),
#                 forestloss_5yr = data.table::frollsum(forestloss_area, 5, align="right"),
#                 forestloss_10yr = data.table::frollsum(forestloss_area, 10, align="right"),
#                 forestloss_cumulative = cumsum(forestloss_area)) %>%
#   left_join(pop)
# 
# # aggregate to the district level (population-weighted mean of commune-level change)
# for_pop = forest %>%
#   dplyr::group_by(areanameen, year) %>%
#   dplyr::mutate(pop_weight = population / sum(population, na.rm=TRUE)) %>%
#   dplyr::summarise(forestloss_3yr_cw = sum(forestloss_3yr * pop_weight, na.rm=TRUE),
#                    forestloss_5yr_cw = sum(forestloss_5yr * pop_weight, na.rm=TRUE),
#                    forestloss_10yr_cw = sum(forestloss_10yr * pop_weight, na.rm=TRUE),
#                    forestloss_cumulative_cw = sum(forestloss_cumulative * pop_weight, na.rm=TRUE),
#                    areaprovin = head(areaprovin, 1))
# 
# # set years outside window at start of time series to NA
# for_pop$forestloss_10yr_cw[ for_pop$year < 2004 ] = NA
# for_pop$forestloss_5yr_cw[ for_pop$year < 1999 ] = NA
# for_pop$forestloss_3yr_cw[ for_pop$year < 1997 ] = NA
# 
# # prov = "Thua Thien - Hue"
# # ggplot(for_pop[ for_pop$areaprovin == prov, ]) +
# #   geom_line(aes(year, forestloss_10yr_cw)) +
# #   facet_wrap(~areanameen)
# 
# 
# # 2. urban expansion rates 
# urban = read.csv("./output/data_processed/landcover/li_urbandynamics/VietnamAll_Urbanisation_Communes_Feb2021.csv", stringsAsFactors = FALSE) %>%
#   dplyr::filter(!is.na(GID_3)) %>%
#   dplyr::arrange(GID_3, year) %>%
#   group_by(GID_3) %>%
#   dplyr::mutate(u2 = replace(urbanexp_km2, is.na(urbanexp_km2), 0)) %>%
#   dplyr::rename("urbanexp" = urbanexp_km2) %>%
#   dplyr::mutate(urbanexp_3yr = data.table::frollsum(urbanexp, 3, align="right"),
#                 urbanexp_5yr = data.table::frollsum(urbanexp, 5, align="right"),
#                 urbanexp_10yr = data.table::frollsum(urbanexp, 10, align="right"),
#                 urbanexp_15yr = data.table::frollsum(urbanexp, 15, align="right"),
#                 urbanexp_cumulative = cumsum(u2)) %>%
#   left_join(pop) %>%
#   dplyr::select( -urbanexp_prop, -u2 )
#               
# # calculate mean at the district level (population-weighted)
# urb_pop = urban %>%
#   dplyr::group_by(areanameen, year) %>%
#   dplyr::mutate(pop_weight = population / sum(population, na.rm=TRUE)) %>%
#   dplyr::summarise(urbanexp_3yr_cw = sum(urbanexp_3yr * pop_weight, na.rm=TRUE),
#                    urbanexp_5yr_cw = sum(urbanexp_5yr * pop_weight, na.rm=TRUE),
#                    urbanexp_10yr_cw = sum(urbanexp_10yr * pop_weight, na.rm=TRUE),
#                    urbanexp_15yr_cw = sum(urbanexp_15yr * pop_weight, na.rm=TRUE),
#                    urbanexp_cumulative_cw = sum(urbanexp_cumulative * pop_weight, na.rm=TRUE),
#                    areaprovin = head(areaprovin, 1))
# 
# # combine urban and forest change together, and save
# commune_dat = left_join(for_pop[ for_pop$year >= 1997, ], urb_pop[ urb_pop$year >= 1997, ])
# write.csv(commune_dat, file="./output/data_processed/landcover/landcoverchange_popweighted_commune_hansenliu.csv", row.names=FALSE)
# 
# 
# 
# prov = "Tay Ninh"
# ggplot(urb_pop[ urb_pop$areaprovin == prov, ]) +
#   geom_line(aes(year, urbanexp_10yr_cw)) +
#   facet_wrap(~areanameen)
# 
# aa = read.csv("./output/data_processed/landcover/landcoverchange_popweighted_commune_hansenliu.csv")
# ggplot(aa[ aa$areaprovin == prov, ]) +
#   geom_line(aes(year, urbanexp_km2_pw)) +
#   facet_wrap(~areanameen)
