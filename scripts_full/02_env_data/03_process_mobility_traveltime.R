


# =========================== estimate land-based travel distance between districts =================================

# using Malaria Atlas Project roads-based friction surface: https://malariaatlas.org/research-project/accessibility_to_cities/
# uses Google Maps and OSM to estimate travel times to cities for nomial year 2015 (but effectively just current roads)
# but can also be used w/ MAP tools to estimate time between 2 locations using friction surface
# do this between every pair of districts in Vietnam: creates national-level connectivity matrix
# https://medium.com/@abertozz/mapping-travel-times-with-malariaatlas-and-friction-surfaces-f4960f584f08

setwd("C:/Users/roryj/Documents/PhD/202007_LSHTM_Dengue/analysis/")

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

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))



# =============================== data ====================================

# polygons for entire of Vietnam (same areaid matches)
# specify unique field identifier (called "areaid") for processing
shp = st_read("./data/shapefiles/d_moss/vietnam_districts_merged.shp")
shp = shp[ ! shp$areanameen %in% c("Truong Sa", "Hoang Sa"), ] # remove offshore districts
shp$areaid = shp$areaid

# population density to use as weighting (aggregate to same res as friction surface, i.e. *10)
# pop = raster("./data/population/viet_population_worldpop/vnm_ppp_2015.tif")
# pop = aggregate(pop, fact=10, fun=sum)
# writeRaster(pop, "./output/data_processed/population/VietPop2015_1kmres.tif", format="GTiff")
pop = raster("./output/data_processed/population/VietPop2015_1kmres.tif")

# friction surface from malariaAtlas project
# friction = malariaAtlas::getRaster(
#   surface = "A global friction surface enumerating land-based travel speed for a nominal year 2015",
#   shp = as_Spatial(shp))
#malariaAtlas::autoplot_MAPraster(friction)
#writeRaster(friction, "./output/data_processed/connectivity/traveltimematrix_map/frictionRaster_MAP_Vietnam.tif", format="GTiff")
friction = raster("./output/data_processed/connectivity/traveltimematrix_map/frictionRaster_MAP_Vietnam.tif")

# ensure same sizes
pop = crop(pop, friction)
pop = resample(pop, friction, fun="ngb")
  


# ===================== 1. pairwise travel times between focal polygons and all others ==========================

# calcTravelTimesMatrix function
#' @param provinces specify provinces of Vietnam to subset (could be easily generalised)
#' @param friction friction surface raster (generated above from MAP)
#' @param pop population raster at same resolution as friction surface (see above)
#' @param method either centroid, popcentre (both fast) or popweight (slower)

calcTravelTimesMatrix = function(provinces, friction, pop, method="centroid"){
  
  # reporting
  if(!method %in% c("centroid", "popcentre", "popweight")){ return("Method must be one of centroid, popcentre or popweight") }
  if(method %in% c("popcentre", "popweight") & paste(dim(pop), collapse=" ") != paste(dim(friction), collapse=" ")){ "Check extent and dimensions of pop and friction rasters"}
  
  # subset spatial extent to relevant provinces
  # if "vietnam" include entire country, otherwise subset to specified province(s)
  # these can be changed depending on shapefile fields
  if(provinces == "Vietnam"){ 
    shpx = shp 
  } else{
    shpx = shp[ shp$areaprovin %in% provinces, ]
  }
  ff = crop(friction, shpx)
  pp = crop(pop, shpx); names(pp) = "population"
  
  # create movement matrix and apply geocorrection
  # https://medium.com/@abertozz/mapping-travel-times-with-malariaatlas-and-friction-surfaces-f4960f584f08
  tmat = gdistance::transition(ff, function(x) 1/mean(x), 8) 
  tmat_gc = gdistance::geoCorrection(tmat)
  
  # create matrix for results
  res_mat = matrix(nrow=nrow(shpx), ncol=nrow(shpx))
  dimnames(res_mat) = list(shpx$areaid, shpx$areaid)
  
  
  # --------------------------------------------------------------------------------------------
  
  # Three different methods:
  # 1. centroid: calculates population weighted mean travel times from centroid of origin polygon to all other (destination) polygons
  # 2. popcentre: calculates population weighted mean travel times from highest population grid cell of origin polygon, to all other(destination) polygons
  # 3. popweight: approximates full population weighted mean travel times with weighting applied in both origin and destination polygons (slower)
    # i.e. the mean travel time between any pair of persons in origin and destination polygons
    # popweight does this selecting a representative sample of point locations within origin polygon (sample weighted by population), calculating travel distance for each
    # and then calculating the mean population-weighted mean travel time to destination across all sample points in origin polygon
  
  
  # ------------- methods based on a single point location in origin polygon -------------------
  
  if(method != "popweight"){
    
    # origin is centroid
    if(method == "centroid"){
      
      # create centroids
      locs = suppressWarnings(as.data.frame(sf::st_coordinates(sf::st_centroid(shpx))))
      names(locs) = c("x", "y")
      locs$areaid = shpx$areaid
      #plot(pp); points(locs$x, locs$y, col="red"); plot(shpx$geometry, add=TRUE)
      
      # ensure centroid doesn't fall in NA space; if so, select other random point within polygon
      locs$test = raster::extract(pp, SpatialPoints(coords = locs[, c("x", "y")]))
      to_reselect = shpx[ shpx$areaid %in% locs$areaid[ is.na(locs$test)], ]
      for(i in 1:nrow(to_reselect)){ 
        px = as.data.frame(spsample(as_Spatial(to_reselect[i, ]), n=1, type = "random"))
        locs$x[ locs$areaid == to_reselect$areaid[i] ] = px$x
        locs$y[ locs$areaid == to_reselect$areaid[i] ] = px$y
      }
      # locs$test = raster::extract(pp, SpatialPoints(coords = locs[, c("x", "y")]))
      # sum(is.na(locs$test))
      
      # manually add centroids for several coastal/island areas
      locs[ locs$areaid == 70178, c("x", "y")] = c(107.468, 21.1471) # Van Don
      locs[ locs$areaid == 70359, c("x", "y")] = c(108.1311, 16.0983) # Lien Chieu
      locs[ locs$areaid == 70671, c("x", "y")] = c(104.6367, 9.8087) # Kien Hai
      locs[ locs$areaid == 70154, c("x", "y")] = c(107.7527, 20.9973) # Kien Hai
    }
    
    # origin is grid cell with greatest population density
    if(method == "popcentre"){
      
      dist_ras = raster::rasterize(shpx, ff, field="areaid"); names(dist_ras) = "areaid"
      rr = as.data.frame(raster::stack(pp, dist_ras), xy=TRUE)
      rr = na.omit(rr)
      getDensLoc = function(x){
        rrx = rr[ rr$areaid == x, ] 
        rrx[ rrx$population == max(rrx$population), c("x", "y", "areaid")] 
      }
      locs = do.call(rbind.data.frame, lapply(unique(rr$areaid), getDensLoc))
      #plot(pp); points(locs$x, locs$y, col="red"); plot(shpx$geometry, add=TRUE)
    }
    
    # for each ith district, estimate population weighted travel time to all other districts
    # create accessibility raster, replace inf values with NA, extract and run function for pop weighting
    for(i in 1:nrow(locs)){
      
      cat(paste(i, "...", sep=""))
      access.raster = gdistance::accCost(tmat_gc, as.matrix(locs[i, c("x", "y")]))
      values(access.raster)[ values(access.raster) == Inf ] = NA
      #plot(access.raster, col=rev(viridis::magma(200))); points(locs[i, 1], locs[i, 2], pch=16, col="darkgreen"); plot(shpx$geometry, add=T)
      
      # extract and run pop weight function
      extr = exactextractr::exact_extract(stack(access.raster, pp), shpx)
      popweightTravTime = function(x){ 
        x = na.omit(x)
        x$popweight = x$population / sum(x$population)
        sum(x$layer * x$popweight * x$coverage_fraction)
        }
      extr = unlist(lapply(extr, popweightTravTime))
      
      # set distance 0 to NAs (as these are offshore)
      extr[ extr == 0 ] = NA
      
      # add into matrix, optionally set diagonal to 0 if calculating mean travel time (i.e. 0 to same district)
      res_mat[ , i ] = extr
      #res_mat[ i, i] = 0
    }
    # ggplot(data =  reshape2::melt(res_mat), aes(x=as.factor(Var1), y=as.factor(Var2), fill=value)) +
    #   geom_tile() + scale_fill_viridis_c(option="magma") +
    #   theme(axis.title=element_blank(), axis.text.x = element_text(angle=90))
    
  }
  
  
  # ---------------- method 'popweight' approximating per-person average travel time between districts ------------------
  
  # approximation of population weighted travel time across both origin and destination
  # selects n locations within district at specified sampling density, with probability of selecting location proportional to population (i.e. pop weighted)
  # calculate distance rasters for each point location
  # extract travel time to all other districts and calculate mean across all origin rasters
 
   if(method == "popweight"){
    
    # ---- 1. select sampling locations
    
    # extract raster to dataframe for all locations
    dist_ras = rasterize(shpx, ff, field="areaid"); names(dist_ras) = "areaid"
    rr = as.data.frame(stack(pp, dist_ras), xy=TRUE)
    rr = na.omit(rr)
    
    # specify density of sampling (1 location per x km2) and calculate num points to be selected per polygon based on area
    # ensure each polygon has at least 1 point
    samp_dens = 50
    shpx$samp_effort = as.numeric(round((st_area(shpx)/10^6)/samp_dens))
    shpx$samp_effort[ shpx$samp_effort == 0 ] = 1
    rr = left_join(rr, shpx[ , c("areaid", "samp_effort")] %>% sf::st_drop_geometry())
    
    # for each polygon, sample n points weighted by pop'n
    pointsSample = function(x){
      rrx = rr[ rr$areaid == x, ]
      rrx$probweight = rrx$population/sum(rrx$population)
      rrx[ sample(1:nrow(rrx), replace=FALSE, size=rrx$samp_effort[1], prob=rrx$probweight), c("x", "y", "areaid") ]
      }
    locs = do.call(rbind.data.frame, lapply(unique(rr$areaid), pointsSample))
    plot(pp); points(locs$x, locs$y, col="red"); plot(shpx$geometry, add=T)
    
    
    # --- 2. create stack of travel time rasters for each location per focal district
    
    # for each polygon calculate travel time to destinations
    for(areaid in unique(locs$areaid)){
      
      # create travel time rasters for all points within origin polygon
      print(paste0("Processing: ", areaid, sep=""))
      locsx = locs[ locs$areaid == areaid, ]
      createAccessRaster = function(x){
        arx = gdistance::accCost(tmat_gc, as.matrix(locsx[x, c("x", "y")]))
        values(arx)[ values(arx) == Inf ] = NA
        names(arx) = paste0("ttx_", x, sep="")
        arx
      }
      access.stack = do.call(raster::stack, lapply(1:nrow(locsx), createAccessRaster))
      
      # extract travel time to all destination polygons, and calculate population weighted mean across all origin locations
      extr = exactextractr::exact_extract(stack(access.stack, pp), shpx)
      popweightTravTime = function(x){
        
        # remove any polygons w/ over 50% NA coverage (corrects odd behaviour around coastline), and omit full NA rows
        over_thresh_na = as.numeric(which((colSums(is.na(x[ , grep("ttx", names(x)), drop=FALSE]))/nrow(x))>=0.50))
        if(length(over_thresh_na) > 0){ x = x[ , -over_thresh_na ] }
        x = na.omit(x)
        x$tt = rowMeans(x[ , grep("ttx", names(x)), drop=FALSE ])
        sum((x$tt * x$population * x$coverage_fraction)/sum(x$population * x$coverage_fraction))
      }
      extr = unlist(lapply(extr, popweightTravTime))
      extr[ extr == 0 ] = NA # set distance 0 to NA (offshore)
      rm(access.stack)
      
      # save to matrix
      mi = which(colnames(res_mat) == areaid)
      res_mat[ , mi ] = extr
      #res_mat[ mi, mi] = 0
    }
    # ggplot(data =  reshape2::melt(res_mat), aes(x=as.factor(Var1), y=as.factor(Var2), fill=value)) +
    #   geom_tile() + scale_fill_viridis_c(option="magma") +
    #   theme(axis.title=element_blank(), axis.text.x = element_text(angle=90))
  }
  
  # -------------------  create travelMatrix object to return ---------------------
  
  result = list(
    provinces = provinces,
    method = method,
    travelmatrix = res_mat
  )
  class(result) = "travelMatrix"
  return(result)
  
} # end of function


# test function and visualise; centroid is a good approximation of the full population weighted approach
# px = "Khanh Hoa"
# h1 = calcTravelTimesMatrix(provinces=px, friction=friction, pop=pop, method="centroid")
# h2 = calcTravelTimesMatrix(provinces=px, friction=friction, pop=pop, method="popcentre")
# h3 = calcTravelTimesMatrix(provinces=px, friction=friction, pop=pop, method="popweight")
# 
# id = 8
# fd = colnames(h1$travelmatrix)[id]
# compx = data.frame(centroid = h1$travelmatrix[ , id],
#                    popcentre = h2$travelmatrix[ , id], 
#                    popweight = h3$travelmatrix[ , id])
# gridExtra::grid.arrange(ggplot(compx) + geom_point(aes(centroid, popcentre)) + geom_abline() + coord_fixed(),
#                         ggplot(compx) + geom_point(aes(centroid, popweight)) + geom_abline() + coord_fixed(),
#                         ggplot(compx) + geom_point(aes(popcentre, popweight)) + geom_abline() + coord_fixed(),
#                         nrow=1, top=fd)


# run for Vietnam using centroid (good approximation of pop weight)
tt_viet = calcTravelTimesMatrix(provinces="Vietnam", friction=friction, pop=pop, method="centroid")
save(tt_viet, file="./output/data_processed/connectivity/VietnamAll_MergedSHP2_PairwiseTravelTime_Centroid_MAP.R")

# run for Vietnam using pop weight (time consuming)
# tt_viet = calcTravelTimesMatrix(provinces="Vietnam", friction=friction, pop=pop, method="popweight")
# save(tt_viet, file="./output/data_processed/connectivity/Vietnam_PairwiseTravelTime_popweight_MAP.R")

# manual adjustment of travel time for offshore districts
# 1. Ly Son (areaid = 70382) - 30 minute ferry from Sa Ky port, Binh Son district (70384), Quang Ngai
# 2. Phu Quy (areaid = 70596) - 3 hour ferry from Phan Thiet, Binh Thuan (70540)
# 3. Phu Quoc (70666) - 1 hour ferry from Ha Tien (70620), Kien Giang
# 4. Kien Hai (70671) - 1 hr from Rach Gia (70652) tourism wharf
# 5. Con Dao (70711) - 3.5 hours from Vung Tau (70615)

# adjust rows and columns
adjustMat = function(mx, from, to, add){
  self_to_self <- mx[ which(rownames(mx)==from),  which(colnames(mx)==from) ] 
  mx[ which(rownames(mx)==from), ] <- mx[ which(rownames(mx)==to), ] + add
  mx[ , which(colnames(mx)==from) ] <- mx[ , which(colnames(mx)==to) ] + add
  mx[ which(rownames(mx)==from),  which(colnames(mx)==from) ]  <- self_to_self
  return(mx)
}

# adjust for offshores
matrix_adj = tt_viet$travelmatrix %>%
  adjustMat(from=70382, to=70384, add=60) %>%
  adjustMat(from=70596, to=70540, add=180) %>%
  adjustMat(from=70666, to=70620, add=60) %>%
  adjustMat(from=70671, to=70652, add=60) %>%
  adjustMat(from=70711, to=70615, add=210) 
  
tt_viet$travelmatrix = matrix_adj
save(tt_viet, file="./output/data_processed/connectivity/VietnamAll_MergedSHP2_PairwiseTravelTime_Centroid_MAP_offshoreadjusted.R")




# ====================== 2. travel time to major city for each district ========================

# travel time to major city (50,000+ inhabitants) surface from malariaAtlas project
# citytrav = malariaAtlas::getRaster(
#   surface = "A global map of travel time to cities to assess inequalities in accessibility in 2015",
#   shp = as_Spatial(shp))
# malariaAtlas::autoplot_MAPraster(citytrav)
#writeRaster(citytrav, "./output/data_processed/connectivity/traveltimecity_map/traveltimetocity_MAP_Vietnam.tif", format="GTiff")
citytrav = raster("./output/data_processed/connectivity/traveltimecity_map/traveltimetocity_MAP_Vietnam.tif")

# combine with population and extract
cc = stack(citytrav, pop)
names(cc) = c("citytrav", "population")

# extract for vietnam using population weighting
extr = exactextractr::exact_extract(cc, shp)
popweightTravTime = function(x){
  if(all(is.na(x$citytrav))){ return(NA) }
  x = na.omit(x)
  x$popweight = x$population / sum(x$population)
  sum(x$popweight * x$citytrav * x$coverage_fraction)
}
shp_cc = shp
shp_cc$traveltime_city =  unlist(lapply(extr, popweightTravTime))
citytrav = shp_cc %>%
  st_drop_geometry() %>%
  dplyr::select(areaid, areanameen, areaprovin, traveltime_city)

# save
#write.csv(citytrav, "./output/data_processed/connectivity/Vietnam_CityTravelTime_MAP.csv", row.names=FALSE)




# ================== 3. travel time to commercial airport for each district ================

# get commercial airport locations
ap = read.csv("./data/flights/worldpop/Monthly-data/Airports_2010.csv", stringsAsFactors = FALSE) %>%
  dplyr::rename("NodeName" = 1, "x"=Lon, "y"=Lat) %>%
  dplyr::filter(!is.na(CountryCode) & CountryCode != "AN") %>%
  dplyr::mutate(Country = countrycode::countrycode(CountryCode, origin="iso2c", destination = "country.name")) %>%
  dplyr::filter(Country == "Vietnam")
  
# plot
#coordinates(ap) = ~x+y
# plot(shp$geometry)
# plot(ap, pch=16, col="red", add=TRUE)

# estimate travel time using friction surface
ff = friction
tmat = gdistance::transition(ff, function(x) 1/mean(x), 8) 
tmat_gc = gdistance::geoCorrection(tmat)
airport_trav = gdistance::accCost(tmat_gc, as.matrix(ap[, c("x", "y")]))
values(airport_trav)[ values(airport_trav) == Inf ] = NA

# stack and calculate population weighted mean
# combine with population and extract
cc = stack(airport_trav, pop)
names(cc) = c("airporttrav", "population")
extr = exactextractr::exact_extract(cc, shp)
popweightTravTime = function(x){
  if(all(is.na(x$airporttrav))){ return(NA) }
  x = na.omit(x)
  x$popweight = x$population / sum(x$population)
  sum(x$popweight * x$airporttrav * x$coverage_fraction)
}
shp_cc = shp
shp_cc$traveltime_airport =  unlist(lapply(extr, popweightTravTime))
airtrav = shp_cc %>%
  st_drop_geometry() %>%
  dplyr::select(areaid, areanameen, areaprovin, traveltime_airport)

# combine with citytrav
trav = left_join(citytrav, airtrav)
write.csv(trav, "./code/viet_dengue_districts/output/covariates/VietnamAll_MergedSHP2_CityAirportTravelTimes_MAP.csv", row.names=FALSE)


# =============== 4. Matrix of Euclidean distances between all pairs of polygons ================

# calculate Euclidean distance between most populated cell of each pair of polygons
# first extract most populated area
shp$idx = 1:nrow(shp)
dist_ras = raster::rasterize(shp, pop, field="idx"); names(dist_ras) = "idx"
rr = as.data.frame(raster::stack(pop, dist_ras), xy=TRUE)
rr = na.omit(rr)
getDensLoc = function(x){
  rrx = rr[ rr$idx == x, ] 
  rrx[ rrx$VietPop2015_1kmres == max(rrx$VietPop2015_1kmres), c("x", "y", "idx")] 
}
locs = do.call(rbind.data.frame, lapply(unique(rr$idx), getDensLoc))
locs = left_join(locs, shp[ , c("areaid", "idx")] %>% st_drop_geometry()) %>%
  dplyr::select(-idx)

# create data frame of source and destination locations
pts_i = locs %>% rename("area_i"=areaid, "x_i"=x, "y_i"=y)
pts_j = locs %>% rename("area_j"=areaid, "x_j"=x, "y_j"=y)
pts = expand.grid(pts_i$area_i, pts_j$area_j) %>%
  dplyr::rename("area_i"=1, "area_j"=2) %>%
  left_join(pts_i) %>%
  left_join(pts_j)

# haversine distance (i.e. great circle) between most populated centroids
havCalc = function(x){
  geosphere::distHaversine(p1 = c(pts$x_i[x], pts$y_i[x]), p2 = c(pts$x_j[x], pts$y_j[x]))
}
dists = sapply(1:nrow(pts), havCalc)
pts$dist_km = dists/1000
write.csv(pts, "./output/data_processed/connectivity/VietnamAll_MergedSHP2_Pairwise_EuclideanDist.csv", row.names=FALSE)
