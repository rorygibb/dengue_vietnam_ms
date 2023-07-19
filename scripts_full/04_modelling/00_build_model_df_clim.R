


# ====================================================================================================
# ======================= BUILDS OBJECTS FOR FITTING AND EVALUATING MODELS ===========================
# ====================================================================================================

# This script is called as source in a modelling script where various fields can be specified

library(dplyr); library(raster); library(rgdal); library(sf)
library(stringr); library(ggplot2); library(lubridate)
library(magrittr); library(INLA); library(spdep)
source("00_plot_themes.R")
source("00_inla_setup_functions_r4.R")

# fields to specify
# project name
# number of bins to group climatic data
# region to subset to (N/S/C)
# plot_graph: plot neighbourhood matrix for inla model?
# province_case_threshold: only keep provinces with threshold of >n cases
# region2: subregion to subset to

if(!exists("projname")){
  projname = "temp"
}

if(!exists("n_clim_bins")){
  n_clim_bins = 40
}

if(!exists("region")){
  region = "all"
}

if(!exists("region2")){
  region2 = NA
} 

if(!exists("plot_graph")){
  plot_graph = TRUE
}

if(!exists("province_case_threshold")){
  province_case_threshold = NA
} else{
  province_case_threshold = province_case_threshold
}



# ============== Set up project file name and output locations ==============

# specify project name: all outputs will be saved in a directory of this name
projname = projname

# create folder structure for saving outputs
save_dir = paste(c("./output/model_outputs/", projname, "/"), collapse="")
if(!dir.exists(save_dir)){ 
  dir.create(save_dir)
  dir.create(paste(save_dir, "model_output/", sep="")) 
  dir.create(paste(save_dir, "errors/", sep="")) 
  dir.create(paste(save_dir, "models/", sep="")) 
  dir.create(paste(save_dir, "fitmetrics/", sep="")) 
}


# ================= Build dengue dataset =================

# districts to be excluded (offshore)
offshore_areas = c(70154, 70339, 70273, 70355, 70698)

# shapefiles: district, province, and vietnam wide
shp = st_read("./data/shapefiles/dengue_districts_shapefile.shp") %>%
  dplyr::filter(!areaid %in% offshore_areas)
shp_prov = st_read("./data/shapefiles/provinces.shp")
st_crs(shp_prov) = st_crs(shp)
shp_prov = st_crop(shp_prov, shp)
shp_vt = st_read("./data/shapefiles/gadm36_VNM_0.shp") %>% st_crop(shp)

# dengue, regions, climate, landuse, connectivity data
dd = read.csv("./output/model_data/ModelData_Dengue_VietAll.csv") %>%
  dplyr::left_join(read.csv("./output/model_data/ModelData_ClimaticRegions.csv")) %>%
  dplyr::left_join(read.csv("./output/model_data/ModelData_SocioEcologicalCovar_VietAll.csv")) %>%
  #dplyr::left_join(read.csv("./output/model_data/ModelData_FlightsMonthly_VietAll.csv")) %>%
  dplyr::left_join(read.csv("./output/model_data/ModelData_ClimateLags_VietAll.csv")) %>%
  # dplyr::left_join(read.csv("./output/model_data/ModelData_ClimateAnomalies_VietAll.csv") %>% 
  #                    dplyr::filter(lubridate::year(as.Date(date))>= 1998)) %>%
  dplyr::mutate(date = as.Date(date))

# exclude all years prior to useable timeseries (because some areas have sligthly shorter timeseries)
exclude_useable = function(x){
  px = dd[ dd$province == x, ]
  px = px[ px$year >= px$year_useable_from[1], ]
}
dd = do.call(rbind.data.frame, lapply(unique(dd$province), exclude_useable))




# =============== group predictors for fitting nonlinear effects ==================

# group predictors for fitting nonlinear effects 
# do this before subsetting to geographical subregions to ensure consistency for prediction/projection

print("Grouping climate predictors")

nbins = n_clim_bins
dx = dd[ , grep("tmean|tmin|tmax|precip|spei", names(dd))]
names(dx) = paste(names(dx), "_g", sep="")

groupCols = function(x){
  x = dx[ , x, drop=FALSE ]
  x[ , 1] = inla.group(x[ , 1], n=nbins)
  x
}
dx = do.call(cbind.data.frame, lapply(1:ncol(dx), groupCols))
dd = cbind(dd, dx)




# ================ subset to specified region(s), if specified =================

if(region != "all"){
  if(!region %in% unique(dd$region3)){
    print("Region not recognised: defaulting to 'all'")
  }
  else{
    dd = dd[ dd$region3 %in% c(region), ]
    shp = shp[ shp$areaid %in% dd$areaid, ]
    shp_prov = shp_prov[ shp_prov$provincena %in% dd$province, ]
  }
}


if(!is.na(region2)){
  if(!region2 %in% unique(dd$region2)){
    print("Region not recognised: defaulting to 'all'")
  }
  else{
    dd = dd[ dd$region2 %in% c(region2), ]
    shp = shp[ shp$areaid %in% dd$areaid, ]
    shp_prov = shp_prov[ shp_prov$provincena %in% dd$province, ]
  }
}



# ============== remove all provinces with cases < specified threshold, if specified ==============

if( ! is.na(province_case_threshold) ){
  
  if(! is.numeric(province_case_threshold)){ 
    
    print("Province case threshold not specifying; defaulting to no threshold")
    
  } else{
      
    provs = dd %>%
      dplyr::group_by(province) %>%
      dplyr::summarise(cases = sum(cases, na.rm=TRUE)) %>%
      dplyr::filter(cases >= province_case_threshold)
    dd = dd[ dd$province %in% provs$province, ]
    
    }
}


# ================ set up spatial neighbourhood matrix for CAR model ===================

# subset shapefile and dengue cases to focal district
shpf = shp[ shp$areaid %in% dd$areaid, ]

# create neighbourhood matrix for CAR spatial
# firstly create lookup refs for polygon ids, create neighbourhood matrix, then add lookups into dataframe
id_ref = data.frame(areaid = shpf$areaid, polyid = 1:nrow(shpf))
district.nb = spdep::poly2nb(sf::as_Spatial(shpf), row.names=id_ref$polyid)
dd = left_join(dd, id_ref)

# manually add neighbours for offshore districts (locale where ferries depart/arrive)
# 1. Ly Son (areaid = 70382) - 30 minute ferry from Sa Ky port, Binh Son district (70384), Quang Ngai
# 2. Phu Quy (areaid = 70596) - 3 hour ferry from Phan Thiet, Binh Thuan (70540)
# 3. Phu Quoc (70666) - 1 hour ferry from Ha Tien (70620), Kien Giang
# 4. Kien Hai (70671) - 1 hr from Rach Gia (70652) tourism wharf
# 5. Con Dao (70711) - 3.5 hours from Vung Tau (70615)
# if(70382 %in% id_ref$areaid){ district.nb[[ id_ref$polyid[ id_ref$areaid == 70382] ]] <- id_ref$polyid[ id_ref$areaid == 70384 ] }
# if(70596 %in% id_ref$areaid){ district.nb[[ id_ref$polyid[ id_ref$areaid == 70596] ]] <- id_ref$polyid[ id_ref$areaid == 70540 ] }
# if(70666 %in% id_ref$areaid){ district.nb[[ id_ref$polyid[ id_ref$areaid == 70666] ]] <- id_ref$polyid[ id_ref$areaid == 70620 ] }
# if(70671 %in% id_ref$areaid){ district.nb[[ id_ref$polyid[ id_ref$areaid == 70671] ]] <- id_ref$polyid[ id_ref$areaid == 70652 ] }
# if(70711 %in% id_ref$areaid){ district.nb[[ id_ref$polyid[ id_ref$areaid == 70711] ]] <- id_ref$polyid[ id_ref$areaid == 70615 ] }

# save neighbourhood matrix with focal district (if not already existing)
nbmatrix_name = paste(save_dir, "adjmatrix_", region, "_",  paste(tolower(projname), collapse="_"), sep="") 
nb2INLA(nbmatrix_name, district.nb)

# plot neighbourhood matrix if specified
if(plot_graph){
  xy = as.data.frame(rgeos::gCentroid(as_Spatial(shpf), byid=TRUE))
  plot(shpf$geometry)
  plot(district.nb, coords = cbind(xy$x, xy$y), col="red", add=T, cex=0.5)
}


# =============== setup covariates, transform, scale and set grouping factors ====================

# defining offset as log population (in hundreds of thousands)
# so model is estimating incidence per 100,000 inhabitants
# along with defining fairly tight priors on intercept precision, this deals with numerical challenges of fitting when estimating incidence per inhabitant
#dd$logpop = log(dd$population_gpw/100000) 
dd$logpop = log(dd$population_census/100000)

# key covariate
dd$logpopdens = log(dd$popdens_census)

# replicate variables for grouping
areaidx = factor(dd$areaid, levels=unique(dd$areaid)[ order(unique(dd$areaid)) ], ordered=TRUE)
dd$areaidx = as.integer(areaidx)
dd$areaidy = as.integer(areaidx)
dd$yearx = as.integer(as.factor(dd$year))
dd$provincex = as.integer(factor(dd$province, levels=unique(dd$province)[ order(unique(dd$province)) ], ordered=TRUE))
dd$provincey = dd$provincex

# different levels of region 
dd$regionx = as.integer(as.factor(dd$region1))
dd$regiony = as.integer(as.factor(dd$region2))
dd$regionz = as.integer(as.factor(dd$region3))
dd$regionf = as.integer(as.factor(dd$region4))

