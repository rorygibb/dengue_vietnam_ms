


# ==================== Population data from census 2009 and 2019 =======================

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

# shapefile
shp = st_read("./code/viet_dengue_districts/data/shapefiles/vietnam_districts_final.shp")
shpm = st_read("./code/viet_dengue_districts/data/shapefiles/vietnam_districts_merged.shp")





# ================= 2019 population =====================

# read population for 2019
pop = read.csv("./code/viet_dengue_districts/data/census_2019/census2019_population.csv") %>%
  dplyr::select(1, 2, 3, 4, 7, 10) %>%
  dplyr::mutate(Urban = replace(Urban, Urban == "-", 0),
                Rural = replace(Rural, Rural == "-", 0)) %>%
  dplyr::mutate(Total = as.numeric(stringr::str_replace_all(Total, stringr::fixed(" "), "")),
                Urban = as.numeric(stringr::str_replace_all(Urban, stringr::fixed(" "), "")),
                Rural = as.numeric(stringr::str_replace_all(Rural, stringr::fixed(" "), ""))) %>%
  dplyr::filter(areaid != 0) %>%
  dplyr::left_join(shp %>% dplyr::select(areaid, areanameen, areaprovin) %>% st_drop_geometry()) %>%
  dplyr::rename_all(tolower) %>%
  dplyr::select(areaprovin, areaid, areanameen, total, urban, rural) %>%
  dplyr::rename("province"=areaprovin) %>%
  dplyr::mutate(year=2019)

# combine for merged polygons
pop_tomerge = pop[ !pop$areaid %in% shpm$areaid, ]
pop_tomerge$areaid = as.character(pop_tomerge$areaid)

merge_lk = shpm[ !shpm$areaid %in% pop$areaid, ] %>%
  st_drop_geometry() %>%
  dplyr::mutate(idx = areaid) %>%
  tidyr::separate_rows(areaid, sep="_")
pop_summarised = left_join(pop_tomerge, merge_lk[ , c("areaid", "idx")]) %>%
  dplyr::select(idx, total, urban, rural) %>%
  dplyr::group_by(idx) %>%
  dplyr::summarise_all("sum", na.rm=TRUE) %>%
  dplyr::rename("areaid" = idx)
pop_summarised = left_join(pop_summarised, 
                           shpm %>% st_drop_geometry() %>% dplyr::select(areaid, areanameen, areaprovin)) %>%
  dplyr::rename("province"=areaprovin) %>%
  dplyr::mutate(year=2019)

# combine
pop = rbind(pop[ !pop$areaid %in% pop_tomerge$areaid, ], pop_summarised)

# # all are in the dengue data
# dd = read.csv("./code/viet_dengue_districts/output/dengue_timeseries/viet_dengue_alldistricts_02022021.csv")
# dd[ !dd$areaid %in% pop$areaid, ]

# save as "pop19"
pop19 = pop

# a = pop19[ pop19$areaid %in% dd$areaid, ]
# b = dd[ dd$areaid %in% a$areaid & dd$year == 2019 & dd$month == 1, c("areaid", "population_gpw") ]
# c = left_join(a, b)
# ggplot(c) + 
#   geom_point(aes(total, population_gpw)) +
#   geom_abline(intercept=0, slope=1)
# c$resid = c$population_gpw - c$total
# c$res_prop = c$resid/c$total
# c[ c$resid > 10000, ]
# c[ c$res_prop > 0.25, ]
# c[ c$province =="Ha Noi", ]




# ================== 2009 data ===================

# population for 2009
pop = read.csv("./code/viet_dengue_districts/data/census_2019/census2009_population_correct.csv") %>%
  dplyr::select(1, 3, 4, 7, 10) %>%
  dplyr::rename("urban"=4, "rural"=5) %>%
  dplyr::mutate(urban = replace(urban, urban == "- ", 0),
                rural = replace(rural, rural == "- ", 0)) %>%
  dplyr::mutate(Total = as.numeric(Total),
                urban = as.numeric(urban), 
                rural = as.numeric(rural)) %>%
  dplyr::filter(areaid != 0 & !is.na(areaid)) %>%
  dplyr::left_join(shp %>% dplyr::select(areaid, areanameen, areaprovin) %>% st_drop_geometry()) %>%
  dplyr::rename_all(tolower) %>%
  dplyr::select(areaprovin, areaid, areanameen, total, urban, rural) %>%
  dplyr::rename("province"=areaprovin) %>%
  dplyr::mutate(year=2009)

# combine for merged polygons
pop_tomerge = pop[ !pop$areaid %in% shpm$areaid, ]
pop_tomerge$areaid = as.character(pop_tomerge$areaid)

merge_lk = shpm[ !shpm$areaid %in% pop$areaid, ] %>%
  st_drop_geometry() %>%
  dplyr::mutate(idx = areaid) %>%
  tidyr::separate_rows(areaid, sep="_")
pop_summarised = left_join(pop_tomerge, merge_lk[ , c("areaid", "idx")]) %>%
  dplyr::select(idx, total, urban, rural) %>%
  dplyr::group_by(idx) %>%
  dplyr::summarise_all("sum", na.rm=TRUE) %>%
  dplyr::rename("areaid" = idx)
pop_summarised = left_join(pop_summarised, 
                           shpm %>% st_drop_geometry() %>% dplyr::select(areaid, areanameen, areaprovin)) %>%
  dplyr::rename("province"=areaprovin) %>%
  dplyr::mutate(year=2009)

# combine
pop = rbind(pop[ !pop$areaid %in% pop_tomerge$areaid, ], pop_summarised)

# save as "pop09"
pop09 = pop




# ================== identify districts that are not present in the 2009 data but are present in 2019 and merged shp ===================

disc = pop19[ !pop19$areaid %in% pop09$areaid, ]


# 1. Ha Noi: Bac Tu Liem - in 2013 Tu Liem (in 2009 census) was split into Nam and Bac Tu Liem
# Important to keep both so partition 2009 population into same ratio as in 2019 population (this is consistent with GPW)

# whole of Tu Liem in 2009
d1 = pop09[ grep("Tu Liem", pop09$areanameen), ]

# both in 2019
d1.2019 = pop19[ grep("Tu Liem", pop19$areanameen), ]

# create Bac/Nam Tu Liem for 2009, keeping same urban/rural ratio as in Tu Liem district
d1.2009 = d1.2019 %>%
  dplyr::mutate(year = 2009)
d1.2009$total[ d1.2009$areanameen == "Bac Tu Liem" ] = 230576
d1.2009$total[ d1.2009$areanameen == "Nam Tu Liem" ] = 186593
d1.2009$urban = round(d1.2009$total * (d1$urban/d1$total))
d1.2009$rural = d1.2009$total - d1.2009$urban

# add back into pop09
pop09 = rbind(
  pop09 %>% dplyr::filter(areanameen != "Nam Tu Liem"),
  d1.2009
) %>%
  dplyr::arrange(province, areaid)


# 2. Tran De in Soc Trang- created from Long Phu and My Xuyen in 2009
# 75,046 from Long Phu and 55,031 from My Xuyen so calculate based on these

d2 =  pop09[ pop09$areanameen %in% c("Tran De", "Long Phu", "My Xuyen"), ]

d2.2019 = pop19[ pop19$areanameen %in% c("Tran De", "Long Phu", "My Xuyen"), ]

# calculate based on partitioning in 2009 and assign same rural/urban split as in 2009 census
d2.2009 = d2.2019 %>% dplyr::mutate(year = 2009)
d2.2009$total[ d2.2009$areanameen == "Long Phu" ] = d2$total[ d2$areanameen == "Long Phu" ] - 75046
d2.2009$total[ d2.2009$areanameen == "My Xuyen" ] = d2$total[ d2$areanameen == "My Xuyen" ] - 55031
d2.2009$total[ d2.2009$areanameen == "Tran De" ] = 75046 + 55031
d2.2009$urban[ d2.2009$areanameen == "Long Phu" ] = round(d2.2009$total[ d2.2009$areanameen == "Long Phu" ] * (d2$urban[ 1 ]/d2$total[1]))
d2.2009$urban[ d2.2009$areanameen == "My Xuyen" ] = round(d2.2009$total[ d2.2009$areanameen == "My Xuyen" ] * (d2$urban[ 2 ]/d2$total[2]))
d2.2009$urban[ d2.2009$areanameen == "Tran De" ] = round(d2.2009$total[ d2.2009$areanameen == "Tran De" ] * mean(d2$urban/d2$total)) # mean of both for Tran De
d2.2009$rural = d2.2009$total - d2.2009$urban

# combine
pop09 = rbind(
  pop09 %>% dplyr::filter(!areaid %in% d2.2009$areaid),
  d2.2009
) %>%
  dplyr::arrange(province, areaid)




# =================== compare 2009 data (based on census) to 2009 data extracted from GPW ===============

# #gpw data
pp = read.csv("./output/data_processed/population/VietPopulation_District_MergedSHP_GPW4_interpolated.csv") %>%
  dplyr::filter(year == 2009)
# pp = left_join(pp, pop09 %>% dplyr::select(areaid, total))
# 
# # view: very strongly correlated in absolute terms
# ggplot(pp) +
#   geom_point(aes(total, population_gpw), size=1.2, col="darkred", alpha=0.75) +
#   geom_abline(i=0, s=1) + theme_bw()
# 
# # histogram of discrepancies - median 2%, mean 5%
# pp$percent_disc = (pp$population_gpw - pp$total)/pp$total
# hist(pp$percent_disc, 200)
# median(abs(pp$percent_disc), na.rm=TRUE)
# mean(abs(pp$percent_disc), na.rm=TRUE)
# pp[ pp$percent_disc>1, ]



# ============== interpolate from 2009 to 2019 for all districts ===============

popn = pop09 %>%
  dplyr::left_join(pop19 %>% dplyr::select(areaid, total, urban, rural) %>% dplyr::rename("total_2019"=total, "urban_2019"=urban, "rural_2019"=rural))

popn$total_peryear = (popn$total_2019 - popn$total)/10
popn$urb_peryear = (popn$urban_2019  - popn$urban)/10

result = popn

# for each year
for(i in 1:10){
  
  popn = popn %>%
    dplyr::mutate(year = year + 1,
                  total = total + total_peryear,
                  urban = urban + urb_peryear, 
                  rural = total - urban)
  
  result = rbind(result, popn)
}

# round
result$total = round(result$total)
result$rural = round(result$rural)
result$urban = round(result$urban)

# result %>%
#   dplyr::filter( year == 2019 ) %>%
#   ggplot() + 
#   geom_point(aes(rural, rural_2019)) + 
#   geom_abline(i=0, s=1)

# comparison of trends in HaNoi from GPW versus census data; census data is a major difference
# result %>%
#   dplyr::filter(province == "Ha Noi") %>%
#   ggplot() + 
#   geom_line(aes(year, total, group=areaid))
# read.csv("./output/data_processed/population/VietPopulation_District_MergedSHP_GPW4_interpolated.csv") %>%
#   dplyr::filter(year %in% 2009:2019) %>%
#   dplyr::filter(areaprovin == "Ha Noi") %>%
#   ggplot() + 
#   geom_line(aes(year, population_gpw, group=areaid))

# cut down to required cols
pop_proj1 = result %>%
  dplyr::select(1:7)



# ================= project forward to 2020 =================

p2020 = data.frame()

for(i in 1:length(unique(pop_proj1$areaid))){
  
  px = pop_proj1[ pop_proj1$areaid == unique(pop_proj1$areaid)[i], ] %>%
    dplyr::arrange(desc(year))
  cr_total = mean(px$total[ 1:5 ] / px$total[ 2:6 ])
  if(is.na(cr_total) | is.nan(cr_total)){ cr_total = 1 }
  cr_urban = mean(px$urban [ 1:5 ] / px$urban[ 2:6 ])
  if(is.na(cr_urban) | is.nan(cr_urban)){ cr_urban = 1 }
  resx = px[ 1, ] %>%
    dplyr::mutate(year = 2020,
                  total = round(total * cr_total),
                  urban = round(urban * cr_urban))
  if(resx$urban > resx$total){ resx$urban = resx$total }
  resx$rural = resx$total - resx$urban
  p2020 = rbind(p2020, resx)
                
}

pop_proj1 = rbind(pop_proj1, p2020) %>%
  dplyr::arrange(areaid, year)




# ================== 

# because census 1999 data are not publically available at sufficient granularity
# use population growth rates from 1998 to 2009 calculated for GPW data at district-level
# these are estimated at district level from comparison of 1999 and 2009 data (i.e. observed trends which in GPW are extrapolated to 2019)
# so should correctly reflect broad trends at district level

pp = read.csv("./output/data_processed/population/VietPopulation_District_MergedSHP_GPW4_interpolated.csv") %>%
  dplyr::filter(year %in% 1998:2009)
pp_dif = data.frame()

for(i in unique(pp$areaid)){
  pp_i = pp[ pp$areaid == i, ]
  pp_i$dif = c(pp_i$population_gpw[ 1:(nrow(pp_i)-1)] / pp_i$population_gpw[ 2:(nrow(pp_i))], NA)
  pp_dif = rbind(pp_dif, pp_i)
}

pop_proj1 = pop_proj1 %>% dplyr::filter(areanameen != "Truong Sa")
pp_dif = pp_dif %>% dplyr::filter(areanameen != "Truong Sa")

# population projection function
backwardProjectPerDistrict = function(areaidx){
  print(areaidx)
  
  distx = pop_proj1 %>%
    dplyr::filter(areaid == areaidx)
  
  projrates_x = pp_dif %>%
    dplyr::filter(areaid == areaidx)
  
  # mean yearly change in proportion pop in urban calculated over following 3 years
  propurb = distx$urban / distx$total
  propurb = mean(propurb[ 1:4 ]/ propurb[2:5])
  if(is.nan(propurb)){ propurb = 0 }
  
  for(i in 2008:1998){
    
    year_i = distx[ distx$year == i + 1, ]
    propurb_i = year_i$urban/year_i$total * propurb # change in propurb
    if(propurb_i > 1){ propurb_i = 1 } # if reaches full
    year_i = year_i %>%
      dplyr::mutate(year = i,
                    total = total * projrates_x$dif[ projrates_x$year == i ],
                    urban = total * propurb_i,
                    rural = total - urban)
    distx = rbind(year_i, distx)
  }
  # ggplot() + geom_line(data = distx, aes(year, total)) + 
  #   geom_point(data = projrates_x, aes(year, population_gpw))
  
  # return
  return(distx)
}

# create pop projection
pop_proj2 = do.call(rbind.data.frame, lapply(unique(pop_proj1$areaid), backwardProjectPerDistrict))



# # =============== test viz ===============
#
# pop_proj2 %>%
#   dplyr::filter(province == "TP. Ho Chi Minh") %>%
#   ggplot() +
#   geom_line(aes(year, total, group=areaid))
# # # 
# pop_proj2 %>%
#   dplyr::filter(province == "Ha Noi") %>%
#   ggplot() +
#   geom_line(aes(year, urban/total, group=areaid))
# # 
# pop_proj2 %>%
#   dplyr::filter(urban > total)


# ================ save ====================

write.csv(pop_proj2, "./code/viet_dengue_districts/output/covariates/Population_Census_2009_2019.csv", row.names=FALSE)

