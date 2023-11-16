

# ====================== Projections of temperature-driven change in dengue risk over long-term ========================

# This script reads in the fitted functions from the full model and long-term ERA5-Land Tmean data (both stored in the GitHib) 
# Analyses present-day spatiotemporal risk patterns and relative risk changes between historical reference period and present day

library(dplyr)
library(sf)
library(ggplot2)
library(vroom)
library(lubridate)
library(raster)

PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

source("./scripts/04_modelling/00_inla_setup_functions_r4.R")
source("./scripts/04_modelling/00_plot_themes.R")

# directory where model objects are stored (large files so external)
models_dir = "C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/output/model_outputs/results_23/"



# ================= shapefile and dengue data =================

# districts to be excluded (offshore)
# only ones with substantial dengue cases are in Kien Giang; could be worth exploring including them but thisis sufficient for now
offshore_areas = c(70154, 70339, 70273, 70355, 70698)

# districts shapefile for Vietnam
shp = st_read("./data/shapefiles/dengue_districts_shapefile.shp") %>%
  dplyr::filter(!areaid %in% offshore_areas)

# provinces shapefile with regions info, and vietnam outline for mapping
shp_prov = st_read("./data/shapefiles/provinces.shp")
st_crs(shp_prov) = st_crs(shp)
shp_prov = st_crop(shp_prov, shp)
rr = read.csv("./data/shapefiles/regions_lookup.csv")
shp_prov = left_join(shp_prov, rr[ , c("provincena", "region")])
shp_vt = st_read("./data/shapefiles/vt_national.shp") %>%
  st_crop(shp)

# dengue, regions, climate, landuse, connectivity data
dd = read.csv("./output/model_data/ModelData_Dengue_VietAll.csv") %>%
  dplyr::left_join(read.csv("./output/model_data/ModelData_ClimaticRegions.csv")) %>%
  dplyr::left_join(read.csv("./output/model_data/ModelData_SocioEcologicalCovar_VietAll.csv")) %>%
  dplyr::left_join(vroom::vroom("./output/model_data/ModelData_ClimateLags_VietAll.csv.gz", col_types = list(areaid=vroom::col_character(), date=vroom::col_character()))) %>%
  dplyr::mutate(date = as.Date(date))

# provincex
dd$provincex = as.integer(as.factor(dd$province))

# polyid indicator
shpf = shp[ shp$areaid %in% dd$areaid, ]
id_ref = data.frame(areaid = shpf$areaid, polyid = 1:nrow(shpf))
dd = left_join(dd, id_ref)

# add polygon area, lat lon, region information
#shp$area_km2 = as.vector(st_area(shp) / 10^6)
shp = cbind(shp, as.data.frame(st_coordinates(st_centroid(shp))) %>% dplyr::rename("longitude"=1, "latitude"=2))
dd = left_join(dd, shp[ , c("areaid", "latitude", "longitude")] %>% st_drop_geometry())
shp = left_join(shp, dd[ !duplicated(dd$areaid), c("areaid", "region1", "region2", "region3") ])




# ================ data and model ====================

# temperature effect from full model
rf = read.csv("./output/model_outputs/fullmodel/fitted_params/fullmodel_fittedclimatefunctions_rw2.csv") %>%
  dplyr::filter(effect == "tmean_1m_g")

# long term tmean from era5-land
era = raster::brick("./data/climate/era5-land/tmean_monthly/tmean_month_viet_19502020.nc")

# scaling quirk using terra and exactextract for NCDF; fix with conversion from raster to terra rather than reading directly (longwinded)
era_terra = terra::rast()
for(i in 1:nlayers(era)){
  rr = raster(era[[1]])
  values(rr) = values(era[[i]])
  rr = terra::rast(rr)
  names(rr) = names(era[[i]])
  era_terra = c(era_terra, rr)
}
temp = exactextractr::exact_extract(era_terra, shp, fun="mean")

temp = cbind(shp[ , "areaid"] %>% st_drop_geometry(), temp) %>%
  reshape2::melt(id.vars = "areaid") %>%
  dplyr::mutate(date = substr(variable, 7, 30),
                date = as.Date(date, format="%Y.%m.%d"),
                tmean = value - 273.15) %>%
  dplyr::select(-variable, -value)

# rollmean to align with same processing of climate data into models
clim = temp %>%
  dplyr::select(areaid, date, tmean) %>%
  dplyr::arrange(areaid, date) %>%
  dplyr::group_by(areaid) %>%
  dplyr::mutate(
    tmean = data.table::frollmean(tmean, 2, align="right")) %>%
  dplyr::ungroup()

# generate 1m lag and round to 3dp
clim = clim %>%
  dplyr::mutate(date = date %m+% months(1)) %>%
  dplyr::rename("tmean_1m" = tmean) %>%
  dplyr::mutate(tmean_1m = round(tmean_1m, 3))




# ----------------- set up nonlinear temperature variable for projection -----------------

# temperature fitted function (using rw2 so evaluated at specific grouping levels)
tvars = data.frame(
  tmean_1m = round(seq(min(clim$tmean_1m, na.rm=TRUE), max(clim$tmean_1m, na.rm=TRUE), by=0.001), 3)
) %>%
  dplyr::left_join(
    rf %>% 
      dplyr::filter(effect == "tmean_1m_g") %>%
      dplyr::select(value, mean) %>% 
      dplyr::rename("tmean_1m"=value) %>%
      dplyr::mutate(tmean_1m = round(tmean_1m, 3))
  )

# interpolate risk between levels of Tmean function to 3dp to align with temperature data
interp = function(df){
  
  # name var
  varname = names(df)[1]
  df = df %>% dplyr::rename("varname"=1)
  
  mm = unique(df$mean); mm = mm[ !is.na(mm)]
  result = data.frame()
  for(i in 2:length(mm)){
    start = mm[i-1]
    end = mm[i]
    in_between = df[ which(df$mean==start) : which(df$mean == end), ]
    in_between$mean = seq(start, end, length.out=nrow(in_between))
    result = rbind(result, in_between)
  }
  result = result %>% distinct()
  
  # anything over the end
  fin = df %>% filter(varname > max(result$varname))
  if(nrow(fin) > 0){
    
    multip_factor = mean(result$mean[ (nrow(result)-100) : nrow(result)] / result$mean[ (nrow(result)-101) : (nrow(result)-1) ])
    for(i in 1:nrow(fin)){
      if(i == 1){ fin$mean[i] = result$mean[ nrow(result) ] * multip_factor 
      } else{
        fin$mean[i] = fin$mean[i-1] * multip_factor
      }
    }
    result = rbind(result, fin)
  }
  
  # anything beforehand
  pre = df %>% filter(varname < min(result$varname))
  if(nrow(pre) > 0){
    
    multip_factor = mean(result$mean[ 1:100 ] / result$mean[ 2:101 ])
    for(i in nrow(pre):1){
      if(i == nrow(pre)){ pre$mean[i] = result$mean[ 1 ] * multip_factor 
      } else{
        pre$mean[i] = pre$mean[i+1] * multip_factor
      }
    }
    result = rbind(pre, result)
  }
  
  names(result)[1] = varname
  return(result)
}

tvars = interp(tvars)

# combine into dataframe
clim = clim %>%
  left_join(tvars %>% dplyr::rename("tmean_effect"=mean))

# fixes for coastal districts
clim = rbind(clim[ clim$areaid != 70361, ], clim[ clim$areaid == 70363, ] %>% dplyr::mutate(areaid = "70361"))
clim = rbind(clim[ clim$areaid != 70382, ], clim[ clim$areaid == 70384, ] %>% dplyr::mutate(areaid = "70382"))
clim = rbind(clim[ clim$areaid != 70596, ], clim[ clim$areaid == 70517, ] %>% dplyr::mutate(areaid = "70596"))
clim = rbind(clim[ clim$areaid != 70671, ], clim[ clim$areaid == 70680, ] %>% dplyr::mutate(areaid = "70671"))
clim = rbind(clim[ clim$areaid != 70711, ], clim[ clim$areaid == 70699, ] %>% dplyr::mutate(areaid = "70711"))



# ================ change in risk between baseline period (1950-1970) and present (2001-2020) ======================

# set historical and present-day periods
clim_for_test = clim %>%
  dplyr::mutate(year = lubridate::year(date)) %>%
  dplyr::filter(year %in% c(1951:1970, 2001:2020)) %>%
  dplyr::mutate(year_grp = ifelse(year < 2000, 1, 2),
                month = lubridate::month(date)) %>%
  dplyr::select(areaid, year, year_grp, date, month, tmean_1m, tmean_effect)

# dengue month
dm = data.frame(month=c(5:12, 1:4), monthdengue = c(1:12), month_name = c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"))
clim_for_test = left_join(clim_for_test, dm)

# monthly linear models; compare between baseline and test period
mnx = clim_for_test %>%
  dplyr::filter(!is.na(tmean_effect)) %>%
  dplyr::group_by(areaid, monthdengue) %>%
  dplyr::summarise(month_name = head(month_name, 1), 
                   tmean_fac_slope = coef(summary(lm(tmean_effect ~ factor(year_grp))))[ 2, 1],
                   tmean_fac_p = coef(summary(lm(tmean_effect ~ factor(year_grp))))[ 2, 4],
                   mean_risk = mean(tmean_effect[ year > 2000 ]))

# mapping mean risk in present day
map1 = shp %>% 
  dplyr::full_join(mnx) %>%
  dplyr::filter(!is.na(monthdengue)) %>%
  left_join(dm) %>%
  dplyr::mutate(month_name = factor(month_name, levels=dm$month_name, ordered=TRUE)) %>%
  ggplot() + 
  geom_sf(aes(fill=exp(mean_risk)), color=NA) + 
  maptheme + 
  facet_wrap(~month_name, nrow=1) + 
  ggtitle("Seasonal temperature-driven dengue risk (monthly mean, 2001-2020)") +
  scale_fill_gradientn(colors=rev(viridisLite::mako(200)), na.value="white", name="Relative\nrisk") +
  geom_sf(data = shp_vt, fill=NA, col="grey20", size=0.25) +
  theme(strip.text = element_text(size=11.5), plot.title=element_text(size=12.5), legend.title = element_text(size=12),
        legend.text = element_text(size=11)) +
  labs(tag = "a") +
  theme(plot.tag = element_text(size=22, face="bold"),
        plot.tag.position = c(.02, .95)) #+

# only show slopes where p < 0.05
lims = (exp(mnx$tmean_fac_slope)-1)*100
lims = c(-max(abs(lims), na.rm=TRUE), max(abs(lims), na.rm=TRUE))

# mapping geography of change
map2 = shp %>% 
  dplyr::full_join(mnx) %>%
  dplyr::filter(!is.na(monthdengue)) %>%
  dplyr::mutate(tmean_fac_slope = replace(tmean_fac_slope, tmean_fac_p > 0.05, NA)) %>%
  left_join(dm) %>%
  dplyr::mutate(month_name = factor(month_name, levels=dm$month_name, ordered=TRUE)) %>%
  ggplot() + 
  geom_sf(aes(fill=(exp(tmean_fac_slope)-1)*100), color=NA) + 
  maptheme + 
  facet_wrap(~month_name, nrow=1) + 
  ggtitle("Change in dengue risk from historical reference (1951-1970) to present (2001-2020)") + 
  geom_sf(data = shp_vt, fill=NA, col="grey20", size=0.25) +
  scale_fill_gradientn(colors=rev(pals::brewer.rdylbu(200)), na.value="white", limits=lims, name="Risk\nchange\n(%)", breaks=c(-50, -25, 0, 25, 50)) +
  #scale_fill_gradient2(high=scales::muted("blue"), low=scales::muted("red"), midpoint=0, name="Annual\nrisk\nchange\n(%)", na.value="grey95") +
  theme(strip.text = element_blank(), plot.title=element_text(size=12.5), legend.title = element_text(size=12), legend.text = element_text(size=11)) +
  labs(tag = "b") +
  theme(plot.tag = element_text(size=22, face="bold"),
        plot.tag.position = c(.02, .95)) #+




# ---------- visualise graphs of risk over time for example cities --------------

locs = shp %>% dplyr::filter(areanameen %in% c("Hoan Kiem", "Nha Trang", "Quan 1", "Buon Ma Thuot", "Binh Thuy"))

clim_for_test = clim %>%
  dplyr::filter(areaid %in% locs$areaid) %>%
  dplyr::mutate(year = lubridate::year(date)) %>%
  dplyr::filter(year %in% c(1951:1970, 2001:2020)) %>%
  dplyr::mutate(year_grp = ifelse(year < 2000, 1, 2),
                month = lubridate::month(date)) %>%
  dplyr::select(areaid, year, year_grp, date, month, tmean_1m, tmean_effect)

dm = data.frame(month=c(5:12, 1:4), monthdengue = c(1:12), month_name = c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"))
clim_for_test = left_join(clim_for_test, dm)

mnx = clim_for_test %>%
  dplyr::filter(!is.na(tmean_effect)) %>%
  dplyr::group_by(areaid, monthdengue, year_grp) %>%
  dplyr::summarise(month_name = head(month_name, 1), 
                   mean_tmean = mean(tmean_1m),
                   se_tmean = sd(tmean_1m) / sqrt(length(tmean_1m)),
                   mean_risk = mean(tmean_effect),
                   sd_risk = sd(tmean_effect),
                   se_risk = sd_risk / sqrt(length(tmean_effect))) %>%
  dplyr::mutate(
    upper = mean_risk + (1.96*se_risk),
    lower = mean_risk - (1.96*se_risk),
    upper_tmean = mean_tmean + (1.96*se_tmean),
    lower_tmean = mean_tmean - (1.96*se_tmean)
  ) %>%
  dplyr::mutate(variable = ifelse(year_grp == 1, "1951-1970", "2001-2020")) %>%
  dplyr::mutate(
    loc = "Buon Ma Thuot (Dak Lak)",
    loc = replace(loc, areaid == "70165", "Ha Noi"),
    loc = replace(loc, areaid == "70468", "Nha Trang (Khanh Hoa)"),
    loc = replace(loc, areaid == "70563", "TP. Ho Chi Minh"),
    loc = replace(loc, areaid == "70648", "Can Tho"),
    loc = factor(loc, levels = c("Ha Noi", 
                                 "Buon Ma Thuot (Dak Lak)", 
                                 "Nha Trang (Khanh Hoa)", 
                                 "TP. Ho Chi Minh", 
                                 "Can Tho"), ordered=TRUE))

# plot
p3 = clim_for_test %>% 
  dplyr::mutate(variable = ifelse(year_grp == 1, "1951-1970", "2001-2020")) %>%
  dplyr::mutate(
    loc = "Buon Ma Thuot (Dak Lak)",
    loc = replace(loc, areaid == "70165", "Ha Noi"),
    loc = replace(loc, areaid == "70468", "Nha Trang (Khanh Hoa)"),
    loc = replace(loc, areaid == "70563", "TP. Ho Chi Minh"),
    loc = replace(loc, areaid == "70648", "Can Tho"),
    loc = factor(loc, levels = c("Ha Noi", 
                                 "Buon Ma Thuot (Dak Lak)", 
                                 "Nha Trang (Khanh Hoa)", 
                                 "TP. Ho Chi Minh", 
                                 "Can Tho"), ordered=TRUE)) %>%
  ggplot() + 
  geom_line(aes(monthdengue, exp(tmean_effect), group=year, color=variable), alpha=0.2, show.legend = FALSE, size=0.2) + 
  geom_point(data = mnx, aes(monthdengue, exp(mean_risk), group=variable, color=variable), position=position_dodge(width=0.5), size=1.3, alpha=1) +  
  geom_line(data = mnx, aes(monthdengue, exp(mean_risk), group=variable, color=variable), position=position_dodge(width=0.5), size=0.5, alpha=1) +
  geom_linerange(data = mnx, aes(monthdengue, ymin=exp(lower), ymax=exp(upper), group=variable, color=variable), position=position_dodge(width=0.5), size=0.6, alpha=1, show.legend = FALSE) +
  facet_wrap(~areaid) +
  theme_classic() +
  lemon::facet_rep_wrap(~loc, nrow=1) +
  scale_x_continuous(breaks=1:12) +
  geom_hline(yintercept=1, lty=2) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=11), 
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=9.25), 
        axis.title.x = element_text(size=11.5),
        axis.title.y = element_text(size=11),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11)) +
  scale_color_viridis_d(begin=0.05, end=0.65, direction=-1, name="Time period") +
  xlab("Dengue month (May to April)") +
  ylab("Dengue risk\n(mean ± 95% CI)") +
  labs(tag = "c") +
  theme(plot.tag = element_text(size=22, face="bold"),
        plot.tag.position = c(.02, .95)) #+
 # coord_cartesian(xlim = c(50, 350), ylim = c(10, 35), clip = "off")

# p3 = gridExtra::grid.arrange(p3)
# p3 = ggpubr::as_ggplot(p3) +
#   cowplot::draw_plot_label(label = c("c"),
#                            fontface = "bold", size = 20,
#                            x = 0.01, y = 0.98)


# ----------- save ------------

# combine with maps and save
library(patchwork)
comb = map1 + map2 + p3 #+ plot_annotation(tag_levels = 'a') + theme(plot.tag = element_text(size = 20, face = "bold"))
comb = comb + plot_layout(guides = "collect", heights=c(1.3, 1.3, 1)) 

ggsave(comb, file="./output/figures/Figure5_TemperatureProjections_LongTerm_full.jpg", device="jpg", units="in", width=18, height=10, dpi=600, scale=0.65)

ggsave(comb, file="./output/figures/Figure5_TemperatureProjections_LongTerm_full.pdf", device="pdf", units="in", width=18, height=10, scale=0.65)




# -------------- patterns of temperature change ----------------

# # monthly linear models; compare between baseline and test period
# clim_for_test = clim %>%
#   dplyr::mutate(year = lubridate::year(date)) %>%
#   dplyr::filter(year %in% c(1951:1970, 2001:2020)) %>%
#   dplyr::mutate(year_grp = ifelse(year < 2000, 1, 2),
#                 month = lubridate::month(date)) %>%
#   dplyr::select(areaid, year, year_grp, date, month, tmean_1m, tmean_effect)
# 
# # dengue month
# dm = data.frame(month=c(5:12, 1:4), monthdengue = c(1:12), month_name = c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"))
# clim_for_test = left_join(clim_for_test, dm)
# 
# mnx = clim_for_test %>%
#   dplyr::filter(!is.na(tmean_effect)) %>%
#   dplyr::group_by(areaid, monthdengue) %>%
#   dplyr::summarise(month_name = head(month_name, 1), 
#                    tmean_fac_slope = coef(summary(lm(tmean_1m ~ factor(year_grp))))[ 2, 1],
#                    tmean_fac_p = coef(summary(lm(tmean_1m ~ factor(year_grp))))[ 2, 4],
#                    mean_risk = mean(tmean_effect[ year > 2000 ]))
# 
# # only show slopes where p < 0.05
# lims = mnx$tmean_fac_slope
# lims = c(-max(abs(lims), na.rm=TRUE), max(abs(lims), na.rm=TRUE))
# 
# # mapping geography of change
# map_temp = shp %>% 
#   dplyr::full_join(mnx) %>%
#   dplyr::filter(!is.na(monthdengue)) %>%
#   dplyr::mutate(tmean_fac_slope = replace(tmean_fac_slope, tmean_fac_p > 0.05, NA)) %>%
#   left_join(dm) %>%
#   dplyr::mutate(month_name = factor(month_name, levels=dm$month_name, ordered=TRUE)) %>%
#   ggplot() + 
#   geom_sf(aes(fill=tmean_fac_slope), color=NA) + 
#   maptheme + 
#   facet_wrap(~month_name, nrow=1) + 
#   ggtitle("Change in Tmean (1-month lag) from historical reference (1951-1970) to present (2001-2020)") + 
#   geom_sf(data = shp_vt, fill=NA, col="grey20", size=0.25) +
#   scale_fill_gradientn(colors=rev(pals::brewer.rdylbu(200)), na.value="white", limits=lims, name="Tmean\nchange\n(C)", breaks=c(-1.5, -1, -0.50, 0, 0.5, 1, 1.5)) +
#   #scale_fill_gradient2(high=scales::muted("blue"), low=scales::muted("red"), midpoint=0, name="Annual\nrisk\nchange\n(%)", na.value="grey95") +
#   theme(strip.text = element_text(size=12), plot.title=element_text(size=12.5), legend.title = element_text(size=12), legend.text = element_text(size=11))

# plot for cities
locs = shp %>% dplyr::filter(areanameen %in% c("Hoan Kiem", "Nha Trang", "Quan 1", "Buon Ma Thuot", "Binh Thuy"))

clim_for_test = clim %>%
  dplyr::filter(areaid %in% locs$areaid) %>%
  dplyr::mutate(year = lubridate::year(date)) %>%
  dplyr::filter(year < 2021) %>%
  dplyr::mutate(decade = paste(substr(year, 1, 3), 0, sep="")) %>%
  dplyr::mutate(year_grp = ifelse(year < 2000, 1, 2),
                month = lubridate::month(date)) %>%
  dplyr::select(areaid, year, decade, year_grp, date, month, tmean_1m, tmean_effect)

dm = data.frame(month=c(5:12, 1:4), monthdengue = c(1:12), month_name = c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"))
clim_for_test = left_join(clim_for_test, dm)

p4 = clim_for_test %>% 
  dplyr::mutate(
    loc = "Buon Ma Thuot (Dak Lak)",
    loc = replace(loc, areaid == "70165", "Ha Noi"),
    loc = replace(loc, areaid == "70468", "Nha Trang (Khanh Hoa)"),
    loc = replace(loc, areaid == "70563", "TP. Ho Chi Minh"),
    loc = replace(loc, areaid == "70648", "Can Tho"),
    loc = factor(loc, levels = c("Ha Noi", 
                                 "Buon Ma Thuot (Dak Lak)", 
                                 "Nha Trang (Khanh Hoa)", 
                                 "TP. Ho Chi Minh", 
                                 "Can Tho"), ordered=TRUE)) %>%
  ggplot() + 
  geom_line(aes(monthdengue, tmean_1m, group=year, color=year), alpha=0.3, size=0.4) + 
  # geom_point(data = mnx, aes(monthdengue, mean_tmean, group=variable, color=variable), position=position_dodge(width=0.5), size=1, alpha=1) +
  # geom_line(data = mnx, aes(monthdengue, mean_tmean, group=variable, color=variable), position=position_dodge(width=0.5), size=0.5, alpha=1) +
  # geom_linerange(data = mnx, aes(monthdengue, ymin=lower_tmean, ymax=upper_tmean, group=variable, color=variable), position=position_dodge(width=0.5), size=0.6, alpha=1, show.legend = FALSE) +
  theme_classic() +
  lemon::facet_rep_wrap(~loc, nrow=3, scales="free_y") +
  scale_x_continuous(breaks=1:12) +
  geom_hline(yintercept = rf$value[ rf$mean == max(rf$mean) ], lty=2) + # dotted line for peak risk
  theme(strip.background = element_blank(),
        strip.text = element_text(size=11), 
        axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=9.25), 
        axis.title.x = element_text(size=11.5),
        axis.title.y = element_text(size=11),
        legend.title = element_text(size=12),
        legend.position = c(0.75, 0.16),
        legend.text = element_text(size=11)) +
  scale_color_viridis_c(begin=0.05, end=0.85, direction=-1, name="Year") +
  xlab("Dengue month (May to April)") +
  ylab("Tmean (1-month lag)") 

ggsave(p4, file="./output/figures/SuppFigure_TmeanProjections_MajorCities.jpg", device="jpg", units="in", width=8, height=11, dpi=600, scale=0.6)







# ---------- Supp Figure: risk over time for high and low change areas --------------

clim_for_test = clim %>%
  dplyr::mutate(year = lubridate::year(date)) %>%
  dplyr::filter(year %in% c(1951:1970, 2001:2020)) %>%
  dplyr::mutate(year_grp = ifelse(year < 2000, 1, 2),
                month = lubridate::month(date)) %>%
  dplyr::select(areaid, year, year_grp, date, month, tmean_1m, tmean_effect)

dm = data.frame(month=c(5:12, 1:4), monthdengue = c(1:12), month_name = c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"))
clim_for_test = left_join(clim_for_test, dm)

high_low = clim_for_test %>%
  dplyr::filter(!is.na(tmean_effect)) %>%
  dplyr::group_by(areaid, monthdengue) %>%
  dplyr::summarise(month_name = head(month_name, 1), 
                   tmean_fac_slope = coef(summary(lm(tmean_effect ~ factor(year_grp))))[ 2, 1],
                   tmean_fac_p = coef(summary(lm(tmean_effect ~ factor(year_grp))))[ 2, 4],
                   mean_risk = mean(tmean_effect[ year > 2000 ])) 
high = high_low %>% dplyr::arrange(desc(tmean_fac_slope)) 
high = unique(high$areaid[1:30])[ 1:10 ]
low = high_low %>% dplyr::arrange(tmean_fac_slope) 
low = unique(low$areaid[1:30])[ 1:10 ]

mnx = clim_for_test %>%
  dplyr::filter(areaid %in% c(high, low)) %>%
  dplyr::filter(!is.na(tmean_effect)) %>%
  dplyr::group_by(areaid, monthdengue, year_grp) %>%
  dplyr::summarise(month_name = head(month_name, 1), 
                   mean_tmean = mean(tmean_1m),
                   se_tmean = sd(tmean_1m) / sqrt(length(tmean_1m)),
                   mean_risk = mean(tmean_effect),
                   sd_risk = sd(tmean_effect),
                   se_risk = sd_risk / sqrt(length(tmean_effect))) %>%
  dplyr::mutate(
    upper = mean_risk + (1.96*se_risk),
    lower = mean_risk - (1.96*se_risk),
    upper_tmean = mean_tmean + (1.96*se_tmean),
    lower_tmean = mean_tmean - (1.96*se_tmean)
    # upper = mean_risk + sd_risk,
    # lower = mean_risk - sd_risk
  ) %>%
  dplyr::mutate(variable = ifelse(year_grp == 1, "1951-1970", "2001-2020")) 

clim_for_test = clim_for_test %>% dplyr::filter(areaid %in% c(high, low))
mnx = left_join(mnx, shp %>% dplyr::select(areaid, areanameen, region1, areaprovin) %>% sf::st_drop_geometry())
mnx$tag = paste(mnx$areanameen, "\n(", mnx$areaprovin, ")\n", mnx$region1, sep="")
clim_for_test = left_join(clim_for_test, shp %>% dplyr::select(areaid, areanameen, region1, areaprovin) %>% sf::st_drop_geometry())
clim_for_test$tag = paste(clim_for_test$areanameen, "\n(", clim_for_test$areaprovin, ")\n", clim_for_test$region1, sep="")


# plot
ps1 = clim_for_test %>% 
  dplyr::filter(areaid %in% high) %>%
  dplyr::mutate(variable = ifelse(year_grp == 1, "1951-1970", "2001-2020")) %>%
  #dplyr::left_join(locs %>% sf::st_drop_geometry() %>% dplyr::select(areaid, areanameen)) %>%
  ggplot() + 
  geom_line(aes(monthdengue, exp(tmean_effect), group=year, color=variable), alpha=0.2, show.legend = FALSE, size=0.2) + 
  geom_point(data = mnx[ mnx$areaid %in% high, ], aes(monthdengue, exp(mean_risk), group=variable, color=variable), position=position_dodge(width=0.5), size=1.3, alpha=1) +  
  geom_line(data = mnx[ mnx$areaid %in% high, ], aes(monthdengue, exp(mean_risk), group=variable, color=variable), position=position_dodge(width=0.5), size=0.5, alpha=1) +
  geom_linerange(data = mnx[ mnx$areaid %in% high, ], aes(monthdengue, ymin=exp(lower), ymax=exp(upper), group=variable, color=variable), position=position_dodge(width=0.5), size=0.6, alpha=1, show.legend = FALSE) +
  theme_classic() +
  lemon::facet_rep_wrap(~tag, nrow=2) +
  scale_x_continuous(breaks=c(1, 3, 6, 9, 12)) +
  geom_hline(yintercept=1, lty=2) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=10.5), 
        axis.text = element_text(size=13), 
        axis.title = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=14), 
        #legend.position="bottom",
        axis.title.x = element_blank(),
        plot.title = element_text(size=18, hjust=0.5)) +
  scale_color_viridis_d(begin=0.05, end=0.65, direction=-1, name="Time period") +
  xlab("Dengue month (May to April)") +
  ylab("Dengue relative risk (mean ± 95% CI)") +
  ggtitle("Districts with greatest relative increase in temperature-driven risk")

# plot
ps2 = clim_for_test %>% 
  dplyr::filter(areaid %in% low) %>%
  dplyr::mutate(variable = ifelse(year_grp == 1, "1951-1970", "2001-2020")) %>%
  #dplyr::left_join(locs %>% sf::st_drop_geometry() %>% dplyr::select(areaid, areanameen)) %>%
  ggplot() + 
  geom_line(aes(monthdengue, exp(tmean_effect), group=year, color=variable), alpha=0.2, show.legend = FALSE, size=0.2) + 
  geom_point(data = mnx[ mnx$areaid %in% low, ], aes(monthdengue, exp(mean_risk), group=variable, color=variable), position=position_dodge(width=0.5), size=1.3, alpha=1) +  
  geom_line(data = mnx[ mnx$areaid %in% low, ], aes(monthdengue, exp(mean_risk), group=variable, color=variable), position=position_dodge(width=0.5), size=0.5, alpha=1) +
  geom_linerange(data = mnx[ mnx$areaid %in% low, ], aes(monthdengue, ymin=exp(lower), ymax=exp(upper), group=variable, color=variable), position=position_dodge(width=0.5), size=0.6, alpha=1, show.legend = FALSE) +
  theme_classic() +
  lemon::facet_rep_wrap(~tag, nrow=2) +
  scale_x_continuous(breaks=c(1, 3, 6, 9, 12)) +
  geom_hline(yintercept=1, lty=2) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=10.5), 
        axis.text = element_text(size=13), 
        axis.title = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=14), 
        #legend.position="bottom",
        plot.title = element_text(size=18, hjust=0.5)) +
  scale_color_viridis_d(begin=0.05, end=0.65, direction=-1, name="Time period") +
  xlab("Dengue month (May to April)") +
  ylab("Dengue relative risk (mean ± 95% CI)") +
  ggtitle("Districts with greatest relative decrease in temperature-driven risk")

library(patchwork)
comb = ps1 + ps2
comb = comb + plot_layout(guides = "collect", heights=c(1, 1), nrow=2) & theme(legend.position = "bottom")
ggsave(comb, file="./output/figures/SuppFigure_TemperatureRiskPredictions_Top10.jpg", device="jpg", units="in", width=18, height=19, dpi=600, scale=0.6)






# ================ Sensitivity test: baseline period 1971-1990 ======================

# set historical and present-day periods
clim_for_test = clim %>%
  dplyr::mutate(year = lubridate::year(date)) %>%
  dplyr::filter(year %in% c(1971:1990, 2001:2020)) %>%
  dplyr::mutate(year_grp = ifelse(year < 2000, 1, 2),
                month = lubridate::month(date)) %>%
  dplyr::select(areaid, year, year_grp, date, month, tmean_1m, tmean_effect)

# dengue month
dm = data.frame(month=c(5:12, 1:4), monthdengue = c(1:12), month_name = c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"))
clim_for_test = left_join(clim_for_test, dm)

# monthly linear models; compare between baseline and test period
mnx = clim_for_test %>%
  dplyr::filter(!is.na(tmean_effect)) %>%
  dplyr::group_by(areaid, monthdengue) %>%
  dplyr::summarise(month_name = head(month_name, 1),
                   tmean_fac_slope = coef(summary(lm(tmean_effect ~ factor(year_grp))))[ 2, 1],
                   tmean_fac_p = coef(summary(lm(tmean_effect ~ factor(year_grp))))[ 2, 4],
                   mean_risk = mean(tmean_effect[ year > 2000 ]))

# only show slopes where p < 0.05
lims = (exp(mnx$tmean_fac_slope)-1)*100
lims = c(-max(abs(lims), na.rm=TRUE), max(abs(lims), na.rm=TRUE))

# mapping geography of change
map2 = shp %>%
  dplyr::full_join(mnx) %>%
  dplyr::filter(!is.na(monthdengue)) %>%
  dplyr::mutate(tmean_fac_slope = replace(tmean_fac_slope, tmean_fac_p > 0.05, NA)) %>%
  left_join(dm) %>%
  dplyr::mutate(month_name = factor(month_name, levels=dm$month_name, ordered=TRUE)) %>%
  ggplot() +
  geom_sf(aes(fill=(exp(tmean_fac_slope)-1)*100), color=NA) +
  maptheme +
  facet_wrap(~month_name, nrow=1) +
  ggtitle("Change in dengue risk from historical reference (1971-1990) to present (2001-2020)") +
  geom_sf(data = shp_vt, fill=NA, col="grey20", size=0.25) +
  scale_fill_gradientn(colors=rev(pals::brewer.rdylbu(200)), na.value="white", limits=lims, name="Risk\nchange\n(%)", breaks=c(-50, -25, 0, 25, 50)) +
  #scale_fill_gradient2(high=scales::muted("blue"), low=scales::muted("red"), midpoint=0, name="Annual\nrisk\nchange\n(%)", na.value="grey95") +
  theme(strip.text = element_blank(), plot.title=element_text(size=12.5), legend.title = element_text(size=12), legend.text = element_text(size=11))



# ---------- visualise graphs of risk over time for some example regions --------------

locs = shp %>% dplyr::filter(areanameen %in% c("Hoan Kiem", "Nha Trang", "Quan 1", "Buon Ma Thuot", "Binh Thuy"))

clim_for_test = clim %>%
  dplyr::filter(areaid %in% locs$areaid) %>%
  dplyr::mutate(year = lubridate::year(date)) %>%
  dplyr::filter(year %in% c(1971:1990, 2001:2020)) %>%
  dplyr::mutate(year_grp = ifelse(year < 2000, 1, 2),
                month = lubridate::month(date)) %>%
  dplyr::select(areaid, year, year_grp, date, month, tmean_1m, tmean_effect)

dm = data.frame(month=c(5:12, 1:4), monthdengue = c(1:12), month_name = c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"))
clim_for_test = left_join(clim_for_test, dm)

mnx = clim_for_test %>%
  dplyr::filter(!is.na(tmean_effect)) %>%
  dplyr::group_by(areaid, monthdengue, year_grp) %>%
  dplyr::summarise(month_name = head(month_name, 1),
                   mean_tmean = mean(tmean_1m),
                   se_tmean = sd(tmean_1m) / sqrt(length(tmean_1m)),
                   mean_risk = mean(tmean_effect),
                   sd_risk = sd(tmean_effect),
                   se_risk = sd_risk / sqrt(length(tmean_effect))) %>%
  dplyr::mutate(
    upper = mean_risk + (1.96*se_risk),
    lower = mean_risk - (1.96*se_risk),
    upper_tmean = mean_tmean + (1.96*se_tmean),
    lower_tmean = mean_tmean - (1.96*se_tmean)
    # upper = mean_risk + sd_risk,
    # lower = mean_risk - sd_risk
  ) %>%
  dplyr::mutate(variable = ifelse(year_grp == 1, "1971-1990", "2001-2020")) %>%
  #dplyr::left_join(locs %>% sf::st_drop_geometry() %>% dplyr::select(areaid, areanameen)) %>%
  dplyr::mutate(
    loc = "Buon Ma Thuot (Dak Lak)",
    loc = replace(loc, areaid == "70165", "Ha Noi"),
    loc = replace(loc, areaid == "70468", "Nha Trang (Khanh Hoa)"),
    loc = replace(loc, areaid == "70563", "TP. Ho Chi Minh"),
    loc = replace(loc, areaid == "70648", "Can Tho"),
    loc = factor(loc, levels = c("Ha Noi",
                                 "Buon Ma Thuot (Dak Lak)",
                                 "Nha Trang (Khanh Hoa)",
                                 "TP. Ho Chi Minh",
                                 "Can Tho"), ordered=TRUE))

# plot
p3 = clim_for_test %>%
  dplyr::mutate(variable = ifelse(year_grp == 1, "1971-1990", "2001-2020")) %>%
  #dplyr::left_join(locs %>% sf::st_drop_geometry() %>% dplyr::select(areaid, areanameen)) %>%
  dplyr::mutate(
    loc = "Buon Ma Thuot (Dak Lak)",
    loc = replace(loc, areaid == "70165", "Ha Noi"),
    loc = replace(loc, areaid == "70468", "Nha Trang (Khanh Hoa)"),
    loc = replace(loc, areaid == "70563", "TP. Ho Chi Minh"),
    loc = replace(loc, areaid == "70648", "Can Tho"),
    loc = factor(loc, levels = c("Ha Noi",
                                 "Buon Ma Thuot (Dak Lak)",
                                 "Nha Trang (Khanh Hoa)",
                                 "TP. Ho Chi Minh",
                                 "Can Tho"), ordered=TRUE)) %>%
  ggplot() +
  geom_line(aes(monthdengue, exp(tmean_effect), group=year, color=variable), alpha=0.2, show.legend = FALSE, size=0.2) +
  geom_point(data = mnx, aes(monthdengue, exp(mean_risk), group=variable, color=variable), position=position_dodge(width=0.5), size=1.3, alpha=1) +
  geom_line(data = mnx, aes(monthdengue, exp(mean_risk), group=variable, color=variable), position=position_dodge(width=0.5), size=0.5, alpha=1) +
  geom_linerange(data = mnx, aes(monthdengue, ymin=exp(lower), ymax=exp(upper), group=variable, color=variable), position=position_dodge(width=0.5), size=0.6, alpha=1, show.legend = FALSE) +
  facet_wrap(~areaid) +
  theme_classic() +
  lemon::facet_rep_wrap(~loc, nrow=1) +
  scale_x_continuous(breaks=1:12) +
  geom_hline(yintercept=1, lty=2) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=9.25),
        axis.title.x = element_text(size=11.5),
        axis.title.y = element_text(size=11),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11)) +
  scale_color_viridis_d(begin=0.05, end=0.65, direction=-1, name="Time period") +
  xlab("Dengue month (May to April)") +
  ylab("Dengue risk (mean ± 95% CI)")


# ----------- save ------------

# # combine with maps and save
# library(patchwork)
# comb = map1 + map2 + p3
# comb = comb + plot_layout(guides = "collect", heights=c(1.2, 1.2, 1))# + plot_annotation(tag_levels = 'a') + theme(plot.tag = element_text(size = 20, face = "bold"))
# ggsave(comb, file="./output/figures/SuppFigure_TemperatureProjections_LongTerm_1981baseline.jpg", device="jpg", units="in", width=18, height=10, dpi=600, scale=0.65)

comb = map2 + p3
comb = comb + plot_layout(guides = "collect", heights=c(1.2, 1))# + plot_annotation(tag_levels = 'a') + theme(plot.tag = element_text(size = 20, face = "bold"))
ggsave(comb, file="./output/figures/SuppFigure_TemperatureProjections_LongTerm_1971baseline.jpg", device="jpg", units="in", width=18, height=7, dpi=600, scale=0.65)


