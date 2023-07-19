

# ================ Examining interactions between water infrastructure and SPEI metrics ================


library(dplyr)
library(sf)
library(ggplot2)
library(vroom)
library(raster)

PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

source("./scripts/04_modelling/00_inla_setup_functions_r4.R")
source("./scripts/04_modelling/00_plot_themes.R")

# directory where model objects are stored (large files so external)
models_dir = "C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/output/model_outputs/results_23/infraspei_south/"



# ================= shapefile and dengue data =================

# districts to be excluded (offshore)
# only ones with substantial dengue cases are in Kien Giang; could be worth exploring including them but thisis sufficient for now
offshore_areas = c(70154, 70339, 70273, 70355, 70698)

# dengue, regions, climate, landuse, connectivity data
dd = read.csv("./output/model_data/ModelData_Dengue_VietAll.csv") %>%
  dplyr::left_join(read.csv("./output/model_data/ModelData_ClimaticRegions.csv")) %>%
  dplyr::left_join(read.csv("./output/model_data/ModelData_SocioEcologicalCovar_VietAll.csv")) %>%
  dplyr::left_join(vroom::vroom("./output/model_data/ModelData_ClimateLags_VietAll.csv.gz", col_types = list(areaid=vroom::col_character(), date=vroom::col_character()))) %>%
  dplyr::mutate(date = as.Date(date)) %>%
  dplyr::filter(region2 %in% c("South", "South Central Coast", "Central Highlands")) %>%
  dplyr::filter(!areaid %in% offshore_areas)

# provincex
dd$provincex = as.integer(as.factor(dd$province))

# districts shapefile for Vietnam
shp = st_read("./data/shapefiles/dengue_districts_shapefile.shp") %>%
  dplyr::filter(areaid %in% dd$areaid)

# provinces shapefile with regions info, and vietnam outline for mapping
shp_prov = st_read("./data/shapefiles/provinces.shp") %>%
  dplyr::filter(provincena %in% dd$province)
st_crs(shp_prov) = st_crs(shp)
shp_prov = st_crop(shp_prov, shp)
rr = read.csv("./data/shapefiles/regions_lookup.csv")
shp_prov = left_join(shp_prov, rr[ , c("provincena", "region")])

# polyid indicator
shpf = shp[ shp$areaid %in% dd$areaid, ]
id_ref = data.frame(areaid = shpf$areaid, polyid = 1:nrow(shpf))
dd = left_join(dd, id_ref)

# add polygon area, lat lon, region information
#shp$area_km2 = as.vector(st_area(shp) / 10^6)
shp = cbind(shp, as.data.frame(st_coordinates(st_centroid(shp))) %>% dplyr::rename("longitude"=1, "latitude"=2))
dd = left_join(dd, shp[ , c("areaid", "latitude", "longitude")] %>% st_drop_geometry())
shp = left_join(shp, dd[ !duplicated(dd$areaid), c("areaid", "region1", "region2", "region3") ])

# three categories (low, normal, high) of piped water and urbanicity
cats = dd %>%
  dplyr::select(areaid, year, water_piped_year, urban_pw) %>%
  distinct()

# terciles of water
tercs = c(0.25, 0.75)
cats$water_c = 1
cats$water_c[ cats$water_piped_year >= tercs[1] ] = 2
cats$water_c[ cats$water_piped_year >= tercs[2] ] = 3
table(cats$water_c)

# urban >= 0.25 or 0.75
cats$urban_c = 1
cats$urban_c[ cats$urban_pw >= 0.25 ] = 2
cats$urban_c[ cats$urban_pw >= 0.75 ] = 3

# combine
dd = dd %>%
  dplyr::left_join(cats[ , c(1, 2, 5, 6)])





# ================= summary information criteria for fitted models ================

# functions for reading in
readfile = function(x){
  foo = read.csv(x)
  if(!"covar" %in% names(foo)){ foo$covar = "" }
  if(!"model_sub" %in% names(foo)){ foo$model_sub = "" }
  if(!"model_filename" %in% names(foo)){ foo$model_filename = "" }
  foo
}

# socioecological models
# tmean better than tmin
ff = list.files(paste(models_dir, "fitmetrics/", sep=""), full.names=TRUE, pattern=".csv")

fx1 = do.call(rbind.data.frame, lapply(ff, readfile)) %>%
  dplyr::mutate(covar = modname) %>%
  dplyr::arrange(waic) %>%
  dplyr::select(-formula, -fx, -modname) %>%
  dplyr::mutate(deltaDIC = dic - dic[ candidate == "full"],
                deltaWAIC = waic - waic[ candidate == "full"],
                deltaLS = logscore - logscore[ candidate == "full"]) %>%
  ungroup() 



# ================= visualise change in model adequacy metrics ===================

mod_names = c("Full model", "SPEI-1 excluded", "SPEI-1 (1m) x Water", "SPEI-1 (1m) x Urban", "SPEI-6 excluded", "SPEI-6 (5m) x Water", "SPEI-6 (5m) x Urban")
md = data.frame(candidate = c("full", 
                              "no_spei1", "spei1water", "spei1urban",
                              "no_spei6", "spei6water", "spei6urban"),
                cand2 = mod_names)

datx = fx1 %>%
  dplyr::select(candidate, deltaDIC, deltaWAIC, deltaLS) %>%
  reshape2::melt(id.vars=1) %>%
  dplyr::filter(!candidate %in% c("baseline")) %>%
  dplyr::left_join(
    data.frame(variable = c("deltaDIC", "deltaWAIC", "deltaLS"), var2 = c("DIC", "WAIC", "log-score"))
  ) %>%
  dplyr::mutate(
    var2 = factor(var2, levels=c("WAIC", "DIC", "log-score"), ordered=TRUE)
  ) %>%
  dplyr::left_join(
    md
  ) %>%
  dplyr::mutate(
    cand2 = factor(cand2, levels=rev(mod_names), ordered=TRUE)
  )

p1 = datx %>%
  dplyr::filter(cand2 != "Full model") %>%
  ggplot() + 
  geom_point(aes(cand2, value), size=3.5, pch=21 ,fill="coral2") + 
  geom_hline(data=datx[ datx$candidate == "full",], aes(yintercept=value), lty=2) +
  facet_wrap(~var2, scales="free_x") +
  coord_flip() +
  theme_bw() +
  theme(strip.text = element_text(size=14),
        strip.background = element_blank(),
        axis.title = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size=11)) +
  ylab("Difference in metric relative to full model") +
  xlab("Model")
ggsave(p1, file="./output/figures/SuppFigure_InfrastructureSPEI_WSmetrics.jpg", device="jpg", units="in", width=9.2, height=2.5, dpi=600)




# ================== visualise fitted functions ====================


# ----------- 1. Piped water * SPEI6 -------------

# read in model and visualise
load(paste(models_dir, "models/", fx1$model_filename[ fx1$candidate == "spei6water" ],  sep=""))

summary(mod_i)

# plot
eff1 = extractRandomINLA(mod_i$summary.random$spei6_5m_g, effect_name="spei", transform=TRUE) 
xlims = c(min(eff1$value), max(eff1$value))
#xlims = range(ddf$spei6_4m)
ylims = range(c(eff1$lower, eff1$upper))

p1 = eff1 %>%
  dplyr::left_join(
    data.frame(group=1:3, level=c("Low (<25%)", "Intermediate", "High (>75%)"))
  ) %>%
  dplyr::mutate(level = factor(level, levels=c("Low (<25%)", "Intermediate", "High (>75%)"), ordered=TRUE)) %>%
  ggplot() + 
  ylim(ylims) +
  xlim(xlims) +
  # geom_rect(aes(xmin=-2, xmax=-1, ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey95", alpha=0.1) +
  # geom_rect(aes(xmin=xlims[1], xmax=-2, ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey90", alpha=0.1) +
  geom_rect(aes(xmin=-2, xmax=2, ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey91", alpha=0.1) +
  geom_rect(aes(xmin=-1, xmax=1, ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey84", alpha=0.1) +
  # geom_rect(aes(xmin=1, xmax=2, ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey95", alpha=0.1) +
  # geom_rect(aes(xmin=2, xmax=xlims[2], ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey90", alpha=0.1) +
  geom_hline(yintercept=c(2, 3), col="grey92", size=0.3) +
  geom_hline(yintercept=1, lty=2, col="grey50", size=0.4) +
  geom_vline(xintercept=0, lty=1, col="grey80") +
  geom_ribbon(aes(value, ymin=lower, ymax=upper, group=level, fill=factor(level)), alpha=0.5) + 
  geom_line(aes(value, median, group=level, col=factor(level)), size=0.8, alpha=1) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
    axis.title = element_text(size=13), axis.text= element_text(size=12),
    axis.title.x = element_blank(),
    legend.position = c(0.8, 0.85),
    legend.title=element_text(face="bold")) + 
  scale_fill_viridis_d(end=0.8, option="viridis", name="Improved\nwater access") + 
  scale_color_viridis_d(end=0.5, option="viridis", name="Improved\nwater access") + 
  xlab("SPEI-6 (5m lag)") + 
  ylab("Dengue relative risk")



# ----------- 2. Urbanisation * SPEI6 -------------

# read in model and visualise
load(paste(models_dir, "models/", fx1$model_filename[ fx1$candidate == "spei6urban" ],  sep=""))

summary(mod_i)

# plot
eff1 = extractRandomINLA(mod_i$summary.random$spei6_5m_g, effect_name="spei", transform=TRUE) 
xlims2 = c(min(eff1$value), max(eff1$value))
ylims2 = range(c(eff1$lower, eff1$upper))

p2 = eff1 %>%
  dplyr::left_join(
    data.frame(group=1:3, level=c("Low (<25%)", "Intermediate", "High (>75)"))
  ) %>%
  dplyr::mutate(level = factor(level, levels=c("Low (<25%)", "Intermediate", "High (>75)"), ordered=TRUE)) %>%
  ggplot() + 
  ylim(ylims2) +
  # geom_rect(aes(xmin=-2, xmax=-1, ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey95", alpha=0.1) +
  # geom_rect(aes(xmin=xlims[1], xmax=-2, ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey90", alpha=0.1) +
  geom_rect(aes(xmin=-2, xmax=2, ymin=ylims2[1], ymax=ylims2[2]), col=NA, fill="grey91", alpha=0.1) +
  geom_rect(aes(xmin=-1, xmax=1, ymin=ylims2[1], ymax=ylims2[2]), col=NA, fill="grey84", alpha=0.1) +
  # geom_rect(aes(xmin=1, xmax=2, ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey95", alpha=0.1) +
  # geom_rect(aes(xmin=2, xmax=xlims[2], ymin=ylims[1], ymax=ylims[2]), col=NA, fill="grey90", alpha=0.1) +
  geom_hline(yintercept=c(2, 3, 4), col="grey92", size=0.3) +
  geom_hline(yintercept=1, lty=2, col="grey50", size=0.4) +
  geom_vline(xintercept=0, lty=1, col="grey80") +
  geom_ribbon(aes(value, ymin=lower, ymax=upper, group=level, fill=factor(level)), alpha=0.5) + 
  geom_line(aes(value, median, group=level, col=factor(level)), size=0.8, alpha=1) +
  theme_classic() +
  theme(
    #legend.position="top", 
    panel.grid.minor = element_blank(),
    axis.title = element_text(size=13), axis.text= element_text(size=12),
    axis.title.x = element_blank(),
    legend.position = c(0.85, 0.8),
    legend.title=element_text(face="bold")) + 
  scale_fill_viridis_d(end=0.8, name="Built-up land") + 
  scale_color_viridis_d(end=0.5, name="Built-up land") + 
  xlab("SPEI-6 (4m lag)") + 
  ylab("Dengue relative risk")




# -------------- 3. histogram of distribution ----------------

# by water
p3a = dd %>% 
  dplyr::mutate(water_c = factor(water_c, levels=4:1, ordered=TRUE)) %>%
  ggplot() + 
  geom_vline(xintercept=0, lty=1, size=1, alpha=0.7, col="grey50") + 
  geom_vline(xintercept=c(-2, -1, 2, 1), lty=1, size=0.5, alpha=0.7, col="grey70") + 
  geom_histogram(aes(x=spei6_5m, y=..count../sum(..count..), group=water_c, fill=water_c), col=NA, alpha=0.7, size=1, bins=45) + 
  theme_classic() + 
  theme(axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=8, color="white"),
        axis.text.y= element_text(size=10),
        axis.text.x= element_text(size=12),
        legend.position="none") + 
  xlab("SPEI-6 (5m lag)") + 
  ylab("Prop. obs") +
  scale_fill_viridis_d(end=0.8, name="Water", direction=-1, guide=guide_legend(reverse=TRUE)) 

# by urban
p3b = dd %>% 
  dplyr::mutate(urban_c = factor(urban_c, levels=4:1, ordered=TRUE)) %>%
  ggplot() + 
  geom_vline(xintercept=0, lty=1, size=1, alpha=0.7, col="grey50") + 
  geom_vline(xintercept=c(-2, -1, 2, 1), lty=1, size=0.5, alpha=0.7, col="grey70") + 
  geom_histogram(aes(x=spei6_5m, y=..count../sum(..count..), group=urban_c, fill=urban_c), col=NA, alpha=0.7, size=1, bins=45) + 
  theme_classic() + 
  theme(axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=6, color="white"),
        axis.text.y= element_text(size=10),
        axis.text.x= element_text(size=12), 
        legend.position="none") + 
  xlab("SPEI-6 (5m lag)") + 
  ylab("Prop. obs") +
  scale_fill_viridis_d(end=0.8, name="Urban", direction=-1, guide=guide_legend(reverse=TRUE)) 



# -------------- create combined plots ---------------

# save
pc1 = gridExtra::grid.arrange(p1, p3a, nrow=2, heights=c(1, 0.3))
ggsave(pc1, file="./output/figures/Figure6_InfraSPEI_noscenarios.jpg", device="jpg", units="in", width=5, height=6.35, dpi=600, scale=0.95)

pc2 = gridExtra::grid.arrange(p2, p3b, nrow=2, heights=c(1, 0.4))
ggsave(pc2, file="./output/figures/SuppFigure_SPEIUrban.jpg", device="jpg", units="in", width=5, height=6.75, dpi=600, scale=0.95)







# ==================== Scenarios of relative risk responses to SPEI in high and low piped water settings =======================

# Example SPEI timeseries from an example location in Mekong Delta (highly endemic, very variable piped water access)
# Show different relative risk functions over time assuming high/low access

# read in model 
load(paste(models_dir, "models/", fx1$model_filename[ fx1$candidate == "spei6water" ],  sep=""))

eff1 = extractRandomINLA(mod_i$summary.random$spei6_5m_g, effect_name="spei", transform=FALSE) 

interpRisk = function(df){
  
  result = data.frame()
  
  df$value = as.numeric(df$value)
  df$value = round(df$value, 2)

  for(i in 2:nrow(df)){
    
    vals = df$value[i-1] + seq(0, 10, by=0.01)
    vals = vals[ vals >= df$value[i-1] & vals <= df$value[i] ]
    
    dfx = data.frame(
      value = vals
    ) %>%
      dplyr::mutate(
        lower = seq(df$lower[i-1], df$lower[i], length.out=length(value)),
        median = seq(df$median[i-1], df$median[i], length.out=length(value)),
        upper = seq(df$upper[i-1], df$upper[i], length.out=length(value))
      )
    
    result = rbind(result, dfx)
  }
  
  result = result %>% distinct()
  
  # append extremes
  lowest = data.frame(
    value = result$value[1] - seq(0, 5, by=0.01),
    lower = result$lower[1],
    median = result$median[1],
    upper = result$upper[1]
  ) %>%
    dplyr::arrange(value)
  
  highest = data.frame(
    value = result$value[nrow(result)] + seq(0, 5, by=0.01),
    lower = result$lower[nrow(result)],
    median = result$median[nrow(result)],
    upper = result$upper[nrow(result)]
  ) %>%
    dplyr::arrange(value)
  
  result = lowest %>%
    rbind(result) %>%
    rbind(highest) %>%
    dplyr::mutate(value = as.numeric(value))
  
  result
}

eff1 = do.call(
  rbind.data.frame,
  list(
    interpRisk(eff1[ eff1$group == 1, ]) %>% dplyr::mutate(water = "Low"),
    interpRisk(eff1[ eff1$group == 2, ]) %>% dplyr::mutate(water = "Intermediate"),
    interpRisk(eff1[ eff1$group == 3, ]) %>% dplyr::mutate(water = "High")
  )
) %>%
  dplyr::rename("spei6_5m"=1)


# SPEI timeseries from Hong Ngu district in Dong Thap province
# chosen because contains a wide range of SPEI values including extremes, in a region with substantial changes in piped water access
# over the time period
proj = dd %>%
  dplyr::mutate(date = as.Date(date)) %>%
  #dplyr::filter(district == "Phu Tan") %>%
  #dplyr::filter(district == "Ha Tien") %>%
  #dplyr::filter(district == "Cho Lach") %>%
  #dplyr::filter(district == "Ben Tre") %>%
  #dplyr::filter(district == "Vung Liem") %>%
  #dplyr::filter(district == "My Tho") %>%
  dplyr::filter(district == "HongNgu") %>%
  dplyr::select(
    province,
    incidence,
    date, 
    district, 
    spei6_5m
  ) %>%
  dplyr::mutate(spei6_5m = round(spei6_5m, 2))

# low
low = left_join(
  proj %>% dplyr::mutate(spei6_5m = as.character(spei6_5m)), 
  eff1[ eff1$water == "Low", ] %>% dplyr::mutate(spei6_5m = as.character(spei6_5m))
)
#low[ , c("lower", "median", "upper")] = low[ , c("lower", "median", "upper")] + -0.07992131

# medium
medium = left_join(
  proj %>% dplyr::mutate(spei6_5m = as.character(spei6_5m)), 
  eff1[ eff1$water == "Intermediate", ] %>% dplyr::mutate(spei6_5m = as.character(spei6_5m))
)
#medium[ , c("lower", "median", "upper")] = medium[ , c("lower", "median", "upper")] + 0.05529191 

# high
high = left_join(
  proj %>% dplyr::mutate(spei6_5m = as.character(spei6_5m)), 
  eff1[ eff1$water == "High", ] %>% dplyr::mutate(spei6_5m = as.character(spei6_5m))
)
#high[ , c("lower", "median", "upper")] = high[ , c("lower", "median", "upper")] + 0.02462948

scen = rbind(low, high) %>%
  rbind(medium)

scaleFUN <- function(x) sprintf("%.1f", x)
sp = scen %>%
  dplyr::mutate(metric = "SPEI-6 (5-month lag)") %>%
  ggplot() + 
  geom_line(aes(date, as.numeric(spei6_5m)), size=0.5) + 
  geom_hline(yintercept=0, lty=2) + 
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        panel.grid.major = element_line(color="grey90", size=0.5),
        panel.grid.minor = element_line(color="grey95", size=0.3),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=11.5),
        axis.title.y = element_text(size=12.8)) + 
  #scale_y_continuous(breaks = seq(-2, 2, by=1), labels=seq(-2, 2, by=1)) +
  scale_y_continuous(labels = scaleFUN) +
  ylab("SPEI-6 (5m lag)")
  
pp = scen %>%
  #dplyr::mutate(water = factor(water, levels=c("Low", "High"), ordered=TRUE)) %>%
  dplyr::mutate(water = factor(water, levels=c("Low", "Intermediate", "High"), ordered=TRUE)) %>%
  ggplot() + 
  geom_ribbon(aes(date, ymin=exp(lower), ymax=exp(upper), group=water, fill=water), alpha=0.2) + 
  geom_line(aes(date, exp(median), group=water, col=water), size=0.5) + 
  theme_classic() +
  #facet_wrap(~water, nrow=2) +
  geom_hline(yintercept=1, lty=2) +
  scale_fill_viridis_d(end=0.8, option="viridis", name="Water access") + 
  scale_color_viridis_d(end=0.6, option="viridis", name="Water access") +
  theme(legend.position=c(0.88, 0.85)) +
  ylab("Relative risk") +
  xlab("Month") +
  theme(
    panel.grid.major = element_line(color="grey90", size=0.5),
    panel.grid.minor = element_line(color="grey95", size=0.3),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12.88),
    axis.title.x = element_text(size=13.2),
    axis.title.y = element_text(size=13.2),
    legend.title = element_blank(),
    legend.text = element_text(size=12)
  )

# combine
scen_plot = gridExtra::grid.arrange(sp, pp, nrow=2, heights=c(0.5, 1))





# ================== combined plot for manuscript ==================

fig_fin = gridExtra::grid.arrange(pc1, scen_plot, nrow=1, widths=c(1, 1.6))

fig_fin = ggpubr::as_ggplot(fig_fin)  +
  cowplot::draw_plot_label(label = c("a", "b"), 
                           fontface = "bold", size = 24, 
                           x = c(0.007, 0.438), y = c(1, 1))
ggsave(fig_fin, file="./output/figures/Figure6_InfraSPEI_scenarios.png", device="png", units="in", width=13, height=6.35, dpi=600, scale=0.9)







