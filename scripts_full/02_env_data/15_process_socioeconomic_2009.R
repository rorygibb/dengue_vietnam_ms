

# ------------------------ MapVietnam socioeconomic variables -----------------

# matches MapVIETNAM (census 2009) socioeconomic data to districts shapefile/naming scheme

# mapVIETNAM is an interactive map of socioeconomic data at the province and district level for Vietnam. 
# Using Census data and periodical data from Ministry of Labor - Invalids and Social Affairs and other related agencies, 
# it offers data on several indicators alongside the poverty rates in that province/district.
# http://www5.worldbank.org/mapvietnam/#en

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



# ------------ load polygons data ------------

# polygons for entire of Vietnam (same areaid matches)
vx = st_read("./data/shapefiles/d_moss/vietnam_districts_final.shp")
vx = vx[ ! vx$areanameen %in% c("Truong Sa", "Hoang Sa"), ] # remove offshore districts

# create lookup
# lookup = vx %>% 
#   st_drop_geometry() %>%
#   select("areaid", "areanameen", "areaprovin") %>%
#   arrange(areaprovin, areanameen)
# write.csv(lookup, "./data/socioeconomic/worldbank_indicators_2009/dmoss_lookup.csv", row.names=FALSE)


# --------------

# read world bank data
wb = read.csv("./data/socioeconomic/worldbank_indicators_2009/data/all_district2_csv.csv", stringsAsFactors = FALSE)

# read lookup chart to match to shapefile
lookup = read.csv("./data/socioeconomic/worldbank_indicators_2009/data/all_district_shapefile_lookup_completedRG.csv", stringsAsFactors = FALSE)

# merge
names(lookup)[1] = "ID"
names(wb)[1] = "ID"
wb = left_join(wb, lookup[ , c("ID", "areaid_dmoss", "areaname_dmoss")] )

# calculate world bank pop estimate
wb$population = as.numeric(unlist(
  lapply(
    strsplit(wb$population_n, ","),
    paste, 
    collapse=""
  )
))

# calculate world bank pop estimate
wb$workingpop = as.numeric(unlist(
  lapply(
    strsplit(wb$workingpop_n, ","),
    paste, 
    collapse=""
  )
))

# calculate world bank household estimate
wb$households = as.numeric(unlist(
  lapply(
    strsplit(wb$households_n, ","),
    paste, 
    collapse=""
  )
))



# --------------- convert numbers to % of households ------------------

# population metrics
wb$poverty_gsowb_headcount_n = as.numeric(unlist(lapply(strsplit(wb$poverty_gsowb_headcount_n, ","), paste, collapse="")))
wb$poverty_popinnationalbottom40percent_n = as.numeric(unlist(lapply(strsplit(wb$poverty_popinnationalbottom40percent_n, ","), paste, collapse="")))
wb$poverty_gsowb_extremepoverty_headcount_n  = as.numeric(unlist(lapply(strsplit(wb$poverty_gsowb_extremepoverty_headcount_n , ","), paste, collapse="")))
wb$mainemployment_agriculture_n   = as.numeric(unlist(lapply(strsplit(wb$mainemployment_agriculture_n, ","), paste, collapse="")))
wb$mainemployment_nonfarmselfemployment_n  = as.numeric(unlist(lapply(strsplit(wb$mainemployment_nonfarmselfemployment_n , ","), paste, collapse="")))
wb$mainemployment_wage_n  = as.numeric(unlist(lapply(strsplit(wb$mainemployment_wage_n, ","), paste, collapse="")))

# households metrics
wb$mainlightsource_electricity_n  = as.numeric(unlist(lapply(strsplit(wb$mainlightsource_electricity_n, ","), paste, collapse="")))
wb$sanitation_flushtoilet_indoor_n  = as.numeric(unlist(lapply(strsplit(wb$sanitation_flushtoilet_indoor_n, ","), paste, collapse="")))
wb$sanitation_flushtoilet_outdoor_n  = as.numeric(unlist(lapply(strsplit(wb$sanitation_flushtoilet_outdoor_n, ","), paste, collapse="")))
wb$sanitation_flushtoilet_any_n  = as.numeric(unlist(lapply(strsplit(wb$sanitation_flushtoilet_any_n, ","), paste, collapse="")))
wb$water_piped_well_n  = as.numeric(unlist(lapply(strsplit(wb$water_piped_well_n, ","), paste, collapse="")))

wb$secondaryschool_overall_1118_n  = as.numeric(unlist(lapply(strsplit(wb$secondaryschool_overall_1118_n, ","), paste, collapse="")))
wb$secondaryschool_lower_n  = as.numeric(unlist(lapply(strsplit(wb$secondaryschool_lower_n, ","), paste, collapse="")))
wb$secondaryschool_upper_n  = as.numeric(unlist(lapply(strsplit(wb$secondaryschool_upper_n, ","), paste, collapse="")))



# ---------------- merge polygons for newly merged ids ---------------

# merge shapefile
shpm = sf::st_read("./data/shapefiles/d_moss/vietnam_districts_merged.shp") %>%
  dplyr::filter(!areanameen %in% c("Truong Sa", "Hoang Sa"))

# to merge
tomerge = wb %>% dplyr::filter(!areaid_dmoss %in% shpm$areaid)

# new areaids
merge_lk = shpm[ !shpm$areaid %in% wb$areaid, ] %>%
  st_drop_geometry() %>%
  dplyr::mutate(idx = areaid) %>%
  tidyr::separate_rows(areaid, sep="_") %>%
  dplyr::mutate(areaid = as.numeric(areaid))

# combine
tomerge = left_join(tomerge, merge_lk[ , c("areaid", "idx")], by=c("areaid_dmoss"="areaid"))

# calculate total ns for ns 
merged1 = tomerge %>%
  dplyr::select(idx, population, households, workingpop, 
                poverty_gsowb_headcount_n, poverty_popinnationalbottom40percent_n, poverty_gsowb_extremepoverty_headcount_n, 
                mainemployment_agriculture_n, mainemployment_nonfarmselfemployment_n, mainemployment_wage_n,
                mainlightsource_electricity_n, sanitation_flushtoilet_indoor_n, sanitation_flushtoilet_outdoor_n, sanitation_flushtoilet_any_n, water_piped_well_n) %>%
  dplyr::group_by(idx) %>%
  dplyr::summarise_all("sum", na.rm=TRUE) %>%
  dplyr::rename("areaid" = idx)

# calculate pop weighted % for %
merged2 = tomerge %>%
  dplyr::select(idx, population, population_under125perday_percent, population_under2perday_percent,
                secondaryschool_overall_1118_percent, secondaryschool_lower_percent, secondaryschool_upper_percent) %>%
  dplyr::group_by(idx) %>%
  dplyr::summarise(population_under125perday_percent = sum(population_under125perday_percent*population)/sum(population),
                   population_under2perday_percent = sum(population_under2perday_percent*population)/sum(population), 
                   secondaryschool_overall_1118_percent = sum(secondaryschool_overall_1118_percent*population)/sum(population), 
                   secondaryschool_lower_percent = sum(secondaryschool_lower_percent*population)/sum(population), 
                   secondaryschool_upper_percent = sum(secondaryschool_upper_percent*population)/sum(population)) %>%
  dplyr::rename("areaid" = idx)

# combine
merged = left_join(merged1, merged2)



# -------------- combine with others --------------

# subset to required vars 
wb = wb %>%
  dplyr::rename("areaid"=areaid_dmoss) 
wb = wb[ , which(names(wb) %in% names(merged))]

# combine with merged
wb = wb[ wb$areaid %in% shpm$areaid, ]
wb = rbind(wb, merged)

# calculate summary vars
soc = wb %>%
  dplyr::mutate(
    
    poverty_headcount = poverty_gsowb_headcount_n / population, 
    poverty_nationalbottom40 = poverty_popinnationalbottom40percent_n / population,
    poverty_extreme = poverty_gsowb_extremepoverty_headcount_n / population,
    poverty_1.25 = population_under125perday_percent,
    poverty_2.00 = population_under2perday_percent,
    
    livelihoodpercent_agriculture = mainemployment_agriculture_n / workingpop,
    livelihoodpercent_nonfarmself = mainemployment_nonfarmselfemployment_n / workingpop,
    livelihoodpercent_wagework = mainemployment_wage_n / workingpop,
    
    electricity = mainlightsource_electricity_n / households,
    sanitation_flushtoilet_indoor = sanitation_flushtoilet_indoor_n / households,
    sanitation_flushtoilet_outdoor = sanitation_flushtoilet_outdoor_n / households,
    sanitation_flushtoilet_any = sanitation_flushtoilet_any_n / households,
    water_piped_well = water_piped_well_n / households,
    
    secondaryschool_any = secondaryschool_overall_1118_percent,
    secondaryschool_lower = secondaryschool_lower_percent,
    secondaryschool_upper = secondaryschool_upper_percent
  ) %>%
  dplyr::select(areaid, population, households, poverty_headcount, poverty_nationalbottom40, poverty_extreme,
                poverty_1.25, poverty_2.00, livelihoodpercent_agriculture, livelihoodpercent_nonfarmself, livelihoodpercent_wagework,
                electricity, sanitation_flushtoilet_indoor, sanitation_flushtoilet_outdoor, sanitation_flushtoilet_any, water_piped_well,
                secondaryschool_any, secondaryschool_lower, secondaryschool_upper)

# save
write.csv(soc, "./code/viet_dengue_districts/output/covariates/Socioeconomic_WorldBank_MAPVietnam_2009.csv", row.names=FALSE)

# shpm %>% 
#   dplyr::left_join(soc) %>%
#   ggplot() + 
#   geom_sf(aes(fill=electricity), col=NA) +
#   maptheme +
#   scale_fill_gradientn(colors=viridisLite::turbo(200)) 





# ======================= PCA on socioeconomic components to identify main axis/axes of variation =========================

# identify main axes of variation through socioeconomic components (similar to other studies using factor analysis to quantify urbanisation)
socm = read.csv("./code/viet_dengue_districts/output/covariates/Socioeconomic_WorldBank_MAPVietnam_2009.csv")
soc = socm
df = soc[ soc$areaid != 70273, ]

# include "none" for flush toilet/education
df$sanitation_flushtoilet_none = (1 - df$sanitation_flushtoilet_any)
df$secondaryschool_none = (1 - df$secondaryschool_any)


# ---------- 1. Socioeconomic PCA: including livelihoods and education ------------

# run for major socioeconomic factors; can see that first component alone accounts for 56% of variance in data; first 2 components account for 78%
# looks broadly like the two components describe (1) urban-to-rural transition; (2) rural development in areas of higher agricultural livelihoods
#pca = prcomp(df[ , c(6:10, 12:13, 15:18) ], center=TRUE, scale=TRUE)
dfx = df[ , c(9, 11, 12, 16, 13:14, 17, 20) ]
names(dfx) = c("Agricultural work", "Wage work", "Electricity", "Piped water", "Indoor toilet", "Outdoor flush toilet", "Secondary school", "No flush toilet")
pca = prcomp(dfx, center=TRUE, scale=TRUE)
summary(pca)
pca

# screeplot shows only 2 have eigenvalue (sd) of >1
# knee point is probably just including 1st PC; arguably worth only keeping this one (explains )
# screeplot(pca, type = "l", npcs = 15, main = "Screeplot of the first 10 PCs")
# abline(h = 1, col="red", lty=5)

# values of the components; strongly +ve correlated to a point (development axis; very broad in middle (i.e. urban/rural transition); slightly -vely correlated at higher values i.e urban shift)
df$soc_pca_1 = -pca$x[ , 1]
df$soc_pca_2 = pca$x[ , 2]



# ---------- 2. Infrastructure PCA -------------

# very similar structuring
dfx = df[ , c(12, 16, 13:14, 20) ]
names(dfx) = c("Electricity", "Piped water", "Indoor toilet", "Outdoor flush toilet", "No flush toilet")
pca = prcomp(dfx, center=TRUE, scale=TRUE)
summary(pca)
pca

# values of the components; strongly +ve correlated to a point (development axis; very broad in middle (i.e. urban/rural transition); slightly -vely correlated at higher values i.e urban shift)
# n.b. invert soc pca 2 for interpretation (increased values correspond to increasing development)
df$infra_pca_1 = -pca$x[ , 1]
df$infra_pca_2 = pca$x[ , 2]

# save
write.csv(df[ , c("areaid", "soc_pca_1", "soc_pca_2", "infra_pca_1", "infra_pca_2") ], "./code/viet_dengue_districts/output/covariates/Socioeconomic_PCAs.csv", row.names=FALSE)



# ------------ Relationship between 2 PCAs ----------

plot(df$soc_pca_1, df$infra_pca_1)
plot(df$soc_pca_2, df$infra_pca_2)

df$socq = INLA::inla.group(df$soc_pca_1, n=10, method="quantile")
df$soc1_qx = as.integer(as.factor(df$socq))
df$soc1_qx[ df$soc1_qx %in% 1:4 ] = 1
df$soc1_qx[ df$soc1_qx %in% 5:8 ] = 2
df$soc1_qx[ df$soc1_qx %in% 9 ] = 3
df$soc1_qx[ df$soc1_qx %in% 10 ] = 4

dfx = df[ , c("areaid", "soc1_qx", "water_piped_well", "sanitation_flushtoilet_indoor", "sanitation_flushtoilet_none", "electricity") ]
names(dfx)[3:6] = c("Piped/well water", "Flush toilet (indoor)", "Flush toilet (none)", "Electricity")
dfx = reshape2::melt(dfx, id.vars=1:2)
dfx$soc1_qx[ dfx$soc1_qx == 1 ] = "1-4"
dfx$soc1_qx[ dfx$soc1_qx == 2 ] = "5-8"
dfx$soc1_qx[ dfx$soc1_qx == 3 ] = "9"
dfx$soc1_qx[ dfx$soc1_qx == 4 ] = "10"
dfx$soc1_qx = factor(dfx$soc1_qx, levels=c("1-4", "5-8", "9", "10"), ordered=TRUE)
px = ggplot(dfx) +
  geom_boxplot(aes(x=factor(soc1_qx), y=value, fill=variable), width=0.5) +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.text=element_text(size=10)) +
  ylab("Proportion of district population") +
  xlab("Decile of socioeconomic urbanisation")
# ggsave(px, file="./output/figures_vietall/PCA_Socioeconomic_decilesbyvariable.png", device="png", dpi=600, units="in", width=8, height=5, scale=0.9)



# ================ plot loadings ================

# calculate coordinates of loadings (correlation between data and that component)
# https://stats.stackexchange.com/questions/276645/arrows-of-underlying-variables-in-pca-biplot-in-r
# https://georgemdallas.wordpress.com/2013/10/30/principal-component-analysis-4-dummies-eigenvectors-eigenvalues-and-dimension-reduction/
x_c = cor(dfx, pca$x[, 1]) * 0.8 * sqrt(nrow(df)-1)
y_c = cor(dfx, pca$x[, 2]) * 0.8 * sqrt(nrow(df)-1)
coords = data.frame(xmin=rep(0, length(x_c)),
                    xmax=as.vector(x_c),
                    ymin=rep(0, length(y_c)),
                    ymax=as.vector(y_c),
                    label=names(dfx))

# plot
pca_plot = ggplot() + 
  geom_point(data=df, aes(infra_pca_1, infra_pca_2), col="grey50", alpha=0.7, size=1.2) +
  geom_segment(data=coords, aes(x=xmin, y=ymin, xend=xmax/6, yend=ymax/6), col="blue", size=1, arrow=arrow(length = unit(0.1, "inches"))) +
  theme_classic() +
  geom_text(data=coords, aes(x=xmax/6.2, y=ymax/5.2, label=label), col="darkred", size=4.5) +
  xlab("PC1 (59% of variance)") + ylab("PC2 (20% of variance)") +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=11))
ggsave(pca_plot, file="./output/figures_vietall/PCA_Infrastructure_loadings.png", device="png", dpi=600, units="in", width=7, height=5)



# ----------------- some prelim analysis of geography and structure of components ---------------

# map components
# provx = c("Dak Lak")
# # shpx = left_join(shpm, df) %>%
# #   dplyr::filter(areaprovin %in% provx)
# shpx = left_join(shpm, df)
# range_pc1 = range(df$soc_pca_1)
# range_pc2 = range(df$soc_pca_2)
# px = ggplot() + 
#   geom_sf(data=shpx, aes(fill=soc_pca_1), col=NA) + 
#   scale_fill_viridis_c(option="magma", direction=-1, name="PC1", limits=range_pc1) + maptheme +
#   ggtitle("PC1 (rural-to-urban)")
# py = ggplot() + 
#   geom_sf(data=shpx, aes(fill=soc_pca_2), col=NA) + 
#   scale_fill_viridis_c(option="magma", direction=-1, name="PC2", limits=range_pc2) + maptheme +
#   ggtitle("PC2 (rural development)")
# gridExtra::grid.arrange(px, py, nrow=1)
# 
# # how are these components correlated to urbanisation i.e. impervious and population density
# urb = read.csv("./code/viet_dengue_districts/output/covariates/LandUseCover_ESA_MergedSHP.csv") %>% dplyr::filter(year == 2009)
# pop = read.csv("./code/viet_dengue_districts/output/covariates/Population_GPW_WP.csv") %>% dplyr::filter(year == 2009)
# dfx = left_join(df, urb)
# dfx = left_join(dfx, pop)
# 
# ggplot(dfx) + geom_point(aes(urban_pw, soc_pca_1))
# ggplot(dfx) + geom_point(aes(log(popdens_pw), soc_pca_1))
# ggplot(dfx) + geom_point(aes(urban_pw, soc_pca_2))
# ggplot(dfx) + geom_point(aes(log(popdens_pw), soc_pca_2))
# 
# cor.test(log(dfx$popdens_pw), dfx$soc_pca_1)
# cor.test(log(dfx$popdens_pw), dfx$soc_pca_2)
# cor.test(dfx$urban_pw, dfx$soc_pca_1)
# cor.test(dfx$urban_pw, dfx$soc_pca_2)
# 
# # pc1 highly correlated to population density; pc2 effectively uncorrelated
# # pc2 peaks at intermediate levels of pc1 and declines at lower values because describes rural-type development
# # quadrants describe, broadly, high poverty and low infrastructure access (bottom left); highly rural areas with good infrastructure (top left)
# # highly developed rural and peri-urban areas (top right), increasingly urbanised areas (bottom right)
# # although quadrants not a clear descriptor
# ggplot(dfx) + geom_point(aes(soc_pca_1, soc_pca_2, col=log(popdens_pw)), size=2) + theme_minimal() + 
#   scale_color_viridis_c(option="viridis") + geom_hline(yintercept=0, lty=2) + geom_vline(xintercept=0, lty=2) +
#   xlab("PC1 (urbanisation)") + ylab("PC2 (rural development)")
# ggplot(dfx) + geom_point(aes(soc_pca_1, soc_pca_2, col=poverty_nationalbottom40 ), size=2) + theme_minimal() + 
#   scale_color_viridis_c(option="viridis") + geom_hline(yintercept=0, lty=2) + geom_vline(xintercept=0, lty=2) +
#   xlab("PC1 (urbanisation)") + ylab("PC2 (rural development)")
# 
# # ggplot(dfx) +  geom_point(aes(soc_pca_2, livelihoodpercent_agriculture)) + theme_minimal()
# # ggplot(dfx) +  geom_point(aes(soc_pca_2, secondaryschool_lower)) + theme_minimal()
# # ggplot(dfx) +  geom_point(aes(soc_pca_2, secondaryschool_upper)) + theme_minimal()
# ggplot(dfx) +  geom_point(aes(soc_pca_1, secondaryschool_upper)) + theme_minimal()
# 
# # compare to regions and provinces
# rr = read.csv("./data/shapefiles/regions_lookup.csv") %>%
#   dplyr::rename("areaprovin"=provincena)
# regx = shpm %>%
#   st_drop_geometry() %>%
#   left_join(rr)
# dfx = left_join(dfx, regx[ , c("areaid", "areaprovin", "region")])
# 
# # some provinces
# kp = c("Ha Noi", "Da Nang", "Can Tho", "Dak Lak", "Khanh Hoa", "Dong Nai", "TP. Ho Chi Minh", "Hai Phong", "Thua Thien - Hue")
# ggplot(dfx[ dfx$areaprovin %in% kp, ]) + 
#   geom_point(aes(soc_pca_1, soc_pca_2, col=factor(areaprovin)), size=2) + theme_minimal() + 
#   geom_hline(yintercept=0, lty=2) + geom_vline(xintercept=0, lty=2) +
#   stat_ellipse(aes(soc_pca_1, soc_pca_2, fill=factor(areaprovin)), geom="polygon", alpha=0.1) +
#   stat_ellipse(aes(soc_pca_1, soc_pca_2, col=factor(areaprovin)), alpha=0.5) +
#   xlab("PC1 (urbanisation)") + ylab("PC2 (rural development)") +
#   facet_wrap(~areaprovin) +
#   theme(legend.position="bottom")
# 
# # regeions
# dfx$region = factor(dfx$region, levels=unique(dfx$region), ordered=TRUE)
# ggplot(dfx[ !is.na(dfx$region), ]) + 
#   geom_point(aes(soc_pca_1, soc_pca_2, col=factor(region)), size=1.5) + theme_minimal() + 
#   geom_hline(yintercept=0, lty=2) + geom_vline(xintercept=0, lty=2) +
#   stat_ellipse(aes(soc_pca_1, soc_pca_2, fill=factor(region)), geom="polygon", alpha=0.2) +
#   #stat_ellipse(aes(soc_pca_1, soc_pca_2, col=factor(region)), alpha=0.5) +
#   xlab("PC1 (urbanisation)") + ylab("PC2 (rural development)") +
#   facet_wrap(~region, nrow=2) +
#   theme(legend.position="none") + 
#   scale_color_viridis_d(begin=0, end=0.8, option="magma") +
#   scale_fill_viridis_d(begin=0, end=0.8, option="magma")
# 
# # regions with city municipalities highlighted: nice perspective on where different cities fall
# dfx$region = factor(dfx$region, levels=unique(dfx$region), ordered=TRUE)
# ggplot(dfx[ !is.na(dfx$region), ]) + 
#   stat_ellipse(data=dfx[ dfx$areaprovin %in% c("Da Nang", "Ha Noi", "Can Tho", "TP. Ho Chi Minh", "Hai Phong"), ], aes(soc_pca_1, soc_pca_2, fill=factor(areaprovin)), geom="polygon", alpha=0.3) +
#   #stat_ellipse(data=dfx[ dfx$areaprovin %in% c("Da Nang", "Ha Noi", "Can Tho", "TP. Ho Chi Minh", "Hai Phong"), ], aes(soc_pca_1, soc_pca_2, col=factor(region)), alpha=0.5) +
#   geom_point(aes(soc_pca_1, soc_pca_2), size=2.4, col="grey30", alpha=0.5) + 
#   geom_point(data=dfx[ dfx$areaprovin %in% c("Da Nang", "Ha Noi", "Can Tho", "TP. Ho Chi Minh", "Hai Phong"), ], aes(soc_pca_1, soc_pca_2), col="darkred", alpha=0.5, size=2.8) +
#   theme_minimal() + 
#   geom_hline(yintercept=0, lty=2) + geom_vline(xintercept=0, lty=2) +
#   xlab("PC1 (urbanisation)") + ylab("PC2 (rural development)") +
#   facet_wrap(~region, nrow=2) +
#   theme(legend.position="bottom",
#         axis.title=element_text(size=13),
#         strip.text=element_text(size=13),
#         legend.title = element_text(size=12),
#         legend.text = element_text(size=12)) + 
#   scale_color_viridis_d(begin=0, end=0.8, option="magma") +
#   scale_fill_viridis_d(begin=0, end=0.8, option="viridis", name="Municipality")
# 
# # can see clear geographical structuring in socioeconomic development, with lower urbanisation and lower rural development in peripheries;
# # whereas economic centres of the country, southeast and red river delta, are highly urbanised
# # although, in the south rural development at a given level of urban-to-rural transition is generally lower (mainly consequence of education status, lower piped water access)
# # PC1: "urbanisation"/"urbal-to-rural transition" (decline in agri livelihoods and increase in nonfarm/wage work, increased electricity, water access, sanitation access and education)
# # would hypothesise 1 to be positively correlated to dengue due to increased population density and built environment development (/and also likely indistinguishable from effect of population density)
# # PC2: "rural development" (correlated to agri livelihoods with increased electricity access, outdoor flush toilet and improved education, but still with lower sanitation access than cities (lower piped water, indoor toilet)
# # would expect a negative correlation with dengue accounting for PC1/urbanisation; i.e. rural areas with lower sanitation systems access
# # might also expect the most rapid urban expansion to be occurring in areas that are middle of PC1 and high PC2, so these could confound any effects of built env
# 
# ggplot(dfx[ dfx$areaprovin %in% kp, ]) + 
#   geom_boxplot(aes(x= areaprovin, y=secondaryschool_lower, group=areaprovin)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle=90))
# ggplot(dfx[ dfx$areaprovin %in% kp, ]) + 
#   geom_boxplot(aes(x= areaprovin, y=sanitation_flushtoilet_none, group=areaprovin)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle=90))
# 
# 
# 
# # -------------- how are PCs correlated to urbanisation dynamics? ----------------
# 
# urbdyn = read.csv("./code/viet_dengue_districts/output/covariates/LandUseDynamics_MergedSHP_Feb2021.csv")
# dfu = left_join(dfx, urbdyn[ , c("areaid", "urbanexp_3yr_d", "urbanexp_5yr_d", "urbanexp_10yr_d")])
# 
# # generally urban expansion rates are fastest in areas with med-high scores on the "urban-to-rural transition" PC
# # and conversely, generally fastest in areas with medium-to-low scores on the "rural development" PC
# # so an important question here is, accounting for socioeconomics, does rapid urban expansion increase dengue outbreaks?
# # and/or does socioeconomic development modify the effect of rapid urban expansion on dengue outbreaks
# # i.e. does rapid expansion in less-developed areas with poorer infrastructure have a greater positive influence on risk
# # than rapid expansion in socioeconomically better-developed areas with better infrastructure and housing quality
# ggplot(dfu) + 
#   geom_point(aes(soc_pca_1, urbanexp_3yr_d, col=region))
# ggplot(dfu) + 
#   geom_point(aes(soc_pca_2, urbanexp_3yr_d, col=region)) 
# ggplot(dfu) + 
#   geom_point(aes(soc_pca_1, urbanexp_10yr_d, col=region))
# ggplot(dfu) + 
#   geom_point(aes(soc_pca_2, urbanexp_10yr_d, col=region))
# 
