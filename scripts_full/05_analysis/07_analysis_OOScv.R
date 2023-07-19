

# ================== Measuring predictive influence of covariates using block cross-validation ===================

# Models were fitted using 5-fold cross validation under 3 block designs.
# Cross validation scripts are stored in the "scripts/modelling" folder and were designed to run on an HPC cluster.
# These scripts output a final file for each block design called "models_completed", which contains summary predictive statistics (MAE/RMSE)
# The first section of this script reads these to produce visualisations.
# The later sections of the scripts read in the full model outputs to calculate more granular summaries;
# those outputs are saved in the folder "output/model_outputs/regional_OOS".

library(dplyr)
library(sf)
library(ggplot2)
library(vroom)
library(raster)

PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

source("./scripts/04_modelling/00_inla_setup_functions_r4.R")
source("./scripts/04_modelling/00_plot_themes.R")



# ------------- MS Figure: visualise overall error metrics --------------

# final output files from cross-validation 
r1 = read.csv("./output/model_outputs/spOOS/models_completed.csv") %>% 
  dplyr::mutate(holdout_type ="Spatial") %>% 
  dplyr::filter(!duplicated(model_identifier)) %>%
  dplyr::filter(unique_id != "dummy")

r2 = read.csv("./output/model_outputs/stempOOS/models_completed.csv") %>% 
  dplyr::mutate(holdout_type = "Spatiotemporal") %>% 
  dplyr::filter(!duplicated(model_identifier)) %>%
  dplyr::filter(unique_id != "dummy") %>%
  dplyr::filter(!unique_id %in% c(733759, 842257)) %>%
  dplyr::filter(n_models_fitted == 5)

r3 = read.csv("./output/model_outputs/tempOOS/models_completed.csv") %>% 
  dplyr::mutate(holdout_type = "Seasonal") %>% 
  dplyr::filter(!duplicated(model_identifier)) %>%
  dplyr::filter(unique_id != "dummy") %>%
  dplyr::filter(!unique_id %in% c(652377, 933064))

# combine
oos_df = do.call(rbind.data.frame, list(r1, r2, r3)) %>%
  dplyr::mutate(mae_oos = as.numeric(mae_oos),
                rmse_oos = as.numeric(rmse_oos))

oos_df = oos_df %>%
  dplyr::left_join(
    data.frame(
      candidate = c("baseline", "flushany", "urbanexp10", "water", "tmean_coolest", "tmean", "traffic", "spei1",
                "urban", "gravity", "full", "spei6"),
      model2 = c("Baseline (all)", "Hygeinic toilet", "Urb. exp (10y)", "Piped water", "Tmean coolest", "Tmean (1m)", 
                 "Traffic", "SPEI-1 (1m)", "Built-up", "Gravity", "Full model", "SPEI-6 (4m)"),
      vartype = c("Baseline (all)", "Infrastructure", "Urbanisation", "Infrastructure", "Climate", "Climate", "Mobility", "Climate",
               "Urbanisation", "Mobility", "Full", "Climate"),
      vartype2 = c("Baseline", "Socio", "Socio", "Socio", "Clim", "Clim", "Socio", "Clim", "Socio", "Socio", "Full", "Clim")
    )
  )

# calculate difference from baseline
oos_df = oos_df %>%
  dplyr::group_by(unique_id) %>%
  dplyr::mutate(delta_mae = mae_oos - mae_oos[ candidate == "full"],
                delta_rmse = rmse_oos - rmse_oos[ candidate == "full"])

sm_summary = oos_df %>%
  dplyr::group_by(holdout_type, model2) %>%
  dplyr::summarise(vartype = head(vartype, 1), 
                   deltamae_mean = mean(delta_mae),
                   deltamae_se = plotrix::std.error(delta_mae),
                   deltarmse_mean = mean(delta_rmse),
                   deltarmse_se = plotrix::std.error(delta_rmse))

# factor order and colours

col_clim = viridis::viridis(200)[40]
col_socio= viridis::viridis(200)[105]

fac_order = sm_summary %>%
  dplyr::filter(holdout_type == "Spatiotemporal") %>%
  dplyr::arrange(deltamae_mean)
#fac_order = c(fac_order$model2[ fac_order$model2 != "Full model" ], "Full model")
fac_order = fac_order$model2[ fac_order$model2 != "Full model" ]

oos_df$holdout_type = factor(oos_df$holdout_type, levels=c("Spatial", "Spatiotemporal", "Seasonal"), ordered=TRUE)
sm_summary$holdout_type = factor(sm_summary$holdout_type, levels=c("Spatial", "Spatiotemporal", "Seasonal"), ordered=TRUE)

p1 = oos_df %>%
  dplyr::filter(model2 != "Full model") %>%
  dplyr::mutate(model2 = factor(model2, levels=fac_order, ordered=TRUE)) %>%
  ggplot() + 
  geom_point(aes(model2, delta_mae, group=unique_id, color=vartype2), size=5, pch=16, alpha=0.15) +
  geom_point(data=sm_summary[ sm_summary$model2 != "Full model", ], aes(model2, deltamae_mean), pch=18, color="black", alpha=1, size=3.25) +
  geom_linerange(data=sm_summary[ sm_summary$model2 != "Full model", ], aes(model2, ymin=deltamae_mean-1.96*deltamae_se, ymax=deltamae_mean+1.96*deltamae_se), color="black", alpha=1, size=0.8) +
  #geom_rect(aes(xmin=xmin1, xmax=xmax1, ymin=ymin1, ymax=0), alpha=0.1, fill="green") + 
  geom_hline(yintercept=0, lty=2) + 
  theme_bw() + 
  ylab(expression(paste(Delta, "MAE (out-of-sample) relative to full model"))) +
  #ylab("Improvement in MAE (out-of-sample)") +
  xlab("Covariate excluded") + 
  theme(strip.background = element_blank(), 
        panel.grid.minor =  element_blank(),
        strip.text=element_text(size=16),
        axis.text.y = element_text(size=12, color="black"), 
        axis.text.x = element_text(size=13), 
        axis.title = element_text(size=15),
        legend.position="none") + 
  facet_wrap(~holdout_type, scales="free_x") +
  coord_flip() +
  scale_color_manual(
    values=c("Baseline"="black", "Clim"=col_clim, "Socio"=col_socio)
  )

p1 = gridExtra::grid.arrange(p1)
p1 = ggpubr::as_ggplot(p1)  +
  cowplot::draw_plot_label(label = c("a", "b", "c"), 
                           fontface = "bold", size = 22.5, 
                           x = c(0.121, 0.415, 0.705), y = c(0.91, 0.91, 0.91))

ggsave(p1, file="./output/figures/Figure4_PredictiveHoldouts.png", device="png", units="in", dpi=150, width=12, height=4)



# 
# # ================ RMSE =================
# 
# oos_df %>%
#   dplyr::filter(model2 != "Full model") %>%
#   dplyr::mutate(model2 = factor(model2, levels=fac_order, ordered=TRUE)) %>%
#   ggplot() + 
#   geom_point(aes(model2, delta_rmse, group=unique_id, color=vartype2), size=5, pch=16, alpha=0.15) +
#   geom_point(data=sm_summary[ sm_summary$model2 != "Full model", ], aes(model2, deltarmse_mean), pch=18, color="black", alpha=1, size=3.25) +
#   geom_linerange(data=sm_summary[ sm_summary$model2 != "Full model", ], aes(model2, ymin=deltarmse_mean-1.96*deltarmse_se, ymax=deltarmse_mean+1.96*deltarmse_se), color="black", alpha=1, size=0.8) +
#   #geom_rect(aes(xmin=xmin1, xmax=xmax1, ymin=ymin1, ymax=0), alpha=0.1, fill="green") + 
#   geom_hline(yintercept=0, lty=2) + 
#   theme_bw() + 
#   ylab(expression(paste(Delta, "RMSE (out-of-sample) relative to full model"))) +
#   #ylab("Improvement in MAE (out-of-sample)") +
#   xlab("Covariate excluded") + 
#   theme(strip.background = element_blank(), 
#         panel.grid.minor =  element_blank(),
#         strip.text=element_text(size=16),
#         axis.text.y = element_text(size=12, color="black"), 
#         axis.text.x = element_text(size=13), 
#         axis.title = element_text(size=15),
#         legend.position="none") + 
#   facet_wrap(~holdout_type, scales="free_x") +
#   coord_flip() +
#   scale_color_manual(
#     values=c("Baseline"="black", "Clim"=col_clim, "Socio"=col_socio)
#   )




# ---------------- Supplement: mapping changes from baseline to full model ---------------

# all combined model names
comb_names = paste("_", oos_df$model_identifier, "_", oos_df$unique_id, "_", sep="")

locs_sp = list.files("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/output/model_outputs/results_23/spOOS_1/models/", full.names=TRUE)
locs_st = list.files("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/output/model_outputs/results_23/stempOOS_1/models/", full.names=TRUE)
locs_temp = list.files("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/output/model_outputs/results_23/tempOOS_1/models/", full.names=TRUE)
locs = c(locs_sp, locs_st, locs_temp)


# calculate MAE and RMSE
calcSpatialOOS = function(x){
  
  print(x)
  loc = unique(locs[ grep(x, locs) ])
  
  csv = do.call(rbind.data.frame, lapply(loc, read.csv)) %>%
    dplyr::mutate(pred = exp(logpop + mean),
                  resid = pred - cases) %>%
    dplyr::group_by(areaid) %>%
    dplyr::summarise(
      model = head(model, 1), holdout_id = head(holdout_id, 1), holdout_type = head(holdout_type, 1),
      mae = mean(abs(resid), na.rm=TRUE), 
      rmse = sqrt( mean( resid^2, na.rm=TRUE )),
      n_mods_included = length(loc)
    )
  return(csv)
}

all_oos = do.call(rbind.data.frame, lapply(comb_names, calcSpatialOOS))

# add baseline mae
all_oos2 = left_join(
  all_oos, 
  all_oos %>% dplyr::filter(model=="baseline") %>% dplyr::select(areaid, holdout_id, holdout_type, mae) %>% dplyr::rename("mae_bl"=mae)
)

# plot and visualise change
dm = all_oos2 %>%
  dplyr::filter(model == "full") %>%
  dplyr::mutate(delta_mae = mae - mae_bl) %>%
  dplyr::group_by(holdout_type, areaid) %>%
  dplyr::summarise(delta_mae = mean(delta_mae, na.rm=TRUE))

# dm %>%
#   ggplot() + 
#   geom_histogram(aes(delta_mae), bins=40) +
#   facet_wrap(~holdout_type)

# broad categories
p1 = shp %>%
  left_join(
    dm
  ) %>%
  dplyr::mutate(
    mm = "< -2",
    mm = replace(mm, delta_mae > -2, "-2 to -1"),
    mm = replace(mm, delta_mae > -1, "-1 to 0"),
    mm = replace(mm, delta_mae > 0, "0 to +1"),
    mm = replace(mm, delta_mae > 1, "+1 to +2"),
    mm = replace(mm, delta_mae > 2, "> +2")
  ) %>%
  dplyr::mutate(
    mm = factor(mm, levels=c("< -2", "-2 to -1", "-1 to 0", "0 to +1", "+1 to +2", "> +2"), ordered=TRUE)
  ) %>%
  dplyr::filter(!is.na(holdout_type)) %>%
  dplyr::mutate(
    holdout_type = replace(holdout_type, holdout_type=="tempOOS", "Seasonal"),
    holdout_type = replace(holdout_type, holdout_type=="spOOS", "Spatial"),
    holdout_type = replace(holdout_type, holdout_type=="stempOOS", "Spatiotemporal"),
    holdout_type = factor(holdout_type, 
                          levels=c("Spatial", "Spatiotemporal", "Seasonal"), 
                          ordered=TRUE)
  ) %>%
  ggplot() + 
  geom_sf(aes(fill=mm), col=NA) + 
  geom_sf(data=shp_prov %>% st_crop(shp[ shp$areaid %in% dm$areaid, ]), fill=NA, color="grey60", alpha=0.5, size=0.25) +
  maptheme + 
  scale_fill_brewer(palette = "BrBG", direction=-1, name=expression(paste(Delta, "MAE"))) + 
  #theme(legend.position=c(0.8, 0.5)) +
  facet_wrap(~holdout_type, ncol=3) +
  theme(legend.text = element_text(size=12), legend.title=element_text(size=13))

ggsave(p1, file="./output/figures/SuppFigure_OOSResultsSpatial.png", device="png", units="in", width=11, height=5, dpi=150)




# ------------ Supp. Figure: regional contributions -------------

# calculate MAE and RMSE
calcRegionalOOS = function(x){
  
  print(x)
  loc = unique(locs[ grep(x, locs) ])
  
  csv = do.call(rbind.data.frame, lapply(loc, read.csv)) %>%
    dplyr::mutate(pred = exp(logpop + mean),
                  resid = pred - cases) %>%
    dplyr::left_join(shp[ , c("areaid", "region2")] %>% st_drop_geometry()) %>%
    dplyr::mutate(region = ifelse(region2 %in% c("Central Highlands", "South", "South Central Coast"), "Southern", "Northern")) %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(
      model = head(model, 1), holdout_id = head(holdout_id, 1), holdout_type = head(holdout_type, 1),
      mae = mean(abs(resid), na.rm=TRUE), 
      rmse = sqrt( mean( resid^2, na.rm=TRUE )),
      n_mods_included = length(loc)
    )
  return(csv)
}

r_oos = do.call(rbind.data.frame, lapply(comb_names, calcRegionalOOS))

oos_df = r_oos %>%
  dplyr::left_join(
    data.frame(
      model = c("baseline", "flushany", "urbanexp10", "water", "tmean_coolest", "tmean", "traffic", "spei1",
                    "urban", "gravity", "full", "spei6"),
      model2 = c("Baseline (all)", "Hygeinic toilet", "Urb. exp (10y)", "Piped water", "Tmean coolest", "Tmean (1m)", 
                 "Traffic", "SPEI-1 (1m)", "Built-up", "Gravity", "Full model", "SPEI-6 (4m)"),
      vartype = c("Baseline (all)", "Infrastructure", "Urbanisation", "Infrastructure", "Climate", "Climate", "Mobility", "Climate",
                  "Urbanisation", "Mobility", "Full", "Climate"),
      vartype2 = c("Baseline", "Socio", "Socio", "Socio", "Clim", "Clim", "Socio", "Clim", "Socio", "Socio", "Full", "Clim")
    )
  )

# calculate difference from baseline
oos_df = oos_df %>%
  dplyr::group_by(region, holdout_id) %>%
  dplyr::mutate(holdout_type = head(holdout_type, 1), 
                delta_mae = mae - mae[ model == "full"],
                delta_rmse = rmse - rmse[ model == "full"])

oos_df$holdout_type[ oos_df$holdout_type == "spOOS" ] = "Spatial"
oos_df$holdout_type[ oos_df$holdout_type == "stempOOS" ] = "Spatiotemporal"
oos_df$holdout_type[ oos_df$holdout_type == "tempOOS" ] = "Seasonal"

sm_summary = oos_df %>%
  dplyr::group_by(region, holdout_type, model2) %>%
  dplyr::summarise(vartype = head(vartype, 1), 
                   deltamae_mean = mean(delta_mae),
                   deltamae_se = plotrix::std.error(delta_mae),
                   deltarmse_mean = mean(delta_rmse),
                   deltarmse_se = plotrix::std.error(delta_rmse))

# factor order and colours

col_clim = viridis::viridis(200)[40]
col_socio= viridis::viridis(200)[105]



oos_df$holdout_type = factor(oos_df$holdout_type, levels=c("Spatial", "Spatiotemporal", "Seasonal"), ordered=TRUE)
sm_summary$holdout_type = factor(sm_summary$holdout_type, levels=c("Spatial", "Spatiotemporal", "Seasonal"), ordered=TRUE)



# plot for north

sm_n = sm_summary %>% dplyr::filter(region == "Northern")
od_n = oos_df %>% dplyr::filter(region == "Northern")

fac_order = sm_n %>%
  dplyr::filter(holdout_type == "Spatiotemporal") %>%
  dplyr::arrange(deltamae_mean)
#fac_order = c(fac_order$model2[ fac_order$model2 != "Full model" ], "Full model")
fac_order = fac_order$model2[ fac_order$model2 != "Full model" ]

p1 = od_n %>%
  dplyr::filter(model2 != "Full model") %>%
  dplyr::mutate(model2 = factor(model2, levels=fac_order, ordered=TRUE)) %>%
  ggplot() + 
  geom_point(aes(model2, delta_mae, group=holdout_id, color=vartype2), size=5, pch=16, alpha=0.15) +
  geom_point(data=sm_n[ sm_n$model2 != "Full model", ], aes(model2, deltamae_mean), pch=18, color="black", alpha=1, size=3.25) +
  geom_linerange(data=sm_n[ sm_n$model2 != "Full model", ], aes(model2, ymin=deltamae_mean-1.96*deltamae_se, ymax=deltamae_mean+1.96*deltamae_se), color="black", alpha=1, size=0.8) +
  #geom_rect(aes(xmin=xmin1, xmax=xmax1, ymin=ymin1, ymax=0), alpha=0.1, fill="green") + 
  geom_hline(yintercept=0, lty=2) + 
  theme_bw() + 
  ylab(expression(paste(Delta, "MAE (out-of-sample) relative to full model"))) +
  #ylab("Improvement in MAE (out-of-sample)") +
  xlab("Covariate excluded") + 
  theme(strip.background = element_blank(), 
        panel.grid.minor =  element_blank(),
        strip.text=element_text(size=15),
        axis.text.y = element_text(size=12, color="black"), 
        axis.text.x = element_text(size=13), 
        axis.title = element_text(size=15),
        plot.title = element_text(size=16, hjust=0.5),
        legend.position="none") + 
  facet_wrap(~holdout_type, scales="free_x") +
  coord_flip() +
  scale_color_manual(
    values=c("Baseline"="black", "Clim"=col_clim, "Socio"=col_socio)
  ) +
  ggtitle("Northern Vietnam\n(North Central, Red River Delta, Northeast, Northwest)")


# plot for south

sm_n = sm_summary %>% dplyr::filter(region == "Southern")
od_n = oos_df %>% dplyr::filter(region == "Southern")

fac_order = sm_n %>%
  dplyr::filter(holdout_type == "Spatiotemporal") %>%
  dplyr::arrange(deltamae_mean)
#fac_order = c(fac_order$model2[ fac_order$model2 != "Full model" ], "Full model")
fac_order = fac_order$model2[ fac_order$model2 != "Full model" ]

p2 = od_n %>%
  dplyr::filter(model2 != "Full model") %>%
  dplyr::mutate(model2 = factor(model2, levels=fac_order, ordered=TRUE)) %>%
  ggplot() + 
  geom_point(aes(model2, delta_mae, group=holdout_id, color=vartype2), size=5, pch=16, alpha=0.15) +
  geom_point(data=sm_n[ sm_n$model2 != "Full model", ], aes(model2, deltamae_mean), pch=18, color="black", alpha=1, size=3.25) +
  geom_linerange(data=sm_n[ sm_n$model2 != "Full model", ], aes(model2, ymin=deltamae_mean-1.96*deltamae_se, ymax=deltamae_mean+1.96*deltamae_se), color="black", alpha=1, size=0.8) +
  #geom_rect(aes(xmin=xmin1, xmax=xmax1, ymin=ymin1, ymax=0), alpha=0.1, fill="green") + 
  geom_hline(yintercept=0, lty=2) + 
  theme_bw() + 
  ylab(expression(paste(Delta, "MAE (out-of-sample) relative to full model"))) +
  #ylab("Improvement in MAE (out-of-sample)") +
  xlab("Covariate excluded") + 
  theme(strip.background = element_blank(), 
        panel.grid.minor =  element_blank(),
        strip.text=element_text(size=15),
        axis.text.y = element_text(size=12, color="black"), 
        axis.text.x = element_text(size=13), 
        axis.title = element_text(size=15),
        plot.title = element_text(size=16, hjust=0.5),
        legend.position="none") + 
  facet_wrap(~holdout_type, scales="free_x") +
  coord_flip() +
  scale_color_manual(
    values=c("Baseline"="black", "Clim"=col_clim, "Socio"=col_socio)
  ) +
  ggtitle("Southern Vietnam\n(Mekong Delta, Southeast, South Central Coast, Central Highlands)")



p1 = gridExtra::grid.arrange(p1, p2, nrow=2)
# p1 = ggpubr::as_ggplot(p1)  +
#   cowplot::draw_plot_label(label = c("a", "b", "c"), 
#                            fontface = "bold", size = 22.5, 
#                            x = c(0.121, 0.415, 0.705), y = c(0.91, 0.91, 0.91))
# 
ggsave(p1, file="./output/figures/SuppFigure_PredictiveHoldouts_Regional.png", device="png", units="in", dpi=600, width=12, height=8)


