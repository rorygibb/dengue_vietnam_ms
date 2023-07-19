
# ================== Results for regional sub-models ================

library(dplyr)
library(sf)
library(ggplot2)
library(raster)

PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

source("./scripts/04_modelling/00_inla_setup_functions_r4.R")
source("./scripts/04_modelling/00_plot_themes.R")

load("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/output/model_outputs/results_23/fullmodel_urbansens/models/05_urbansens_nb_model_1.R")
mod_05 = mod_i

load("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/output/model_outputs/results_23/fullmodel_urbansens/models/07_urbansens_nb_model_1.R")
mod_07 = mod_i

load("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/output/model_outputs/results_23/fullmodel_urbansens/models/09_urbansens_nb_model_1.R")
mod_09 = mod_i



# =============== Extract ranefs and fixefs =================

ranefs = data.frame()
byms = data.frame()
fixefs = data.frame()

for(uu in c(0.5, 0.7, 0.9)){
  
  if(uu == 0.5){ mod_i = mod_05 }
  if(uu == 0.7){ mod_i = mod_07 }
  if(uu == 0.9){ mod_i = mod_09 }
  
  # ranefs
  rr = mod_i$summary.random
  rf = lapply(rr[3:length(rr)], extractRandomINLA, effect_name = "x", transform=TRUE)
  for(i in 1:length(rf)){
    rf[[i]]$effect = names(rr[ 3:length(rr)][i] )
  }
  rf = do.call(rbind.data.frame, rf) %>%
    dplyr::mutate(urbsens = uu)
  
  # byms
  byx = extractRandomINLA(mod_i$summary.random$polyid, effect_name="bym", model_is_bym = TRUE) %>%
    dplyr::mutate(urbsens = uu)
  
  # fixefs
  fx = extractFixedINLA(mod_i, model_name=region, transform=FALSE) %>%
    dplyr::mutate(urbsens = uu)
  
  # combine
  ranefs = rbind(ranefs, rf)
  fixefs = rbind(fixefs, fx)
  byms = rbind(byms, byx)
}

# view
ggplot(ranefs) +
  geom_line(aes(value, mean, group=urbsens, col=factor(urbsens))) +
  geom_ribbon(aes(value, ymin=lower, ymax=upper, group=urbsens, fill=factor(urbsens)), alpha=0.2) +
  theme_bw() +
  geom_hline(yintercept=1, lty=2) +
  facet_wrap(~effect, scales="free")

fixefs %>%
  dplyr::filter(param != "Intercept") %>%
  ggplot() +
  geom_point(aes(param, mean, group=urbsens, col=factor(urbsens)), position=position_dodge(width=0.5), size=2) +
  geom_linerange(aes(param, ymin=lower, ymax=upper, group=urbsens, col=factor(urbsens)), position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0, lty=2) +
  theme_bw()





# extract ranefs and create standardised plots
rf = ranefs
effs = unique(rf$effect)
plots = vector("list", length(effs))
plotRF = function(x){
  
  rr=rf[ rf$effect == effs[x], ]
  
  rr = left_join(
    rr,
    data.frame(
      effect = effs,
      effname = c("Hygeinic toilet access", "Piped water access", "Tmean (1-month lag)", "SPEI-1 (1-month lag)", "SPEI-6 (5-month lag)"),
      xlab = c("Prop. households", "Prop. households", "Tmean (°C)", "SPEI-1", "SPEI-6"),
      type = c("Infrastructure", "Infrastructure", "Climate", "Climate", "Climate")
    )
  )
  
  if(rr$effect[1] == "tmean_1m_g"){
    rr$lower[ rr$value < 12.5 ] = NA
    rr$upper[ rr$value < 12.5 ] = NA
  }
  
  px = ggplot(rr) + 
    geom_ribbon(aes(value, ymin=lower, ymax=upper, group=urbsens), fill="grey80", alpha=0.3) +
    geom_line(aes(value, mean, group=urbsens, col=factor(urbsens)), size=0.9) +
    geom_hline(yintercept=1, lty=2) +
    facet_wrap(~effname) +
    theme_classic() +
    scale_color_viridis_d(end=0.75, name="Urban holdout\nthreshold") + 
    ylab("Relative risk") + xlab(rr$xlab[1]) +
    theme(axis.text = element_text(size=11),
          axis.title = element_text(size=12), 
          strip.background = element_blank(),
          strip.text = element_text(size=13),
          legend.position="none") 
  
  # add vline for 0 in SPEI demarcating wet/dry
  if(grepl("spei", effs[x])){
    px = px + geom_vline(xintercept=0, lty=1, size=0.25, col="grey70", alpha=0.4)
  }
  
  plots[[x]] = px
}
prf = lapply(1:length(effs), plotRF)


# fixed effects plot
pfx = fixefs %>%
  dplyr::filter(param != "Intercept") %>%
  dplyr::left_join(
    data.frame(
      param=c("gravityf_log", "urban_s", "urbanexp10_log", "traffic_kmperinhab_log", "tmean_coolestmonth_s"),
      paramname = c("Gravity\n(log)", "Built-up\nland", "Urban\nexpansion\nrate", "Traffic\nper inhab\n(log)", "Tmean\n(coolest\nmonth)"),
      type = c("Mobility", "Urbanisation", "Urbanisation", "Mobility", "Climate")
    )
  ) %>%
  dplyr::mutate(paramname = factor(paramname, levels=c("Tmean\n(coolest\nmonth)", "Gravity\n(log)", "Traffic\nper inhab\n(log)",  "Built-up\nland", "Urban\nexpansion\nrate"), ordered=TRUE)) %>%
  ggplot() +
  geom_point(aes(paramname, mean, group=urbsens, col=factor(urbsens)), position=position_dodge(width=0.5), size=2) +
  geom_linerange(aes(paramname, ymin=lower, ymax=upper, group=urbsens, col=factor(urbsens)), position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0, lty=2) +
  theme_classic() +
  scale_color_viridis_d(end=0.75, name="Urban holdout\nthreshold") + 
  ylab("Linear fixed effect\n(posterior median + 95% CI)") + xlab("") + 
  xlab("") + 
  #scale_color_viridis_d(begin=0.2, end=0.6) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0, size=11),
        axis.text.y = element_text(size=11),
        axis.title.y = element_text(size=12), 
        legend.position=c(0.8, 0.8), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))

# combine and save
pc1 = gridExtra::grid.arrange(pfx, prf[[1]], prf[[2]] + ylab(""), nrow=1, widths=c(1.1, 0.9, 0.9)) 
pc2 = gridExtra::grid.arrange(prf[[3]], prf[[4]] + ylab(""), prf[[5]] + ylab(""), nrow=1)
pcomb =  gridExtra::grid.arrange(pc1, pc2, nrow=2)                       
ggsave(pcomb, file="./output/figures/SuppFigure_UrbanHoldoutSensitivityTest.png", dpi=600, device="png", units="in", width=11.5, height=7, scale=0.9)

