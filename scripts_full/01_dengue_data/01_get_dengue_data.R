

# ====================== Access dengue surveillance data ======================

# working directory
PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

# dependencies
pacman::p_load("dplyr", "raster", "rgdal", "sf", "ecmwfr", "stringr", "ggplot2", "lubridate", "magrittr", "vroom")

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))




# --------------- dengue dataset -------------------

# Vietnam dengue timeseries 1998-2021
# generated from raw surveillance Excel sheets in separate private repository:
# https://github.com/rorygibb/viet_dengue_data/tree/main/output/dengue_timeseries

# N.B. for publication, data are provided for 4 provinces (Ha Noi, Khanh Hoa, Dak Lak and Dong Nai)
# contact details for requesting the full dataset are provided in the MS

read.csv("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/analysis/code/viet_dengue_data/output/dengue_timeseries/dengue_districts_19982021.csv") %>%
  dplyr::filter(yeardengue %in% 1998:2020) %>%
  dplyr::filter(province %in% c("Ha Noi", "Khanh Hoa", "Dak Lak", "Dong Nai")) %>%
  write.csv("./data/dengue/dengue_districts_19982020.csv", row.names = FALSE)

