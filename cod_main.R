# Load necessary packages
library(readxl)
library(ggpubr)
library(plyr)
library(tidyverse)
library(gridExtra)
library(writexl)

# Set the path for saving results
res_path = 'resultados/'

# Read data from Excel file
datos = read_excel('data.xlsx')

# Save a copy of the original data
datos_o = datos

# Standardize numeric variables in the data
datos = datos %>%
  mutate_at(vars(Matsuda:TyG_WHtR), scale)

# Define a sequence of potential values for parameter "p"
(pot = seq(0.98, 1.5, 0.01))

# Source scripts for additional functions
# f_mm() for loading min-max normalization
source('repo/function_minmax.R')

# f_vn() for loading observation window
source('repo/function_window.R')

# f_po() for loading minimizing distance
source('repo/function_distance.R')


# ANALISIS ####
var_names = colnames(datos)[-c(1, 2)]

threshold_list = list()
tot_pts = list()

# Loop over variable names (var_names)
for (vni in var_names) {
  # Print a header for the current variable name
  cat('\n\n\t', which(var_names == vni), '-', vni, '\n')
  
  # Select relevant columns from the data
  dfs = datos %>%
    dplyr::select(G = Status,
                  X = Matsuda,
                  Y = all_of(vni))
  
  # Find the observation window for X and Y
  pvxy = with(dfs, f_vn(X, Y, G))
  pvx_m = pvxy$med[1]
  pvy_m = pvxy$med[2]
  
  # Calculate mean values for X and Y
  prx = mean(dfs$X)
  pry = mean(dfs$Y)
  
  # Initialize a list to store results for each data point
  res_i = lapply(seq(nrow(dfs)), function(i) {
    if (i %% 10 == 0) {
      cat(i, '- ')
    }
    dfs_f = dfs[-i, ]
    f_po(
      x = as.vector(dfs_f$X),
      y = as.vector(dfs_f$Y),
      g = as.vector(dfs_f$G),
      mrx = prx,
      mry = pry,
      mvx = pvx_m,
      mvy = pvy_m
    )
  })
  
  # Initialize vectors to store breakpoint values and J statistics
  bpts = NULL
  js = NULL
  
  # Loop over results and calculate J statistic for each breakpoint
  for (ri in res_i) {
    bpt = ri[2]
    # Create a contingency table based on conditions
    tbl = table(y = ifelse(dfs$Y < bpt, 'NIR', 'IR'),
                x = ifelse(dfs$X > pvx_m, 'NIR', 'IR'))
    # Calculate sensitivity, specificity, and J statistic
    VP = tbl[1, 1]
    VN = tbl[2, 2]
    FP = tbl[1, 2]
    FN = tbl[2, 1]
    sen = VP / (VP + FN)
    esp = VN / (VN + FP)
    j = sen + esp - 1
    js = append(js, j)
    bpts = append(bpts, bpt)
  }
  
  # Save results for each variable name
  tot_pts[[vni]] = res_i
  
  # Calculate the threshold that maximizes the J statistic
  bpt = mean(bpts[js == max(js)])
  
  # Scale the threshold to the original data scale
  med_o = attr(dfs$Y, 'scaled:center')
  des_o = attr(dfs$Y, 'scaled:scale')
  bpt_o = bpt * des_o + med_o
  
  # Create a contingency table based on the selected threshold
  tbl = table(y = ifelse(dfs$Y < bpt, 'NIR', 'IR'),
              x = ifelse(dfs$X > pvx_m, 'NIR', 'IR'))
  
  # Calculate sensitivity and specificity for the selected threshold
  VP = tbl[1, 1]
  VN = tbl[2, 2]
  FP = tbl[1, 2]
  FN = tbl[2, 1]
  sen = VP / (VP + FN)
  esp = VN / (VN + FP)
  
  # Save threshold information for each variable name
  threshold_list[[vni]] = c(
    bpt = bpt,
    bpt_o = bpt_o,
    sen = sen,
    esp = esp
  )
}


df_threshold = ldply(threshold_list, .id = 'var_name') %>%
  dplyr::rename(
    cutoff_s = bpt,
    cutoff_o = bpt_o,
    Sensivity = sen,
    Specificity = esp
  )

df_threshold