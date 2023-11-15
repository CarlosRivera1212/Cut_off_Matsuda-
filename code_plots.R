# PLOTS AND FIGURES

# Extract variable names from the data
var_names = colnames(datos)[-c(1, 2)]

# Source variable name labels
source('var_name_label.R')

# Initialize lists to store plots
plt_pts_list = list()
plt_vio_list = list()

# Loop over variable names
for (vni in var_names) {
  # Print a header for the current variable name
  cat('\n\n\t', which(var_names == vni), '-', vni, '\n')
  
  # Select relevant columns from the original data
  dfs = datos_o %>%
    dplyr::select(G = Status,
                  X = Matsuda,
                  Y = all_of(vni))
  
  # Find the observation window for X and Y
  pvxy = with(dfs, f_vn(X, Y, G))
  pvx_mn = pvxy$int[1]
  pvx_mx = pvxy$int[2]
  pvx_m = pvxy$med[1]
  
  pvy_mn = pvxy$int[3]
  pvy_mx = pvxy$int[4]
  pvy_m = pvxy$med[2]
  
  # Round the threshold value to 2 decimal places
  thr_i = round(threshold_list[[vni]]['bpt_o'], 2)
  sen_i = round(threshold_list[[vni]]['sen'], 2)
  esp_i = round(threshold_list[[vni]]['esp'], 2)
  
  # Create a string for sensitivity and specificity information
  se_i = paste0('Sensitivity: ', sen_i,
                '\n',
                'Specificity: ', esp_i)
  
  # Create a violin plot
  plt_vio_list[[vni]] = ggplot(dfs) +
    geom_violin(aes(x = G, y = Y, fill = G)) +
    geom_boxplot(aes(x = G, y = Y, fill = G), width = 0.3) +
    geom_hline(yintercept = thr_i, linetype = 'dashed') +
    geom_label(
      x = 1.5,
      y = thr_i,
      label = round(thr_i, 3),
      fill = 'white'
    ) +
    geom_label(
      x = Inf,
      y = Inf,
      label = se_i,
      vjust = 1,
      hjust = 1,
      label.r = unit(0, 'cm')
    ) +
    annotate(
      'text',
      x = -Inf,
      y = Inf,
      label = var_letras_violin[[vni]],
      hjust = -0.5,
      vjust = 1.5,
      size = 6
    ) +
    labs(x = 'Status',
         y = var_names_label[[vni]]) +
    scale_fill_brewer(palette = 'Set2') +
    theme_bw() +
    theme(legend.position = 'none')
  
  # Create a scatter plot
  plt_pts_list[[vni]] = ggplot(dfs) +
    aes(X, Y) +
    geom_point(aes(color = G)) +
    scale_color_brewer(palette = 'Set2') +
    geom_vline(
      xintercept = c(pvx_mn, pvx_mx, pvx_m),
      linetype = c('solid', 'solid', 'dashed'),
      color = '#999999'
    ) +
    geom_hline(
      yintercept = c(pvy_mn, pvy_mx, pvy_m),
      linetype = c('solid', 'solid', 'dashed'),
      color = '#999999'
    ) +
    geom_hline(yintercept = thr_i) +
    geom_label(x = pvx_m,
               y = min(dfs$Y),
               label = round(pvx_m, 3)) +
    geom_label(x = max(dfs$X) * 0.9,
               y = thr_i,
               label = thr_i) +
    annotate(
      'text',
      x = Inf,
      y = Inf,
      label = var_letras_scatter[[vni]],
      hjust = 4,
      vjust = 1.5,
      size = 6
    ) +
    labs(x = var_names_label$Matsuda,
         y = var_names_label[[vni]],
         color = 'Status') +
    theme_bw()
}
