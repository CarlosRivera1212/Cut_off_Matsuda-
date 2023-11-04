library(readxl)
library(ggpubr)
library(plyr)
library(tidyverse)
library(gridExtra)
library(writexl)

setwd("D:/DRIVE_101323/My Drive/PC/UNAL/D/articulos/gluc")
res_path = 'resultados/'
datos = read_excel('INDICES_LEPTIN_imputado_101923.xlsx')

datos_o = datos

datos = datos %>% 
  mutate_at(vars(Matsuda:TyG_WHtR), scale)


(pot = seq(0.98, 1.5, 0.01))

f_mm = function(x){ (x - min(x)) / (max(x) - min(x)) }
f_vn = function(x,y,g){
  px_mn = max(x[g == 'IR'])
  px_mx = min(x[g == 'NIR'])
  px_m = (px_mn + px_mx)/2
  
  pre_y_mn = rev(tapply(y, g, min))
  pre_y_mx = tapply(y, g, max)
  py_mn = pre_y_mn[which.min(pre_y_mx-pre_y_mn)]
  py_mx = pre_y_mx[which.min(pre_y_mx-pre_y_mn)]
  py_m = (py_mn + py_mx)/2
  
  return(list(
    med = c(px_m, py_m),
    int = c(px_mn, px_mx, py_mn, py_mx)
  ))
}

f_po = function(x,y,g, mrx, mry, mvx, mvy){
  r = cor(x,y, method = 'spearman') ** 2
  
  d = as.matrix(dist(cbind(f_mm(x), f_mm(y))))
  # d[d>max(d)/2] = 0
  di = 1/d
  # di = (1/d)**(-1)
  di[is.infinite(di)] = 0
  
  dif = NULL
  for (i in pot){
    dip = di ** i
    W = dip/apply(dip, 1, sum)
    W[is.na(W)] = 0
    
    # xm = mean( (W %*% x) * r + (1-r)*mrx )
    # ym = mean( (W %*% y) * r + (1-r)*mry )
    xm = mean( (W %*% x) * r + (1-r)*x )
    ym = mean( (W %*% y) * r + (1-r)*y )
    
    dif_i = sqrt((xm-mvx)**2 + (ym-mvy)**2)
    dif = append(dif, dif_i)
  }
  
  pot_sel = pot[which.min(dif)]
  dip = di ** pot_sel
  W = dip/apply(dip, 1, sum)
  W[is.na(W)] = 0
  
  # Wx = (W %*% x) * r + (1-r)*mrx
  # Wy = (W %*% y) * r + (1-r)*mry
  Wx = (W %*% x) * r + (1-r)*x
  Wy = (W %*% y) * r + (1-r)*y
  
  return(c(mwx = mean(Wx), mwy = mean(Wy)))
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# ANALISIS ####
var_names = colnames(datos)[-c(1,2)]

threshold_list = list()
tot_pts = list()

for(vni in var_names){
  cat('\n\n\t', which(var_names == vni), '-', vni, '\n')
  
  dfs = datos %>%
    dplyr::select(G = Status,
                  X = Matsuda,
                  Y = all_of(vni))
  
  pvxy = with(dfs, f_vn(X,Y,G))
  pvx_m = pvxy$med[1]
  pvy_m = pvxy$med[2]
  
  prx = mean(dfs$X)
  pry = mean(dfs$Y)
  
  res_i = lapply(seq(nrow(dfs)), function(i){
    if(i %% 10 == 0){ cat(i,'- ') }
    dfs_f = dfs[-i,]
    f_po(x = as.vector(dfs_f$X),
         y = as.vector(dfs_f$Y),
         g = as.vector(dfs_f$G),
         mrx = prx, mry = pry,
         mvx = pvx_m, mvy = pvy_m)
  })
  
  bpts = NULL
  js = NULL
  for (ri in res_i){
    bpt = ri[2]
    if (vni %in% c('QUICKI')) {
      tbl = table(y = ifelse(dfs$Y<bpt, 'IR', 'NIR'),
                  x = ifelse(dfs$X>pvx_m, 'NIR', 'IR'))
    } else {
      tbl = table(y = ifelse(dfs$Y<bpt, 'NIR', 'IR'),
                  x = ifelse(dfs$X>pvx_m, 'NIR', 'IR'))
    }
    VP = tbl[1,1]
    VN = tbl[2,2]
    FP = tbl[1,2]
    FN = tbl[2,1]
    
    sen = VP / (VP+FN)
    esp = VN / (VN+FP)
    j = sen+esp-1
    js = append(js, j)
    bpts = append(bpts, bpt)
  }
  
  tot_pts[[vni]] = res_i
  bpt = mean(bpts[js == max(js)])
  
  med_o = attr(dfs$Y, 'scaled:center')
  des_o = attr(dfs$Y, 'scaled:scale')
  bpt_o = bpt*des_o + med_o
  
  if (vni %in% c('QUICKI')) {
    tbl = table(y = ifelse(dfs$Y<bpt, 'IR', 'NIR'),
                x = ifelse(dfs$X>pvx_m, 'NIR', 'IR'))
  } else {
    tbl = table(y = ifelse(dfs$Y<bpt, 'NIR', 'IR'),
                x = ifelse(dfs$X>pvx_m, 'NIR', 'IR'))
  }
  
  VP = tbl[1,1]
  VN = tbl[2,2]
  FP = tbl[1,2]
  FN = tbl[2,1]
  
  sen = VP / (VP+FN)
  esp = VN / (VN+FP)
  
  threshold_list[[vni]] = c(bpt = bpt, bpt_o = bpt_o, sen = sen, esp = esp)
}


saveRDS(threshold_list, paste0(res_path,'threshold_102223.rds'))
threshold_list = readRDS('resultados/threshold_102223.rds')

df_threshold = ldply(threshold_list, .id = 'var_name') %>%
  dplyr::rename(cutoff_s = bpt,
                cutoff_o = bpt_o,
                Sensivity = sen,
                Specificity = esp)
df_threshold %>% 
  mutate(j=Sensivity+Specificity-1)
write_xlsx(df_threshold, 'resultados/df_threshold_102223.xlsx')

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# GRAFICOS ####
df_threshold = read_excel('resultados/df_threshold_102223.xlsx')
var_names = colnames(datos)[-c(1,2)]

# load('resultados/var_names_label.rda')
source('var_name_label.R')
var_names_label
plt_pts_list = list()
plt_vio_list = list()

for(vni in var_names){
  cat('\n\n\t', which(var_names == vni), '-', vni, '\n')
  
  dfs = datos_o %>% 
    dplyr::select(G = Status,
           X = Matsuda,
           Y = all_of(vni))
  
  pvxy = with(dfs, f_vn(X,Y,G))
  pvx_mn = pvxy$int[1]
  pvx_mx = pvxy$int[2]
  pvx_m = pvxy$med[1]
  
  pvy_mn = pvxy$int[3]
  pvy_mx = pvxy$int[4]
  pvy_m = pvxy$med[2]
  
  thr_i = round(threshold_list[[vni]]['bpt_o'], 2)
  sen_i = round(threshold_list[[vni]]['sen'], 2)
  esp_i = round(threshold_list[[vni]]['esp'], 2)
  
  se_i = paste0(
    'Sensitivity: ', sen_i,
    '\n',
    'Specificity: ', esp_i
  )
  
  plt_vio_list[[vni]] = ggplot(dfs)+
    geom_violin(aes(x=G, y=Y, fill = G))+
    geom_boxplot(aes(x=G, y=Y, fill = G), width=0.3)+
    geom_hline(yintercept = thr_i, linetype='dashed')+
    geom_label(x = 1.5, y = thr_i, label=round(thr_i,3), fill='white')+
    geom_label(x=Inf, y=Inf, label=se_i, vjust=1, hjust=1,
               label.r = unit(0, 'cm'))+
    labs(x = 'Status',
         y = var_names_label[[vni]])+
    scale_fill_brewer(palette = 'Set2')+
    theme_bw()+
    theme(legend.position = 'none')
  
  # plt_vio_list[[vni]]
  
  plt_pts_list[[vni]] = ggplot(dfs)+
    aes(X, Y)+
    geom_point(aes(color = G))+
    scale_color_brewer(palette = 'Set2')+
    
    geom_vline(xintercept = c(pvx_mn, pvx_mx, pvx_m),
               linetype = c('solid', 'solid', 'dashed'),
               color = '#999999')+
    geom_hline(yintercept = c(pvy_mn, pvy_mx, pvy_m),
               linetype = c('solid', 'solid', 'dashed'),
               color = '#999999')+
    
    geom_hline(yintercept = thr_i)+
    
    geom_label(x = pvx_m, y = min(dfs$Y), label=round(pvx_m,3))+
    geom_label(x = max(dfs$X)*0.9, y = thr_i, label=thr_i)+
    
    labs(x = 'Matsuda',
         y = var_names_label[[vni]],
         color = 'Status')+
    theme_bw()
}

for(ni in names(plt_pts_list)){
  print(ni)
  plt_i = plt_pts_list[[ni]]
  ggsave(paste0(ni,'.tiff'), plt_i, 'tiff',
         paste0(res_path, 'graf_102223/scatter/'),
         1, 16, 9, 'cm', 300)
}

for(ni in names(plt_vio_list)){
  print(ni)
  plt_i = plt_vio_list[[ni]]
  ggsave(paste0(ni,'.tiff'), plt_i, 'tiff',
         paste0(res_path, 'graf_102223/vio/'),
         1, 16, 9, 'cm', 300)
}
