# packages needed
library(sf); library(foreach); library(rfishbase); library(data.table); library(dplyr); library(mco); library(sp); library(exactextractr); library(igraph)
# raster package is also needed

# retrieve paths to input files
source('R/master_paths_local.R')

# functions for CI calculation
source('R/functions_connectivity.R')

# retrieved from server
library(caRamel)
op <- readRDS('~/surfdrive/tmp/optifish_proc_20220420/caramel_ic_vol_sed_ci_gen48079_pop100.rds')
plot_caramel(op)

# overlay plots
plot_pareto(op$objectives,maximized = rep(T,4),objnames = c('InCap','Vol','Sed','CI'))


# 
dec <- round(op$parameters,0) %>% as.data.frame()
ob <- op$objectives %>% as.data.frame()
# 
# # sort based on IC <<< WHY? ask Rafa
dec <- dec[sort(ob$V1,index.return=T)$ix,]
ob <- ob[sort(ob$V1,index.return=T)$ix,]
# 
# 
# INCLUSION PROBABILITY
# probability of each dam being included
dams <- read_sf('data/Dams Mekong MRC and PRC.gpkg') %>%
  # filter(Status %in% c('E','C')) %>%
  st_transform(4326)
dams$incl <- apply(dec,2,mean) %>% as.numeric
# 
plot(dams[,'incl'])

# load river network and watershed for plotting

# # watershed
# hb_data <- foreach(i = c('as'),.combine = 'rbind') %do% read_sf(paste0(hb_directory,'/global_lev12/hybas_',i,'_lev12_v1c.shp'))
# # add basin area
# main_bas_area <- hb_data %>%
#   as_tibble() %>%
#   select(HYBAS_ID,MAIN_BAS,SUB_AREA) %>%
#   group_by(MAIN_BAS) %>%
#   summarize(MAIN_BAS_AREA = sum(SUB_AREA))
# hb_data <- inner_join(hb_data,main_bas_area,by='MAIN_BAS') %>%
#   filter(MAIN_BAS %in% '4120017020')
# write_sf(hb_data,'proc/basins_mekong.gpkg')



ws <- hb_data %>% group_by(MAIN_BAS) %>% summarise(do_union = T)



riv <- read_sf('data/RiverATLAS/RiverATLAS_v10_as.shp')
riv_mekong <- riv %>% filter(HYBAS_L12 %in% hb_data$HYBAS_ID)
write_sf(riv_mekong,'proc/rivers_mekong.gpkg')


r <- riv %>% filter(DIS_AV_CMS > 10) %>%
  filter(HYBAS_L12 %in% hb_data$HYBAS_ID) %>%
  mutate(merge_col = 1) %>%
  group_by(merge_col) %>% summarise(do_union = T)


plot(st_geometry(ws))
plot(st_geometry(r),add=T)
plot(dams[,'incl'],pch = 19,add=T)

# 
# # ws <- read_sf('')
# library(ggplot2)
# ggplot(ob) + geom_point(aes(x = V1, y = V2, color = V3)) + xlab('IC') + ylab('sediments')
# 

# next steps
# 1 - check how to improve the speed performance of sediment part in fitness function, probably need to pre-map the dams matrices and simply multiply them in the fitness function - check how to best multiply matrices in R
# --- checked, this is the best that can be done, the bottleneck is the size of the matrix ~6k by 6k

# 2 - check whether volume or flow deviation (storage/streamflow) is the best way to substitute IC
# 3 - do the pre-filtering on dams based on whether they occurr on the tributary of the hb unit (see google docs)
# 4 - run on Alice


