# packages needed
library(sf); library(foreach); library(rfishbase); library(data.table); library(dplyr); library(mco); library(sp); library(exactextractr); library(igraph)
# raster package is also needed

# retrieve paths to input files
source('R/master_paths_local.R')

# functions for CI calculation
source('R/functions_connectivity.R')

# retrieved from server
optim <- readRDS('~/surfdrive/tmp/optim_proc_20220202/optimize_mekong_ic_vol_sed_ci_gen200_pop100.rds')


plot(optim)
plot(-optim$value[,1],-optim$value[,2])
plot(-optim$value[,3],-optim$value[,4])
plot(-optim$value[,2],-optim$value[,4])
plot(-optim$value[,1],-optim$value[,3])
# 
dec <- round(optim$par,0) %>% as.data.frame()
ob <- -optim$value %>% as.data.frame()
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
plot(st_geometry(dams))
plot(dams[,'incl'])
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


