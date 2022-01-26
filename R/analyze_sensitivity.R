# retrieved from server
optim <- readRDS('~/surfdrive/tmp/optim_proc_20220126/optimize_mekong_vol_sed_ci_gen200_pop100.rds')


plot(optim)
plot(-optim$value[,1],-optim$value[,2])

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
# ggplot(ob) + geom_point(aes(x = V1, y = V2, color = V3)) + xlab('IC') + ylab('sediments')
# 

# next steps
# 1 - check how to improve the speed performance of sediment part in fitness function, probably need to pre-map the dams matrices and simply multiply them in the fitness function - check how to best multiply matrices in R
# --- checked, this is the best that can be done, the bottleneck is the size of the matrix ~6k by 6k

# 2 - check whether volume or flow deviation (storage/streamflow) is the best way to substitute IC
# 3 - do the pre-filtering on dams based on whether they occurr on the tributary of the hb unit (see google docs)
# 4 - run on Alice


