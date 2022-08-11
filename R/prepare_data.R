
# ##############################################################################
# CALCULATION BLOCK ############################################################

# HydroBASINS data -------------------------------------------------------------
cat('Preparing data..\n')
# read hydrobasins data
hb_data <- foreach(i = c('as'),.combine = 'rbind') %do% read_sf(paste0(hb_directory,'/global_lev12/hybas_',i,'_lev12_v1c.shp'))

# add basin area
main_bas_area <- hb_data %>%
  as_tibble() %>%
  select(HYBAS_ID,MAIN_BAS,SUB_AREA) %>%
  group_by(MAIN_BAS) %>%
  summarize(MAIN_BAS_AREA = sum(SUB_AREA))
hb_data <- inner_join(hb_data,main_bas_area,by='MAIN_BAS') %>%
  filter(MAIN_BAS %in% main_bas_id)

# ------------------------------------------------------------------------------

# Species data -----------------------------------------------------------------

# from pre/process_species_data.R
sp_data <- read.csv('proc/sp_data_hybas12_mekong.csv') %>% as.data.table()

# ------------------------------------------------------------------------------

# Merge interbasins ------------------------------------------------------------
cat('\n# -------------------------------------------------------------------\n')
cat('Mapping interbasins..\n')
st <- Sys.time()

# load dams data and intersect with hybas
dams <- read_sf(dams_file) %>%
  # filter(Status %in% c('E','C')) %>%
  st_transform(4326)
dams <- st_intersection(dams,hb_data %>% dplyr::select(HYBAS_ID)) %>% as_tibble() %>% dplyr::select(-geom)


# select dams hybas id for the current basin
hb_dams <- merge(hb_data,dams,by='HYBAS_ID')
dams_no_cur <- 0
# store dams no first for data frame
dams_no_cur <- nrow(hb_dams)
# clean duplicates
hb_dams <- hb_dams[!duplicated(hb_dams[,1:2]),]
# sort them, higher PFAFSTETTER number = more upstream
hb_dams <- hb_dams[order(hb_dams$PFAF_ID,decreasing = T),]

# find upstrea IDs of each dam
# st <- Sys.time()
master_upstream_list <- find_upstream_ids(t=as.data.frame(hb_data[,c('HYBAS_ID','NEXT_DOWN')])[,-3], #removed the geom column
                                          IDs=unique(hb_dams$HYBAS_ID))
# Sys.time() - st

# include the outlet (which is connected to all hb units)
master_upstream_list <- append(master_upstream_list,list(outlet = hb_data$HYBAS_ID))
names(master_upstream_list)[length(master_upstream_list)] <- outlet

# create INTERBASINS groups
inter_list <- list()
for(i in 1:length(master_upstream_list)){
  # store the hybas_id of the upstream basins
  inter_list[[i]] <- master_upstream_list[[i]]
  # need to exclude basins that are in the upstream groups
  if(i > 1) inter_list[[i]] <- inter_list[[i]][!inter_list[[i]] %in% do.call('c',inter_list[1:(i-1)])]
}
names(inter_list) <- names(master_upstream_list)

# # check all hb_basins are in there
# lapply(inter_list,length) %>% unlist %>% sum
# nrow(hb_data) # yes!

# identify next down for each interbasin
inter_basins <- data.frame(INTER_HYBAS_ID = names(inter_list), INTER_NEXT_DOWN = NA, INTER_ID = 1:length(inter_list), INTER_NEXT = NA)
for(i in 1:nrow(inter_basins)){
  nid <- hb_data$NEXT_DOWN[hb_data$HYBAS_ID == inter_basins$INTER_HYBAS_ID[i]]
  if(nid != 0){
    # search for group containing nexd down basin
    nid_group <- foreach(j = 1:length(inter_list),.combine = 'c') %do% (nid %in% inter_list[[j]]) %>% which
    inter_basins$INTER_NEXT_DOWN[i] <- inter_basins$INTER_HYBAS_ID[nid_group]
  }else{
    inter_basins$INTER_NEXT_DOWN[i] <- 0
    nid_group <- 0
  }
  inter_basins$INTER_NEXT[i] <- nid_group
}

# create correspondence table with INTER-BASINS, hybas ID and up-down basins
inter_basin_corr <- foreach(i = 1:length(names(inter_list)),.combine = 'rbind') %do% {
  data.frame(INTER_HYBAS_ID = names(inter_list)[i], HYBAS_ID = inter_list[[i]])
} %>% left_join(.,inter_basins)

# create polygon shapefile of interbasins
inter_basin_sf <- left_join(hb_data,inter_basin_corr) %>%
  select(-HYBAS_ID, -NEXT_DOWN, -NEXT_SINK) %>%
  group_by(INTER_ID,INTER_NEXT) %>%
  summarise(SUB_AREA = sum(SUB_AREA))

# recalculate a mater upstream list based on inter_basin
# find upstrea IDs of each dam
# st <- Sys.time()
master_upstream_list <- find_upstream_ids(t=as.data.frame(inter_basin_sf[,c('INTER_ID','INTER_NEXT')])[,-3], #removed the geom column
                                          IDs=inter_basin_sf$INTER_ID)
# Sys.time() - st

# remap the species to inter_basin and store L for each inter_basin
# calculate L
alpha <- 0.55
sp_data$L <- sp_data$SUB_AREA**alpha
sp_data_inter <- sp_data %>%
  left_join(inter_basin_corr) %>%
  group_by(INTER_ID, binomial) %>%
  summarise(L_tot = sum(L),diad = unique(diad)) %>% ungroup() %>% as.data.frame

# in the fitness algorithm merge groups that have 0s with the next downstream group

et <- Sys.time() - st
cat('\nCompleted in',round(et,2), attr(et,'units'),'\n')
# ------------------------------------------------------------------------------

# Passability ------------------------------------------------------------------

# split diad and pota species
sp_data_inter_p <- sp_data_inter %>% filter(diad == 'f')
sp_data_inter_d <- sp_data_inter %>% filter(diad == 't')

# caluclate Length matrix and sum it for potamodromous species
L <- lapply(split(sp_data_inter_p,sp_data_inter_p$binomial), function(x){
  bas_sp <- x %>%
    select(INTER_ID,l = L_tot) %>%
    right_join(inter_basins %>% select(INTER_ID,INTER_NEXT),by = "INTER_ID")
  bas_sp$l[is.na(bas_sp$l)] <- 0
  bas_sp <- bas_sp[order(bas_sp$INTER_ID),]
  
  l <- bas_sp$l
  I <- matrix(1, nrow = length(l), ncol = length(l))
  L <- (I*l * t(I*l)) / sum(l)**2 
  return(L)
}) %>% Reduce('+',.)

# calculate l vector for diadromous species
ld <- lapply(split(sp_data_inter_d,sp_data_inter_d$binomial), function(x){
  bas_sp <- x %>%
    select(INTER_ID,l = L_tot) %>%
    right_join(inter_basins %>% select(INTER_ID,INTER_NEXT),by = "INTER_ID")
  bas_sp$l[is.na(bas_sp$l)] <- 0
  bas_sp <- bas_sp[order(bas_sp$INTER_ID),]
  
  l <- bas_sp$l/sum(bas_sp$l)
  return(l)
}) %>% Reduce('+',.)

# graph for passability
bas_tot <- inter_basins %>% select(INTER_ID,INTER_NEXT)
df <- data.frame(from = bas_tot$INTER_ID,to = bas_tot$INTER_NEXT)
outlet <- as.character(df$from[which(df$to == 0)])

# graph
df[which(df$to == 0),]$to <- outlet %>% as.numeric()
g_master_pass <- graph_from_data_frame(d=df)
# plot(g_master_pass)



# Sedimentation data -----------------------------------------------------------
cat('\n# -------------------------------------------------------------------\n')
cat('Sedimentation data..\n')
st <- Sys.time()

# load 'geomorphic provinces' large areas with similar sediment YIELD [t/km2/yr]
GP<-read_sf(mekong_geomorphic_prov)
# hydrobasins
HYBAS <- inter_basin_sf %>% st_transform(st_crs(read_sf(hybas_rafa)))
# two ways of doing this 

# create a raster and extract sediment yields in each HB using exact_extract. 
# Advantage: Handles situations where a hydrobasins cuts accross multiple geomoeprhic provinces. 
GPr<-raster::raster(crs = raster::crs(GP), vals = 0, 
                    resolution = c(5000, 5000), ext = raster::extent(GP)) %>%
  raster::rasterize(as_Spatial(GP), ., field = "SedYield")

HYBAS$YS_1<-exact_extract(GPr,HYBAS,fun = "mean") 
HYBAS$QS_1<-HYBAS$YS_1*HYBAS$SUB_AREA

# outlet
df <- data.frame(from = HYBAS$INTER_ID,to = HYBAS$INTER_NEXT)
outlet <- as.character(df$from[which(df$to == 0)])

df[which(df$to == 0),]$to <- outlet %>% as.numeric()
g <- graph_from_data_frame(d=df)
I <- get.adjacency(g)

# calculate all connected nodes, for each row place a 1 where downstream nodes are connected
for(i in row.names(I)) {
  I[i,as.numeric(shortest_paths(graph=g,from = i,to=as.character(outlet))$vpath[[1]])] <- 1
}


dams <- read_sf(dams_file) %>%
  # filter(Status %in% c('E','C')) %>%
  st_transform(4326)
sf_use_s2(FALSE)
dams <- st_intersection(dams,inter_basin_sf %>% dplyr::select(INTER_ID)) %>% as_tibble() %>% dplyr::select(-geom)

# load dams dataset and calculate release efficiency
dams$CI <- pull(dams,name_col_V)*10**6/(pull(dams,name_col_Q)*60*60*24*365.25)
dams$TE_sus <- 1-0.05/sqrt(dams$CI) # trapping efficiency
dams$TE_sus[dams$TE_sus < 0] <- 0
dams$TE_sus[is.na(dams$TE_sus)] <- 0
dams$TE <- (1-f_bed)*dams$TE_sus + f_bed*TE_bed 
dams$RE <- 1-dams$TE # release efficiency

# calculate the source load per hydrobasin
qs <- data.frame(QS = HYBAS$QS_1)
row.names(qs) <- as.character(HYBAS$INTER_ID)
# remove NAs
qs[is.na(qs)] <- 0

# map outlet
out <- which(row.names(I) == outlet)

# premap dams release efficiency to network
RE_M_array <- foreach(j = 1:length(dams$INTER_ID)) %do% {
  M <- I
  id <- dams$INTER_ID[j] # first three rows can be defined outside the fitness function
  up_ <- row.names(I)[I[,as.character(id)] == 1]
  dw_ <- colnames(I)[I[as.character(id),] == 1]
  # multiply RE for all nodes downstream the dam
  M[up_,dw_] <- M[up_,dw_] * dams$RE[j] # fixed from dams$RE[dams$HYBAS_ID == id] to dams$RE[j] otherwise non-unique solutions identified (multiple dams belong to one HB unit)
  return(M)
}
names(RE_M_array) <- dams$INTER_ID

et <- Sys.time() - st
cat('\nCompleted in',round(et,2), attr(et,'units'),'\n')
# ------------------------------------------------------------------------------


# FITNESS FUNCTION --------------------------------------------------------------------------

# simplify dams table for fitness function
dams_data <- dams[,c('INTER_ID','Status','DamHeight',name_col_IC, name_col_V)]

# # assign passability
dams_data$pass <- 0
# dams_data$pass[dams_data$DamHeight < 20] <- 0.4
# dams_data$pass[dams_data$DamHeight < 10] <- 0.6
# dams_data$pass[dams_data$DamHeight < 1] <- 0.8
# 
# dams_data$pass <- scales::rescale(dams_data$DamHeight, to = c(1,0), from = c(0,40))
# dams_data$pass[dams_data$pass < 0] <- 0

dams_e <- dams_data[dams_data$Status != 'P',]
dams_f <- dams_data[dams_data$Status == 'P',]

# no species
n_sp <- round(sum(L) + sum(ld),0)

# no dams
n = nrow(dams_f)
if(all_dams) n = n + nrow(dams_e)

