LOCAL = FALSE

# packages needed
library(sf); library(foreach); library(rfishbase); library(data.table); 
library(dplyr); library(mco); library(sp); library(exactextractr); library(igraph)
# raster package is also needed

# retrieve paths to input files
source('R/master_paths.R')
if(LOCAL) source('R/master_paths_local.R')


# functions for CI calculation
source('R/functions_connectivity.R')

# HYBAS ID of Mekong
main_bas_id <- outlet <- 4120017020 #this corresponds to the outlet!

# ------------------------------------------------------------------------------
# HydroBASINS data -------------------------------------------------------------

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
# sp_data <- read.csv('proc/sp_data_mekong.csv') %>% as.data.table()

# redo these steps as in level 08 above <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
cat('\nRetrieving diadromous species from fishbase..')

if(!file.exists('proc/fishbase_data.csv')){
  # # load fishbase metadata
  fishbase <- species(fields = c('Species','AnaCat')) %>% # get species table for all species
    rename(binomial = Species)
  #
  fishbase$AnaCat <- as.factor(fishbase$AnaCat)
  levels(fishbase$AnaCat) <- c(rep('Diad.',23),'Non.','Ocea.','Ocea.','Pota.','Pota.')
  
  write.csv(fishbase,'proc/fishbase_data.csv',row.names = F)
  
}else{
  fishbase <- read.csv('proc/fishbase_data.csv') %>% as_tibble  
}

sp_data <- readRDS(sp_ranges_file) %>%
  inner_join(.,hb_data %>% as_tibble() %>% select(HYBAS_ID,MAIN_BAS,SUB_AREA,MAIN_BAS_AREA),by="HYBAS_ID") %>%
  as.data.table(.)

# assign diadromous-non diadromous category
sp_data$diad <- 'f'
sp_data$diad[sp_data$binomial %in% fishbase$binomial[fishbase$AnaCat == 'Diad.']] <- 't'


# ------------------------------------------------------------------------------
# Merge interbasins ------------------------------------------------------------

# load dams data and intersect with hybas
dams <- read_sf(dams_file) %>%
  # filter(Status %in% c('E','C')) %>%
  st_transform(4326)
sf_use_s2(FALSE)
dams <- st_intersection(dams,hb_data %>% dplyr::select(HYBAS_ID)) %>% as_tibble() %>% dplyr::select(-geom)


# select dams hybas id for the current basin
dcur <- merge(hb_data,dams,by='HYBAS_ID')
dams_no_cur <- 0
# store dams no first for data frame
dams_no_cur <- nrow(dcur)
# clean duplicates
dcur <- dcur[!duplicated(dcur[,1:2]),]
# sort them, higher PFAFSTETTER number = more upstream
dcur <- dcur[order(dcur$PFAF_ID,decreasing = T),]

# find upstrea IDs of each dam
st <- Sys.time()
master_upstream_list <- find_upstream_ids(t=as.data.frame(hb_data[,c('HYBAS_ID','NEXT_DOWN')])[,-3], #removed the geom column
                                          IDs=unique(dcur$HYBAS_ID))
Sys.time() - st

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
st <- Sys.time()
master_upstream_list <- find_upstream_ids(t=as.data.frame(inter_basin_sf[,c('INTER_ID','INTER_NEXT')])[,-3], #removed the geom column
                                          IDs=inter_basin_sf$INTER_ID)
Sys.time() - st

# remap the species to inter_basin and store L for each inter_basin
# calculate L
alpha <- 0.55
sp_data$L <- sp_data$SUB_AREA**alpha
sp_data_inter <- sp_data %>%
  left_join(inter_basin_corr) %>%
  group_by(INTER_ID, binomial) %>%
  summarise(L_tot = sum(L),diad = unique(diad)) %>% ungroup()


# assign inter ids to dcur
dcur <- dcur %>%
  left_join(inter_basins %>% mutate(INTER_HYBAS_ID = as.numeric(INTER_HYBAS_ID)),by = c(HYBAS_ID = 'INTER_HYBAS_ID'))

# in the fitness algorithm merge groups that have 0s with the next downstream group


# ------------------------------------------------------------------------------
# Sedimentation data -----------------------------------------------------------

# load 'geomorphic provinces' large areas with similar sediment YIELD [t/km2/yr]
GP<-read_sf(mekong_geomorphic_prov)
# hydrobasins
HYBAS <- inter_basin_sf %>% st_transform(st_crs(read_sf(hybas_rafa)))
# two ways of doing this 

# create a raster and extract sediment yields in each HB using exat extract. 
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
dams$CI <- dams$GrossStora*10**6/(dams$MeanQ_m3s*60*60*24*365.25)
dams$TE_sus <- 1-0.05/sqrt(dams$CI) # trapping efficiency
dams$TE_sus[dams$TE_sus < 0] <- 0
dams$TE_sus[is.na(dams$TE_sus)] <- 0
TE_bed <- 1
f_bed <- 0.1
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

# FITNESS FUNCTION --------------------------------------------------------------------------

fitness <- function(x, sedim = T, nc = 1){ # make it modular to switch on and off different modules
  # library(dplyr)
  # nc<-24 #set no. cores
  decision <- round(x,0)
  ids <- dams$INTER_ID[decision == 1] %>% unique
  
  # installed capacity ----------------------------------------
  totIC <- sum(dams$InstalledC*decision)
  
  # volume ----------------------------------------
  totVL <- sum(dams$GrossStora*decision)
  
  # sedimentation ---------------------------------------------
  # st <- Sys.time()
  if(sedim == T){
    
    # select dams in current portfolio
    dams_sed <- dams[decision == 1,]
    
    if(nrow(dams_sed) > 0){
      # calculate total RE matrix (= multiply single dam RE matrices)
      MQS = Reduce('*',RE_M_array[as.character(dams_sed$INTER_ID)])
      
      # multiply RE by load to get final load at mouth
      totQS <- as.numeric(qs[row.names(MQS),'QS'] %*% MQS)[out]
    }else{
      totQS <- sum(qs$QS)
    }
    
  }
  # Sys.time() - st
  
  # fragmentation ----------------------------------------------
  dc <- dcur[dcur$INTER_ID %in% ids,]
  
  if(nrow(dc) > 0){
    
    # connectivity index
    upstream_list <- master_upstream_list[names(master_upstream_list) %in% dc$INTER_ID]
    
    # create ID groups from most upstream to downstream
    # exclude IDs in downstream groups already present in upstream groups
    groups_cur <- list()
    for(i in 1:nrow(dc)){
      # store the hybas_id of the upstream basins
      groups_cur[[i]] <- upstream_list[[as.character(dc$INTER_ID[i])]]
      # need to exclude basins that are in the upstream groups
      if(i > 1) groups_cur[[i]] <- groups_cur[[i]][!groups_cur[[i]] %in% do.call('c',groups_cur[1:(i-1)])]
    }
    names(groups_cur) <- 1:length(groups_cur)
    groups_cur <- groups_cur[sapply(groups_cur,function(x) length(x) != 0)]
    
    groups_cur_df <- mapply(function(x,y) data.frame(INTER_ID = x, group_cur = y),x=groups_cur,y=1:length(groups_cur), SIMPLIFY = F) %>% do.call('rbind',.)
    
  }else{
    groups_cur_df <- data.frame(INTER_ID = inter_basins$INTER_ID, group_cur = 0)
  }
  
  # map grouping of interbasins based on selected dams
  sbas_sp = left_join(sp_data_inter,groups_cur_df, by = 'INTER_ID')
  sbas_sp$group_cur[is.na(sbas_sp$group_cur)] <- 0
  
  sps <- unique(as.character(sbas_sp$binomial))
  
  if(nc > 1){
    st <- Sys.time()
    CI <- parallel::mclapply(split(sps,cut(1:length(sps),nc,labels = F)), function(spp){
      lapply(spp, function(sp){
        occ <- sbas_sp[sbas_sp$binomial == sp,]
        
        L <- occ$L_tot
        
        if(occ$diad[1] == 't'){
          # total area
          A = sum(L)
          sa_cur = sum(L[occ$group_cur == min(occ$group_cur)])
          # cat <- 'diadromous'
        }else{
          # total area
          A = sum(L)**2
          # sum pf areas for patches from different groups
          sa_cur <- sapply(unique(occ$group_cur), function(i) sum(L[occ$group_cur == i])**2) %>% sum
          # cat <- 'potamodromous'
        }
        
        # output
        return((sa_cur/A)*100)
        
        
      }) %>% do.call('c',.)},mc.cores = nc) %>% unlist
    Sys.time() - st
    
  }else{
    # st <- Sys.time()
    CI <- lapply(sps, function(sp){
      occ <- sbas_sp[sbas_sp$binomial == sp,]
      
      L <- occ$L_tot
      
      if(occ$diad[1] == 't'){
        # total area
        A = sum(L)
        sa_cur = sum(L[occ$group_cur == min(occ$group_cur)])
        # cat <- 'diadromous'
      }else{
        # total area
        A = sum(L)**2
        # sum pf areas for patches from different groups
        sa_cur <- sapply(unique(occ$group_cur), function(i) sum(L[occ$group_cur == i])**2) %>% sum
        # cat <- 'potamodromous'
      }
      
      # output
      return((sa_cur/A)*100)
      
    }) %>% do.call('c',.)
    # Sys.time() - st
    
  }
  
  
  
  # Return output --------------------------------------------------------------
  if(sedim == T){
    return(c(-totIC,-totVL,-totQS,-mean(CI,na.rm=T)))
  }else{
    return(c(-totIC,-totVL,-mean(CI,na.rm=T)))
  }
}

# Run the algorithm -------------------------------------------------------------------------

# st <- Sys.time()
# fitness(sample(c(0,1),nrow(dams),replace = T),nc=8)
# Sys.time() - st

st <- Sys.time()
fitness(rep(0,nrow(dams)), sedim = T, nc=1)
Sys.time() - st


st <- Sys.time()
fitness(rep(1,nrow(dams)), sedim = T, nc=1)
Sys.time() - st
# should yield (based on previous script with all hybas12 units):
# > fitness(rep(1,nrow(dams)), sedim = T, nc=8)
# [1] -5.879448e+04 -1.750559e+05 -1.065452e+07 -4.782631e+01


n <- nrow(dams)
st <- Sys.time()
optim <- nsga2(fn = fitness,
               nc = 1, # passed to fitness
               sedim = T, # passed to fitness
               idim = n, 
               odim = 4, 
               generations = 1000,
               popsize = 40, 
               mprob = 0.2, 
               cprob = 0.8,
               lower.bounds = rep(0, n), upper.bounds = rep(1, n))
Sys.time() - st

saveRDS(optim,'proc/optimize_lev12_mekong_ic_vol_sed_ci_gen1000_pop40.rds')

# # GA does not work on multi-objectives, only single objective, but is parallelized
# library(GA)
# st <- Sys.time()
# a=ga(
#   type = 'real-valued',
#   fitness = fitness,
#   nc = 1, # passed to fitness
#   sedim = T, # passed to fitness
#   # idim = n, 
#   # odim = 4, 
#   maxiter = 4,
#   popSize = 40, 
#   pmutation = 0.2, 
#   pcrossover = 0.8,
#   lower = rep(0, n), upper = rep(1, n))
# Sys.time() - st
# 
# st <- Sys.time()
# a=ga(
#   type = 'real-valued',
#   fitness = fitness,
#   nc = 1, # passed to fitness
#   sedim = T, # passed to fitness
#   # idim = n, 
#   # odim = 4, 
#   maxiter = 4,
#   popSize = 40, 
#   pmutation = 0.2, 
#   pcrossover = 0.8,
#   lower = rep(0, n), upper = rep(1, n), parallel = T)
# Sys.time() - st

# could try NSGA3 implementation
# https://cran.r-project.org/web/packages/nsga3/nsga3.pdf
# but requires quite some changing of the objective function structure etc

