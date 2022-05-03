
# ##############################################################################
# SETTINGS BLOCK ###############################################################

LOCAL = FALSE
# g <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# caramel optimization setups --------------------------------------------------
pop_run <- 100
tot_run <- 6000000
init_pop_run <- 100
arch_run <- 100
# ------------------------------------------------------------------------------

# fitness function setup -------------------------------------------------------
sedimentation = T
fragmentation = T
water = T
energy = T
# ------------------------------------------------------------------------------

# dams table details -----------------------------------------------------------
name_col_IC <- 'InstalledC'
name_col_V <- 'GrossStora'
name_col_Q <- 'MeanQ_m3s' # to calculate release efficiency only
# ------------------------------------------------------------------------------

# coeffs for sedimentation calculations ----------------------------------------
TE_bed <- 1
f_bed <- 0.1
# ------------------------------------------------------------------------------

# MAIN BASIN ID ----------------------------------------------------------------
main_bas_id <- outlet <- 4120017020 #this corresponds to the outlet!
# ------------------------------------------------------------------------------

# packages & ancillary functions -----------------------------------------------
library(sf); library(foreach); library(rfishbase); 
library(data.table); library(dplyr); library(mco); 
library(sp); library(exactextractr); library(igraph)

# settings for sf
sf_use_s2(FALSE)

# read paths to input files
source('R/master_paths.R')
if(LOCAL) source('R/master_paths_local.R')

# functions for CI calculation
source('R/functions_connectivity.R')
# ------------------------------------------------------------------------------

# ##############################################################################
# CALCULATION BLOCK ############################################################

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
sp_data <- read.csv('proc/sp_data_hybas12_mekong.csv') %>% as.data.table()

# ------------------------------------------------------------------------------

# Merge interbasins ------------------------------------------------------------

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
st <- Sys.time()
master_upstream_list <- find_upstream_ids(t=as.data.frame(hb_data[,c('HYBAS_ID','NEXT_DOWN')])[,-3], #removed the geom column
                                          IDs=unique(hb_dams$HYBAS_ID))
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
  summarise(L_tot = sum(L),diad = unique(diad)) %>% ungroup() %>% as.data.frame

# in the fitness algorithm merge groups that have 0s with the next downstream group

# ------------------------------------------------------------------------------

# Sedimentation data -----------------------------------------------------------

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

# FITNESS FUNCTION --------------------------------------------------------------------------

# simplify dams table for fitness function
dams_data <- dams[,c('INTER_ID',name_col_IC, name_col_V)]

# required data from global env:
# dams_data: table with INTER_ID, InstalledC, GrossStora [data frame]
# RE_M_array: premapped dams release efficiency to the interbasins network [array of matrices]
# qs: source load per interbasin [one column data frame]
# out: ID of outlet interbasin
# inter_basins: table with all interbasins as specified by INTER_ID
# master_upstream_list: list with ids of upstream interbasins for each interbasin
# sp_data_inter: table with species (and diadromy) in each INTER_ID and pre-mapped total length per INTER_ID: INTER_ID,binomial,L_tot,diad

fitness <- 
  function(index_caramel){ 
    
    # select dams based on index
    decision <- round(x[index_caramel,],0)
    dams_inter_ids <- unique(dams_data$INTER_ID[decision == 1])
    
    # result vector
    result_fitness <- numeric()
    
    # installed capacity ----------------------------------------
    if(energy) result_fitness <- c(result_fitness,sum(pull(dams_data,name_col_IC)*decision))
    
    # volume ----------------------------------------
    if(water) result_fitness <- c(result_fitness,sum(pull(dams_data,name_col_V)*decision))
    
    # sedimentation ---------------------------------------------
    if(sedimentation){
      if(length(dams_inter_ids) > 0){
        # calculate total RE matrix (= multiply single dam RE matrices)
        MQS = Reduce('*',RE_M_array[as.character(dams_data$INTER_ID[decision == 1])])
        
        # multiply RE by load to get final load at mouth
        totQS <- as.numeric(qs[row.names(MQS),'QS'] %*% MQS)[out]
      }else{
        totQS <- sum(qs$QS)
      }
      result_fitness <- c(result_fitness,totQS)
    }
    
    
    # fragmentation ----------------------------------------------
    
    # inter_ <- inter_dams[inter_dams$INTER_ID %in% dams_inter_ids,]
    
    if(fragmentation){
      if(length(dams_inter_ids) > 0){
        
        # connectivity index
        upstream_list <- master_upstream_list[names(master_upstream_list) %in% dams_inter_ids]
        
        # create ID groups from most upstream to downstream
        # exclude IDs in downstream groups already present in upstream groups
        groups <- list()
        for(i in 1:length(dams_inter_ids)){
          # store the hybas_id of the upstream basins
          groups[[i]] <- upstream_list[[as.character(dams_inter_ids[i])]]
          # need to exclude basins that are in the upstream groups
          if(i > 1) groups[[i]] <- groups[[i]][!groups[[i]] %in% do.call('c',groups[1:(i-1)])]
        }
        names(groups) <- 1:length(groups)
        groups <- groups[sapply(groups,function(x) length(x) != 0)]
        
        # create a table with group ID assigned to each inter ID
        groups_df <- do.call('rbind',
                             mapply(function(x,y) data.frame(INTER_ID = x, group = y),x=groups,y=1:length(groups), SIMPLIFY = F)
        )
        
      }else{
        # if there are no dams, create a table where each inter ID is assigned to the same group
        groups_df <- data.frame(INTER_ID = inter_basins$INTER_ID, group = 0)
      }
      
      # map grouping of interbasins based on selected dams
      sbas_sp = left_join(sp_data_inter,groups_df, by = 'INTER_ID')
      sbas_sp$group[is.na(sbas_sp$group)] <- 0
      
      # get all the species names
      spp <- unique(as.character(sbas_sp$binomial))
      
      # calculate CI for each species
      totCI <- 
        sum(
          do.call('c',
                  lapply(spp, function(sp){
                    occ <- sbas_sp[sbas_sp$binomial == sp,]
                    
                    L <- occ$L_tot
                    
                    if(occ$diad[1] == 't'){
                      # total area
                      A = sum(L)
                      sa_cur = sum(L[occ$group == min(occ$group)])
                      # cat <- 'diadromous'
                    }else{
                      # total area
                      A = sum(L)**2
                      # sum pf areas for patches from different groups
                      sa_cur <- sum(
                        sapply(unique(occ$group), 
                               function(i) sum(L[occ$group == i])**2)
                      )
                      # cat <- 'potamodromous'
                    }
                    
                    # output
                    return((sa_cur/A)*100)
                    
                  })
          ),
          na.rm=T)
      
      result_fitness <- c(result_fitness,totCI) 
    }
    
    
    # Return output --------------------------------------------------------------
    return( result_fitness )
    
  }

# Run the algorithm -------------------------------------------------------------------------

library(caRamel)
n = nrow(dams)

# need init function
InitFitness <- function(cl,numcores){    
  # packages
  parLapply( cl, 1:numcores, function(xx){
    require('dplyr')
    require('Matrix')
  })
  # variables needed:
  # dams_data: table with INTER_ID, InstalledC, GrossStora [data frame]
  # RE_M_array: premapped dams release efficiency to the interbasins network [array of matrices]
  # qs: source load per interbasin [one column data frame]
  # out: ID of outlet interbasin
  # inter_basins: table with all interbasins as specified by INTER_ID
  # master_upstream_list: list with ids of upstream interbasins for each interbasin
  # sp_data_inter: table with species (and diadromy) in each INTER_ID and pre-mapped total length per INTER_ID: INTER_ID,binomial,L_tot,diad
  clusterExport(
    cl=cl, 
    varlist=c("dams_data",
              "RE_M_array","qs","out",
              "master_upstream_list","inter_basins",
              "sp_data_inter",
              "name_col_IC", "name_col_V",
              "sedimentation","fragmentation","water","energy")
  )
} 

# define number of objectives based on settings
nobjs <- sum(sedimentation,fragmentation,water,energy)
op <- caRamel(
  nobj = nobjs,
  nvar = n,
  minmax = rep(TRUE,nobjs),
  bounds = matrix(c(rep(0,n),rep(1,n)), ncol = 2, nrow = n),
  func = fitness,
  repart_gene = c(pop_run/4,pop_run/4,pop_run/4,pop_run/4),
  funcinit = InitFitness,
  popsize = init_pop_run,
  archsize = arch_run,
  maxrun = (tot_run+init_pop_run),
  prec = matrix(1.e-5, nrow = 1, ncol = nobjs),
  carallel = TRUE,
  graph = FALSE,
  sensitivity = FALSE
)

plot_caramel(op)
plot_pareto(op$objectives,maximized = rep(T,nobjs),objnames = c('InCap','Vol','Sed','CI')[c(energy,water,sedimentation,fragmentation)])

# save the right variables in the output name
save_str <- c('ic','vol','sed','ci')[c(energy,water,sedimentation,fragmentation)]

saveRDS(op,paste0('proc/caramel_',paste(save_str,collapse = '_'),'_gen',nrow(op$save_crit),'_pop',pop_run,'.rds'))
