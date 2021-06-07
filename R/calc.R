# packages needed
library(sf); library(foreach); library(rfishbase); library(data.table); library(dplyr); library(mco)

# sp_ranges_file <- '~/surfdrive/Documents/projects/fishsuit/proc/species_ranges_raw_on_hybas12.rds' #local
sp_ranges_file <- '~/data/data/species_ranges_raw_on_hybas12.rds'
# hb_directory <- 'data/HydroBASINS/global_lev12/' #local
hb_directory <- '~/data/data/HydroBASINS/global_lev12/'
# dams_file <- 'data/Dams Mekong MRC and PRC.gpkg' #local
dams_file <- '~/data/data/Mekong dams/Dams Mekong MRC and PRC.gpkg'

# HydroBASINS data ------------------------------------------------------------------------------------------
# read hydrobasins data
hb_data <- foreach(i = c('as'),.combine = 'rbind') %do% read_sf(paste0(hb_directory,'hybas_',i,'_lev12_v1c.shp'))
# add basin area
# main_bas_area <- do.call('rbind',lapply(split(hb_data_frame,hb_data_frame$MAIN_BAS),function(x) data.frame(MAIN_BAS = unique(x$MAIN_BAS),MAIN_BAS_AREA = sum(x$SUB_AREA))))
cat('\nCompiling main basin area..')

main_bas_area <- hb_data %>%
  as_tibble() %>%
  select(HYBAS_ID,MAIN_BAS,SUB_AREA) %>%
  group_by(MAIN_BAS) %>%
  summarize(MAIN_BAS_AREA = sum(SUB_AREA))

hb_data <- inner_join(hb_data,main_bas_area,by='MAIN_BAS') %>%
  filter(MAIN_BAS %in% c(4120017020))#,4120023810,4120023060))

# select diadromous and non-diadromous species------------------------------------------------------------

cat('\nRetrieving diadromous species from fishbase..')

# # load fishbase metadata
fishbase <- species(fields = c('Species','AnaCat')) %>% # get species table for all species
  rename(binomial = Species)
#
fishbase$AnaCat <- as.factor(fishbase$AnaCat)
levels(fishbase$AnaCat) <- c(rep('Diad.',6),'Non.','Ocea.','Ocea.','Pota.','Pota.')
# table(fishbase$AnaCat)
#
# # Species range data --------------------------------------------------------------------------------------
sp_data <- readRDS(sp_ranges_file) %>%
  inner_join(.,hb_data %>% as_tibble() %>% select(HYBAS_ID,MAIN_BAS,SUB_AREA,MAIN_BAS_AREA),by="HYBAS_ID") %>%
  as.data.table(.)

# assign diadromous-non diadromous category
sp_data$diad <- 'f'
sp_data$diad[sp_data$binomial %in% fishbase$binomial[fishbase$AnaCat == 'Diad.']] <- 't'

# Dams data ---------------------------------------------------------------------------------------------
# reference dams here but keep all records and corresponding HYBAS
# then later in the function filter dams based on decision vector and select unique HYBAS IDs

dams <- read_sf(dams_file) %>%
  filter(Status %in% c('E','C')) %>%
  st_transform(4326)
dams <- st_intersection(dams,hb_data %>% select(HYBAS_ID)) %>% as_tibble() %>% select(-geom)

source('R/functions_connectivity.R')
# FUNCTION THAT CALCULATES CI PER MAIN BASIN ------------------------------------------------------------------

# Danube 2120008490

# create fictitious sp data
# sp_data <- hb_data %>%
#   as.data.table() %>%
#   filter(MAIN_BAS == 4120017020) %>%
#   mutate(binomial = 'fragile', diad = 'f') %>%
#   select('HYBAS_ID', 'binomial', 'MAIN_BAS', 'SUB_AREA', 'MAIN_BAS_AREA', 'diad')

# basin_connectivity <- function(main_bas_id){
main_bas_id <- 4120017020

sbas <- hb_data%>%
  filter(MAIN_BAS == main_bas_id)

# select dams for the current basin
dcur <- merge(sbas,dams,by='HYBAS_ID')
dams_no_cur <- 0

# if(nrow(dcur) > 0){

# divide the basins in groups
# store dams no first for data frame <<<<<<<<<<<
dams_no_cur <- nrow(dcur)
# clean duplicates
dcur <- dcur[!duplicated(dcur[,1:2]),]
# sort them, higher PFAFSTETTER number = more upstream
dcur <- dcur[order(dcur$PFAF_ID,decreasing = T),]

# find upstrea IDs of each dam 
st <- Sys.time()
master_upstream_list <- find_upstream_ids(t=as.data.frame(sbas[,c('HYBAS_ID','NEXT_DOWN')])[,-3], #removed the geom column
                                          IDs=dcur$HYBAS_ID) # dfut has both cur and fut
Sys.time() - st

# function that allocates the row to a different group based on hybas_id
# if no group is found, returns 0
assign_group <- function(hybas_id,group){
  a <- 0
  for(j in 1:length(group)){
    if(hybas_id %in% group[[j]]){
      a <- j
      break
    }
  }
  return(a) 
}

fitness <- function(x){
  nc<-24 #set no. cores
  decision <- round(x,0)
  
  # installed capacity
  totIC <- sum(dams$InstalledC*decision)
  
  ids <- dams$HYBAS_ID[decision == 1] %>% unique
  dc <- dcur[dcur$HYBAS_ID %in% ids,]
  
  # connectivity index
  upstream_list <- master_upstream_list[names(master_upstream_list) %in% dc$HYBAS_ID]
  
  # create ID groups from most upstream to downstream
  # exclude IDs in downstream groups already present in upstream groups
  # if(nrow(dc) > 0){
  groups_cur <- list()
  for(i in 1:nrow(dc)){
    # store the hybas_id of the upstream basins
    groups_cur[[i]] <- upstream_list[[as.character(dc$HYBAS_ID[i])]]
    # need to exclude basins that are in the upstream groups
    if(i > 1) groups_cur[[i]] <- groups_cur[[i]][!groups_cur[[i]] %in% do.call('c',groups_cur[1:(i-1)])]
  }
  names(groups_cur) <- dc$HYBAS_ID
  groups_cur <- groups_cur[sapply(groups_cur,function(x) length(x) != 0)]
  
  
  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  # load species data here and eventually loop through
  sbas_sp <- sp_data %>%
    filter(MAIN_BAS == main_bas_id)
  
  sps <- unique(as.character(sbas_sp$binomial))
  
  # st <- Sys.time()
  CI <- parallel::mclapply(split(sps,cut(1:length(sps),nc,labels = F)), function(sp){
    
    occ <- sbas_sp[sbas_sp$binomial == sp,]
    occ$group_cur <- sapply(occ$HYBAS_ID,function(x) assign_group(x,groups_cur)) #<<<time consuming
    
    alpha <- 0.55
    L <- occ$SUB_AREA**alpha
    
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

    
  },mc.cores = nc) %>% do.call('c',.)
  # Sys.time() - st
  
  return(c(-totIC,-mean(CI,na.rm=T)))
}

st <- Sys.time()
fitness(sample(c(0,1),56,replace = T))
Sys.time() - st

n <- nrow(dams)
st <- Sys.time()
optim <- nsga2(fn = fitness, idim = n, odim = 2, generations = 1000,
               mprob = 0.2, popsize = 40, cprob = 0.8,
               lower.bounds = rep(0, n), upper.bounds = rep(1, n))
Sys.time() - st
# plot(optim)
# plot(-optim$value[,1],-optim$value[,2])

saveRDS(optim,'proc/optim_mekong.rds')

# dec <- round(optim$par,0) %>% as.data.frame()
# ob <- -optim$value %>% as.data.frame()
# 
# # sort based on IC
# dec <- dec[sort(ob$V1,index.return=T)$ix,]
# ob <- ob[sort(ob$V1,index.return=T)$ix,]
# 
# 
# # INCLUSION PROBABILITY
# # probability of each dam being included
# dams <- read_sf('data/Dams Mekong MRC and PRC.gpkg') %>%
#   filter(Status %in% c('E','C')) %>%
#   st_transform(4326)
# dams$incl <- apply(dec,2,mean) %>% as.numeric
# 
# plot(st_geometry(dams))
# plot(dams[,'incl'])

# ws <- read_sf('')

