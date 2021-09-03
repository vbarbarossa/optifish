# packages needed
library(sf); library(foreach); library(rfishbase); library(data.table); library(dplyr); library(mco)

sp_ranges_file <- '~/surfdrive/Documents/projects/fishsuit/proc/species_ranges_raw_on_hybas12.rds' #local
# sp_ranges_file <- '~/data/data/species_ranges_raw_on_hybas12.rds'
hb_directory <- 'data/HydroBASINS/global_lev12/' #local
# hb_directory <- '~/data/data/HydroBASINS/global_lev12/'
dams_file <- 'data/Dams Mekong MRC and PRC.gpkg' #local
# dams_file <- '~/data/data/Mekong dams/Dams Mekong MRC and PRC.gpkg'

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
levels(fishbase$AnaCat) <- c(rep('Diad.',23),'Non.','Ocea.','Ocea.','Pota.','Pota.')
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
  # filter(Status %in% c('E','C')) %>%
  st_transform(4326)
sf_use_s2(FALSE)
dams <- st_intersection(dams,hb_data %>% dplyr::select(HYBAS_ID)) %>% as_tibble() %>% dplyr::select(-geom)

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

# for sedimentation
library(sf)
library(sp)
library(raster)
library(ggplot2)
library(exactextractr)
library(dplyr)

# load 'geomorphic provinces' large areas with similar sediment YIELD [t/km2/yr]
# GP<-read_sf('./R/data/Mekong geomorphic provinces/Mekong_Geomorphic_provinces.gpkg')
GP<-read_sf('~/surfdrive/Shared/optifish/GIS/Mekong dams/Geomorphic provinces/Mekong_Geomorphic_provinces.gpkg')

# hydrobasins
# HYBAS <- read_sf('./R/data/HYBAS_LEVEL_12_MEKONG/Hybas Level 12 Mekong prj.shp')
HYBAS <- read_sf('~/surfdrive/Shared/optifish/GIS/Mekong dams/Hybas as used by Rafa/Hybas Level 12 Mekong prj.shp')
# two ways of doing this 

# 1: create a raster and extract sediment yields in each HB using exat extract. Advantage: Handles situations where a hydrobasins cuts accross multiple geomoeprhic provinces. 
GPr<-raster(crs = crs(GP), vals = 0, resolution = c(5000, 5000), ext = extent(GP)) %>%
  rasterize(as_Spatial(GP), ., field = "SedYield")

HYBAS$YS_1<-exact_extract(GPr,HYBAS,fun = "mean") 
HYBAS$QS_1<-HYBAS$YS_1*HYBAS$SUB_AREA


# 2 using spatial overlay - somewhat easier and faster, but can't handle if a HYBAS intersects multiple geomorphic provinces
HYBAS$YS_2<-over(as_Spatial(HYBAS),as_Spatial(GP))["SedYield"] %>% data.matrix(.)
HYBAS$QS_2<-HYBAS$YS_2*HYBAS$SUB_AREA

library(igraph)
df <- data.frame(from = HYBAS$HYBAS_ID,to = HYBAS$NEXT_DOWN)

outlet <- as.character(df$from[which(df$to == 0)])

df[which(df$to == 0),]$to <- outlet %>% as.numeric()
g <- graph_from_data_frame(d=df)
I <- get.adjacency(g)

# 
# m <- as.matrix(I)
for(i in row.names(I)) {
  I[i,as.numeric(shortest_paths(graph=g,from = i,to=as.character(outlet))$vpath[[1]])] <- 1
}

# load dams dataset
dams$CI <- dams$GrossStora*10**6/(dams$MeanQ_m3s*60*60*24*365.25)
dams$TE_sus <- 1-0.05/sqrt(dams$CI) # trapping efficiency
dams$TE_sus[dams$TE_sus < 0] <- 0
dams$TE_sus[is.na(dams$TE_sus)] <- 0
TE_bed <- 1
f_bed <- 0.1
dams$TE <- (1-f_bed)*dams$TE_sus + f_bed*TE_bed 
dams$RE <- 1-dams$TE # release efficiency

qs <- data.frame(QS = HYBAS$QS_1)
row.names(qs) <- as.character(HYBAS$HYBAS_ID)
# remove NAs
qs[is.na(qs)] <- 0

# map outlet
outlet <- "4120017020"
out <- which(row.names(I) == outlet)

# premap dams release efficiency to network
damsMQS <- foreach(j = 1:length(dams$HYBAS_ID)) %do% {
  M <- I
  id <- dams$HYBAS_ID[j] # first three rows can be defined outside the fitness function
  up_ <- row.names(I)[I[,as.character(id)] == 1]
  dw_ <- colnames(I)[I[as.character(id),] == 1]
  M[up_,dw_] <- M[up_,dw_] * dams$RE[dams$HYBAS_ID == id]
  return(M)
}
names(damsMQS) <- dams$HYBAS_ID

# even slower to convert it to stars
# library(stars)
# rl <- lapply(damsMQS,function(x){
#   r <- raster(as.matrix(x))
#   crs(r) <- 4326
#   r <- st_as_stars(r)
#   return(r)
# })
# 
# rs = do.call('c',rl)
# rs = st_redimension(rs)


fitness <- function(x, sedim = T, nc = 24){ # make it modular to switch on and off different modules
  # nc<-24 #set no. cores
  decision <- round(x,0)
  ids <- dams$HYBAS_ID[decision == 1] %>% unique
  
  # installed capacity ----------------------------------------
  totIC <- sum(dams$InstalledC*decision)
  
  # sedimentation ---------------------------------------------
  if(sedim == T){
    
    # reoder the QS based on matrix rows
    dams_sed <- dams[decision == 1,]
    
    # MQS = I
    # # [1:n,n:most_downstream] >> [cols of I, rows of I]
    # # <<<<<<<<<< BOTTLENECK
    # st <- Sys.time()
    # for(j in 1:length(dams_sed$HYBAS_ID)){
    #   id <- dams_sed$HYBAS_ID[j] # first three rows can be defined outside the fitness function
    #   up_ <- row.names(I)[I[,as.character(id)] == 1]
    #   dw_ <- colnames(I)[I[as.character(id),] == 1]
    #   MQS[up_,dw_] <- MQS[up_,dw_] * dams_sed$RE[dams_sed$HYBAS_ID == id]
    # }
    # Sys.time() - st
    
    # same speed than above
    MQS = Reduce('*',damsMQS[as.character(dams_sed$HYBAS_ID)])
    
    # data.frame(HYBAS_ID = as.numeric(colnames(MQS)),
    #              QS_sum = as.numeric(qs[row.names(MQS),'QS'] %*% MQS))
    totQS <- as.numeric(qs[row.names(MQS),'QS'] %*% MQS)[out]
  }
  
  
  # fragmentation ----------------------------------------------
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
  
  if(sedim == T){
    return(c(-totIC,-totQS,-mean(CI,na.rm=T))) 
  }else{
    return(c(-totIC,-mean(CI,na.rm=T)))
  }
}

st <- Sys.time()
fitness(sample(c(0,1),nrow(dams),replace = T),nc=8)
Sys.time() - st

# max computation time ~7.5 sec
st <- Sys.time()
fitness(rep(1,nrow(dams)), sedim = T, nc=8)
Sys.time() - st

# without sedimentation part ~0.6 sec
st <- Sys.time()
fitness(rep(1,nrow(dams)), sedim = F, nc=8)
Sys.time() - st

n <- nrow(dams)
st <- Sys.time()
optim <- nsga2(fn = fitness, idim = n, odim = 3, generations = 1,
               mprob = 0.2, popsize = 12, cprob = 0.8,
               lower.bounds = rep(0, n), upper.bounds = rep(1, n))
Sys.time() - st
# plot(optim)
# plot(-optim$value[,1],-optim$value[,2])

saveRDS(optim,'proc/optim_mekong_pop80_gen5000.rds')

dec <- round(optim$par,0) %>% as.data.frame()
ob <- -optim$value %>% as.data.frame()

# sort based on IC
dec <- dec[sort(ob$V1,index.return=T)$ix,]
ob <- ob[sort(ob$V1,index.return=T)$ix,]


# INCLUSION PROBABILITY
# probability of each dam being included
dams <- read_sf('data/Dams Mekong MRC and PRC.gpkg') %>%
  # filter(Status %in% c('E','C')) %>%
  st_transform(4326)
dams$incl <- apply(dec,2,mean) %>% as.numeric

plot(st_geometry(dams))
plot(dams[,'incl'])

# ws <- read_sf('')
ggplot(ob) + geom_point(aes(x = V1, y = V2, color = V3)) + xlab('IC') + ylab('sediments')


# next steps
# 1 - check how to improve the speed performance of sediment part in fitness function, probably need to pre-map the dams matrices and simply multiply them in the fitness function - check how to best multiply matrices in R
# 2 - check whether volume or flow deviation (storage/streamflow) is the best way to substitute IC
# 3 - do the pre-filtering on dams based on whether they occurr on the tributary of the hb unit (see google docs)
# 4 - run on Alice


