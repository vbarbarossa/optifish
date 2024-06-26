# source('R/MASTER.R')
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # see here: https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location-in-rstudio

# packages needed
library(sf); library(foreach);

library(data.table); library(dplyr); 
# library(vroom)

# HydroBASINS data ------------------------------------------------------------------------------------------
# read hydrobasins data
hb_data <- # read_sf('data/HydroBASINS/global_lev12/hybas_as_lev12_v1c.shp')
  read_sf('data/HYBAS_LEVEL_12_MEKONG/Hybas Level 12 Mekong prj.shp')

hb_data <-st_transform(x=hb_data,crs=4326)

# add basin area
# main_bas_area <- do.call('rbind',lapply(split(hb_data_frame,hb_data_frame$MAIN_BAS),function(x) data.frame(MAIN_BAS = unique(x$MAIN_BAS),MAIN_BAS_AREA = sum(x$SUB_AREA))))
cat('\nCompiling main basin area..')

main_bas_area <- hb_data %>%
  as_tibble() %>%
  select(HYBAS_ID,MAIN_BAS,SUB_AREA) %>%
  group_by(MAIN_BAS) %>%
  summarize(MAIN_BAS_AREA = sum(SUB_AREA))

hb_data <- inner_join(hb_data,main_bas_area,by='MAIN_BAS') %>%
  filter(MAIN_BAS %in% c(4120017020))

# select diadromous and non-diadromous species------------------------------------------------------------

cat('\nRetrieving diadromous species from fishbase..')

# Dams data ---------------------------------------------------------------------------------------------
# reference dams here but keep all records and corresponding HYBAS
# then later in the function filter dams based on decision vector and select unique HYBAS IDs

dams <- read_sf('data/Dam data/Dams Mekong MRC and PRC.gpkg') %>%
  filter(Status %in% c('E','C')) %>%
  st_transform(4326)
dams <- st_intersection(dams,hb_data %>% select(HYBAS_ID)) %>% as_tibble() %>% select(-geom)

source('functions_connectivity.R')
# FUNCTION THAT CALCULATES CI PER MAIN BASIN ------------------------------------------------------------------

# Danube 2120008490

# create fictitious sp data
sp_data <- hb_data %>%
  as.data.table() %>%
  filter(MAIN_BAS == 4120017020) %>%
  mutate(binomial = 'fragile', diad = 'f') %>%
  select('HYBAS_ID', 'binomial', 'MAIN_BAS', 'SUB_AREA', 'MAIN_BAS_AREA', 'diad')

# basin_connectivity <- function(main_bas_id){
main_bas_id <- 4120017020

sbas <- hb_data%>%
  filter(MAIN_BAS == main_bas_id)

# select dams for the current basin
dcur <- merge(sbas,dams,by='HYBAS_ID')
# dfut <- merge(sbas,dams,by='HYBAS_ID')

dams_no_cur <- 0
# dams_no_fut <- 0

# if(nrow(dcur) > 0){

# divide the basins in groups
# store dams no first for data frame <<<<<<<<<<<
dams_no_cur <- nrow(dcur)
# dams_no_fut <- nrow(dfut)

# clean duplicates
dcur <- dcur[!duplicated(dcur[,1:2]),]
# dfut <- rbind(dcur[,1:9],dfut[,1:9])
# dfut <- dfut[!duplicated(dfut[,1:2]),]

# sort them, higher PFAFSTETTER number = more upstream
dcur <- dcur[order(dcur$PFAF_ID,decreasing = T),]
# dfut <- dfut[order(dfut$PFAF_ID,decreasing = T),]

# find upstrea IDs of each dam 
# <<<<<<<<<<<<<<<<<<<<<<<<< TIME CONSUMING BOTTLENECK <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# <<<<<<<<<<<<<<<<<<<<<<<<< can map the dams connectivity once and then update it based on the scenarios?
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

st <- Sys.time()
master_upstream_list <- find_upstream_ids(t=as.data.frame(sbas[,c('HYBAS_ID','NEXT_DOWN')])[,-3], #removed the geom column
                                   IDs=dcur$HYBAS_ID) # dfut has both cur and fut
Sys.time() - st

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# when running changed configuration using the same upstream list, need to filter out removed dams


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
  
  tab <- foreach(sp = unique(as.character(sbas_sp$binomial)),.combine = 'rbind') %do% {
    
    occ <- sbas_sp[sbas_sp$binomial == sp,]
    occ$group_cur <- sapply(occ$HYBAS_ID,function(x) assign_group(x,groups_cur)) #<<<time consuming
    # occ$group_fut <- sapply(occ$HYBAS_ID,function(x) assign_group(x,groups_fut)) #<<<time consuming
    
    foreach(alpha = c(0.55),.combine = 'rbind') %do% {
      
      L <- occ$SUB_AREA**alpha
      
      if(occ$diad[1] == 't'){
        # total area
        A = sum(L)
        sa_cur = sum(L[occ$group_cur == min(occ$group_cur)])
        # sa_fut = sum(L[occ$group_fut == min(occ$group_fut)])
        cat <- 'diadromous'
      }else{
        # total area
        A = sum(L)**2
        # sum pf areas for patches from different groups
        sa_cur = 0
        # sa_fut = 0
        for(i in unique(occ$group_cur)) sa_cur = sa_cur + sum(L[occ$group_cur == i])**2
        # for(i in unique(occ$group_fut)) sa_fut = sa_fut + sum(L[occ$group_fut == i])**2
        cat <- 'potamodromous'
        
      }
      
      # output
      # main basin id, species name, no.patches, sum area patches, connectivity current, connectivity future
      
      data.frame(
        MAIN_BAS = unique(occ$MAIN_BAS),
        MAIN_BAS_AREA = occ$MAIN_BAS_AREA[1],
        binomial = sp,
        patches.no = nrow(occ),
        patches.cum.area = sum(occ$SUB_AREA),
        dams.cur.no = dams_no_cur,
        # dams.fut.no = dams_no_fut,
        dams.cur.no.bas = nrow(dcur),
        # dams.fut.no.bas = nrow(dfut) - nrow(dcur),
        category = cat,
        connectivity.cur = (sa_cur/A)*100,
        # connectivity.fut = (sa_fut/A)*100,
        alpha = alpha
      )
      
    }
  }
  
  # here plug in connectivity
  # totCI <- a number based on the config of the dams
  
  return(c(-totIC,-tab$connectivity.cur))
  
}

st <- Sys.time()
fitness(sample(c(0,1),56,replace = T))
Sys.time() - st



n <- nrow(dams)
library(mco)
st <- Sys.time()
optim <- nsga2(fn = fitness, idim = n, odim = 2, generations = 20,
               mprob = 0.2, popsize = 24, cprob = 0.8,
               lower.bounds = rep(0, n), upper.bounds = rep(1, n))
Sys.time() - st
# plot(optim)
plot(-optim$value[,1],-optim$value[,2])


dec <- round(optim$par,0) %>% as.data.frame()
ob <- -optim$value %>% as.data.frame()

# sort based on IC
dec <- dec[sort(ob$V1,index.return=T)$ix,]
ob <- ob[sort(ob$V1,index.return=T)$ix,]


# INCLUSION PROBABILITY
# probability of each dam being included
dams <- read_sf('data/Dams Mekong MRC and PRC.gpkg') %>%
  filter(Status %in% c('E','C')) %>%
  st_transform(4326)
dams$incl <- apply(dec,2,mean) %>% as.numeric

plot(st_geometry(dams))
plot(dams[,'incl'])

# ws <- read_sf('')

#### Test impacts of hyper parameters on model performance and runtime----
library(hypervolume)
library(scales)


generations<-list()
popsize<-list()
runtime<-list()
dec<-list()
ob<-list()
hv<-list()

i<-0 #counter 

for (ngen in c(100,1000,10000,100000)){ # number of generation
  
    for (ps in seq(from=100, to = 1000, by = 100)){ # initial pop size
      
    i<-i+1
    print(paste('run ', i))
    
      
      timeStart <- Sys.time()
      optim <- nsga2(fn = fitness, idim = n, odim = 2, generations = 20,
                     mprob = 0.2, popsize = 24, cprob = 0.8,
                     lower.bounds = rep(0, n), upper.bounds = rep(1, n))
      timeEnd<-Sys.time() 
      
      # save: 
      generations[[i]]<-ngen # number of generations
      popsize[[i]]<-ps #population size
      runtime[[i]]<-as.numeric(difftime(timeEnd, timeStart, units='mins')) # runtime
      dec[[i]]<-round(optim$par,0) # decisionms
      ob[[i]]<- -optim$value %>% as.data.frame() # objective values 
      
      # hyper volumes: normalize first....
      hv[[i]]<-hypervolume(method='box', 
                           cbind(rescale(ob[[i]]$V1,to =c(0,1)),
                             rescale(ob[[i]]$V2,to =c(0,1))))@Volume
      
      
    }
  
}

# backup 
ob_orig<-ob

hyperP_orig<-data.frame(popsize=as.numeric(popsize),
                        generations=as.numeric(generations),
                        runtime=as.numeric(runtime),
                        hypervolume=as.numeric(hv))

# plotting ----

# Make dataframe collection results for all trials 
hyperP_results<-data.frame(popsize=as.numeric(popsize),
                              generations=as.numeric(generations),
                              runtime=as.numeric(runtime),
                              hypervolume=as.numeric(hv))

# create unique factor labels for each run
hyperP_results$labels<-as.factor(paste(
  paste( 'popsize', as.character(hyperP_results$popsize)),
  paste(as.character(hyperP_results$generations),'generations')))

# create an order for the labels (i.e., from small n generation to large, same for popsize. Required to later sort the paretofronts on ggplot)
hyperP_results$labelorder<-with(hyperP_results, order(generations,popsize))


library(ggplot2)
library(forcats)
library(reshape2)
library(patchwork)
library(pals)


#Plot sensitivity of hyper parameters
  
  p1<-ggplot(data=hyperP_results, aes(x=popsize,y=runtime,color=as.factor(generations)))+
    geom_line()
  p2<-ggplot(data=hyperP_results, aes(x=popsize,y=hypervolume,color=as.factor(generations)))+
    geom_line()
  
  p3<-p1+p2

  ggsave('NSGAII benchmarking.png',p3)

# Create data from the list of PO portfolios of each rum
df<-melt(ob,id.vars = c('V1','V2'))
df$popsize<-hyperP_results$popsize[df$L1]
df$generations<-hyperP_results$generations[df$L1]

# assign labels for plotting and resort the labels (wchih are factors) according to popsize and generations. Otherwise, colors are associated to the set of PO point from eahc run in a random way. 
df$labels<-hyperP_results$labels[df$L1]
df$labelorder<-hyperP_results$labelorder[df$L1]

# factor ordering
 df%>%mutate(labels = fct_reorder(labels, labelorder))

# Plot pareto optimal tradeoffs for each run.  
p4<-ggplot(data=df, aes(x=V1, y=V2,color=labels))+
  geom_point()+scale_colour_manual(name = "Species Names", values = rev(plasma(length(unique(df$labels)))))

ggsave('NSGAII benchmarking PO solutions.png',p4)
