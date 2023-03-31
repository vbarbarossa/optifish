# packages needed
library(sf); library(foreach); library(rfishbase); 
library(data.table); library(dplyr); library(exactextractr); 
library(igraph); library(ggplot2); library(rnaturalearth)
# raster package is also needed

# settings for sf
sf_use_s2(FALSE)

# retrieve paths to input files
source('R/master_paths_local.R')

# functions for CI calculation
source('R/functions_connectivity.R')

# PARETO FRONTS ----------------------------------------------------------------

# calculate current dams based on fitness
all_dams = T
name_col_IC <- 'InstalledC'
name_col_V <- 'GrossStora'
name_col_Q <- 'MeanQ_m3s' # flow of the river -> to calculate trap efficiency 
TE_bed <- 1 # trapping efficiency for coarse bed material
f_bed <- 0.1 # fraction of coarse bed material in the total sediments
main_bas_id <- outlet <- 4120017020 #this corresponds to the outlet!
source('R/prepare_data.R')

# simplify dams table for fitness function
dams_data <- dams[,c('Code','INTER_ID','Status','DamHeight','ReferenceY',name_col_IC, name_col_V)]
dams_data$ReferenceY[dams_data$ReferenceY == 0] <- min(dams_data$ReferenceY[dams_data$ReferenceY != 0])
# # assign passability
dams_data$pass <- 0

# lots of P dams have reference year = 2007 or higher
dams_data$ReferenceY[dams_data$Status == 'P'] <- 2050
dams_e <- dams_data[dams_data$Status != 'P',]
dams_f <- dams_data[dams_data$Status == 'P',]

# no species
n_sp <- round(sum(L) + sum(ld),0)

# no dams
n = nrow(dams_f)
n = n + nrow(dams_e)

sedimentation = F
fragmentation = T
energy = T
water = F

calc_objs <- 
  function(x){ 
    
    # select dams based on index
    decision <- x
    
    if(all_dams){
      dams_tot <- rbind(dams_e,dams_f)[decision == 1,] 
    }else{
      dams_tot <- rbind(dams_e,dams_f[decision == 1,]) 
    }
    
    dams_inter_ids <- sort(unique(dams_tot$INTER_ID))
    
    # result vector
    result_fitness <- numeric()
    
    # installed capacity ----------------------------------------
    if(energy) result_fitness <- c(result_fitness,sum(pull(dams_tot,name_col_IC)))
    
    # volume ----------------------------------------
    if(water) result_fitness <- c(result_fitness,sum(pull(dams_tot,name_col_V)))
    
    # sedimentation ---------------------------------------------
    if(sedimentation){
      if(length(dams_inter_ids) > 0){
        # calculate total RE matrix (= multiply single dam RE matrices)
        MQS = Reduce('*',RE_M_array[as.character(dams_tot$INTER_ID)])
        
        # multiply RE by load to get final load at mouth
        totQS <- as.numeric(qs[row.names(MQS),'QS'] %*% MQS)[out]
      }else{
        totQS <- sum(qs$QS)
      }
      result_fitness <- c(result_fitness,totQS)
    }
    
    
    # fragmentation ----------------------------------------------
    
    if(fragmentation){
      
      # calculate passability of all dams
      d <- left_join(
        bas_tot,
        # multiply dams occupying the same basin
        summarize(group_by(dams_tot,INTER_ID),pass = prod(pass)),
        by = "INTER_ID"
      )
      d$pass[is.na(d$pass)] <- 1 # assign 1 to dams not included in this set
      
      # calculate the Cij matrix
      Cij <- exp(distances(g_master_pass, mode = 'all', weights = (-1*log(d$pass))) * -1)
      
      totCI <- (sum(Cij * L,na.rm=T) + sum(Cij[,1]*ld))/n_sp
      
      result_fitness <- c(result_fitness,totCI) 
    }
    
    
    # Return output --------------------------------------------------------------
    return( -result_fitness )
    
  }

years_range <- min(dams_e$ReferenceY):max(dams_e$ReferenceY)
op_pres <- foreach(i = years_range,.combine = 'rbind') %do% {
  -calc_objs(c(as.integer(dams_e$ReferenceY <= i),rep(0,nrow(dams_f))))
} %>% as.data.frame()
row.names(op_pres) <- 1:nrow(op_pres)
colnames(op_pres) <- c('ic','ci')
op_pres$scenario <- 'present'
op_pres$year <- c(years_range)


gen_size <- 1000
op <- foreach(i = c('pristine','pristine_pass_bin','pristine_pass_step','future','removal','mitigation_bin','mitigation_step'),.combine = 'rbind') %do% {
  
  t1 <- readRDS(paste0('proc/nsga2_',i,'_ic_ci_gen',gen_size,'_pop100.rds'))
  t1 <- -t1@fitness %>% as.data.frame()
  t1$scenario <- i
  return(t1)
}
colnames(op)[1:2] <- c('ic','ci')
op$year <- ''
op <- rbind(op_pres,op)

op$scenario <- factor(op$scenario, levels = c('pristine','pristine_pass_bin','pristine_pass_step','present','future','removal','mitigation_bin','mitigation_step'))
op$ci <- op$ci*100
p <- ggplot(data = op,
            aes(x = ic, y = ci, color = scenario, label = year)) +
  geom_point(alpha = 1) +
  geom_step(direction = 'hv') +
  # ggrepel::geom_text_repel(force = 1, direction = 'both', max.iter = 100000) +
  # geom_text(hjust=1, vjust=1) +
  scale_color_manual(values = RColorBrewer::brewer.pal(8,'Dark2')[c(1:length(levels(op$scenario)))]) +
  labs(color = '') +
  xlab('Installed Capacity [MW]') +
  ylab('Connectivity Index [%]') +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  coord_cartesian(expand = T) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'right',
        strip.background = element_blank(), strip.placement = 'outside',
        text = element_text(size = 13),
        strip.text = element_text(size = 13)
  )
p


ggsave('figs/pareto_scenarios_all.jpg',p,width = 150, height = 120, units = 'mm',dpi = 600)

plot_scenarios <- function(vars = c('pristine','present','future'),
                           save_str = 'figs/pareto_scenarios_pres_fut.jpg',
                           legend_pos = 'none'){
  op_p <- op %>% filter(scenario %in% vars)
  p <- ggplot(data = op_p,
              aes(x = ic, y = ci, color = scenario, label = year)) +
    geom_point(alpha = 1) +
    geom_step(direction = 'hv') +
    # ggrepel::geom_text_repel(force = 1, direction = 'both', max.iter = 100000) +
    # geom_text(hjust=1, vjust=1) +
    scale_color_manual(values = RColorBrewer::brewer.pal(8,'Dark2')[c(1:length(levels(op_p$scenario)))]) +
    labs(color = '') +
    xlab('Installed Capacity [MW]') +
    ylab('Connectivity Index [%]') +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    coord_cartesian(expand = T) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = legend_pos,
          strip.background = element_blank(), strip.placement = 'outside',
          text = element_text(size = 13),
          strip.text = element_text(size = 13)
    )
  p
  
  ggsave(save_str,p,width = 120, height = 120, units = 'mm',dpi = 600)
  
}

plot_scenarios(vars = c('pristine','present','future'),
               save_str = 'figs/pareto_scenarios_pres_fut.jpg')

plot_scenarios(vars = c('pristine','present','future','removal'),
               save_str = 'figs/pareto_scenarios_pres_fut_rem.jpg')

plot_scenarios(vars = c('pristine_pass_step','present','future','removal','mitigation_step'),
               save_str = 'figs/pareto_scenarios_pres_fut_rem_mit_step.jpg')



op <- readRDS('~/surfdrive/tmp/optifish_proc_20220830/nsga2_allDams_ic_ci_gen50000_pop100.rds')
dec <- op@population %>% as.data.frame() # decision
ob <- -op@fitness %>% as.data.frame() # objective
# # sort based on IC <<< WHY? ask Rafa
dec <- dec[sort(ob$V1,index.return=T)$ix,]
ob <- ob[sort(ob$V1,index.return=T)$ix,]

# look at the inclusion probability of portfolios up to total existing IC
dec[ob$V1 <= op_pres$ic[nrow(op_pres)],] %>% apply(.,2,mean) %>% as.numeric


# INCLUSION PROBABILITY MAP ----------------------------------------------------

# dams data
dams <- read_sf('data/Dams Mekong MRC and PRC.gpkg') %>%
  st_transform(4326)
# inclusion probability based on optimization run
dams$incl <- apply(dec,2,mean) %>% as.numeric





# basins layer
bas <- read_sf('proc/basins_mekong.gpkg') %>% group_by(MAIN_BAS) %>% summarise(do_union = T)

# rivers layer
riv <- read_sf('proc/rivers_mekong.gpkg') %>% 
  filter(DIS_AV_CMS > 10)
riv$Log_Q_avg <- log10(riv$DIS_AV_CMS)

riv1 <- riv[riv$Log_Q_avg < 1.5,]
riv2 <- riv[riv$Log_Q_avg >=1.5 & riv$Log_Q_avg <2,]
riv3 <- riv[riv$Log_Q_avg >=2 & riv$Log_Q_avg <3,]
riv4 <- riv[riv$Log_Q_avg >=3 & riv$Log_Q_avg <4,]
riv5 <- riv[riv$Log_Q_avg >=4,]

dams$log_IC <- log10(dams$InstalledC)
dams$incl <- dams$incl*100

dams$status <- ifelse(dams$Status == 'P','future','existing')

p_base <- ggplot() +
  geom_sf(data = bas, color = 'grey20',fill='grey90') + #no border
  geom_sf(data = riv1,size = 0.1,color = 'grey40') +
  geom_sf(data = riv2,size = 0.3,color = 'grey40') +
  geom_sf(data = riv3,size = 0.6,color = 'grey40') +
  geom_sf(data = riv4,size = 0.9,color = 'grey40') +
  geom_sf(data = riv5,size = 1.2,color = 'grey40') +
  geom_sf(data = dams, aes(size = InstalledC, fill = status), shape = 21,alpha = 0.8) +
  scale_fill_manual(values = c('white','red')) +
  # scale_fill_continuous(type = 'viridis', limits = c(0,100), breaks = c(0, 25,50,75,100),labels = c(0,25,50,75,100)) +
  labs(fill = '', size = 'Inst. Cap. [MW]') +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    panel.grid = element_line(colour = 'transparent'),
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    legend.position = 'right'
  )
p_base

# inclusion probability
p_incl <- ggplot() +
  geom_sf(data = bas, color = 'grey20',fill='grey90') + #no border
  geom_sf(data = riv1,size = 0.1,color = 'grey40') +
  geom_sf(data = riv2,size = 0.3,color = 'grey40') +
  geom_sf(data = riv3,size = 0.6,color = 'grey40') +
  geom_sf(data = riv4,size = 0.9,color = 'grey40') +
  geom_sf(data = riv5,size = 1.2,color = 'grey40') +
  geom_sf(data = dams, aes(size = InstalledC,fill = incl), shape = 21,alpha = 0.8) +
  scale_fill_viridis_c(option = 'C', limits = c(0,100),
                       breaks = c(0, 25,50,75,100),labels = c(0,25,50,75,100),
                       begin = 0, end = 1) +
  # scale_fill_continuous(type = 'viridis', limits = c(0,100), breaks = c(0, 25,50,75,100),labels = c(0,25,50,75,100)) +
  labs(fill = 'Inclusion prob. [%]', size = 'Inst. Cap. [MW]') +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    panel.grid = element_line(colour = 'transparent'),
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    legend.position = 'right'
  )
p_incl

ggsave('figs/inclusion_map_mekong_ic_ci.jpg',p_incl,width = 150, height = 200, units = 'mm',dpi = 600)


# BASINS LAYER MEKONG ----------------------------------------------------------
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
# write_sf(hb_data,'proc/basins_mekong.gpkg',delete_layer = T)

# RIVERS LAYER MEKONG ----------------------------------------------------------
# riv <- read_sf('data/RiverATLAS/RiverATLAS_v10_as.shp')
# riv_mekong <- riv %>% filter(HYBAS_L12 %in% hb_data$HYBAS_ID)
# write_sf(riv_mekong,'proc/rivers_mekong.gpkg')

# INCLUSION PROBABILITY MAP ----------------------------------------------------

# dams data
dams <- read_sf('data/Dams Mekong MRC and PRC.gpkg') %>%
  st_transform(4326)
# inclusion probability based on optimization run
dams$incl <- apply(dec,2,mean) %>% as.numeric

# basins layer
bas <- read_sf('proc/basins_mekong.gpkg') %>% group_by(MAIN_BAS) %>% summarise(do_union = T)

# rivers layer
riv <- read_sf('proc/rivers_mekong.gpkg') %>% 
  filter(DIS_AV_CMS > 10)
riv$Log_Q_avg <- log10(riv$DIS_AV_CMS)

riv1 <- riv[riv$Log_Q_avg < 1.5,]
riv2 <- riv[riv$Log_Q_avg >=1.5 & riv$Log_Q_avg <2,]
riv3 <- riv[riv$Log_Q_avg >=2 & riv$Log_Q_avg <3,]
riv4 <- riv[riv$Log_Q_avg >=3 & riv$Log_Q_avg <4,]
riv5 <- riv[riv$Log_Q_avg >=4,]

dams$log_IC <- log10(dams$InstalledC)
dams$incl <- dams$incl*100

dams$status <- ifelse(dams$Status == 'P','future','existing')

p_base <- ggplot() +
  geom_sf(data = bas, color = 'grey20',fill='grey90') + #no border
  geom_sf(data = riv1,size = 0.1,color = 'grey40') +
  geom_sf(data = riv2,size = 0.3,color = 'grey40') +
  geom_sf(data = riv3,size = 0.6,color = 'grey40') +
  geom_sf(data = riv4,size = 0.9,color = 'grey40') +
  geom_sf(data = riv5,size = 1.2,color = 'grey40') +
  geom_sf(data = dams, aes(size = InstalledC, fill = status), shape = 21,alpha = 0.8) +
  scale_fill_manual(values = c('white','red')) +
  # scale_fill_continuous(type = 'viridis', limits = c(0,100), breaks = c(0, 25,50,75,100),labels = c(0,25,50,75,100)) +
  labs(fill = '', size = 'Inst. Cap. [MW]') +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    panel.grid = element_line(colour = 'transparent'),
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    legend.position = 'right'
  )
p_base

# richness data
hb_data <- foreach(i = c('as'),.combine = 'rbind') %do% read_sf(paste0(hb_directory,'/global_lev12/hybas_',i,'_lev12_v1c.shp'))
main_bas_area <- hb_data %>%
  as_tibble() %>%
  select(HYBAS_ID,MAIN_BAS,SUB_AREA) %>%
  group_by(MAIN_BAS) %>%
  summarize(MAIN_BAS_AREA = sum(SUB_AREA))
hb_data <- inner_join(hb_data,main_bas_area,by='MAIN_BAS') %>%
  filter(MAIN_BAS %in% main_bas_id)
# retrieve threat status and commerical relevance
sr <- read.csv('proc/sp_data_hybas12_mekong.csv') %>% 
  as_tibble() %>%
  group_by(HYBAS_ID) %>%
  summarise(SR = length(unique(binomial))) %>%
  left_join(hb_data,.)

sp <- read.csv('proc/sp_data_hybas12_mekong.csv') %>% 
  as_tibble() %>%
  right_join(hb_data,., by = 'HYBAS_ID') %>%
  group_by(binomial) %>%
  summarise(do_union = T)
sp$area <- st_area(sp)
write_sf(sp,'proc/species_ranges_mekong.gpkg')

plot(st_geometry(sp[sp$binomial == 'Pangasianodon gigas',])) #Giant catfish
plot(st_geometry(sp[sp$binomial == 'Amblyceps foratum',])) #
plot(st_geometry(sp[sp$binomial == 'Chitala blanci',])) #Royal knifefish

p_fish <- ggplot() +
  geom_sf(data = bas, color = 'grey20',fill='grey90') + #no border
  geom_sf(data = riv1,size = 0.1,color = 'grey40') +
  geom_sf(data = riv2,size = 0.3,color = 'grey40') +
  geom_sf(data = riv3,size = 0.6,color = 'grey40') +
  geom_sf(data = riv4,size = 0.9,color = 'grey40') +
  geom_sf(data = riv5,size = 1.2,color = 'grey40') +
  geom_sf(data = sp[sp$binomial == 'Barbodes rhombeus',], fill='blue', color = 'transparent', alpha = 0.5) +
  geom_sf(data = sp[sp$binomial == 'Pangasianodon gigas',], fill='red', color = 'transparent', alpha = 0.5) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    panel.grid = element_line(colour = 'transparent'),
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    legend.position = 'right'
  )
p_fish


p_sr <- ggplot() +
  geom_sf(data = bas, color = 'grey20',fill='grey90') + #no border
  geom_sf(data = sr, aes(fill = SR), color = 'transparent', alpha = 1) +
  scale_fill_viridis_c(option = 'D',
                       # breaks = c(0, 25,50,75,100),labels = c(0,25,50,75,100),
                       begin = 0, end = 1) +
  labs(fill = 'No. species') +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    panel.grid = element_line(colour = 'transparent'),
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    legend.position = 'right'
  )
p_sr

p_incl <- ggplot() +
  geom_sf(data = bas, color = 'grey20',fill='grey90') + #no border
  geom_sf(data = riv1,size = 0.1,color = 'grey40') +
  geom_sf(data = riv2,size = 0.3,color = 'grey40') +
  geom_sf(data = riv3,size = 0.6,color = 'grey40') +
  geom_sf(data = riv4,size = 0.9,color = 'grey40') +
  geom_sf(data = riv5,size = 1.2,color = 'grey40') +
  geom_sf(data = dams, aes(size = InstalledC,fill = incl), shape = 21,alpha = 0.8) +
  scale_fill_viridis_c(option = 'C', limits = c(0,100),
                       breaks = c(0, 25,50,75,100),labels = c(0,25,50,75,100),
                       begin = 0, end = 1) +
  # scale_fill_continuous(type = 'viridis', limits = c(0,100), breaks = c(0, 25,50,75,100),labels = c(0,25,50,75,100)) +
  labs(fill = 'Inclusion prob. [%]', size = 'Inst. Cap. [MW]') +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    panel.grid = element_line(colour = 'transparent'),
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    legend.position = 'right'
  )
p_incl

ggsave('figs/sr_map_mekong.jpg',p_sr,width = 150, height = 200, units = 'mm',dpi = 600)
ggsave('figs/dams_map_mekong.jpg',p_base,width = 150, height = 200, units = 'mm',dpi = 600)
ggsave('figs/inclusion_map_mekong.jpg',p_incl,width = 150, height = 200, units = 'mm',dpi = 600)
ggsave('figs/fish_example_map_mekong.jpg',p_fish,width = 150, height = 200, units = 'mm',dpi = 600)


# BASIN LOCATION MAP -----------------------------------------------------------
# plot basin location on globe
crs_custom <- "+proj=ortho +lat_0=20 +lon_0=100"
land <- ne_countries(scale = 'small', type = "countries", returnclass = "sf") %>%
  st_transform(crs_custom)

bas_box <- st_as_sfc(st_bbox(bas) + c(-5,-5,5,5)) 

boundary <- st_bbox(ne_countries(scale = 'small', type = "countries", returnclass = "sf")) %>%
  st_as_sfc() %>%
  st_transform(crs_custom)

p_wrld <- ggplot()  +
  geom_sf(data = land, fill="grey", color = NA) +
  geom_sf(data = bas_box, fill = NA, color = "red") +
  geom_sf(data = bas, fill = "red", color = NA) +
  theme_minimal()
p_wrld

ggsave('figs/location_mekong.jpg',p_wrld,width = 50, height = 50, units = 'mm',dpi = 600)



# calculate simplified CI for each dam configuration
source('R/objectives_separate.R')

# for simplified_ci
hb_data <- foreach(i = c('as'),.combine = 'rbind') %do% read_sf(paste0(hb_directory,'/global_lev12/hybas_',i,'_lev12_v1c.shp'))
main_bas_area <- hb_data %>%
  as_tibble() %>%
  select(HYBAS_ID,MAIN_BAS,SUB_AREA) %>%
  group_by(MAIN_BAS) %>%
  summarize(MAIN_BAS_AREA = sum(SUB_AREA))
hb_data <- inner_join(hb_data,main_bas_area,by='MAIN_BAS') %>%
  filter(MAIN_BAS %in% main_bas_id)

tob$ci <- tob$ci/length(unique(read.csv('proc/sp_data_hybas12_mekong.csv')$binomial))

# ggplot(tob) +
#   geom_point(aes(x = ic, y = ci_m)) +
#   geom_point(aes(x = ic, y = ci_s), color = 'red')



token <- 'd361026f05b472e57b0ffe1fa5c9a768aaf3d8391abbb464293e9efe2bbbf733'
library(rredlist)
iucn_code <- foreach(ts = c("DD", "LC", "NT", "VU", "EN","CR", "EW", "EX", "LRlc", "LRnt", "LRcd"),.combine = 'rbind') %do%{
  t <- rl_sp_category(ts,key = token)$result %>% 
    as_tibble() %>%
    mutate(code = ts) %>%
    select(binomial = scientific_name,code)
  return(t)
} %>% arrange(binomial)

# retrieve threat status and commerical relevance
sp_data <- read.csv('proc/sp_data_hybas12_mekong.csv') %>% 
  as_tibble() %>%
  left_join(rfishbase::species() %>% select(Species, Importance), by = c('binomial' = 'Species')) %>%
  left_join(iucn_code)

# separate diadromous-nondiadromous
# threatened species according to IUCN
# commercially relevant species
sp_data_simple <- hb_data %>% as_tibble() %>% select(HYBAS_ID, MAIN_BAS, SUB_AREA, MAIN_BAS_AREA) %>% mutate(binomial = 'sp1', diad = 'f') %>% as.data.table()
tob$ci_simple <- apply(dec,1,function(x) ci(x,sp_data_simple))
sp_data_commercial <- sp_data %>% filter(Importance %in% c('commercial','highly commercial','minor commercial','subsistence fisheries'))
tob$ci_commercial <- apply(dec,1,function(x) ci(x,sp_data_commercial)/length(unique(sp_data_commercial$binomial)))
sp_data_threatened <- sp_data %>% filter(code %in% c('VU','EN','CR'))
tob$ci_threatened <- apply(dec,1,function(x) ci(x,sp_data_threatened)/length(unique(sp_data_threatened$binomial)))
sp_data_diad <- sp_data %>% filter(diad == 't')
sp_data_nondiad <- sp_data %>% filter(diad == 'f')
tob$ci_diadromous <- apply(dec,1,function(x) ci(x,sp_data_diad)/length(unique(sp_data_diad$binomial)))
tob$ci_nondiadromous <- apply(dec,1,function(x) ci(x,sp_data_nondiad)/length(unique(sp_data_nondiad$binomial)))

tob_l = tidyr::pivot_longer(tob, cols = starts_with('ci'), names_prefix = 'ci_')
tob_l$name[tob_l$name == 'ci'] <- 'all'

# all, commercial, threatened
p1 <- ggplot(data = tob_l %>% filter(name %in% c('all','commercial','threatened')),
             aes(x = ic, y = value, color = name)) +
  geom_point(alpha = 1) +
  # geom_smooth(se = F) +
  scale_color_manual(values = RColorBrewer::brewer.pal(9,'Set1')[c(9,5,1)]) +
  labs(color = '') +
  xlab('Installed Capacity [MW]') +
  ylab('Connectivity Index [%]') +
  coord_cartesian(ylim = c(0,100), expand = F) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'top')

# all, simple
p2 <- ggplot(data = tob_l %>% filter(name %in% c('all','simple')),
             aes(x = ic, y = value, color = name)) +
  geom_point(alpha = 1) +
  # geom_smooth(se = F) +
  scale_color_manual(values = RColorBrewer::brewer.pal(9,'Set1')[c(9,2)]) +
  labs(color = '') +
  xlab('Installed Capacity [MW]') +
  ylab('Connectivity Index [%]') +
  coord_cartesian(ylim = c(0,100), expand = F) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'top')

# all, diadromous, nondiadromous
p3 <- ggplot(data = tob_l %>% filter(name %in% c('all','diadromous','nondiadromous')),
             aes(x = ic, y = value, color = name)) +
  geom_point(alpha = 1) +
  # geom_smooth(se = F) +
  scale_color_manual(values = RColorBrewer::brewer.pal(9,'Set1')[c(9,3:4)]) +
  labs(color = '') +
  xlab('Installed Capacity [MW]') +
  ylab('Connectivity Index [%]') +
  coord_cartesian(ylim = c(0,100), expand = F) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'top')

ggsave('figs/pareto_alt1.jpg',p1,width = 100, height = 105, units = 'mm',dpi = 600)
ggsave('figs/pareto_alt2.jpg',p2,width = 100, height = 105, units = 'mm',dpi = 600)
ggsave('figs/pareto_alt3.jpg',p3,width = 100, height = 105, units = 'mm',dpi = 600)


# # with 2 y-axes
# coeff <- 2000000
# 
# ggplot(tob, aes(x=ic)) +
#   
#   geom_point( aes(y=ci), size=2, color='grey60') + 
#   geom_point( aes(y=sed / coeff), size=2, color='red') +
#   
#   scale_y_continuous(
#     
#     # Features of the first axis
#     name = "CI",
#     
#     # Add a second axis and specify its features
#     sec.axis = sec_axis(~.*coeff, name="Sed")
#   ) + 
#   
#   theme_bw() +
#   
#   theme(
#     axis.title.y = element_text(color = 'grey60', size=13),
#     axis.title.y.right = element_text(color = 'red', size=13)
#   )

# redo figures only based on 2 objects
op <- readRDS('~/surfdrive/tmp/optifish_proc_20220516/caramel_ic_ci_gen80130_pop100.rds')

dec <- round(op$parameters,0) %>% as.data.frame() # decision
ob <- op$objectives %>% as.data.frame() # objective
# # sort based on IC <<< WHY? ask Rafa
dec <- dec[sort(ob$V1,index.return=T)$ix,]
ob <- ob[sort(ob$V1,index.return=T)$ix,]


# PARETO FRONTS ----------------------------------------------------------------

# get the pareto front tables
tob <- ob
colnames(tob) <- c('ic','ci')

ggplot(tob) +
  geom_point(aes(x = ic, y = ci))

tob$ci <- tob$ci/length(unique(read.csv('proc/sp_data_hybas12_mekong.csv')$binomial))


tob_l = tidyr::pivot_longer(tob, cols = c('ci'))
tob_l$name[tob_l$name == 'ci'] <- 'Connectivity Index [%]'

p <- ggplot(data = tob_l,
            aes(x = ic, y = value, color = name)) +
  geom_point(alpha = 1) +
  # geom_smooth(se = F) +
  scale_color_manual(values = RColorBrewer::brewer.pal(9,'Set1')[c(3,7)]) +
  labs(color = '') +
  xlab('Installed Capacity [MW]') +
  ylab('') +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  # scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  coord_cartesian(expand = T) +
  facet_grid(name~., scales = 'free_y', switch = 'y',
             # labeller = labeller('Connectivity Index' = ci, 'Sedimentation' = sed)
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'none',
        strip.background = element_blank(), strip.placement = 'outside',
        text = element_text(size = 13),
        strip.text = element_text(size = 13)
  )
p

ggsave('figs/pareto_ci_vs_IC.jpg',p,width = 110, height = 100, units = 'mm',dpi = 600)


sp_data_simple <- hb_data %>% as_tibble() %>% select(HYBAS_ID, MAIN_BAS, SUB_AREA, MAIN_BAS_AREA) %>% mutate(binomial = 'sp1', diad = 'f') %>% as.data.table()
tob$ci_simple <- apply(dec,1,function(x) ci(x,sp_data_simple))
sp_data_commercial <- sp_data %>% filter(Importance %in% c('commercial','highly commercial','minor commercial','subsistence fisheries'))
tob$ci_commercial <- apply(dec,1,function(x) ci(x,sp_data_commercial)/length(unique(sp_data_commercial$binomial)))
sp_data_threatened <- sp_data %>% filter(code %in% c('VU','EN','CR'))
tob$ci_threatened <- apply(dec,1,function(x) ci(x,sp_data_threatened)/length(unique(sp_data_threatened$binomial)))
sp_data_diad <- sp_data %>% filter(diad == 't')
sp_data_nondiad <- sp_data %>% filter(diad == 'f')
tob$ci_diadromous <- apply(dec,1,function(x) ci(x,sp_data_diad)/length(unique(sp_data_diad$binomial)))
tob$ci_nondiadromous <- apply(dec,1,function(x) ci(x,sp_data_nondiad)/length(unique(sp_data_nondiad$binomial)))

tob_l = tidyr::pivot_longer(tob, cols = starts_with('ci'), names_prefix = 'ci_')
tob_l$name[tob_l$name == 'ci'] <- 'all'

# all, commercial, threatened
p1 <- ggplot(data = tob_l %>% filter(name %in% c('all','commercial','threatened')),
             aes(x = ic, y = value, color = name)) +
  geom_point(alpha = 1) +
  # geom_smooth(se = F) +
  scale_color_manual(values = RColorBrewer::brewer.pal(9,'Set1')[c(9,5,1)]) +
  labs(color = '') +
  xlab('Installed Capacity [MW]') +
  ylab('Connectivity Index [%]') +
  # coord_cartesian(ylim = c(0,100), expand = T) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'top')

# all, simple
p2 <- ggplot(data = tob_l %>% filter(name %in% c('all','simple')),
             aes(x = ic, y = value, color = name)) +
  geom_point(alpha = 1) +
  # geom_smooth(se = F) +
  scale_color_manual(values = RColorBrewer::brewer.pal(9,'Set1')[c(9,2)]) +
  labs(color = '') +
  xlab('Installed Capacity [MW]') +
  ylab('Connectivity Index [%]') +
  # coord_cartesian(ylim = c(0,100), expand = T) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'top')

# all, diadromous, nondiadromous
p3 <- ggplot(data = tob_l %>% filter(name %in% c('all','diadromous','nondiadromous')),
             aes(x = ic, y = value, color = name)) +
  geom_point(alpha = 1) +
  # geom_smooth(se = F) +
  scale_color_manual(values = RColorBrewer::brewer.pal(9,'Set1')[c(9,3:4)]) +
  labs(color = '') +
  xlab('Installed Capacity [MW]') +
  ylab('Connectivity Index [%]') +
  # coord_cartesian(ylim = c(0,100), expand = T) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'top')

ggsave('figs/pareto_ic_vs_ci_alt1.jpg',p1,width = 100, height = 105, units = 'mm',dpi = 600)
ggsave('figs/pareto_ic_vs_ci_alt2.jpg',p2,width = 100, height = 105, units = 'mm',dpi = 600)
ggsave('figs/pareto_ic_vs_ci_alt3.jpg',p3,width = 100, height = 105, units = 'mm',dpi = 600)
