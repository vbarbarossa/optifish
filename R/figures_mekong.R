# packages needed
library(sf); library(foreach); library(rfishbase); library(data.table); 
library(dplyr); library(mco); library(sp); library(exactextractr); library(igraph);
library(caRamel); library(ggplot2); library(rnaturalearth)
# raster package is also needed

# settings for sf
sf_use_s2(FALSE)

# retrieve paths to input files
source('R/master_paths_local.R')

# functions for CI calculation
source('R/functions_connectivity.R')

# OPTIMIZATION RESULTS ---------------------------------------------------------
# retrieved from server
op <- readRDS('~/surfdrive/tmp/optifish_proc_20220420/caramel_ic_vol_sed_ci_gen48079_pop100.rds')

dec <- round(op$parameters,0) %>% as.data.frame() # decision
ob <- op$objectives %>% as.data.frame() # objective
# # sort based on IC <<< WHY? ask Rafa
dec <- dec[sort(ob$V1,index.return=T)$ix,]
ob <- ob[sort(ob$V1,index.return=T)$ix,]

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

library(ggplot2)
p_incl <- ggplot() +
  geom_sf(data = bas, color = 'grey20',fill='grey90') + #no border
  geom_sf(data = riv1,size = 0.1,color = 'grey40') +
  geom_sf(data = riv2,size = 0.3,color = 'grey40') +
  geom_sf(data = riv3,size = 0.6,color = 'grey40') +
  geom_sf(data = riv4,size = 0.9,color = 'grey40') +
  geom_sf(data = riv5,size = 1.2,color = 'grey40') +
  geom_sf(data = dams, aes(size = InstalledC,fill = incl), shape = 21,alpha = 0.8) +
  scale_fill_viridis_c(option = 'C',
                       # breaks = c(0, 25,50,75,100),labels = c(0,25,50,75,100),
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

ggsave('figs/inclusion_map_mekong.jpg',p_incl,width = 150, height = 200, units = 'mm',dpi = 600)


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

# PARETO FRONT MAPS ------------------------------------------------------------

# get the pareto front tables
tob <- ob
colnames(tob) <- c('ic','vol','sed','ci')

ggplot(tob) +
  geom_point(aes(x = ic, y = ci))

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

ggplot(tob) +
  geom_point(aes(x = ic, y = ci_m)) +
  geom_point(aes(x = ic, y = ci_s), color = 'red')



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
  geom_point(alpha = 0.2) +
  geom_smooth(se = F) +
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
  geom_point(alpha = 0.2) +
  geom_smooth(se = F) +
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
  geom_point(alpha = 0.2) +
  geom_smooth(se = F) +
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


tob_l = tidyr::pivot_longer(tob, cols = c('ci','sed'))
tob_l$name[tob_l$name == 'ci'] <- 'Connectivity Index [%]'
tob_l$name[tob_l$name == 'sed'] <- 'Sediment supply to delta [t/year]'

p <- ggplot(data = tob_l,
            aes(x = ic, y = value, color = name)) +
  geom_point(alpha = 0.2) +
  geom_smooth(se = F) +
  scale_color_manual(values = RColorBrewer::brewer.pal(9,'Set1')[c(3,7)]) +
  labs(color = '') +
  xlab('Installed Capacity [MW]') +
  ylab('') +
  coord_cartesian(expand = F) +
  facet_wrap('name', scales = 'free_y', ncol = 1,
             labeller = labeller('Connectivity Index' = ci, 'Sedimentation' = sed)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 13),
        strip.text = element_text(size = 13)
        )
p

ggsave('figs/pareto_ci_sed_vs_IC.jpg',p,width = 100, height = 180, units = 'mm',dpi = 600)

