# packages needed
library(sf); library(foreach); library(rfishbase); library(dplyr)
# raster package is also needed

sf_use_s2(FALSE)

# retrieve paths to input files
source('R/master_paths_local.R')

# functions for CI calculation
source('R/functions_connectivity.R')

# HYBAS ID of Mekong
main_bas_id <- 4120017020 #this corresponds to the outlet!

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

# determine centroids with sf
points <- st_centroid(hb_data)

# ------------------------------------------------------------------------------
# Reference species to hybas level 08 ------------------------------------------

sp <- read_sf('~/surfdrive/data/fish_databases/fishsuit/species_ranges_merged.gpkg')

# reference to hydrobasins level 12
lst <- st_contains(sp,points,sparse = T)
# lst is a sparse matrix where each entry is a row of sp and contains a list of hybas12 points falling within that species polygos

# make database where each entry is a hybas ID and 
# loop through the species
# for each species, create a table with hybasID and species id_no

# should update this part using dplyr as in extract_customRanges2hybas12.R
tab <- lapply(seq_along(lst),function(i){
  hb <- points$HYBAS_ID[lst[[i]]]
  if(length(hb) > 0){
    return(
      data.frame(HYBAS_ID = hb,
                 binomial = sp$binomial[i])
    )
  }
}
) %>% do.call('rbind',.) %>% distinct()

saveRDS(tab,'proc/species_ranges_merged_on_hybas12_mekong.rds')

# ------------------------------------------------------------------------------
# Fishbase metadata ------------------------------------------------------------

# retrieve fishbase metadata for diadromous-non diadromous species
if(!file.exists('proc/fishbase_data.csv')){
  # load fishbase metadata
  fishbase <- species(fields = c('Species','AnaCat')) %>% # get species table for all species
    rename(binomial = Species)
  fishbase$AnaCat <- as.factor(fishbase$AnaCat)
  levels(fishbase$AnaCat) <- c(rep('Diad.',23),'Non.','Ocea.','Ocea.','Pota.','Pota.')
  write.csv(fishbase,'proc/fishbase_data.csv',row.names = F)
  
}else{
  fishbase <- read.csv('proc/fishbase_data.csv') %>% as_tibble  
}

# ------------------------------------------------------------------------------
# Final species dataset --------------------------------------------------------

# Species range data
sp_data <- tab %>%
  inner_join(.,hb_data %>% as_tibble() %>% 
               select(HYBAS_ID,MAIN_BAS,SUB_AREA,MAIN_BAS_AREA),by="HYBAS_ID") %>%
  as.data.table(.)

# assign diadromous-non diadromous category
sp_data$diad <- 'f'
sp_data$diad[sp_data$binomial %in% fishbase$binomial[fishbase$AnaCat == 'Diad.']] <- 't'

write.csv(sp_data,'proc/sp_data_hybas12_mekong.csv',row.names = F)
