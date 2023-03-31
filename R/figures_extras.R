library(dplyr); library(sf)

# SR and interbasin ID on HBs for plotting
sr_data <- read.csv('proc/sp_data_hybas12_mekong.csv') %>%
  group_by(HYBAS_ID) %>%
  summarise(SR_tot = length(unique(binomial)),
            SR_diad = length(unique(binomial[diad != 'f'])),
            SR_nondiad = length(unique(binomial[diad == 'f'])))

bas <- read_sf('proc/basins_mekong.gpkg')

# get interbasin ID
source(textConnection(readLines('R/optimize_rmoo_scenarios.R')[1:63]))

bas_sr_inter <- left_join(bas,inter_basin_corr %>% select(HYBAS_ID, INTER_ID, INTER_NEXT)) %>%
  left_join(sr_data)

write_sf(bas_sr_inter,'proc/basins_mekong_spRichness_interBasins.gpkg')
