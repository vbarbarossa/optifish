# packages needed
library(sf); library(foreach); library(rfishbase); 
library(data.table); library(dplyr); library(exactextractr); 
library(igraph); library(ggplot2); library(rnaturalearth)
# raster package is also needed

# specify runs set
for(runs_set in c('','40sp','threatened')){
  
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
  
  if(runs_set != ''){
    source(paste0('R/prepare_data_LMB_',runs_set,'.R'))
  }else{
    source('R/prepare_data_LMB.R')
  }
  
  # simplify dams table for fitness function
  dams_data <- left_join(
    dams[,c('Code','INTER_ID','Status','DamHeight',name_col_IC, name_col_V)],
    # use manually adjusted operation year, called ReferenceY for backwards compatibility
    read.csv('data/dams_adjusted_year.csv') %>% select(Code,ReferenceY = 'COD_adj', Status_adj)
  )
  
  # check
  dams_data %>% group_by(Status) %>% summarise(minY = min(ReferenceY), maxY = max(ReferenceY))
  dams_data %>% group_by(Status_adj) %>% summarise(minY = min(ReferenceY), maxY = max(ReferenceY))
  
  # adjust status:
  dams_status_adjusted <- read.csv('data/dams_adjusted_year.csv')
  
  # simplify dams table for fitness function
  dams_data <- left_join(
    dams[,c('Code','INTER_ID','DamHeight',name_col_IC, name_col_V)],
    dams_status_adjusted %>% select(Code, Status = 'Status_adj',ReferenceY = 'COD_adj')
  )
  
  # # assign passability
  dams_data$pass <- 0
  
  dams_e <- dams_data[dams_data$Status != 'P',]
  dams_f <- dams_data[dams_data$Status == 'P',]
  
  # no species
  n_sp <- round(sum(L) + sum(ld),0)
  
  # no dams
  n = nrow(dams_f) + nrow(dams_e)
  
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
  
  
  dams_e <- dams_data
  years_range <- sort(dams_e$ReferenceY)
  op_pres <- foreach(i = years_range,.combine = 'rbind') %do% {
    -calc_objs(c(as.integer(dams_e$ReferenceY <= i),rep(0,nrow(dams_f))))
  } %>% as.data.frame()
  row.names(op_pres) <- 1:nrow(op_pres)
  colnames(op_pres) <- c('ic','ci')
  op_pres$scenario <- 'present'
  op_pres$year <- c(years_range)
  
  
  gen_size <- 1000
  
  # ,'pristine_pass_bin','pristine_pass_step','future','removal','mitigation_bin','mitigation_step'
  op <- foreach(i = c('pristine','future','removal_bin30','mitigation_bin30'),.combine = 'rbind') %do% {
    
    t1 <- readRDS(paste0('proc/LMB',runs_set,'_nsga2_',i,'_ic_ci_gen',gen_size,'_pop100.rds'))
    t1 <- -t1@fitness %>% as.data.frame()
    t1$scenario <- i
    return(t1)
  }
  
  colnames(op)[1:2] <- c('ic','ci')
  op$year <- ''
  op <- rbind(op_pres,op)
  
  op$scenario[op$scenario == 'mitigation_bin30'] <- 'mitigation_bin'
  op$scenario[op$scenario == 'removal_bin30'] <- 'removal'
  
  op$scenario[op$scenario == 'future'] <- 'strategic'
  op$scenario[which(as.numeric(op$year) >= 2021)] <- 'future'
  op$scenario <- factor(op$scenario, levels = c('pristine','present','future','strategic','removal','mitigation_bin'))
  op$ci <- op$ci*100
  
  op$runs_set <- runs_set
  
  if(runs_set != ''){
    write.csv(op,paste0('proc/pareto_tab_',runs_set,'.csv'),row.names = F)
  }else{
    write.csv(op,paste0('proc/pareto_tab_all.csv'),row.names = F)
  }
  
  rm(list=ls())
}



# create one table
op <- bind_rows(read.csv('proc/pareto_tab_all.csv'),
                read.csv('proc/pareto_tab_40sp.csv'),
                read.csv('proc/pareto_tab_threatened.csv'))
op$runs_set[is.na(op$runs_set)] <- 'all'
op_plot <- op
op_plot$scenario <- factor(op_plot$scenario, levels = c('pristine','present','future','strategic','mitigation_bin','removal'))
levels(op_plot$scenario) <- c('benchmark','existing','planned','strategic','passage','removal') #,'removal','mitigation'
# op_plot$vpanel <- 1
# op_plot$vpanel[op_plot$scenario %in% c()]


op_plot$runs_set <- factor(op_plot$runs_set, levels = c('all','threatened','40sp'))
levels(op_plot$runs_set) <- c('All species','IUCN threatened','MRC important')

op_plot$runs_set <- as.character(op_plot$runs_set)
op$ic <- op$ic/1000
fills <- c(c('grey40','black','white'),
           RColorBrewer::brewer.pal(3,'YlGnBu'))

op_plot$row <- 1
op_plot$row[op_plot$scenario %in% c('strategic','passage','removal')] <- 2

op_plot$ic <- op_plot$ic/1000

p_pareto <- ggplot() +
  # geom_step(data = op_plot %>% filter(scenario == c('existing')),direction = 'hv', show.legend = F, linetype = 2) +
  geom_point(data = op_plot %>% filter(runs_set == 'All species') %>% select(ic,ci,scenario),
             aes(x = ic, y = ci,fill = scenario),shape = 21,alpha = 0.4, size = 2, color = 'grey80') +
  geom_point(data = op_plot %>% filter(runs_set != 'All species'),
             aes(x = ic, y = ci,fill = scenario), shape = 24,alpha = 1, size = 2,  color = 'black') +
  scale_fill_manual(values = fills) +
  # scale_shape_manual(values = shapes) +
  xlab('Installed Capacity [GW]') +
  ylab('Connectivity Index [%]') +
  # scale_x_continuous(labels = function(x) format(x, big.interval=T ,scientific = TRUE)) +
  coord_cartesian(expand = T) +
  theme_bw() +
  guides(fill = 'none') +
  facet_wrap('runs_set',ncol = 2) +
  theme(panel.grid = element_blank(), 
        # panel.border = element_blank(),
        # axis.line = element_line(),
        legend.position = c(0.2,0.15), 
        legend.title = element_blank(),
        legend.direction = 'vertical',
        legend.background = element_blank(),
        strip.background = element_blank(), 
        strip.placement = 'outside',
        text = element_text(size = 14)
  )
p_pareto

ggsave('figs/paretos_species_focus.pdf',p_pareto,
       width = 180, height = 140, units = 'mm', scale = 1.2)

