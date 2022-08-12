
# ##############################################################################
# SETTINGS BLOCK ###############################################################
# github: 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

g <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
LOCAL = FALSE

# caramel optimization setups --------------------------------------------------
pop_run <- 100
tot_run <- c(5000000,5000000)[g]
init_pop_run <- 100
arch_run <- 100
# ------------------------------------------------------------------------------

# fitness function setup -------------------------------------------------------
sedimentation = F
fragmentation = T
energy = T
water = F
all_dams = c(T,F)[g]
# ------------------------------------------------------------------------------

# dams table details -----------------------------------------------------------
name_col_IC <- 'InstalledC'
name_col_V <- 'GrossStora'
name_col_Q <- 'MeanQ_m3s' # flow of the river -> to calculate trap efficiency 
# ------------------------------------------------------------------------------

# coeffs for sedimentation calculations ----------------------------------------
TE_bed <- 1 # trapping efficiency for coarse bed material
f_bed <- 0.1 # fraction of coarse bed material in the total sediments
# ------------------------------------------------------------------------------

# MAIN BASIN ID ----------------------------------------------------------------
main_bas_id <- outlet <- 4120017020 #this corresponds to the outlet!
# ------------------------------------------------------------------------------

# packages & ancillary functions -----------------------------------------------
# install.packages('sf')
# install.packages('foreach')
# install.packages('data.table')
# install.packages('dplyr')
# install.packages('sp')
# remotes::install_github("ATFutures/dodgr")
# remotes::install_github("ropensci/rfishbase")
# remotes::install_github("isciences/exactextractr")
# remotes::install_github('fzao/caRamel',dependencies=T)
library(sf); library(foreach); library(rfishbase); 
library(data.table); library(dplyr); library(exactextractr); 
library(igraph); library(caRamel)

# settings for sf
sf_use_s2(FALSE)

# read paths to input files
source('R/master_paths.R')
if(LOCAL) source('R/master_paths_local.R')

# functions for CI calculation
source('R/functions_connectivity.R')

# prepare the data for the optimization
source('R/prepare_data.R')
# ------------------------------------------------------------------------------

# FITNESS FUNCTION --------------------------------------------------------------------------

# simplify dams table for fitness function
dams_data <- dams[,c('INTER_ID','Status','DamHeight',name_col_IC, name_col_V)]

# # assign passability
dams_data$pass <- 0
# dams_data$pass[dams_data$DamHeight < 20] <- 0.4
# dams_data$pass[dams_data$DamHeight < 10] <- 0.6
# dams_data$pass[dams_data$DamHeight < 1] <- 0.8
# 
# dams_data$pass <- scales::rescale(dams_data$DamHeight, to = c(1,0), from = c(0,40))
# dams_data$pass[dams_data$pass < 0] <- 0

dams_e <- dams_data[dams_data$Status != 'P',]
dams_f <- dams_data[dams_data$Status == 'P',]

# no species
n_sp <- round(sum(L) + sum(ld),0)

# no dams
n = nrow(dams_f)
if(all_dams) n = n + nrow(dams_e)


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
    return( result_fitness )
    
  }


# Run the algorithm -------------------------------------------------------------------------

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
    varlist=c("dams_e","dams_f",
              "RE_M_array","qs","out",
              "master_upstream_list","inter_basins",
              "sp_data_inter",
              "name_col_IC", "name_col_V",
              "sedimentation","fragmentation","water","energy","all_dams",
              "g_master_pass","bas_tot","n_sp","L","ld")
  )
} 

cat('\n# -------------------------------------------------------------------\n')
st <- Sys.time()
cat('Starting optimization on', as.character(st), '\n\n')

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
  prec = matrix(c(58,0.001), nrow = 1, ncol = nobjs),
  carallel = TRUE,
  graph = FALSE,
  sensitivity = FALSE,
  verbose = T
)

# plot_caramel(op)
plot_pareto(op$objectives,maximized = rep(T,nobjs),objnames = c('InCap','Vol','Sed','CI')[c(energy,water,sedimentation,fragmentation)])

# save the right variables in the output name
save_str <- c('ic','vol','sed','ci')[c(energy,water,sedimentation,fragmentation)]
fut_str <- 'allDams'
if(!all_dams) fut_str <- 'futureDams'

cat('\nSuccessful: ',op$success)
cat('\nSaving..')
saveRDS(op,paste0('proc/caramel_',fut_str,'_',paste(save_str,collapse = '_'),'_gen',nrow(op$save_crit),'_pop',pop_run,'.rds'))

et <- Sys.time() - st
cat('\nCompleted in',round(et,2), attr(et,'units'),'\n')
cat('Ended on', as.character(Sys.time()),'\n')
cat('# -E-vissero-tutti-felici-e-contenti----------------------------------#\n')
