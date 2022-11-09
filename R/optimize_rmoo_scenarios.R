
# ##############################################################################
# SETTINGS BLOCK ###############################################################
# github: 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

# g <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
LOCAL = TRUE

# rmoo optimization setups --------------------------------------------------
pop_size <- 100
gen_size <- 1000 #already a good 
MONITOR <- FALSE
# ------------------------------------------------------------------------------

# fitness function setup -------------------------------------------------------
sedimentation = F
fragmentation = T
energy = T
water = F
all_dams = F
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
library(igraph);

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

# FITNESS FUNCTION -------------------------------------------------------------

# fitness function
calc_objs <- 
  function(x){ 
    
    # select dams based on index
    decision <- x
    
    # extract set of dams included
    # this way the order of the decision vector is the same for any scenario
    dams_tot <- dams_data
    dams_tot[fut_index,'incl'] <- dams_tot[fut_index,'incl']*decision
    dams_tot <- dams_tot[dams_tot$incl == 1,]
    
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
# ------------------------------------------------------------------------------

# DAMS DATA FILTERING ----------------------------------------------------------

# simplify dams table for fitness function
dams_data <- dams[,c('Code','INTER_ID','Status','DamHeight',name_col_IC, name_col_V)]

assign_passability <- function(height){
  return(
    sapply(height, function(x){
      res <- 0
      if(x < 20) res <- 0.4
      if(x < 10) res <- 0.6
      if(x < 1) res <- 0.8
      return(res) 
    }
    )
  )
}
# dams_data$pass[dams_data$DamHeight < 20] <- 0.4
# dams_data$pass[dams_data$DamHeight < 10] <- 0.6
# dams_data$pass[dams_data$DamHeight < 1] <- 0.8
# 
# dams_data$pass <- scales::rescale(dams_data$DamHeight, to = c(1,0), from = c(0,40))
# dams_data$pass[dams_data$pass < 0] <- 0

# defauls inclusion index, needed later in fitness function
dams_data$incl <- 1

# convert to data frame
dams_data <- as.data.frame(dams_data)

# PRISTINE w/ PASSABILITY SCENARIO ---------------------------------------------
scen <- 'pristine_pass'
cat('\n',scen,' RUN --------------------------------------------------------\n')

# assign passability
dams_data$pass <- assign_passability(dams_data$DamHeight)

# dams_data$status <- 'current'
# dams_data$status[dams_data$Status == 'P'] <- 'future'
dams_data$status <- 'future' # set all dams to future

# no species
n_sp <- round(sum(L) + sum(ld),0)
# no dams
n = nrow(dams_data[dams_data$status == 'future',])
# index of future dams
fut_index <- dams_data$status == 'future'

st <- Sys.time()
cat('Starting optimization on', as.character(st), '\n\n')

# try the function
# calc_objs(rep(1,n))
# calc_objs(rep(0,n))

# define number of objectives based on settings
nobjs <- sum(sedimentation,fragmentation,water,energy)

# save the right variables in the output name
save_str <- c('ic','vol','sed','ci')[c(energy,water,sedimentation,fragmentation)]

# run the optimization
st <- Sys.time()
op <- rmoo::nsga2(
  type = 'binary', nBits = n, fitness = calc_objs,
  popSize = pop_size,
  nObj = nobjs,
  pcrossover = 0.8,
  pmutation = 0.2,
  maxiter = gen_size,
  suggestions = rbind(rep(1,n),rep(0,n)),
  summary = FALSE,
  monitor = MONITOR,
  names = save_str
)
Sys.time() - st


# plot(-op@fitness[,1],-op@fitness[,2])

cat('\nSaving..')
saveRDS(op,paste0('proc/nsga2_',scen,'_',paste(save_str,collapse = '_'),'_gen',gen_size,'_pop',pop_size,'.rds'))

et <- Sys.time() - st
cat('\nCompleted in',round(et,2), attr(et,'units'),'\n')

# PRISTINE SCENARIO ------------------------------------------------------------
scen <- 'pristine'
cat('\n',scen,' RUN --------------------------------------------------------\n')
# assign passability
dams_data$pass <- 0

# dams_data$status <- 'current'
# dams_data$status[dams_data$Status == 'P'] <- 'future'
dams_data$status <- 'future' # set all dams to future

# no species
n_sp <- round(sum(L) + sum(ld),0)
# no dams
n = nrow(dams_data[dams_data$status == 'future',])
# index of future dams
fut_index <- dams_data$status == 'future'

st <- Sys.time()
cat('Starting optimization on', as.character(st), '\n\n')

# try the function
# calc_objs(rep(1,n))
# calc_objs(rep(0,n))

# define number of objectives based on settings
nobjs <- sum(sedimentation,fragmentation,water,energy)

# save the right variables in the output name
save_str <- c('ic','vol','sed','ci')[c(energy,water,sedimentation,fragmentation)]

# run the optimization
st <- Sys.time()
op <- rmoo::nsga2(
  type = 'binary', nBits = n, fitness = calc_objs,
  popSize = pop_size,
  nObj = nobjs,
  pcrossover = 0.8,
  pmutation = 0.2,
  maxiter = gen_size,
  suggestions = rbind(rep(1,n),rep(0,n)),
  summary = FALSE,
  monitor = MONITOR,
  names = save_str
)
Sys.time() - st


# plot(-op@fitness[,1],-op@fitness[,2])

cat('\nSaving..')
saveRDS(op,paste0('proc/nsga2_',scen,'_',paste(save_str,collapse = '_'),'_gen',gen_size,'_pop',pop_size,'.rds'))

et <- Sys.time() - st
cat('\nCompleted in',round(et,2), attr(et,'units'),'\n')


# FUTURE SCENARIO ------------------------------------------------------------
scen <- 'future'
cat('\n',scen,' RUN --------------------------------------------------------\n')
# assign passability
dams_data$pass <- 0

dams_data$status <- 'current'
dams_data$status[dams_data$Status == 'P'] <- 'future'

# no species
n_sp <- round(sum(L) + sum(ld),0)
# no dams
n = nrow(dams_data[dams_data$status == 'future',])
# index of future dams
fut_index <- dams_data$status == 'future'

st <- Sys.time()
cat('Starting optimization on', as.character(st), '\n\n')

# try the function
# calc_objs(rep(1,n))
# calc_objs(rep(0,n))

# define number of objectives based on settings
nobjs <- sum(sedimentation,fragmentation,water,energy)

# save the right variables in the output name
save_str <- c('ic','vol','sed','ci')[c(energy,water,sedimentation,fragmentation)]

# run the optimization
st <- Sys.time()
op <- rmoo::nsga2(
  type = 'binary', nBits = n, fitness = calc_objs,
  popSize = pop_size,
  nObj = nobjs,
  pcrossover = 0.8,
  pmutation = 0.2,
  maxiter = gen_size,
  suggestions = rbind(rep(1,n),rep(0,n)),
  summary = FALSE,
  monitor = MONITOR,
  names = save_str
)
Sys.time() - st


# plot(-op@fitness[,1],-op@fitness[,2])

cat('\nSaving..')
saveRDS(op,paste0('proc/nsga2_',scen,'_',paste(save_str,collapse = '_'),'_gen',gen_size,'_pop',pop_size,'.rds'))

et <- Sys.time() - st
cat('\nCompleted in',round(et,2), attr(et,'units'),'\n')


# REMOVAL SCENARIO -------------------------------------------------------------
scen <- 'removal'
cat('\n',scen,' RUN --------------------------------------------------------\n')
# assign passability
dams_data$pass <- 0

# remove dams that are not in pristine runs at equal IC
ICmax <- sum(pull(dams_data[dams_data$Status != 'P',],name_col_IC))

# pareto pristine
pp <- readRDS(paste0('proc/nsga2_pristine_',paste(save_str,collapse = '_'),'_gen',gen_size,'_pop',pop_size,'.rds'))
# find row with most similar IC to ICmax
ICmax_row <- which(abs(-pp@fitness[,1] - ICmax) == min(abs(-pp@fitness[,1] - ICmax)))[1]

# set dams
dams_data$status <- 'current'
dams_data$status[dams_data$Status == 'P'] <- 'future'
dams_data$status[pp@population[ICmax_row,] == 0] <- 'future'

# no species
n_sp <- round(sum(L) + sum(ld),0)
# no dams
n = nrow(dams_data[dams_data$status == 'future',])
# index of future dams
fut_index <- dams_data$status == 'future'

st <- Sys.time()
cat('Starting optimization on', as.character(st), '\n\n')

# try the function
# calc_objs(rep(1,n))
# calc_objs(rep(0,n))

# define number of objectives based on settings
nobjs <- sum(sedimentation,fragmentation,water,energy)

# save the right variables in the output name
save_str <- c('ic','vol','sed','ci')[c(energy,water,sedimentation,fragmentation)]

# run the optimization
st <- Sys.time()
op <- rmoo::nsga2(
  type = 'binary', nBits = n, fitness = calc_objs,
  popSize = pop_size,
  nObj = nobjs,
  pcrossover = 0.8,
  pmutation = 0.2,
  maxiter = gen_size,
  suggestions = rbind(rep(1,n),rep(0,n)),
  summary = FALSE,
  monitor = MONITOR,
  names = save_str
)
Sys.time() - st


# plot(-op@fitness[,1],-op@fitness[,2])

cat('\nSaving..')
saveRDS(op,paste0('proc/nsga2_',scen,'_',paste(save_str,collapse = '_'),'_gen',gen_size,'_pop',pop_size,'.rds'))

et <- Sys.time() - st
cat('\nCompleted in',round(et,2), attr(et,'units'),'\n')

# MITIGATION SCENARIO -------------------------------------------------------------
scen <- 'mitigation'
cat('\n',scen,' RUN --------------------------------------------------------\n')
# assign passability
dams_data$pass <- 0

# set passability for dams that are not in pristine runs at equal IC
ICmax <- sum(pull(dams_data[dams_data$Status != 'P',],name_col_IC))

# pareto pristine
pp <- readRDS(paste0('proc/nsga2_pristine_',paste(save_str,collapse = '_'),'_gen',gen_size,'_pop',pop_size,'.rds'))
# find row with most similar IC to ICmax
ICmax_row <- which(abs(-pp@fitness[,1] - ICmax) == min(abs(-pp@fitness[,1] - ICmax)))[1]

# set dams
dams_data$status <- 'current'
dams_data$status[dams_data$Status == 'P'] <- 'future'

# set passability for dams to mitigate
dams_data$pass[(pp@population[ICmax_row,] == 0) & (dams_data$status == 'current')] <- assign_passability(height = dams_data$DamHeight[(pp@population[ICmax_row,] == 0) & (dams_data$status == 'current')])

# set passability also for future dams??? <<<<<<<<<<<<<
# dams_data$pass[dams_data$status == 'future'] <- assign_passability(dams_data$DamHeight[dams_data$status == 'future'])


# no species
n_sp <- round(sum(L) + sum(ld),0)
# no dams
n = nrow(dams_data[dams_data$status == 'future',])
# index of future dams
fut_index <- dams_data$status == 'future'

st <- Sys.time()
cat('Starting optimization on', as.character(st), '\n\n')

# try the function
# calc_objs(rep(1,n))
# calc_objs(rep(0,n))

# define number of objectives based on settings
nobjs <- sum(sedimentation,fragmentation,water,energy)

# save the right variables in the output name
save_str <- c('ic','vol','sed','ci')[c(energy,water,sedimentation,fragmentation)]

# run the optimization
st <- Sys.time()
op <- rmoo::nsga2(
  type = 'binary', nBits = n, fitness = calc_objs,
  popSize = pop_size,
  nObj = nobjs,
  pcrossover = 0.8,
  pmutation = 0.2,
  maxiter = gen_size,
  suggestions = rbind(rep(1,n),rep(0,n)),
  summary = FALSE,
  monitor = MONITOR,
  names = save_str
)
Sys.time() - st


cat('\nSaving..')
saveRDS(op,paste0('proc/nsga2_',scen,'_',paste(save_str,collapse = '_'),'_gen',gen_size,'_pop',pop_size,'.rds'))

et <- Sys.time() - st
cat('\nCompleted in',round(et,2), attr(et,'units'),'\n')


# MITIGATION2 SCENARIO ---------------------------------------------------------
scen <- 'mitigation2'
cat('\n',scen,' RUN --------------------------------------------------------\n')
# assign passability
dams_data$pass <- 0

# set passability for dams that are not in pristine runs at equal IC
ICmax <- sum(pull(dams_data[dams_data$Status != 'P',],name_col_IC))

# pareto pristine
pp <- readRDS(paste0('proc/nsga2_pristine_',paste(save_str,collapse = '_'),'_gen',gen_size,'_pop',pop_size,'.rds'))
# find row with most similar IC to ICmax
ICmax_row <- which(abs(-pp@fitness[,1] - ICmax) == min(abs(-pp@fitness[,1] - ICmax)))[1]

# set dams
dams_data$status <- 'current'
dams_data$status[dams_data$Status == 'P'] <- 'future'

# set passability for dams to mitigate
dams_data$pass[(pp@population[ICmax_row,] == 0) & (dams_data$status == 'current')] <- assign_passability(height = dams_data$DamHeight[(pp@population[ICmax_row,] == 0) & (dams_data$status == 'current')])

# set passability also for future dams??? <<<<<<<<<<<<<
dams_data$pass[dams_data$status == 'future'] <- assign_passability(dams_data$DamHeight[dams_data$status == 'future'])


# no species
n_sp <- round(sum(L) + sum(ld),0)
# no dams
n = nrow(dams_data[dams_data$status == 'future',])
# index of future dams
fut_index <- dams_data$status == 'future'

st <- Sys.time()
cat('Starting optimization on', as.character(st), '\n\n')

# try the function
# calc_objs(rep(1,n))
# calc_objs(rep(0,n))

# define number of objectives based on settings
nobjs <- sum(sedimentation,fragmentation,water,energy)

# save the right variables in the output name
save_str <- c('ic','vol','sed','ci')[c(energy,water,sedimentation,fragmentation)]

# run the optimization
st <- Sys.time()
op <- rmoo::nsga2(
  type = 'binary', nBits = n, fitness = calc_objs,
  popSize = pop_size,
  nObj = nobjs,
  pcrossover = 0.8,
  pmutation = 0.2,
  maxiter = gen_size,
  suggestions = rbind(rep(1,n),rep(0,n)),
  summary = FALSE,
  monitor = MONITOR,
  names = save_str
)
Sys.time() - st


cat('\nSaving..')
saveRDS(op,paste0('proc/nsga2_',scen,'_',paste(save_str,collapse = '_'),'_gen',gen_size,'_pop',pop_size,'.rds'))

et <- Sys.time() - st
cat('\nCompleted in',round(et,2), attr(et,'units'),'\n')




cat('# -E-vissero-tutti-felici-e-contenti----------------------------------#\n')
