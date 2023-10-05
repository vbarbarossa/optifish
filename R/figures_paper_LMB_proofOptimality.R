# packages needed
library(sf); library(foreach); library(rfishbase); 
library(data.table); library(dplyr); library(exactextractr); 
library(igraph); library(ggplot2); library(rnaturalearth)
# raster package is also needed

# specify runs set
tab <- foreach(gens = c(50,100,500,1000,5000,10000),.combine = 'rbind') %do% {
  
  # ,'pristine_pass_bin','pristine_pass_step','future','removal','mitigation_bin','mitigation_step'
  op <- readRDS(paste0('proc/LMB_PROOFOPT_nsga2_pristine_ic_ci_gen',gens,'_pop100.rds'))
  op <- -op@fitness %>% as.data.frame()
  
  colnames(op)[1:2] <- c('ic','ci')
  op$gen_size <- gens
  return(op)
}

tab$gen_size <- as.factor(tab$gen_size)

tab$ic <- tab$ic/1000
tab$ci <- tab$ci*100

p_pareto <- ggplot() +
  # geom_step(data = op_plot %>% filter(scenario == c('existing')),direction = 'hv', show.legend = F, linetype = 2) +
  geom_step(data = tab,
             aes(x = ic, y = ci, color = gen_size), size = 1) +
  
  # scale_fill_manual(values = fills) +
  # scale_shape_manual(values = shapes) +
  xlab('Installed Capacity [GW]') +
  ylab('Connectivity Index [%]') +
  # scale_x_continuous(labels = function(x) format(x, big.interval=T ,scientific = TRUE)) +
  coord_cartesian(expand = T) +
  labs(color = 'No. iterations') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        # panel.border = element_blank(),
        # axis.line = element_line(),
        # legend.position = c(0.2,0.15), 
        # legend.title = element_blank(),
        legend.direction = 'vertical',
        legend.background = element_blank(),
        strip.background = element_blank(), 
        strip.placement = 'outside',
        text = element_text(size = 14)
  )
p_pareto

ggsave('figs/paretos_prrofOptimality.pdf',p_pareto,
       width = 120, height = 100, units = 'mm', scale = 1.3)

