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

# p1<-ggplot(data=HYBAS,aes(fill=QS_1))+geom_sf(color = NA)+ggtitle('Rasterized Method')
# p2<-ggplot(data=HYBAS,aes(fill=QS_2))+geom_sf(color = NA)+ggtitle("Vectorized Method")
# p3<-ggplot(data=HYBAS,aes(fill=QS_1-QS_2))+geom_sf(color = NA)+ggtitle("Difference")
# 
# library(patchwork)
# p1+p2+p3

# 3 map connectivity matrix

# create correspondence ID


# # create edges matrix based on HBasins
# edges <- as.matrix(HYBAS[,c('HYBAS_ID','NEXT_DOWN')] %>% as.data.frame() %>% dplyr::select(-geometry))
# # edges <- matrix(c(1:5,3,0,4,))
# 
# # create empty adjecency matrix
# codes <- unique(c(edges))
# A <- matrix(0, nrow=length(codes), ncol=length(codes), dimnames=list(codes, codes))
# 
# # fill it
# A[edges] <- 1


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

sed_trans <- function(I){
  
  # reoder the QS based on matrix rows
  qs <- data.frame(QS = HYBAS$QS_1)
  row.names(qs) <- as.character(HYBAS$HYBAS_ID)
  # remove NAs
  qs[is.na(qs)] <- 0
  
  MQS = I
  # [1:n,n:most_downstream] >> [cols of I, rows of I]
  for(j in 1:length(dams$HYBAS_ID)){
    id <- dams$HYBAS_ID[j]
    up_ <- row.names(I)[I[,as.character(id)] == 1]
    dw_ <- colnames(I)[I[as.character(id),] == 1]
    MQS[up_,dw_] <- MQS[up_,dw_] * dams$RE[dams$HYBAS_ID == id]  
  }
  
  return(
    data.frame(HYBAS_ID = as.numeric(colnames(MQS)),
               QS_sum = as.numeric(qs[row.names(M),'QS'] %*% M))
  )
  
}

st <- Sys.time()
a = sed_trans(I=I)
Sys.time() - st



