library(sf)
library(sp)
library(raster)
library(ggplot2)
library(exactextractr)


# load 'geomorphic provinces' large areas with similar sediment YIELD [t/km2/yr]
GP<-read_sf('./R/data/Mekong geomorphic provinces/Mekong_Geomorphic_provinces.gpkg')

# hydrobasins

HYBAS <- 
  read_sf('./R/data/HYBAS_LEVEL_12_MEKONG/Hybas Level 12 Mekong prj.shp')

# two ways of doing this 

# 1: create a raster and extract sediment yields in each HB using exat extract. Advantage: Handles situations where a HYBAS intersects multiple geomorphic provinces 
GPr<-raster(crs = crs(GP), vals = 0, resolution = c(1000, 1000), ext = extent(GP)) %>%
  rasterize(as_Spatial(GP), ., field = "SedYield")

HYBAS$YS_1<-exact_extract(GPr,HYBAS,fun = "mean") 
HYBAS$QS_1<-HYBAS$YS_1*HYBAS$SUB_AREA


# 2 using spatial overlay - somewhat easier and faster, but can't handle if a HYBAS intersects multiple geomorphic provinces
HYBAS$YS_2<-over(as_Spatial(HYBAS),as_Spatial(GP))["SedYield"] %>% data.matrix(.)
HYBAS$QS_2<-HYBAS$YS_2*HYBAS$SUB_AREA

p1<-ggplot(data=HYBAS,aes(fill=QS_1))+geom_sf(color = NA)+ggtitle('Rasterized Method')
p2<-ggplot(data=HYBAS,aes(fill=QS_2))+geom_sf(color = NA)+ggtitle("Vectorized Method")
p3<-ggplot(data=HYBAS,aes(fill=QS_1-QS_2))+geom_sf(color = NA)+ggtitle("Difference")

library('patchwork')
test<-1
p1+p2+p3

sum(HYBAS$QS_1,na.rm=T)
sum(HYBAS$QS_2,na.rm=T)

# 3 map connectivity matrix

# create correspondence ID


# create edges matrix based on HBasins
edges <- as.matrix(HYBAS[,c('HYBAS_ID','NEXT_DOWN')] %>% as.data.frame() %>% dplyr::select(-geometry))
edges <- matrix(c(1:5,3,0,4,))

# create empty adjecency matrix
codes <- unique(c(edges))
A <- matrix(0, nrow=length(codes), ncol=length(codes), dimnames=list(codes, codes))

# fill it
A[edges] <- 1


library(igraph)
df <- data.frame(from = HYBAS$HYBAS_ID,to = HYBAS$NEXT_DOWN)

outlet <- as.character(df$from[which(df$to == 0)])

df <- df[-which(df$to == 0),]
g <- graph_from_data_frame(d=df)
I <- get.adjacency(g)

# 
# m <- as.matrix(I)
for(i in row.names(I)) {
  I[i,as.numeric(shortest_paths(graph=g,from = i,to=as.character(outlet))$vpath[[1]])] <- 1
}

sed_trans <- function(m=I){
  
  # reoder the QS based on matrix rows
  qs <- data.frame(QS = HYBAS$QS_1)
  row.names(qs) <- as.character(HYBAS$HYBAS_ID)
  
  # order based on row names
  # multiply the matrix
  # sum
  
}
