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
  
# 1: create a raster and extract sediment yields in each HB using exat extract. Advantage: Handles situations where a hydrobasins cuts accross multiple geomoeprhic provinces. 
  GPr<-raster(crs = crs(GP), vals = 0, resolution = c(5000, 5000), ext = extent(GP)) %>%
    rasterize(as_Spatial(GP), ., field = "SedYield")
  
  HYBAS$YS_1<-exact_extract(GPr,HYBAS,fun = "mean") 
  HYBAS$QS_1<-HYBAS$YS_1*HYBAS$SUB_AREA

  
# 2 using spatial overlay - somewhat easier and faster, but can't handle if a HYBAS intersects multiple geomorphic provinces
  HYBAS$YS_2<-over(as_Spatial(HB),as_Spatial(GP))["SedYield"] %>% data.matrix(.)
  HYBAS$QS_2<-HYBAS$YS_2*HYBAS$SUB_AREA
  
p1<-ggplot(data=HYBAS,aes(fill=QS_1))+geom_sf(color = NA)+ggtitle('Rasterized Method')
p2<-ggplot(data=HYBAS,aes(fill=QS_2))+geom_sf(color = NA)+ggtitle("Vectorized Method")
p3<-ggplot(data=HYBAS,aes(fill=QS_1-QS_2))+geom_sf(color = NA)+ggtitle("Difference")

library('patchwork')

p1+p2+p3