#### Make slope aspect and horizon maps for CCR

#### packages
library(tidyverse)
library(terra)
library(sf)
library(whitebox)


##set up whitebox
wbt_version()
wbt_init()

#### set file paths
map_dir <- "./CCR_NEW/Spatial_Data_comp/Create_Maps/"
dem <- paste0(map_dir, "Outputs/ccr_dem.tif")


#### Slope
wbt_slope(dem = dem, output = file.path(map_dir, "Outputs/ccr_slope.tif"), units = 'degrees')

plot(rast(paste0(map_dir, "Outputs/ccr_slope.tif")))


#### aspect is standard (0 == 360 == NORTH )
wbt_aspect(dem = dem, output = file.path(map_dir, "Outputs/ccr_aspect.tif"))

plot(rast(paste0(map_dir, "Outputs/ccr_aspect.tif")))

##color map
# asp_cols <- colorRampPalette(c("red", "yellow", "green", "cyan", "blue", "magenta", "red"))(72)   # 5-degree bins
# plot(rast(paste0(map_dir, "Outputs/ccr_aspect.tif")), main = "Aspect", col=asp_cols)



#### Horizons

##Get default east and west horizons
#east
wbt_horizon_angle(dem = dem, output = file.path(map_dir, "Outputs/temp_horizon_east.tif"), azimuth = 270)

plot(rast(paste0(map_dir, "Outputs/temp_horizon_east.tif")))

#west
wbt_horizon_angle(dem = dem, output = file.path(map_dir, "Outputs/temp_horizon_west.tif"), azimuth = 90)

plot(rast(paste0(map_dir, "Outputs/temp_horizon_west.tif")))


##match RHESSys format
horizon_east_sin = sin(rast(file.path(map_dir, "Outputs/temp_horizon_east.tif")) * pi/180)
horizon_west_sin = sin(rast(file.path(map_dir, "Outputs/temp_horizon_west.tif")) * pi/180)

plot(horizon_east_sin)
plot(horizon_west_sin)

writeRaster(horizon_east_sin, file.path(map_dir, "Outputs/ccr_e_horizon.tif"), overwrite = T)
writeRaster(horizon_west_sin, file.path(map_dir, "Outputs/ccr_w_horizon.tif"), overwrite = T)



