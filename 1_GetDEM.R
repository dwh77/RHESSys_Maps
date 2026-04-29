#### GET DEM
## maps need to make maps for CCR basins

#### packages
library(tidyverse)
library(terra)
# library(sf)
library(FedData) #getting DEM
# library(nhdplusTools) #getting flowlines from NHD
# library(mapview) #make interactive map
# library(mapedit) #add points on map to get pour point locations


#### Set directory ----
getwd()
map_dir <- "./Downloaded_Data/"



### Get DEM ----
# READ: ONCE YOU'VE DONE THIS SECTION ONCE YOU HAVE THE DATA LOCALLY AND CAN READ IN FROM THERE

#### set basin boundary Make box for watershed area
## For HPB square
# basin_bounds = c( 37.3800, 37.3540, -79.9715, -80.0071) # N S E W
# ## for CCR square
# basin_bounds = c( 37.4600, 37.3300, -79.8900, -80.1000) # N S E W
#
# basin_ext = rast(crs = "EPSG:4326")
# ext(basin_ext) = basin_bounds[c(4,3,2,1)]
# values(basin_ext) = 1
#
# plot(basin_ext)

# #writeRaster(x = basin_ext, filename = "rhutils_trials/spatial_source/basin_ext.tif",overwrite =T)

#this both reads in a DEM locally and wrties a tif file to the extraction.dir
# NED <- FedData::get_ned(template = basin_ext,
#                label="CCR",
#                extraction.dir = paste0(map_dir, "FED_DEM/"),
#                force.redo = T)
#
# plot(NED)

# #once acquired can read in locally
# NED <- rast(paste0(map_dir, "./FED_DEM/CCR_NED_1.tif"))
# plot(NED)



#### reproject DEM to UTM17N  ----
NED <- rast(paste0(map_dir, "./FED_DEM/CCR_NED_1.tif"))
plot(NED)

# 1. Project DEM to UTM17N (using watershed as CRS/grid reference)
dem_proj <- project(NED, "EPSG:26917")
plot(dem_proj)

# writeRaster(x = dem_proj, filename = "Downloaded_data/FED_DEM/CCR_NED_1_nad83utm17N.tif",overwrite =T)


#### Rastersize basins to DEM from shpfile ----
####Read in basins shapefiles and rasterize (these shapefiles are original download from streamstats)
##Read in shp file
ccr_shp <- vect(paste0(map_dir, "CCR_watershed_shpfile_StreamStats/ccr_ws/layers/globalwatershed.shp"))
plot(ccr_shp)

#reproject shp file
ccr_shp_proj <- project(ccr_shp, crs(dem_proj))
plot(ccr_shp_proj)

# Verify they overlap before rasterizing
plot(dem_proj)
plot(ccr_shp_proj, add = TRUE, col = NA, border = "red", lwd = 2)

# Rasterize and crop to area of watershed
ccr_shp_rast <- rasterize(ccr_shp_proj, dem_proj, field = 1)
plot(ccr_shp_rast)

ccr_mask <- crop(ccr_shp_rast, ccr_shp_proj) |> mask(ccr_shp_proj)

plot(ccr_mask)

#save raster of watershed area
writeRaster(x = ccr_mask, filename = "Downloaded_data/Watershed_rasters/ccr_ws.tif",overwrite =T)



#### TRIM DEM to watershed ----
 # SKIP FOR NOW: CAN JUST TRIM DEM IN GRASS WITH BASIN RASTER WE JUST MADE


# 2. Mask and crop in one shot
#    mask() accepts a SpatRaster as mask directly
dem_final <- crop(dem_proj, ccr_mask) |>
  mask(ccr_mask)


##check final plots
plot(dem_final)
plot(ccr_mask)

#write rasters for these
# writeRaster(dem_final, paste0(map_dir, "Outputs/ccr_dem.tif"), overwrite = TRUE)









