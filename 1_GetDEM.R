#### GET DEM
## maps need to make maps for CCR basins
## has code for getting NHD flowlines at end if that's of interest

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



#---------------------------- SKIPPING THIS FOR NOW: BUT THIS ALLOWS YOU TO GET NHD FLOWLINES 
#### Get NHD flowlines ----

# #get basin in right format
# basin_bounds = c( 37.4600, 37.3300, -79.8900, -80.1000) # N S E W
# 
# basin_ext = rast(crs = "EPSG:4326")
# ext(basin_ext) = basin_bounds[c(4,3,2,1)]
# values(basin_ext) = 1
# 
# plot(basin_ext)
# 
# basin_sf <- st_as_sf(as.polygons(ext(basin_ext), crs = crs(basin_ext)))
# 
# 
# ## This was trial with basic get_nhdplus but it misses all the intermittent stream classes
# # # Download NHDPlus flowlines clipped to that area
# # flowlines <- nhdplusTools::get_nhdplus(AOI = basin_sf, realization = "flowline")
# #
# # plot(flowlines$streamorde)
# #
# # # Rasterize flowlines to match your reference raster exactly
# # flowlines_reproj <- st_transform(flowlines, crs(basin_sf))
# # plot(flowlines_reproj)
# #
# # streams_rast     <- rasterize(vect(flowlines_reproj), basin_ext, field = 1)
# #
# # plot(streams_rast)
# #
# # plot(basin_ext)
# # plot(NED)
# # plot(streams_rast, add = T, col = "orange")
# #
# # ## looks like I'm not getting all the ephemeral streams... see whats up
# # table(flowlines$ftype)
# # unique(flowlines$ftype)
# # unique(flowlines$fcode)
# 
# 
# ## Get expanded nhd plus data that includes intermittent streams
# 
# # #Step 1: download once (only need to do this once)
# # nhdplusTools::download_nhdplushr(
# #   nhd_dir = paste0(map_dir, "NHD_flowlines_download/"),
# #   hu_list = "0301"
# # )
# 
# # Step 2: then read from local files
# flowlines_hr <- nhdplusTools::get_nhdplushr(
#   hr_dir = paste0(map_dir, "NHD_flowlines_download/03"),
#   layers = "NHDFlowline"
# )$NHDFlowline
# 
# #check projections
# crs(flowlines_hr)
# # crs(basin_sf)
# 
# # ##Update this with DEM thats in UTM
# dem_final <- rast(paste0(map_dir, "Outputs/ccr_dem.tif"))
# plot(dem_final)
# crs(dem_final)
# 
# basin_sf <- st_as_sf(as.polygons(ext(dem_final), crs = crs(dem_final)))
# plot(basin_sf)
# 
# #match projection for flowlines to basin
# flowlines_hr_reproj <- st_transform(flowlines_hr, crs(dem_final))
# 
# #filter flowlines to CCR box
# flowlines_fin <- flowlines_hr_reproj |>  sf::st_filter(basin_sf)
# 
# #plot flowlines and check class stypes
# plot(flowlines_fin)
# 
# table(flowlines_fin$FTYPE)
# unique(flowlines_fin$FTYPE)
# unique(flowlines_fin$FCODE)
# 
# #write flowlines fin
# nhd_shp_path <- paste0(map_dir, "Outputs/nhd_flowlines_utm17n.shp")
# st_write(flowlines_fin, nhd_shp_path, delete_dsn = TRUE)
# 
# #turn flowlines into a rater
# streams_rast <- rasterize(vect(flowlines_fin), dem_final, field = 1)
# streams_rast <- mask(streams_rast, dem_final)
# 
# #overlay flowliens on basin and DEM
# plot(dem_final)
# plot(streams_rast, add = T, col = "orange")
# plot(streams_rast, col = "black")
# 
# #save flowlines raster
# writeRaster(streams_rast, paste0(map_dir, "Outputs/NHD_ccr_streams.tif"), overwrite = TRUE)
# 
# 
# 
# ## look at intermittent v perennial
# # Reclassify FType to readable labels
# # https://files.hawaii.gov/dbedt/op/gis/data/NHD%20Complete%20FCode%20Attribute%20Value%20List.pdf
# flowlines_fin$stream_type <- dplyr::case_when(
#   flowlines_fin$FCODE == 46006 ~ "Perennial",
#   flowlines_fin$FCODE == 46003 ~ "Intermittent",
#   flowlines_fin$FCODE == 55800  ~ "Artifical Path",
#   flowlines_fin$FCODE == 55800  ~ "Connector",
#   TRUE ~ "Other"
# )
# 
# # Plot
# ggplot(flowlines_fin) +
#   geom_sf(aes(color = stream_type)) +
#   scale_color_manual(values = c("Ephemeral" = "tan",
#                                 "Intermittent" = "steelblue",
#                                 "Perennial" = "navy",
#                                 "Other" = "grey50"))
# 
# 
# 
# 
# #### make interactive plot to pull out pour points ----
# mapview(dem_final, alpha = 0.6) + mapview(streams_rast, col.regions = "blue")
# 
# # Click your pour points on the map, then click "Done"
# pour_points <- mapedit::drawFeatures(
#   mapview(dem_final, alpha = 0.6) + mapview(streams_rast, col.regions = "blue")
# )
# 
# # Extract the coordinates as a dataframe
# coords <- st_coordinates(pour_points)
# print(coords)
# 
# # Optionally save as a shapefile for use in the wbt snapping step later
# st_write(pour_points, paste0(map_dir, "Outputs/pourpoints/pour_points24.shp"), delete_dsn = TRUE)






