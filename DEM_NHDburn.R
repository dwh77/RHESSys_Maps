#### Get NHD flowlines and burn to DEM

#### packages
library(tidyverse)
library(terra)
library(sf)
library(nhdplusTools)

getwd()
#### Read in DEM and CCR outline from 1_GetDEM.R
dem <- rast("./Downloaded_Data/FED_DEM/CCR_NED_1_nad83utm17N.tif")
crs(dem)
plot(dem)
  
ccrWS <- rast("Downloaded_data/Watershed_rasters/ccr_ws.tif")
crs(ccrWS)
plot(ccrWS)

plot(dem)
plot(ccrWS, add= T)

#### Get NHD flowlines ----
## Set up basin area to download
basin_bounds = c( 37.4600, 37.3300, -79.8900, -80.1000) # N S E W

basin_ext = rast(crs = "EPSG:4326")
ext(basin_ext) = basin_bounds[c(4,3,2,1)]
values(basin_ext) = 1

plot(basin_ext)

basin_sf <- sf::st_as_sf(as.polygons(ext(basin_ext), crs = crs(basin_ext)))


# Download NHDPlus flowlines clipped to that area
flowlines <- nhdplusTools::get_nhdplus(AOI = basin_sf, realization = "flowline")

plot(flowlines$streamorde)

# Rasterize flowlines to match your reference raster exactly
flowlines_reproj <- st_transform(flowlines, crs(ccrWS))
plot(flowlines_reproj)
plot(flowlines_reproj$streamorde)
plot(flowlines_reproj["streamorde"])

z <- crop(flowlines_reproj, ccrWS)

#export shapefile
st_write(flowlines_reproj, "Downloaded_Data/NHD_flowlines_download/nhd_flowlines_utm17n.shp", delete_dsn = TRUE)


#make raster
streams_rast     <- rasterize(vect(flowlines_reproj), ccrWS, field = 1)

plot(streams_rast)
crs(streams_rast)

#look at streams to basin
plot(dem)
plot(streams_rast, add = T)

# ## looks like I'm not getting all the ephemeral streams... see whats up
# table(flowlines$ftype)
# unique(flowlines$ftype)
# unique(flowlines$fcode)



### DONT USE: This is to get expanded nhd plus data that includes additional intermittent streams
# #Step 1: download once (only need to do this once to get data locally)
# nhdplusTools::download_nhdplushr(
#   nhd_dir = "Downloaded_Data/NHD_flowlines_download/",
#   hu_list = "0301"
# )
# 
# # Step 2: then read from local files
# flowlines_hr <- nhdplusTools::get_nhdplushr(
#   hr_dir = "Downloaded_Data/NHD_flowlines_download/03",
#   layers = "NHDFlowline"
# )$NHDFlowline
# 
# 
# #check projections
# crs(flowlines_hr)
# 
# # ##Update this with DEM thats in UTM
# basin_sf <- st_as_sf(as.polygons(ext(ccrWS), crs = crs(ccrWS)))
# plot(basin_sf)
# 
# #match projection for flowlines to basin
# flowlines_hr_reproj <- st_transform(flowlines_hr, crs(ccrWS))
# 
# #filter flowlines to CCR box
# flowlines_fin <- flowlines_hr_reproj |>  sf::st_filter(basin_sf)
# 
# #plot flowlines and check class stypes
# plot(flowlines_fin)
# crs(flowlines_fin)
# 
# table(flowlines_fin$FTYPE)
# unique(flowlines_fin$FTYPE)
# unique(flowlines_fin$FCODE)
# 
# #write flowlines fin
# nhd_shp_path <- "Downloaded_Data/NHD_flowlines_download/nhd_flowlines_utm17n.shp"
# st_write(flowlines_fin, nhd_shp_path, delete_dsn = TRUE)
# 
# #turn flowlines into a rater
# streams_rast <- rasterize(vect(flowlines_fin), ccrWS, field = 1)
# plot(streams_rast)
# streams_rast <- mask(streams_rast, ccrWS)
# plot(streams_rast)
# 
# #overlay flowliens on basin and DEM
# plot(dem)
# plot(streams_rast, add = T, col = "orange")
# plot(streams_rast, col = "black")
# 
# #save flowlines raster
# writeRaster(streams_rast, "Downloaded_Data/NHD_flowlines_download/NHD_ccr_streams.tif", overwrite = TRUE)
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
#                                 "Other" = "red"))




#### make interactive plot to pull out pour points ----
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






#### Clean up DEM and burn NHD streams ----

#### Clean up DEM 
plot(dem)

dem_masked <- dem |> crop(ccrWS) |> mask(ccrWS)
plot(dem_masked)

## see what weird shadow eleavtion near dam is
plot(ifel(dem_masked < 356.6, dem_masked, 0), main = "Elev < 356.6")
plot(ifel(dem_masked < 355, dem_masked, 0), main = "Elev < 355")
plot(ifel(dem_masked < 354, dem_masked, 0), main = "Elev < 354")

#just going to call everything less than known full pond as 356.6 to keep consistent
dem_fullpond <- ifel(dem <= 356.6, 356.6, dem)
plot(dem_fullpond)

writeRaster(dem_fullpond, "./CCR_files/CCR_nhdburn/spatial_data/ccr_DEM_fullpond.tif", overwrite = TRUE)


#### Burn NHD streams 
# Must be identical resolution, extent, and CRS before burning
streams_aligned <- project(streams_rast, dem_fullpond, method = "near")
plot(streams_aligned, col = "black")

#Skip the gap filling for now
# # Binary: 1 = stream, NA = not stream
# streams_bin <- ifel(streams_aligned > 0, 1, NA)
# plot(streams_bin)
# 
# # Convert NA → 0 so focal works cleanly
# streams_0 <- ifel(is.na(streams_bin), 0, streams_bin)
# plot(streams_0)
# 
# # 3x3 focal max — fills single-cell gaps
# streams_filled <- focal(streams_0, w = 3, fun = "max", na.policy = "omit")
# plot(streams_filled)
# 
# # If breaks are larger (diagonal jumps), try a 5x5 window
# # streams_filled <- focal(streams_0, w = 5, fun = "max")
# 
# # Convert back to 1/NA for burning
# streams_closed <- ifel(streams_filled > 0, 1, NA)
# 
# plot(streams_closed, main = "Streams with gaps closed")


#### burn in streams 
streams_closed <- streams_aligned 
# This lowers the DEM elevation along stream cells by `burn_depth` meters,
burn_depth <- 3  # meters — increase if streams still get ignored; can set to 100 to clearly see streams and check for gaps

dem_burned <- dem_fullpond
dem_burned[streams_closed == 1] <- dem[streams_closed == 1] - burn_depth

plot(dem_burned, main = paste0("Burned DEM (-", burn_depth, "m on streams)"))

writeRaster(dem_burned, "./CCR_files/CCR_nhdburn/spatial_data/ccr_DEM_burned2.tif", overwrite = TRUE)












