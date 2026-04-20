library(terra)
library(sf)
library(tmap) #for nicer plots

#Note bathy comes from 2022 bathymetry EDI publication and the TIN was converted to a raster in ArcGIS Pro

##reproject bathy to DEM
getwd()
dem <- rast("./Downloaded_Data/FED_DEM/CCR_NED_1_nad83utm17N.tif")
plot(dem)

bathy <- rast("./Downloaded_Data/CCR_bathy/CCR_bathy_raster.tif")
plot(bathy)

#make 0's NAs
bathy <- ifel(bathy == 0, NA, bathy)
plot(bathy)
crs(bathy)

##match projection to dem
bathy_proj <- project(bathy, crs(dem))
plot(bathy_proj)
crs(bathy_proj)

##match sample size
bathy_proj_resamp <- resample(bathy_proj, dem, method = "bilinear")
plot(bathy_proj_resamp)



#### Append bathy to DEM ----
#note that DEM is in meters elevation and bathy is already in elevation so no need to covert raster units
summary(dem)
summary(bathy_proj_resamp)


# Option B: DEM has values (e.g. water surface) over the reservoir
# Use the bathy extent as a mask to blank out the reservoir in the DEM first
reservoir_mask <- !is.na(bathy_proj_resamp)
plot(reservoir_mask)

dem_noRes <- mask(dem, reservoir_mask, maskvalue = TRUE)
plot(dem_noRes)

combined <- merge(dem_noRes, bathy_proj_resamp)
plot(combined)

#check this really worked
plot(dem) #still hard to tell, try diff color pallete

tm_shape(dem)+
  tm_raster(style = "cont", palette = "Spectral", legend.show = T)+
  tm_scale_bar()+ tm_title("DEM only")

tm_shape(combined)+
  tm_raster(style = "cont", palette = "Spectral", legend.show = T)+
  tm_scale_bar()+ tm_title("DEM plus bathy")


writeRaster(combined, "./Downloaded_Data/CCR_bathy/CCR_bathy_DEM_bind.tif", overwrite = TRUE)


