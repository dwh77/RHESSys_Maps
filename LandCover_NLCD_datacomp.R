## get vegetation from rhutils 1.4b

# devtools::install_github("ropensci/FedData")
library(FedData)
library(terra)

#function to put labels on map
textlocTL = function(srcmap) {
  pct = 0.95
  exts = ext(srcmap)
  x = exts[2] - ((exts[2] - exts[1])*pct)
  y = exts[3] + ((exts[4] - exts[3])*pct)
  return(c(x,y))
}



### set basin boundary Make box for watershed area
## for CCR square
basin_bounds = c( 37.4600, 37.3300, -79.8900, -80.1000) # N S E W
basin_ext = rast(crs = "EPSG:4326")
ext(basin_ext) = basin_bounds[c(4,3,2,1)]
values(basin_ext) = 1
plot(basin_ext)

#use these two lines below if getting a fetch_memory error
library(httr)
set_config(config(http_version = 1))  # 1 = HTTP/1.1

NLCD <- get_nlcd(template = basin_ext, 
                 label = 'CCRbox', 
                 year = 2021, 
                 force.redo = T, 
                 extraction.dir = paste0("Downloaded_Data/NLCD/"))


##plot nlcd
plot(NLCD)
crs(NLCD)
NLCD_reproj <- project(NLCD, "EPSG:26917")
plot(NLCD_reproj)
crs(NLCD_reproj)


##mask to CCR
mask_map = rast("Downloaded_Data/Watershed_rasters/ccr_watershed1.tif")
plot(mask_map)
ccr_mask <- project(mask_map, "EPSG:26917")

NLCD_proj = project(as(NLCD_reproj, "SpatRaster"), ccr_mask, method = "near")
plot(NLCD_proj)
ccr_vect <- as.polygons(ccr_mask)  # or st_as_sf() if using sf
plot(ccr_vect, add = TRUE, col = NA, border = "black", lwd = 3)
NLCD_crop = crop(NLCD_proj, ccr_mask)
NLCD_mask = mask(NLCD_crop, ccr_mask)
plot(NLCD_mask)

#### get land cover info table
atttab = levels(NLCD_mask)[[1]]
atttabclip = atttab[atttab$ID %in% unique(values(NLCD_mask, na.rm=T,mat=F)),]
row.names(atttabclip) = NULL
histtable = data.frame(ID = names(summary(as.factor(values(NLCD_mask, na.rm=T)))), Count = summary(as.factor(values(NLCD_mask, na.rm=T))))
histmerge = merge(histtable, atttabclip, by = "ID")

textlab = paste(capture.output(print(histmerge, row.names = F)), collapse = "\n")

plot(NLCD_mask, type = "classes", main = "NLCD baseline classes")
tl = textlocTL(NLCD_mask)
text(x = tl[1], y = tl[2], labels =textlab , col = "black", cex = 1, adj = c(0,1))

hist(NLCD_mask)
summary(as.factor(values(NLCD_mask)))
histmerge

#  ==================== RECLASS =====================
#updated trying to match numbers in stratum param library: https://github.com/RHESSys/ParameterLibrary/tree/master/Stratum
# 1 = non veg
# 2 = grass
# 3 = deciduous
# 4 = evergreen
#I'm calling mixed as deciduous for now

NLCD_mask_reclass = terra::classify(NLCD_mask,
                                    data.frame(c(11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95),
                                               c(1,  3,  3,  3,   3,  3,  3,  4,  3,  3,  3,  3,  3,  3,  3) )) #skipping grass for now
                                               #c(1,  1,  1,  1,   1,  1,  3,  4,  3,  2,  2,  2,  2,  2,  2) )) #previous try

##FROM RHUTILS SCRIPT
# NLCD_mask_reclass = terra::classify(NLCD_mask, data.frame(c(11, 21, 22, 31, 42, 52, 71, 90, 95),
#                                                           c(31, 31, 31, 31, 7,  50, 50,  7, 50) ))



# Plot to compare NLCD to RHESSys numbers
plot(NLCD_mask, type = "classes", main ="NLCD Land Cover")
plot(NLCD_mask_reclass, type = "classes", main = "RHESSys numbers")

# writeRaster(NLCD_mask_reclass, "Downloaded_Data/NLCD/NLCD_rhessys_veg_cover.tif", overwrite = T)

