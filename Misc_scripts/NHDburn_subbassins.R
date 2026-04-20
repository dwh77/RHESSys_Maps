library(terra)
library(whitebox)
library(sf)

wbt_init()

#### Paths ----
map_dir  <- "./CCR_NEW/Spatial_Data_comp/Create_Maps/"
out_dir  <- paste0(map_dir, "Hydro_Maps/")
dir.create(out_dir, showWarnings = FALSE)

dem_path    <- paste0(map_dir, "Outputs/ccr_dem.tif")
basin_file    <- paste0(map_dir, "Outputs/ccr_watershed.tif")
nhd_shp_path <- paste0(map_dir, "Outputs/nhd_flowlines_utm17n.shp")

stream_path <- paste0(map_dir, "Outputs/NHD_ccr_streams.tif")

dem     <- rast(dem_path)
streams <- rast(stream_path)



# ── Step 1: Align NHD stream raster to DEM grid and fill gaps ────────────────────────────
# Must be identical resolution, extent, and CRS before burning
streams_aligned <- project(streams, dem, method = "near")
plot(streams_aligned, col = "black")

# # Write to disk
# streams_bin_file <- paste0(out_dir, "nhd_streams_aligned.tif")
# writeRaster(streams_bin, streams_bin_file, overwrite = TRUE)

# Binary: 1 = stream, NA = not stream
streams_bin <- ifel(streams_aligned > 0, 1, NA)
plot(streams_bin)

# Convert NA → 0 so focal works cleanly
streams_0 <- ifel(is.na(streams_bin), 0, streams_bin)
plot(streams_0)

# 3x3 focal max — fills single-cell gaps
streams_filled <- focal(streams_0, w = 3, fun = "max", na.policy = "omit")
plot(streams_filled)

# If breaks are larger (diagonal jumps), try a 5x5 window
# streams_filled <- focal(streams_0, w = 5, fun = "max")

# Convert back to 1/NA for burning
streams_closed <- ifel(streams_filled > 0, 1, NA)

plot(streams_closed, main = "Streams with gaps closed")

# Write and proceed with manual burn
streams_closed_file <- paste0(out_dir, "nhd_streams_closed.tif")
writeRaster(streams_closed, streams_closed_file, overwrite = TRUE)









# ── Step 1: wbt_fill_burn() — the right tool for this job ──────────────────
# Takes: DEM raster + vector stream lines
# Does:  rasterizes streams, thins to single-cell width, burns them in,
#        then fills all depressions. One clean step.
# Note:  streams must be a VECTOR file (.shp), not a raster

burned_filled_file <- paste0(out_dir, "dem_burned_filled.tif")

wbt_fill_burn(
  dem     = dem_path,
  streams = nhd_shp_path,      # vector lines, not raster
  output  = burned_filled_file
)

plot(rast(burned_filled_file), main = "Burned + Filled DEM")

# Sanity check — difference from original should show stream cells lowered
dem_orig <- rast(dem_path)
diff <- dem_orig - dem_bf
diff[diff == 0] <- NA
plot(diff, main = "Where DEM was modified (positive = lowered)")
# You should see thin lines along stream channels








# ── Step 2: Burn streams into DEM ──────────────────────────────────────────
# This lowers the DEM elevation along stream cells by `burn_depth` meters,
# forcing D8 flow to follow the NHD network.
burn_depth <- 5  # meters — increase if streams still get ignored; can set to 100 to clearly see streams and check for gaps

dem_burned_file <- paste0(out_dir, "dem_burned.tif")

# Alternative burn approach (more explicit control):
# Manually lower stream cells in terra, then write out
dem_burned <- dem
dem_burned[streams_closed == 1] <- dem[streams_closed == 1] - burn_depth
writeRaster(dem_burned, dem_burned_file, overwrite = TRUE)

plot(dem_burned, main = paste0("Burned DEM (-", burn_depth, "m on streams)"))


























































# ── Step 3: Breach depressions on burned DEM ───────────────────────────────
breached_file <- paste0(out_dir, "dem_breached.tif")

wbt_breach_depressions(dem = dem_burned_file, output = breached_file)
# wbt_breach_depressions_least_cost(
#   dem    = dem_burned_file,
#   output = breached_file,
#   dist   = 10,
#   fill   = F
# )

plot(rast(breached_file))


# ── Step 4: D8 Flow Direction ───────────────────────────────────────────────
d8_file <- paste0(out_dir, "flow_dir_d8.tif")
wbt_d8_pointer(
  dem    = breached_file,
  output = d8_file
)
plot(rast(d8_file))


# ── Step 5: Flow Accumulation ───────────────────────────────────────────────
acc_file <- paste0(out_dir, "flow_accum.tif")
wbt_d8_flow_accumulation(
  input    = breached_file,
  output   = acc_file,
  out_type = "cells"
)

plot(rast(acc_file), main = "Flow Accumulation")


# ── Step 6: Extract stream network ─────────────────────────────────────────
# Because you burned the NHD network in, you have two options here:

# Option A: Use your NHD raster directly as the stream map (cleanest)
#           Do this if you trust the NHD network completely.
stream_file <- streams_closed_file

# Option B: Re-extract from accumulation (picks up streams forced by burning)
#           Use a LOW threshold — the burned cells will dominate
# stream_file <- paste0(out_dir, "streams_derived.tif")
# wbt_extract_streams(
#   flow_accum = acc_file,
#   output     = stream_file,
#   threshold  = 5000  # lower threshold since burning guides routing
# )

plot(rast(stream_file), main = "streams", col = "black")



# ── Picking up after Option A stream gap-closing ────────────────────────────
# At this point you have:
#   streams_closed_file  — gap-closed NHD stream raster
#   dem_burned_file      — DEM with streams burned in
#   breached_file        — breached/filled DEM
#   d8_file              — D8 flow direction
#   acc_file             — flow accumulation

# ── Step 1: Load and prepare your pour points shapefile ─────────────────────
# Your shapefile needs a column with unique integer IDs for each point.
# These IDs will become the subbasin values in the output raster.

pour_pts <- st_read(paste0(map_dir, "Outputs/pourpoints/pour_points24.shp"))

# Check what columns you have
print(names(pour_pts))
print(pour_pts)

# Make sure CRS matches DEM
pour_pts <- st_transform(pour_pts, crs = crs(dem, proj = TRUE))
crs(pour_pts)

# Ensure you have a clean integer ID column — rename/create as needed.
# wbt_watershed() uses the 'FID' or first integer field it finds.
# Safest to be explicit: add a column called 'id' with sequential integers.
pour_pts$id <- seq_len(nrow(pour_pts))

# Write back to disk
pour_pts_shp <- paste0(out_dir, "pour_points_clean.shp")
st_write(pour_pts, pour_pts_shp, delete_dsn = TRUE)

# ── Step 2: Snap pour points to stream network ──────────────────────────────
# Critical — points must land exactly on a stream cell.
# snap_dist is in map units (meters for UTM).
# Use a distance large enough to reach the stream but not jump to wrong reach.
snapped_shp <- paste0(out_dir, "pour_points_snapped.shp")

wbt_jenson_snap_pour_points(
  pour_pts  = pour_pts_shp,
  streams   = streams_closed_file,
  output    = snapped_shp,
  snap_dist = 5   # meters — adjust based on how close your points are to streams
)

# Verify snapping looks right
snapped_pts <- st_read(snapped_shp)
print(snapped_pts)

# Visual check — do snapped points sit on stream cells?
streams_r <- rast(streams_closed_file)
plot(streams_r, col = "blue", legend = FALSE, main = "Snapped pour points on streams")
plot(st_geometry(pour_pts),   add = TRUE, col = "red",   pch = 16, cex = 1.2)  # original
plot(st_geometry(snapped_pts), add = TRUE, col = "green", pch = 3,  cex = 1.5)  # snapped

# ── Step 3: Delineate subbasins from pour points ────────────────────────────
# wbt_watershed() with multiple pour points creates one subbasin per point.
# Each subbasin gets the integer ID value from the pour point.
# Areas not draining to any pour point get NoData.

subbasin_file <- paste0(out_dir, "subbasins_custom.tif")
wbt_watershed(
  d8_pntr  = d8_file,
  pour_pts = snapped_shp,
  output   = subbasin_file
)
subbasins <- rast(subbasin_file)
plot(subbasins, main = "Custom Subbasins from Pour Points")

# ── Step 4: Check coverage ───────────────────────────────────────────────────
# Are there any unassigned areas inside the watershed?
# This happens if a pour point was missed or snapping failed.
basin_mask <- rast(basin_file)  # your overall watershed basin from before
ws_mask    <- ifel(basin_mask > 0, 1, NA)
plot(ws_mask)

subbasins_clipped <- crop(subbasins, ws_mask) |> mask(ws_mask)
plot(subbasins_clipped)
# Cells in watershed but not assigned to any subbasin
unassigned <- ifel(is.na(subbasins_clipped) & !is.na(ws_mask), 1, NA)
plot(unassigned, main = "Unassigned cells (should be empty)")
freq(unassigned)  # should be 0 or very small

# # ── Step 5: Clip stream network to basin and finalize ───────────────────────
# streams_final   <- crop(rast(streams_closed_file), ws_mask) |> mask(ws_mask)
# flow_dir_final  <- crop(rast(d8_file),             ws_mask) |> mask(ws_mask)
# flow_acc_final  <- crop(rast(acc_file),             ws_mask) |> mask(ws_mask)
#
# # ── Step 6: Write all outputs ────────────────────────────────────────────────
# writeRaster(subbasins_clipped, paste0(out_dir, "subbasins_final.tif"), overwrite = TRUE)
# writeRaster(streams_final,     paste0(out_dir, "streams_final.tif"),   overwrite = TRUE)
# writeRaster(flow_dir_final,    paste0(out_dir, "flow_dir_final.tif"),  overwrite = TRUE)
# writeRaster(flow_acc_final,    paste0(out_dir, "flow_acc_final.tif"),  overwrite = TRUE)
#
# # ── Final check plot ─────────────────────────────────────────────────────────
# par(mfrow = c(1, 2))
# plot(subbasins_clipped, main = "Subbasins")
# plot(streams_final, col = "blue", legend = FALSE, add = FALSE, main = "Streams")
# par(mfrow = c(1, 1))








