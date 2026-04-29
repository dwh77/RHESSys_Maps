# =============================================================================
# CCR Watershed Patch Map — 500 Patches, Squarer Patches Edition
# =============================================================================
#
# DESIGN
# ------
# Stage 1 — Reservoir patches
#   Cells with DEM <= 356.6 m within each hillslope start as one patch.
#   Any reservoir patch larger than `res_split_threshold` cells is split
#   into up to 3 sub-patches via k-means on (x, y) coordinates — this
#   produces along-length splits since the reservoir is roughly 1D.
#
# Stage 2 — Allocate remaining budget proportional to non-reservoir area.
#
# Stage 3 — Build squarer non-reservoir patches via 2D binning.
#   Per hillslope, factor n_allot into n_elev x n_pos based on the
#   hillslope's aspect ratio (relief vs horizontal extent). Then:
#     - elev_band = quantile bin of DEM (1..n_elev)
#     - pos_band  = quantile bin of along-axis coord (1..n_pos)
#   Each (elev_band, pos_band) pair is one patch.
#
# WHY 2D BINNING
# --------------
# Pure elevation bands produce long thin strips along contour lines. Adding a
# perpendicular position dimension breaks those strips into blockier shapes.
# The factoring tries to make each block roughly square in (elev, pos) space,
# which on the ground produces roughly square patches.
#
# INPUTS
#   dem        <- rast("./spatial_data/ccr_dem.tif")
#   hillslopes <- rast("./spatial_data/ccr_halfbasin1K_filled.tif")
#   nlcd       <- rast("../../Downloaded_Data/NLCD/NLCD_rhessys_veg_cover_mixed.tif")
#
# OUTPUT
#   ccr_patch_map_500.tif
# =============================================================================

library(terra)

# ── Settings ──────────────────────────────────────────────────────────────────

target_patches      <- 500
reservoir_elev      <- 356.7
res_split_threshold <- 80      # reservoir patches > this size get split
res_max_subpatches  <- 3
elev_weight         <- 1.5     # elevation influence in non-reservoir k-means
out_path       <- "CCR_files/CCR_nhdburn/ccr_patch_map_500.tif"

# ── Load and align ────────────────────────────────────────────────────────────

dem        <- rast("CCR_files/CCR_nhdburn/spatial_data/ccr_dem.tif")
hillslopes <- rast("CCR_files/CCR_nhdburn/spatial_data/ccr_halfbasin1K_filled.tif")
nlcd       <- rast("Downloaded_Data/NLCD/NLCD_rhessys_veg_cover_mixed.tif")

hillslopes <- resample(hillslopes, dem, method = "near")
hillslopes <- mask(hillslopes, dem)
hillslopes <- as.int(hillslopes)

xy <- xyFromCell(dem, seq_len(ncell(dem)))
df <- data.frame(
  cell = seq_len(ncell(dem)),
  x    = xy[, 1],
  y    = xy[, 2],
  dem  = values(dem,        mat = FALSE),
  hill = values(hillslopes, mat = FALSE)
)
df <- df[!is.na(df$dem) & !is.na(df$hill), ]
df$is_res <- df$dem <= reservoir_elev

message("Total valid cells: ",   nrow(df))
message("Reservoir cells: ",     sum(df$is_res))
message("Hillslopes: ",          length(unique(df$hill)))

df$patch <- NA_integer_

# =============================================================================
# STAGE 1 — Reservoir patches
# =============================================================================

message("\n--- Stage 1: Reservoir patches ---")

res_hills <- sort(unique(df$hill[df$is_res]))
next_id   <- 0L

for (h in res_hills) {
  rows    <- which(df$hill == h & df$is_res)
  n_cells <- length(rows)
  
  if (n_cells <= res_split_threshold) {
    next_id <- next_id + 1L
    df$patch[rows] <- next_id
  } else {
    k <- min(res_max_subpatches, max(2, floor(n_cells / res_split_threshold)))
    coords <- df[rows, c("x", "y")]
    set.seed(1)
    km <- kmeans(coords, centers = k, nstart = 5, iter.max = 50)
    df$patch[rows] <- next_id + km$cluster
    next_id <- next_id + k
    message("  Hillslope ", h, ": ", n_cells, " reservoir cells -> ", k, " patches")
  }
}

n_res <- next_id
message("Total reservoir patches: ", n_res)

# =============================================================================
# STAGE 2 — Allocate remaining budget
# =============================================================================

remaining_budget <- target_patches - n_res
message("\n--- Stage 2: Allocating ", remaining_budget, " non-reservoir patches ---")

nonres_counts <- table(df$hill[!df$is_res])
hill_ids_nr   <- as.integer(names(nonres_counts))
hill_areas    <- as.integer(nonres_counts)
total_nonres  <- sum(hill_areas)

hill_alloc <- round(remaining_budget * hill_areas / total_nonres)
hill_alloc[hill_alloc < 1] <- 1

while (sum(hill_alloc) > remaining_budget) {
  biggest <- which.max(hill_alloc)
  hill_alloc[biggest] <- hill_alloc[biggest] - 1
}
while (sum(hill_alloc) < remaining_budget) {
  underserved <- which.max(hill_areas / hill_alloc)
  hill_alloc[underserved] <- hill_alloc[underserved] + 1
}

names(hill_alloc) <- hill_ids_nr
message("Allocation summary across hillslopes:")
print(summary(hill_alloc))

# =============================================================================
# STAGE 3 — Non-reservoir k-means per hillslope
# =============================================================================

message("\n--- Stage 3: Elevation-weighted k-means per hillslope ---")

for (i in seq_along(hill_ids_nr)) {
  
  h       <- hill_ids_nr[i]
  n_allot <- hill_alloc[i]
  
  rows <- which(df$hill == h & !df$is_res)
  if (length(rows) == 0) next
  
  if (n_allot == 1 || length(rows) <= n_allot) {
    # Trivial: one patch (or fewer cells than allocated patches)
    next_id <- next_id + 1L
    df$patch[rows] <- next_id
    next
  }
  
  # Build feature matrix: standardized (x, y, elev * weight)
  # Standardization uses each variable's own sd within this hillslope
  xs <- df$x[rows]
  ys <- df$y[rows]
  es <- df$dem[rows]
  
  xs_s <- if (sd(xs) > 0) (xs - mean(xs)) / sd(xs) else xs * 0
  ys_s <- if (sd(ys) > 0) (ys - mean(ys)) / sd(ys) else ys * 0
  es_s <- if (sd(es) > 0) (es - mean(es)) / sd(es) else es * 0
  
  feats <- cbind(xs_s, ys_s, es_s * elev_weight)
  
  set.seed(i)  # reproducible per hillslope
  km <- tryCatch(
    kmeans(feats, centers = n_allot, nstart = 10, iter.max = 100),
    error = function(e) NULL
  )
  
  if (is.null(km)) {
    # Fall back: collapse to single patch if k-means fails
    next_id <- next_id + 1L
    df$patch[rows] <- next_id
    next
  }
  
  df$patch[rows] <- next_id + km$cluster
  next_id <- next_id + n_allot
}

# =============================================================================
# Compact and write
# =============================================================================

message("\n--- Finalizing ---")

unique_patches <- sort(unique(df$patch[!is.na(df$patch)]))
id_map <- setNames(seq_along(unique_patches), unique_patches)
df$patch <- id_map[as.character(df$patch)]

n_total <- max(df$patch, na.rm = TRUE)
message("Total patches: ", n_total, "  (target: ", target_patches, ")")

patch_map <- dem
values(patch_map) <- NA_integer_
patch_map[df$cell] <- df$patch
patch_map <- as.int(patch_map)

writeRaster(patch_map, out_path, overwrite = TRUE, datatype = "INT4S")
message("Patch map written to: ", out_path)

# =============================================================================
# Visualization
# =============================================================================
plot(patch_map)
unique(patch_map)

n <- max(values(patch_map), na.rm = TRUE)
set.seed(42)
plot(patch_map, col = sample(rainbow(n)), type = "classes",
     legend = FALSE, main = paste(n, "patches"))

# Draw patch outlines on top of the DEM
plot(dem, col = grey.colors(50), legend = FALSE, main = "Patch boundaries")
boundaries <- boundaries(patch_map, classes = TRUE)
plot(boundaries, add = TRUE, col = c(NA, "red"), legend = FALSE)

# Build a 3-class background: NA / reservoir / non-reservoir
bg <- dem
bg[] <- NA
bg[dem <= 356.7] <- 1   # reservoir
bg[dem >  356.7] <- 2   # non-reservoir

plot(bg,
     col = c("steelblue", "grey85"),
     type = "classes",
     levels = c("Reservoir (≤ 356.6 m)", "Upland"),
     main = "Reservoir extent + patch boundaries")

# Overlay patch boundaries on top
b <- boundaries(patch_map, classes = TRUE)
b[b == 0] <- NA   # drop interior cells so only edges plot
plot(b, add = TRUE, col = "red", legend = FALSE)


# =============================================================================
# VALIDATION
# =============================================================================

# message("\n--- Validation ---")
# 
# hill_per_patch <- tapply(df$hill, df$patch, function(x) length(unique(x)))
# n_cross <- sum(hill_per_patch > 1)
# message("Patches crossing hillslope boundary: ", n_cross, "  (should be 0)")
# 
# patch_sizes <- table(df$patch)
# message("\nPatch size summary (cells):")
# print(summary(as.integer(patch_sizes)))
# 
# patch_extents <- aggregate(cbind(x, y) ~ patch, data = df,
#                            FUN = function(v) diff(range(v)))
# patch_extents$aspect <- pmax(patch_extents$x, patch_extents$y) /
#   pmax(pmin(patch_extents$x, patch_extents$y), 1)
# message("\nPatch aspect ratio (1 = square, higher = elongated):")
# print(summary(patch_extents$aspect))
# 
# message("\nDone. Output: ", out_path)


