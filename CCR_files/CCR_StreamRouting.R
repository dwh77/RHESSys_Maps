###stream routing for CCR

library(terra)

#### Functions ----

#### Make stream table
generate_streamtable <- function(
    stream_rast,
    dem_rast,
    patch_rast,
    zone_rast,
    subbasin_rast,
    hill_rast,
    output_file,
    outlet_id,
    ManningsN      = 0.035,
    streamTopWidth = 2.0,
    streamBotWidth = 1.0,
    streamDepth    = 0.5
) {
  stream   <- rast(stream_rast)
  dem      <- rast(dem_rast)
  patch    <- rast(patch_rast)
  zone     <- rast(zone_rast)
  subbasin <- rast(subbasin_rast)
  hill     <- rast(hill_rast)
  
  stream_vals <- values(stream)
  dem_vals    <- values(dem)
  patch_vals  <- values(patch)
  zone_vals   <- values(zone)
  hill_vals   <- values(hill)
  
  reach_ids <- sort(unique(na.omit(stream_vals)))
  n_reaches <- length(reach_ids)
  message("Found ", n_reaches, " stream reaches.")
  
  # ── Build full adjacency graph ───────────────────────────────────────────────
  message("Building reach adjacency graph...")
  reach_neighbors <- lapply(reach_ids, function(rid) {
    reach_cells <- which(stream_vals == rid)
    adj_cells <- unique(unlist(lapply(reach_cells, function(cell) {
      nbrs <- adjacent(stream, cell, directions = 8)[1, ]
      nbrs[!is.na(nbrs)]
    })))
    sort(unique(na.omit(stream_vals[adj_cells][stream_vals[adj_cells] != rid])))
  })
  names(reach_neighbors) <- as.character(reach_ids)
  
  # ── Trace topology upstream from known outlet ────────────────────────────────
  # Recursive traversal guarantees: exactly one outlet, correct many-to-one
  # confluences, no elevation ambiguity
  message("Tracing network upstream from outlet reach ", outlet_id, "...")
  
  downstream_of <- setNames(rep(NA_real_, n_reaches), as.character(reach_ids))
  visited <- c()
  
  assign_downstream <- function(rid, came_from = NA) {
    if (rid %in% visited) return()
    visited <<- c(visited, rid)
    if (!is.na(came_from)) {
      downstream_of[as.character(rid)] <<- came_from
    }
    unvisited_neighbors <- reach_neighbors[[as.character(rid)]]
    unvisited_neighbors <- unvisited_neighbors[!unvisited_neighbors %in% visited]
    for (n in unvisited_neighbors) {
      assign_downstream(n, came_from = rid)
    }
  }
  
  assign_downstream(outlet_id)
  
  # ── Derive upstream lookup from downstream assignments ───────────────────────
  upstream_of <- lapply(as.character(reach_ids), function(rid) {
    reach_ids[which(downstream_of[as.character(reach_ids)] == as.numeric(rid))]
  })
  names(upstream_of) <- as.character(reach_ids)
  
  # ── Topology report ──────────────────────────────────────────────────────────
  disconnected <- reach_ids[is.na(downstream_of[as.character(reach_ids)]) & 
                              reach_ids != outlet_id]
  message("\n=== Topology Summary ===")
  message("Outlet reach: ", outlet_id)
  message("Disconnected reaches (should be none): ", 
          ifelse(length(disconnected) == 0, "none", paste(disconnected, collapse = ", ")))
  if (length(disconnected) > 0) {
    warning("Disconnected reaches found: ", paste(disconnected, collapse = ", "),
            "\nThese reaches are not connected to the outlet in the stream raster.")
  }
  for (rid in reach_ids) {
    ups  <- upstream_of[[as.character(rid)]]
    down <- downstream_of[as.character(rid)]
    message("  Reach ", rid,
            " | up: ", ifelse(length(ups) == 0, "none", paste(ups, collapse = ",")),
            " | down: ", ifelse(is.na(down), "OUTLET", down))
  }
  
  # ── Adjacent patches ─────────────────────────────────────────────────────────
  message("\nFinding adjacent patches and hillslopes...")
  get_adjacent_patches <- function(rid) {
    reach_cells <- which(stream_vals == rid)
    nbr_cells <- unique(unlist(lapply(reach_cells, function(cell) {
      nbrs <- adjacent(stream, cell, directions = 8)[1, ]
      nbrs[!is.na(nbrs)]
    })))
    nbr_cells <- nbr_cells[is.na(stream_vals[nbr_cells])]
    if (length(nbr_cells) == 0) return(data.frame())
    adj_df <- data.frame(
      patch_id = patch_vals[nbr_cells],
      zone_id  = zone_vals[nbr_cells],
      hill_id  = hill_vals[nbr_cells]
    )
    unique(na.omit(adj_df))
  }
  
  # ── Reach geometry ───────────────────────────────────────────────────────────
  get_reach_geometry <- function(rid) {
    cells <- which(stream_vals == rid)
    if (length(cells) == 0) {
      warning("Reach ", rid, " has 0 cells")
      return(list(length = NA, slope = NA))
    }
    elev_vals  <- dem_vals[cells]
    res_m      <- res(dem)[1]
    elev_range <- diff(range(elev_vals, na.rm = TRUE))
    length_m   <- max(length(cells) * res_m, res_m)
    slope      <- max(elev_range / length_m, 0.0001)
    list(length = round(length_m, 2), slope = round(slope, 6))
  }
  
  # ── Write streamtable ────────────────────────────────────────────────────────
  message("\nWriting streamtable to: ", output_file)
  con <- file(output_file, "w")
  writeLines(as.character(n_reaches), con)
  
  for (rid in reach_ids) {
    geom  <- get_reach_geometry(rid)
    adj   <- get_adjacent_patches(rid)
    n_adj <- nrow(adj)
    ups   <- upstream_of[[as.character(rid)]]
    downs <- downstream_of[as.character(rid)]
    downs <- downs[!is.na(downs)]
    
    line <- paste(
      rid,
      streamBotWidth, streamTopWidth, streamDepth,
      geom$slope, ManningsN, geom$length,
      n_adj
    )
    if (n_adj > 0) {
      triplets <- paste(
        apply(adj, 1, function(row)
          paste(row["patch_id"], row["zone_id"], row["hill_id"])),
        collapse = " "
      )
      line <- paste(line, triplets)
    }
    line <- paste(line, length(ups), paste(ups, collapse = " "))
    line <- paste(line, length(downs), paste(downs, collapse = " "))
    writeLines(trimws(line), con)
  }
  
  close(con)
  message("Done. Streamtable has ", n_reaches, " reaches.")
}



#### Function to read in streamtable
read_streamtable <- function(file) {
  lines <- readLines(file)
  n_reaches <- as.integer(lines[1])
  
  rows <- lapply(lines[-1], function(line) {
    vals <- as.numeric(strsplit(trimws(line), "\\s+")[[1]])
    idx <- 1
    
    reach_id   <- vals[idx]; idx <- idx + 1
    bot_width  <- vals[idx]; idx <- idx + 1
    top_width  <- vals[idx]; idx <- idx + 1
    max_height <- vals[idx]; idx <- idx + 1
    slope      <- vals[idx]; idx <- idx + 1
    manning    <- vals[idx]; idx <- idx + 1
    length_m   <- vals[idx]; idx <- idx + 1
    n_adj      <- vals[idx]; idx <- idx + 1
    
    # Read all patch/zone/hill triplets
    patch_ids <- zone_ids <- hill_ids <- c()
    if (n_adj > 0) {
      for (i in 1:n_adj) {
        patch_ids <- c(patch_ids, vals[idx]);     idx <- idx + 1
        zone_ids  <- c(zone_ids,  vals[idx]);     idx <- idx + 1
        hill_ids  <- c(hill_ids,  vals[idx]);     idx <- idx + 1
      }
    }
    
    n_up     <- vals[idx]; idx <- idx + 1
    up_ids   <- if (n_up > 0) vals[idx:(idx + n_up - 1)] else NA
    idx      <- idx + n_up
    
    n_down   <- vals[idx]; idx <- idx + 1
    down_ids <- if (n_down > 0) vals[idx:(idx + n_down - 1)] else NA
    
    data.frame(
      reach_id       = reach_id,
      bot_width      = bot_width,
      top_width      = top_width,
      max_height     = max_height,
      slope          = slope,
      manning        = manning,
      length_m       = length_m,
      n_adj_patches  = n_adj,
      patch_ids      = paste(patch_ids, collapse = ","),
      zone_ids       = paste(zone_ids,  collapse = ","),
      hill_ids       = paste(hill_ids,  collapse = ","),
      n_upstream     = n_up,
      upstream_ids   = paste(up_ids,   collapse = ","),
      n_downstream   = n_down,
      downstream_ids = paste(down_ids, collapse = ",")
    )
  })
  
  do.call(rbind, rows)
} #end function



#### Read in basin and stream map to ID outlet ----
## sanity check on basin IDs and stream IDs
basin_r <- rast("CCR_files/CCR_nhdburn/spatial_data/ccr_basin1K_filled.tif")
stream_r <- rast("CCR_files/CCR_nhdburn/spatial_data/ccr_stream1K.tif")


## look at basin with number labels
#basins with labels
# Get centroids of each basin polygon
basin_polys <- as.polygons(basin_r)
basin_cents <- centroids(basin_polys)

# Extract coordinates and basin ID values
cents_df <- as.data.frame(basin_cents, geom = "XY")

# Plot the raster then overlay labels
plot(basin_r, type = "classes", legend = FALSE)
text(basin_cents, labels = basin_cents$`basin1K@PERMANENT`, cex = 0.6)

library(mapview)
library(sf)

# Convert to sf for mapview
basin_sf <- st_as_sf(basin_polys)

stream_polys <- as.polygons(stream_r)
stream_sf <- st_as_sf(stream_polys)

# Plot interactively - zoom, pan, click polygons to see attributes
#can remove the plus to just do one by one
mapview(basin_sf, zcol = "basin1K@PERMANENT", label = "basin1K@PERMANENT") +
  mapview(stream_sf, zcol = "stream1K@PERMANENT", label = "stream1K@PERMANENT")




#### Make stream table ----
generate_streamtable(
  stream_rast   = "CCR_files/CCR_nhdburn/spatial_data/ccr_stream1K.tif",
  dem_rast      = "CCR_files/CCR_nhdburn/spatial_data/ccr_dem.tif",
  patch_rast    = "CCR_files/CCR_nhdburn/spatial_data/ccr_patch_map_kmeans1000.tif",
  zone_rast     = "CCR_files/CCR_nhdburn/spatial_data/ccr_patch_map_kmeans1000.tif",
  subbasin_rast = "CCR_files/CCR_nhdburn/spatial_data/ccr_basin1K_filled.tif",
  hill_rast     = "CCR_files/CCR_nhdburn/spatial_data/ccr_basin1K_filled.tif",
  output_file   = "CCR_files/CCR_nhdburn/ccr_mixed.stream",
  outlet_id     = 2,   # <-- your known outlet
  ManningsN      = 0.035,
  streamTopWidth = 2.0,
  streamBotWidth = 1.0,
  streamDepth    = 0.5
)






####read in and check file ----
st <- read_streamtable("CCR_files/CCR_nhdburn/ccr_mixed1K.stream")
# sthpb <- read_streamtable("C:/Users/dwh18/OneDrive/Desktop/R_Projects/RHESSys_Tutorial/HPB_files_NewMaps/worldfiles/stream.hpb")




##plot reaches with no upstream reaches
st_noup_df <- st |> filter(n_upstream == 0)
st_noup <- (st_noup_df$reach_id)

stream_headwater <- ifel(stream_r %in% st_noup, 1000, 1)
plot(stream_headwater, type="classes")

##Plot reach that has no downstream reaches
st_nodown_df <- st |> filter(n_downstream == 0)
st_nodown_df <- (st_nodown_df$reach_id)

stream_outlet <- ifel(stream_r %in% st_nodown_df, 1000, 1)
plot(stream_outlet, type="classes")

#stream outlet
stream_outlet <- ifel(stream_r == 2, 1, 1000)
plot(stream_outlet, type="classes", , main = "reach 2 is 1")

#whats upstream
stream_outlet <- ifel(stream_r == 4, 1, 1000)
plot(stream_outlet, type="classes", main = "reach 4 is 1")


## notes on stream # to basin: can use viewer from above to check
#2 is the outlet
#62 is HPB: 64 drains into this; so gets main HPB piece
# 36 is SMB
# 28 is CCS 



