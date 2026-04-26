###stream routing for CCR

library(terra)

generate_streamtable <- function(
    stream_rast,
    dem_rast,
    patch_rast,
    zone_rast,
    subbasin_rast,
    hill_rast,
    output_file,
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
  
  reach_ids <- sort(unique(na.omit(values(stream))))
  n_reaches <- length(reach_ids)
  message("Found ", n_reaches, " stream reaches.")
  
  # ── Mean elevation per reach ────────────────────────────────────────────────
  reach_elev <- sapply(reach_ids, function(rid) {
    cells <- which(values(stream) == rid)
    mean(values(dem)[cells], na.rm = TRUE)
  })
  names(reach_elev) <- reach_ids
  
  # ── Reach topology ──────────────────────────────────────────────────────────
  get_downstream_reach <- function(rid) {
    reach_cells <- which(values(stream) == rid)
    
    adj_stream_cells <- unique(unlist(lapply(reach_cells, function(cell) {
      nbrs <- adjacent(stream, cell, directions = 8)[1, ]
      nbrs <- nbrs[!is.na(nbrs)]
      nbrs[!is.na(values(stream)[nbrs]) & values(stream)[nbrs] != rid]
    })))
    
    if (length(adj_stream_cells) == 0) return(NA)
    
    adj_ids <- unique(na.omit(values(stream)[adj_stream_cells]))
    if (length(adj_ids) == 0) return(NA)
    
    my_elev <- reach_elev[as.character(rid)]
    downstream_candidates <- adj_ids[reach_elev[as.character(adj_ids)] < my_elev]
    
    if (length(downstream_candidates) == 0) return(NA)
    
    best <- downstream_candidates[which.max(reach_elev[as.character(downstream_candidates)])]
    return(best)
  }
  
  message("Building reach topology...")
  downstream_of <- setNames(
    sapply(reach_ids, get_downstream_reach),
    reach_ids
  )
  
  upstream_of <- lapply(reach_ids, function(rid) {
    reach_ids[!is.na(downstream_of) & downstream_of == rid]
  })
  names(upstream_of) <- reach_ids
  
  # ── Adjacent patches ────────────────────────────────────────────────────────
  message("Finding adjacent patches and hillslopes...")
  
  get_adjacent_patches <- function(rid) {
    reach_cells <- which(values(stream) == rid)
    
    nbr_cells <- unique(unlist(lapply(reach_cells, function(cell) {
      nbrs <- adjacent(stream, cell, directions = 8)[1, ]
      nbrs <- nbrs[!is.na(nbrs)]
      nbrs[is.na(values(stream)[nbrs])]
    })))
    
    if (length(nbr_cells) == 0) return(data.frame())
    
    adj_df <- data.frame(
      patch_id = values(patch)[nbr_cells],
      zone_id  = values(zone)[nbr_cells],
      hill_id  = values(hill)[nbr_cells]
    )
    adj_df <- na.omit(adj_df)
    adj_df <- unique(adj_df)
    
    return(adj_df)
  }
  
  # ── Reach geometry — inputs already in meters ───────────────────────────────
  get_reach_geometry <- function(rid) {
    cells <- which(as.vector(values(stream)) == rid)
    
    if (length(cells) == 0) {
      warning("Reach ", rid, " has 0 cells — skipping")
      return(list(length = NA, slope = NA))
    }
    
    elev_vals  <- as.vector(values(dem))[cells]
    res_m      <- res(dem)[1]  # already in meters
    
    elev_range <- diff(range(elev_vals, na.rm = TRUE))
    length_m   <- max(length(cells) * res_m, res_m)
    slope      <- max(elev_range / length_m, 0.0001)
    
    list(length = round(length_m, 2), slope = round(slope, 6))
  }
  
  # ── Write streamtable ───────────────────────────────────────────────────────
  message("Writing streamtable to: ", output_file)
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
          paste(row["patch_id"], row["zone_id"], row["hill_id"])
        ),
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



# ── Run it ───────────────────────────────────────────────────────────────────
generate_streamtable(
  stream_rast   = "CCR_files/spatial_data/ccr_stream1K.tif",
  dem_rast      = "CCR_files/spatial_data/ccr_dem.tif",
  patch_rast    = "CCR_files/spatial_data/ccr_patch_map_kmeans3000.tif",
  zone_rast     = "CCR_files/spatial_data/ccr_patch_map_kmeans3000.tif",
  subbasin_rast = "CCR_files/spatial_data/ccr_basin1K_filled.tif",
  hill_rast     = "CCR_files/spatial_data/ccr_basin1K_filled.tif",
  output_file   = "CCR_files/CCR_NLCDveg/ccr_mixed.stream",
  ManningsN      = 0.035,
  streamTopWidth = 2.0,
  streamBotWidth = 1.0,
  streamDepth    = 0.5
)



### function to readin file to check
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



####read in and check file
st <- read_streamtable("CCR_files/CCR_NLCDveg/ccr_mixed.stream")
# st <- read_streamtable("C:/Users/dwh18/OneDrive/Desktop/R_Projects/RHESSys_Tutorial/HPB_files_NewMaps/worldfiles/stream.hpb")


## sanity check on basin IDs and stream IDs
basin_r <- rast("CCR_files/spatial_data/ccr_basin1K_filled.tif")
plot(basin_r)

stream_r <- rast("CCR_files/spatial_data/ccr_stream1K.tif")

plot(stream_r, type = "classes")

##plot reaches with no upstream reaches
st_noup_df <- st |> filter(n_upstream == 0)
st_noup <- (st_noup_df$reach_id)

stream_headwater <- ifel(stream_r %in% st_noup, 1000, 1)
plot(stream_headwater, type="classes")

##Plot reach that has no downstream reaches
stream_outlet <- ifel(stream_r == 14, 1, 1000)
plot(stream_outlet, type="classes", , main = "reach 14 is 1")

#this looks like little piece that sticks above outflow.... 
#it must be pushing to the deepest depth in the bathy..

#whats upstream
stream_outlet <- ifel(stream_r == 16, 1, 1000)
plot(stream_outlet, type="classes", main = "reach 16 is 1")

#hmm arm to right; 16 has 2, 18, and 20 w
stream_outlet <- ifel(stream_r %in% c(2,18,16,20,14), 1, 1000)
plot(stream_outlet, type="classes")

#run trhough those above
stream_outlet <- ifel(stream_r ==2, 1, 1000)
plot(stream_outlet, type="classes", main = "reach 2 is 1")

#2 is the one that should be the outlet 
stream_outlet <- ifel(stream_r ==4, 1, 1000)
plot(stream_outlet, type="classes", main = "reach 4 is 1")

