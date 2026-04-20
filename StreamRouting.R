#### Trying stream routing generation in R
#https://github.com/RHESSys/RHESSys/wiki/In-Stream-Routing

library(terra)

generate_streamtable <- function(
    stream_rast,
    dem_rast,
    patch_rast,
    zone_rast,       # added back - pass patch map for your setup
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
  zone     <- rast(zone_rast)      # read zone raster directly
  subbasin <- rast(subbasin_rast)
  hill     <- rast(hill_rast)

  reach_ids <- sort(unique(na.omit(values(stream))))
  n_reaches <- length(reach_ids)
  message("Found ", n_reaches, " stream reaches.")

  # Mean elevation per reach — used to determine upstream/downstream direction
  reach_elev <- sapply(reach_ids, function(rid) {
    cells <- which(values(stream) == rid)
    mean(values(dem)[cells], na.rm = TRUE)
  })
  names(reach_elev) <- reach_ids

  # Build reach topology using mean elevation comparison.
  # For each reach, find all spatially adjacent reaches; those with lower mean
  # elevation are downstream candidates. Pick the highest among them (nearest).
  # If none are lower, this reach is the outlet.
  get_downstream_reach <- function(rid) {
    reach_cells <- which(values(stream) == rid)

    # All stream cells directly adjacent (8-connected) to this reach
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

    if (length(downstream_candidates) == 0) return(NA)  # this is the outlet

    # Among candidates, take the one with the highest elevation (most direct)
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

  # Adjacent patches - now reads zone directly from zone raster
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
      zone_id  = values(zone)[nbr_cells],   # pulled directly from zone raster
      hill_id  = values(hill)[nbr_cells]
    )
    adj_df <- na.omit(adj_df)
    adj_df <- unique(adj_df)

    return(adj_df)
  }

#was for degree to meter conversion
  # get_reach_geometry <- function(rid) {
  #   cells <- which(as.vector(values(stream)) == rid)
  #
  #   if (length(cells) == 0) {
  #     warning("Reach ", rid, " has 0 cells — skipping")
  #     return(list(length = NA, slope = NA))
  #   }
  #
  #   elev_vals   <- as.vector(values(dem))[cells]
  #   lat_center  <- mean(yFromCell(dem, cells), na.rm = TRUE)
  #   res_m_x     <- res(dem)[1] * 111320 * cos(lat_center * pi / 180)
  #   res_m_y     <- res(dem)[2] * 111320
  #   res_m       <- sqrt(res_m_x * res_m_y)
  #
  #   elev_range  <- diff(range(elev_vals, na.rm = TRUE))
  #   length_m    <- max(length(cells) * res_m, res_m)   # floor at 1 pixel
  #   slope       <- max(elev_range / length_m, 0.0001)
  #
  #   list(length = round(length_m, 2), slope = round(slope, 6))
  # }

  get_reach_geometry <- function(rid) {
    cells <- which(values(stream) == rid)
    elev_vals <- values(dem)[cells]
    elev_range <- diff(range(elev_vals, na.rm = TRUE))
    res_m    <- res(dem)[1]
    length_m <- length(cells) * res_m
    slope    <- max(elev_range / length_m, 0.0001)
    list(length = round(length_m, 2), slope = round(slope, 6))
  }

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


#### use function
generate_streamtable(
  stream_rast   = "./HPB_files_NewMaps/spatial_data/NAD83UTM17N_epsg26917/hpb_stream750.tif",
  dem_rast      = "./HPB_files_NewMaps/spatial_data/NAD83UTM17N_epsg26917/hpb_dem.tif",
  patch_rast    = "./HPB_files_NewMaps/spatial_data/NAD83UTM17N_epsg26917/hpb_patch_kmeans300.tif",
  zone_rast     = "./HPB_files_NewMaps/spatial_data/NAD83UTM17N_epsg26917/hpb_patch_kmeans300.tif",  # zone = patch
  subbasin_rast = "./HPB_files_NewMaps/spatial_data/NAD83UTM17N_epsg26917/hpb_basin750_filled_15apr26.tif",
  hill_rast     = "./HPB_files_NewMaps/spatial_data/NAD83UTM17N_epsg26917/hpb_basin750_filled_15apr26.tif",  # hill = subbasin
  output_file   = "./HPB_files_NewMaps/worldfiles/stream.hpb",
  ManningsN      = 0.035,
  streamTopWidth = 3.0,
  streamBotWidth = 1.5,
  streamDepth    = 0.6
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
st <- read_streamtable("./HPB_files_NewMaps/worldfiles/stream.hpb")





