#### Make map of CCR and streams that can be used in ASLO talk 

library(terra)
library(tidyverse)
library(tidyterra)
library(sf)
library(maptiles)       # for satellite basemap
library(ggspatial)      # for scale bar & north arrow
library(ggnewscale)

# ── Load data ────────────────────────────────────────────────────────────────
ws        <- rast("./CCR_files/CCR/spatial_data/ccr_ws.tif")
streams1K <- rast("./CCR_files/CCR/spatial_data/ccr_stream1K.tif")
dem       <- rast("./CCR_files/CCR/spatial_data/ccr_dem.tif")

# ── Reproject ─────────────────────────────────────────────────────────────────
target_crs   <- "EPSG:3857"
ws_proj      <- project(ws,        target_crs)
streams_proj <- project(streams1K, target_crs)
dem_proj     <- project(dem,       target_crs)

library(tmap)
tm_shape(dem)+
  tm_raster(style = "cont", palette = "-Spectral", legend.show = F)
  #tm_scale_bar()+ tm_title("DEM only")

# ── Watershed boundary ────────────────────────────────────────────────────────
ws_poly <- ws_proj |>
  as.polygons() |>
  st_as_sf() |>
  st_union()

# ── Reservoir mask ────────────────────────────────────────────────────────────
reservoir    <- ifel(dem_proj <= 356.7, 1, NA)
plot(reservoir)

# ── Stream mask ───────────────────────────────────────────────────────────────
streams_mask <- ifel(streams_proj > 0, 1, NA)
plot(streams_mask)

# ── Convert masks to factor rasters (needed for scale_fill_manual) ────────────
streams_mask <- as.factor(streams_mask)
reservoir    <- as.factor(reservoir)

# ── Satellite basemap ─────────────────────────────────────────────────────────
# bbox_sf <- st_as_sf(as.polygons(ext(ws_proj), crs = target_crs))
# basemap <- get_tiles(bbox_sf, provider = "Esri.WorldImagery", zoom = 13, crop = TRUE)

# ── Expand basemap bbox before fetching tiles ─────────────────────────────────
# Expand the extent by 20% on each side for breathing room
ws_ext     <- ext(ws_proj)
x_pad      <- (ws_ext$xmax - ws_ext$xmin) * 0.05
y_pad      <- (ws_ext$ymax - ws_ext$ymin) * 0.05
bbox_exp   <- ext(ws_ext$xmin - x_pad, ws_ext$xmax + x_pad,
                  ws_ext$ymin - y_pad, ws_ext$ymax + y_pad)
bbox_sf    <- st_as_sf(as.polygons(bbox_exp, crs = target_crs))
basemap    <- get_tiles(bbox_sf, provider = "Esri.WorldImagery", zoom = 13, crop = TRUE)


# ── Colours ───────────────────────────────────────────────────────────────────
stream_col    <- "white"  #"#6BC4C8"
reservoir_col <- "#2166AC"
boundary_col  <- "black"

# ── Site locations ─────────────────────────────────────────────────────────────
ccr_dam <- data.frame(name = "CCR Dam", lon = -79.958,  lat = 37.3697)
hpb     <- data.frame(name = "HPB",     lon = -79.9733, lat = 37.3646)
ccs     <- data.frame(name = "CCS",     lon = -79.9937, lat = 37.3958)
smb     <- data.frame(name = "SMB",     lon = -79.9771, lat = 37.4169)

sites <- rbind(ccr_dam, hpb, ccs, smb) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(3857)


# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot() +
  
  # 1. Satellite basemap
  layer_spatial(basemap) +
  
  # 2. Streams
  geom_spatraster(data = streams_mask, na.rm = TRUE) +
  scale_fill_manual(
    values      = c("1" = stream_col),
    na.value    = "transparent",
    na.translate = FALSE,
    name        = NULL,
    labels      = "Stream network"
  ) +
  
  # 3. Reservoir on top
  new_scale_fill() +
  geom_spatraster(data = reservoir, na.rm = TRUE) +
  scale_fill_manual(
    values      = c("1" = reservoir_col),
    na.value    = "transparent",
    na.translate = FALSE,
    name        = NULL,
    labels      = "Reservoir"
  ) +
  
  # 4. Watershed boundary
  geom_sf(data = ws_poly, fill = "transparent",
          colour = boundary_col, linewidth = 0.8,
          show.legend = FALSE) +
  
  # 5. Site points
  new_scale_fill() +
  geom_sf(data   = sites,
          aes(fill = name, shape = name),
          colour = "black",
          size   = 3,
          stroke = 0.5) +
  scale_fill_manual(
    name   = NULL,
    values = c(
      "CCR Dam" = "#FF4444",
      "HPB"     = "#4FC3F7",
      "CCS"     = "#81C784",
      "SMB"     = "#FFB74D"
    )
  ) +
  scale_shape_manual(
    name   = NULL,
    values = c(
      "CCR Dam" = 24,
      "HPB"     = 22,
      "CCS"     = 22,
      "SMB"     = 22
    )
  ) +
  
  # ── Annotations ─────────────────────────────────────────────────────────────
  annotation_scale(
    location  = "br", width_hint = 0.2,
    text_col  = "white", line_col = "white",
    bar_cols  = c("white", "grey40")
  ) +
  annotation_north_arrow(
    location = "br", pad_y = unit(0.65, "cm"),
    style    = north_arrow_fancy_orienteering(
      fill     = c("white", "grey30"),
      text_col = "white"
    )
  ) +
  
  # ── Extent ──────────────────────────────────────────────────────────────────
  coord_sf(
    xlim   = c(ws_ext$xmin - x_pad, ws_ext$xmax + x_pad),
    ylim   = c(ws_ext$ymin - y_pad, ws_ext$ymax + y_pad),
    expand = FALSE,
    datum  = sf::st_crs(4326)
  ) +
  
  # ── Theme ───────────────────────────────────────────────────────────────────
  theme_void() +
  theme(
    plot.background      = element_rect(fill = "white", colour = NA),
    axis.text            = element_text(colour = "black", size = 7),
    axis.text.x          = element_text(margin = margin(t = 4)),
    axis.text.y          = element_text(margin = margin(r = 4)),
    axis.ticks           = element_line(colour = "black", linewidth = 0.3),
    axis.ticks.length    = unit(0.15, "cm"),
    legend.position      = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = alpha("black", 0.55), colour = NA),
    legend.text          = element_text(colour = "white", size = 9),
    legend.key.size      = unit(0.4, "cm"),
    legend.key           = element_blank(),      # removes the gray box
    legend.margin        = margin(6, 8, 6, 8),
    legend.spacing.y     = unit(0.1, "cm")
  )

p

# ggsave("ccr_watershed_map.png", p,
#        width = 10, height = 8, dpi = 300, bg = "black")




############### US MAP ----------------------
library(maps)

# ── Site location ─────────────────────────────────────────────────────────────
site <- data.frame(lon = -79.97572, lat = 37.39507)

# ── US map data ───────────────────────────────────────────────────────────────
us_states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))

# ── Plot ──────────────────────────────────────────────────────────────────────
us_map <- ggplot() +
  
  geom_sf(data = us_states,
          fill     = "grey50",
          colour   = "grey20",
          linewidth = 0.3) +
  
  geom_point(data = site,
             aes(x = lon, y = lat),
             colour = "#FF4444",
             size   = 5,
             shape  = 21,
             fill   = "#FF4444") +
  
  # Optional: subtle label
  # geom_text(data = site,
  #           aes(x = lon, y = lat),
  #           label    = "CCR",
  #           colour   = "white",
  #           size     = 3,
  #           nudge_x  = 2.5,
  #           nudge_y  = 0.8) +
  
  coord_sf(
    xlim = c(-125, -66),
    ylim = c(24, 50),
    expand = FALSE
  ) +
  
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", colour = "white", linewidth = 0.5),
    plot.margin     = margin(4, 4, 4, 4)
  )

us_map

# ggsave("ccr_us_map.png", us_map,
#        width = 6, height = 4, dpi = 300, bg = "grey10")
