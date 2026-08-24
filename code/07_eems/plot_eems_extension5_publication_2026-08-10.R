#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

project_lib <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/.Rlib"
.libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(rEEMSplots)
  library(sf)
})

effective_size <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 4 || var(x) == 0) return(NA_real_)
  rho <- as.numeric(acf(
    x,
    plot = FALSE,
    lag.max = min(n - 1L, floor(10 * log10(n)))
  )$acf)[-1]
  if (!length(rho)) return(n)
  pair_count <- floor(length(rho) / 2)
  pair_sums <- rho[seq_len(pair_count) * 2 - 1] + rho[seq_len(pair_count) * 2]
  first_negative <- which(pair_sums < 0)[1]
  if (!is.na(first_negative)) pair_sums <- pair_sums[seq_len(first_negative - 1)]
  tau <- 1 + 2 * sum(pair_sums)
  n / max(tau, 1)
}

base_dir <- "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
down_dir <- file.path(base_dir, "downstream_analysis_2026-02-23")
run_dir <- file.path(down_dir, "eems_analysis_ess200_extension5_2026-08-05")
fig_dir <- file.path(run_dir, "figures")
table_dir <- file.path(run_dir, "tables")
data_prefix <- file.path(run_dir, "data", "trow_eems")
coords_file <- "/Users/gspellman/Trowbridgii_analyses/Sample_geographic_coordinants.txt"
popmap_file <- file.path(base_dir, "popmap.tsv")
state_shp <- path.expand("~/.local/share/cartopy/shapefiles/natural_earth/cultural/ne_10m_admin_1_states_provinces_lakes.shp")
land_shp <- path.expand("~/.local/share/cartopy/shapefiles/natural_earth/physical/ne_10m_land.shp")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

mcmc_paths <- file.path(
  run_dir,
  "mcmc",
  paste0("chain", 1:4, "_fixed250_extension5")
)
if (!all(file.exists(file.path(mcmc_paths, "eemsrun.txt")))) {
  stop("One or more EEMS chains are incomplete")
}

sample_order <- readLines(paste0(data_prefix, ".order"), warn = FALSE)
sample_order <- sample_order[nzchar(sample_order)]
coord <- read.table(paste0(data_prefix, ".coord"), header = FALSE)
colnames(coord) <- c("Longitude", "Latitude")
if (length(sample_order) != 98L || "MVZ216210" %in% sample_order) {
  stop("The expected 98-sample dataset excluding MVZ216210 was not found")
}

# Reconstruct the terrestrial prediction domain because the historical .outer
# file concatenated polygon rings. This mask affects display only, not posterior
# samples or the fixed EEMS deme graph.
sample_points <- st_as_sf(
  data.frame(Longitude = coord$Longitude, Latitude = coord$Latitude),
  coords = c("Longitude", "Latitude"), crs = 4326
)
sample_points_mercator <- st_transform(sample_points, 3857)
sample_domain <- st_buffer(
  st_convex_hull(st_union(st_geometry(sample_points_mercator))),
  dist = 100000
)
land <- st_make_valid(suppressWarnings(st_read(land_shp, quiet = TRUE)))
land <- st_transform(land, 3857)
land <- suppressWarnings(st_crop(land, st_bbox(st_buffer(sample_domain, 50000))))
habitat <- st_make_valid(suppressWarnings(
  st_intersection(st_union(st_geometry(land)), sample_domain)
))
habitat_parts <- suppressWarnings(st_cast(habitat, "POLYGON"))
inside_counts <- vapply(seq_along(habitat_parts), function(i) {
  sum(lengths(st_intersects(st_geometry(sample_points_mercator), habitat_parts[i])) > 0)
}, integer(1))
habitat_plot <- habitat_parts[which.max(inside_counts)]
habitat_plot <- st_transform(
  st_simplify(habitat_plot, dTolerance = 1000, preserveTopology = TRUE), 4326
)
if (max(inside_counts) != length(sample_order)) {
  stop("Reconstructed terrestrial prediction domain does not contain every sample")
}
outer_xy <- st_coordinates(habitat_plot)
if ("L1" %in% colnames(outer_xy)) {
  outer_xy <- outer_xy[outer_xy[, "L1"] == 1, , drop = FALSE]
}
outer_xy <- round(outer_xy[, c("X", "Y"), drop = FALSE], 6)

# rEEMSplots reads outer.txt only to draw the habitat. Preserve the engine copy
# and substitute the valid display ring in each completed chain directory.
for (p in mcmc_paths) {
  engine_outer <- file.path(p, "outer.txt")
  backup_outer <- file.path(p, "outer_engine_full_resolution.txt")
  if (!file.exists(backup_outer)) file.copy(engine_outer, backup_outer)
  write.table(outer_xy, engine_outer, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

popmap <- read.delim(popmap_file, header = FALSE, col.names = c("Sample", "Population"))
pop_lookup <- setNames(popmap$Population, popmap$Sample)

pop_colors <- c(
  North = "#4169E1",
  North_Coast = "#4CBB17",
  Sierra_1 = "#800080",
  Sierra_2 = "#00BFC4",
  Sierra_3 = "#E41A1C",
  South_Coast = "#000000"
)
sample_colors <- unname(pop_colors[pop_lookup[sample_order]])

admin1 <- suppressWarnings(st_read(state_shp, quiet = TRUE))
admin1 <- st_make_valid(admin1)
admin1 <- admin1[admin1$adm0_a3 %in% c("USA", "CAN"), ]
bbox <- st_bbox(c(
  xmin = min(coord$Longitude) - 1,
  xmax = max(coord$Longitude) + 1,
  ymin = min(coord$Latitude) - 1,
  ymax = max(coord$Latitude) + 1
), crs = st_crs(4326))
admin_crop <- suppressWarnings(st_crop(admin1, bbox))

overlay_map <- function(plot_title) {
  plot(st_geometry(admin_crop), add = TRUE, border = "#444444", col = NA, lwd = 0.65)
  points(
    coord$Longitude,
    coord$Latitude,
    pch = 21,
    bg = sample_colors,
    col = "white",
    cex = 0.72,
    lwd = 0.45
  )
  lon_ticks <- seq(-124, -118, by = 2)
  lat_ticks <- seq(36, 48, by = 2)
  axis(1, at = lon_ticks, labels = paste0(abs(lon_ticks), "\u00b0W"), cex.axis = 0.95)
  axis(2, at = lat_ticks, labels = paste0(lat_ticks, "\u00b0N"), cex.axis = 0.95)
  mtext("Longitude", side = 1, line = 2.4, cex = 1.0)
  mtext("Latitude", side = 2, line = 2.6, cex = 1.0)
  title(main = plot_title, cex.main = 1.15, line = 1.0)
  box(lwd = 0.8)
}

plot_prefix <- file.path(fig_dir, "EEMS_extension5_posterior_combined")

native_plot_marker <- paste0(plot_prefix, "-pilogl01.pdf")
generate_native_diagnostics <- FALSE
if (generate_native_diagnostics && !file.exists(native_plot_marker)) {
  eems.plots(
    mcmcpath = mcmc_paths,
    plotpath = plot_prefix,
    longlat = TRUE,
    plot.width = 7.2,
    plot.height = 9.2,
    out.png = FALSE,
    add.grid = FALSE,
    add.outline = TRUE,
    col.outline = "#222222",
    lwd.outline = 1.2,
    add.demes = FALSE,
    col.demes = "#202020",
    pch.demes = 19,
    min.cex.demes = 0.35,
    max.cex.demes = 0.85,
    add.abline = TRUE,
    add.r.squared = TRUE,
    add.title = FALSE,
    m.plot.xy = overlay_map("Estimated effective migration surface"),
    q.plot.xy = overlay_map("Estimated effective diversity surface")
  )
}

ess_rows <- list()
for (i in seq_along(mcmc_paths)) {
  d <- read.table(file.path(mcmc_paths[i], "mcmcpilogl.txt"), header = FALSE)
  colnames(d) <- c("LogPrior", "LogLikelihood")
  for (parameter in names(d)) {
    ess_rows[[length(ess_rows) + 1]] <- data.frame(
      Chain = i,
      Parameter = parameter,
      ESS = effective_size(d[[parameter]])
    )
  }
}
ess_table <- do.call(rbind, ess_rows)
write.table(
  ess_table,
  file.path(table_dir, "eems_chain_ess.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

combined <- lapply(mcmc_paths, function(p) {
  d <- read.table(file.path(p, "mcmcpilogl.txt"), header = FALSE)
  colnames(d) <- c("LogPrior", "LogLikelihood")
  d
})
ll_means <- vapply(combined, function(x) mean(x$LogLikelihood), numeric(1))
ll_sds <- vapply(combined, function(x) sd(x$LogLikelihood), numeric(1))
diagnostics <- data.frame(
  Chain = seq_along(mcmc_paths),
  Mean_log_likelihood = ll_means,
  SD_log_likelihood = ll_sds,
  Min_trace_ESS = vapply(seq_along(mcmc_paths), function(i) min(ess_table$ESS[ess_table$Chain == i]), numeric(1))
)
write.table(
  diagnostics,
  file.path(table_dir, "eems_chain_diagnostics.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Build clean publication maps directly from the posterior-averaged EEMS grid.
average_surface <- function(is_migration) {
  read_dimns <- getFromNamespace("read.dimns", "rEEMSplots")
  read_voronoi <- getFromNamespace("read.voronoi", "rEEMSplots")
  transform_rates <- getFromNamespace("transform.rates", "rEEMSplots")
  points_to_raster <- getFromNamespace("points_to_raster", "rEEMSplots")

  dimns <- read_dimns(mcmc_paths[1], longlat = TRUE)
  z_sum <- matrix(0, dimns$nxmrks, dimns$nymrks)
  n_draws <- 0
  for (p in mcmc_paths) {
    v <- read_voronoi(p, longlat = TRUE, is.mrates = is_migration, log_scale = TRUE)
    transformed <- transform_rates(
      dimns, v$tiles, v$rates, v$xseed, v$yseed, zero_mean = TRUE
    )
    z_sum <- z_sum + transformed$Zvals
    n_draws <- n_draws + transformed$niters
  }
  points_to_raster(z_sum / n_draws, dimns, list(proj.in = NULL))
}

migration_cache <- file.path(table_dir, "eems_extension5_migration_surface.rds")
diversity_cache <- file.path(table_dir, "eems_extension5_diversity_surface.rds")
if (file.exists(migration_cache) && file.exists(diversity_cache)) {
  migration_surface <- readRDS(migration_cache)
  diversity_surface <- readRDS(diversity_cache)
} else {
  migration_surface <- average_surface(TRUE)
  diversity_surface <- average_surface(FALSE)
  saveRDS(migration_surface, migration_cache)
  saveRDS(diversity_surface, diversity_cache)
}
surface_colors <- colorRampPalette(
  c("#8C2D04", "#FD8D3C", "#FEE8C8", "#FFFFFF", "#C7E9F1", "#41B6C4", "#045A8D")
)(20)

draw_surface_map <- function(surface, output_stem, title_text, scale_label, minimum_limit) {
  surface <- raster::mask(surface, methods::as(habitat_plot, "Spatial"))
  z <- raster::values(surface)
  data_min <- min(z, na.rm = TRUE)
  data_max <- max(z, na.rm = TRUE)
  # Use the complete posterior surface range. Symmetric limits retain zero as
  # the neutral midpoint without clipping or saturating extreme cells.
  z_limit <- max(minimum_limit, abs(c(data_min, data_max)))
  breaks <- seq(-z_limit, z_limit, length.out = length(surface_colors) + 1)
  ticks <- seq(-z_limit, z_limit, length.out = 5)
  ticks[abs(ticks) < .Machine$double.eps^0.5] <- 0
  tick_labels <- formatC(ticks, format = "f", digits = 3)
  map_bbox <- st_bbox(habitat_plot)
  aspect <- 1 / cos(mean(c(map_bbox["ymin"], map_bbox["ymax"])) * pi / 180)
  x_grid <- raster::xFromCol(surface, seq_len(raster::ncol(surface)))
  row_order <- rev(seq_len(raster::nrow(surface)))
  y_grid <- raster::yFromRow(surface, row_order)
  z_grid <- t(matrix(
    raster::values(surface),
    ncol = raster::ncol(surface),
    byrow = TRUE
  )[row_order, , drop = FALSE])

  render <- function(device) {
    if (device == "pdf") {
      pdf(paste0(output_stem, ".pdf"), width = 7.2, height = 9.2, useDingbats = FALSE)
    } else {
      png(
        paste0(output_stem, ".png"),
        width = 7.2,
        height = 9.2,
        units = "in",
        res = 600,
        type = "quartz"
      )
    }
    layout(matrix(c(1, 2), nrow = 1), widths = c(5.0, 1.2))
    par(mar = c(5.2, 5.4, 3.4, 0.7), mgp = c(3.0, 0.8, 0), tcl = -0.25)
    image(
      x_grid,
      y_grid,
      z_grid,
      col = surface_colors,
      breaks = breaks,
      axes = FALSE,
      asp = aspect,
      xlim = c(map_bbox["xmin"], map_bbox["xmax"]),
      ylim = c(map_bbox["ymin"], map_bbox["ymax"]),
      xlab = "",
      ylab = "",
      useRaster = TRUE
    )
    title(main = title_text, cex.main = 1.2, line = 1.0)
    plot(st_geometry(admin_crop), add = TRUE, border = "#4A4A4A", col = NA, lwd = 0.65)
    plot(st_geometry(habitat_plot), add = TRUE, border = "#161616", col = NA, lwd = 1.1)
    points(
      coord$Longitude, coord$Latitude,
      pch = 21, bg = sample_colors, col = "white", cex = 0.78, lwd = 0.45
    )
    lon_ticks <- seq(-124, -118, by = 2)
    lat_ticks <- seq(36, 48, by = 2)
    axis(1, at = lon_ticks, labels = paste0(abs(lon_ticks), "\u00b0W"), cex.axis = 0.95)
    axis(2, at = lat_ticks, labels = paste0(lat_ticks, "\u00b0N"), cex.axis = 0.95, las = 1)
    mtext("Longitude", side = 1, line = 3.0, cex = 1.0)
    mtext("Latitude", side = 2, line = 3.7, cex = 1.0)
    box(lwd = 0.8)

    par(mar = c(5.2, 0.1, 4.8, 3.2))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(-z_limit, z_limit), xaxs = "i", yaxs = "i")
    rect(
      0.05, breaks[-length(breaks)], 0.62, breaks[-1],
      col = surface_colors, border = NA
    )
    axis(
      4, at = ticks, labels = tick_labels, pos = 0.62,
      las = 1, cex.axis = 0.9, lwd = 0.9, lwd.ticks = 0.9
    )
    mtext(scale_label, side = 3, line = 1.0, cex = 1.0)
    dev.off()
  }
  render("pdf")
  render("png")
  data.frame(
    Surface = title_text,
    Data_minimum = data_min,
    Data_maximum = data_max,
    Color_scale_minimum = -z_limit,
    Color_scale_maximum = z_limit
  )
}

migration_scale <- draw_surface_map(
  migration_surface,
  file.path(fig_dir, "EEMS_extension5_effective_migration_publication"),
  "Estimated effective migration surface",
  expression(log[10](m / bar(m))),
  minimum_limit = 1
)
diversity_scale <- draw_surface_map(
  diversity_surface,
  file.path(fig_dir, "EEMS_extension5_effective_diversity_publication"),
  "Estimated effective diversity surface",
  expression(log[10](q / bar(q))),
  minimum_limit = 0.1
)
write.table(
  rbind(migration_scale, diversity_scale),
  file.path(table_dir, "eems_extension5_publication_color_scale_ranges.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

message("EEMS posterior maps and diagnostics written to: ", run_dir)
