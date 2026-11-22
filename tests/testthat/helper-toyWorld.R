## A toy LandR world for the tests: nothing is downloaded.
##
## 10 x 10 pixels of 240 m in EPSG:3978. Cell 100 (bottom right) is inactive (NA everywhere).
## Columns 1-5 are ecoregion "1_11" and pixelGroup 1; columns 6-10 are "1_12" and pixelGroup 2.
## Every pixelGroup has Pice_mar and Pinu_ban, age 80. Betu_pap has no maxB in "1_12".
## CanLaD is 30 m, so 8 x 8 CanLaD pixels per pixel.

toySpp <- c("Pice_mar", "Pinu_ban", "Betu_pap")

toyRTM <- function() {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2400, ymin = 0, ymax = 2400,
                   crs = "EPSG:3978", vals = 1)
  r[100] <- NA
  r
}

toyEcoregionMap <- function() {
  r <- toyRTM()
  col <- terra::colFromCell(r, seq_len(terra::ncell(r)))
  r[] <- ifelse(col <= 5, 1L, 2L)
  r[100] <- NA
  levels(r) <- data.frame(ID = 1:2, ecoregionGroup = c("1_11", "1_12"))
  r
}

toyPixelGroupMap <- function() {
  r <- toyRTM()
  col <- terra::colFromCell(r, seq_len(terra::ncell(r)))
  r[] <- ifelse(col <= 5, 1L, 2L)
  r[100] <- NA
  r
}

toyCohortData <- function() {
  data.table::data.table(
    pixelGroup = rep(1:2, each = 2), ecoregionGroup = factor(rep(c("1_11", "1_12"), each = 2)),
    speciesCode = factor(rep(toySpp[1:2], 2), levels = toySpp), age = 80L, B = 1000L,
    mortality = 0L, aNPPAct = 0L
  )
}

toySpeciesEcoregion <- function() {
  se <- data.table::CJ(year = 0, speciesCode = factor(toySpp, levels = toySpp),
                       ecoregionGroup = factor(c("1_11", "1_12")))
  se[, `:=`(maxB = 1000L, maxANPP = 100L)]
  se[speciesCode == "Betu_pap" & ecoregionGroup == "1_12", maxB := 0L]
  se
}

toySppEquiv <- function() data.table::data.table(LandR = toySpp, SCANFI = toySpp)

## SCANFI percent cover. By default Pice_mar 50 everywhere; Pinu_ban 4 (below the threshold of 5)
## and Betu_pap 30 everywhere. `...` are cell = c(Pice_mar, Pinu_ban, Betu_pap) overrides.
toySpeciesLayers <- function(overrides = list()) {
  r <- toyRTM()
  lyrs <- lapply(c(50, 4, 30), function(v) { x <- r; x[!is.na(x)] <- v; x })
  for (cell in names(overrides)) for (i in 1:3) lyrs[[i]][as.integer(cell)] <- overrides[[cell]][i]
  out <- terra::rast(lyrs)
  names(out) <- toySpp
  out
}

## One year's CanLaD raster (30 m) in `dir`, with `cells` = list(cellNumber = list(class, fraction))
## of the toy 240 m pixels: that fraction of the pixel's 30 m pixels (whole rows) gets `class`.
toyCanLaD <- function(dir, year, cells, file = "canlad_annual_%d_v1_1_20260508.tif") {
  rtm <- toyRTM()
  fine <- terra::rast(nrows = 80, ncols = 80, xmin = 0, xmax = 2400, ymin = 0, ymax = 2400,
                      crs = "EPSG:3978", vals = 0L)
  for (cell in as.integer(names(cells))) {
    rc <- terra::rowColFromCell(rtm, cell)
    nRows <- round(cells[[as.character(cell)]][[2]] * 8)
    if (nRows == 0) next
    rows <- (rc[1, 1] - 1) * 8 + seq_len(nRows)
    fine[rows, (rc[1, 2] - 1) * 8 + 1:8] <- cells[[as.character(cell)]][[1]]
  }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  terra::writeRaster(fine, file.path(dir, sprintf(file, year)), overwrite = TRUE, datatype = "INT1U")
  invisible(dir)
}

toyModuleRoot <- function() normalizePath(testthat::test_path("..", ".."), winslash = "/")

toyPaths <- function() {
  root <- file.path(tempdir(), "historicHarvestToy")
  paths <- list(cachePath = file.path(root, "cache"), inputPath = file.path(root, "inputs"),
                modulePath = dirname(toyModuleRoot()), outputPath = file.path(root, "outputs"))
  for (p in paths[c("cachePath", "inputPath", "outputPath")])
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  paths
}

toyObjects <- function() {
  list(rasterToMatch = toyRTM(), sppEquiv = toySppEquiv(), cohortData = toyCohortData(),
       pixelGroupMap = toyPixelGroupMap(), ecoregionMap = toyEcoregionMap(),
       speciesEcoregion = toySpeciesEcoregion())
}

## simInit() only, at `year`, with canladUrl = NA so that nothing is downloaded
toySimInit <- function(year, params = list()) {
  withr::local_options(list(spades.moduleCodeChecks = FALSE, spades.useRequire = FALSE,
                            reproducible.verbose = -2, reproducible.useCloud = FALSE),
                       .local_envir = parent.frame())
  p <- utils::modifyList(list(.useCache = FALSE, canladUrl = NA_character_), params)
  suppressMessages(SpaDES.core::simInit(
    times = list(start = year, end = year), modules = "historicHarvest", objects = toyObjects(),
    params = list(historicHarvest = p), paths = toyPaths()
  ))
}
