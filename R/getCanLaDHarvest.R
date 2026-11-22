## Years of CanLaD v1.1 annual rasters
canladYears <- 1985L:2025L

#' Where harvest happened in one year, on the simulation's grid
#'
#' Reads one year's CanLaD annual disturbance raster (30 m; see the module's `canladUrl`),
#' crops it to `rasterToMatch` with [reproducible::prepInputs()], keeps the harvest classes,
#' and projects the result onto `rasterToMatch`. A `rasterToMatch` pixel is cut when at least
#' `minFraction` of the CanLaD pixels in it are harvest; sampling one 30 m pixel per
#' simulation pixel instead would mark isolated pixels at random.
#'
#' @param year Integer. The year (1985 to 2025 for CanLaD v1.1).
#' @param rasterToMatch `SpatRaster`. The grid and mask to return.
#' @param url Folder holding the annual rasters, or `NULL` to use a file already in
#'   `destinationPath`.
#' @param file File name of one year's raster, with `%d` where the year goes.
#' @param destinationPath Where the raster is (or is downloaded to).
#' @param harvestClasses Integer. CanLaD classes that count as harvest.
#' @param minFraction Numeric 0-1. See above.
#'
#' @return A `SpatRaster` on `rasterToMatch`: 1 cut, 0 not cut, `NA` where `rasterToMatch` is `NA`.
#'
#' @export
getCanLaDHarvest <- function(year, rasterToMatch, url = NULL, file = "canlad_annual_%d_v1_1_20260508.tif",
                             destinationPath, harvestClasses = 2L, minFraction = 0.5) {
  file <- sprintf(file, as.integer(year))
  canlad <- reproducible::prepInputs(
    targetFile = file, url = if (is.null(url)) NULL else paste(url, file, sep = "/"),
    destinationPath = destinationPath, cropTo = rasterToMatch, projectTo = NA, maskTo = NA,
    writeTo = NULL
  )
  harvest <- terra::classify(canlad, cbind(as.integer(harvestClasses), 1L), others = 0L)
  fraction <- terra::project(harvest, rasterToMatch, method = "average")
  cut <- terra::subst(terra::ifel(fraction >= minFraction, 1L, 0L), NA, 0L) ## no CanLaD data: not cut
  terra::mask(cut, rasterToMatch)
}
