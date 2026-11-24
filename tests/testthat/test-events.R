test_that("the harvest event writes the objects processPartialDist() reads, for a year with harvest", {
  dir <- file.path(tempdir(), "canlad5")
  toyCanLaD(dir, 2001, list("12" = list(2, 1), "18" = list(2, 1)))
  sim <- toySimInit(2001)
  file.copy(file.path(dir, sprintf(toyDefault("canladFile"), 2001)), SpaDES.core::inputPath(sim), overwrite = TRUE)
  sim <- suppressMessages(SpaDES.core::spades(sim))
  expect_s4_class(sim$rstCurrentHarvest, "SpatRaster")
  expect_identical(which(as.vector(terra::values(sim$rstCurrentHarvest)) == 1), c(12L, 18L))
  expect_identical(sim$partialDistLoc, sim$rstCurrentHarvest)
  tab <- sim$partialDistTable
  expect_identical(sort(as.character(tab$speciesCode)), c("Pice_mar", "Pinu_ban"))
  expect_true(all(tab$distYear == 2001))
  ## every age: processPartialDist()'s window
  expect_true(all(tab$ageMin == 0L) && all(tab$ageMax >= 500L))
})

test_that("a year without harvest leaves rstCurrentHarvest NULL, so a triggerObj skips it", {
  dir <- file.path(tempdir(), "canlad6")
  toyCanLaD(dir, 2002, list())
  sim <- toySimInit(2002)
  file.copy(file.path(dir, sprintf(toyDefault("canladFile"), 2002)), SpaDES.core::inputPath(sim), overwrite = TRUE)
  sim <- suppressMessages(SpaDES.core::spades(sim))
  expect_null(sim$rstCurrentHarvest)
  expect_null(sim$partialDistTable)
})

test_that("the module does not touch rstCurrentBurn, so cut pixels get no fire regeneration", {
  md <- SpaDES.core::moduleMetadata(module = "historicHarvest", path = dirname(toyModuleRoot()))
  outputs <- stats::setNames(md$outputObjects$objectClass, md$outputObjects$objectName)
  expect_identical(outputs[order(names(outputs))],
                   c(partialDistLoc = "SpatRaster", partialDistTable = "data.table",
                     rstCurrentHarvest = "SpatRaster"))
  expect_false("rstCurrentBurn" %in% c(md$inputObjects$objectName, md$outputObjects$objectName))
})

test_that("metadata: timeunit year, the CanLaD defaults", {
  md <- SpaDES.core::moduleMetadata(module = "historicHarvest", path = dirname(toyModuleRoot()))
  expect_identical(md$name, "historicHarvest")
  expect_identical(md$timeunit, "year")
  inputs <- md$inputObjects$objectName
  expect_true(all(c("rasterToMatch", "sppEquiv", "cohortData", "pixelGroupMap", "ecoregionMap",
                    "speciesEcoregion") %in% inputs))
  pars <- stats::setNames(md$parameters$paramClass, md$parameters$paramName)
  expect_true(all(c("canladUrl", "canladFile", "harvestClasses", "harvestMinFraction",
                    "minCoverThreshold", "sppEquivCol", "initialB") %in% names(pars)))
})

test_that(".runInitialTime delays the first harvest event, so the start year is not harvested", {
  dir <- file.path(tempdir(), "canlad7")
  toyCanLaD(dir, 2004, list("12" = list(2, 1)))
  withr::local_options(list(spades.moduleCodeChecks = FALSE, spades.useRequire = FALSE,
                            reproducible.verbose = -2, reproducible.useCloud = FALSE))
  sim <- suppressMessages(SpaDES.core::simInit(
    times = list(start = 2003, end = 2004), modules = "historicHarvest", objects = toyObjects(),
    params = list(historicHarvest = list(.useCache = FALSE, canladUrl = NA_character_, .runInitialTime = 2004)),
    paths = toyPaths()))
  file.copy(file.path(dir, sprintf(toyDefault("canladFile"), 2004)), SpaDES.core::inputPath(sim), overwrite = TRUE)
  sim <- suppressMessages(SpaDES.core::spades(sim))
  done <- SpaDES.core::completed(sim)
  expect_identical(as.numeric(done[eventType == "harvest", eventTime]), 2004)
  expect_identical(which(as.vector(terra::values(sim$rstCurrentHarvest)) == 1), 12L)
})
