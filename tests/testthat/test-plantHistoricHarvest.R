test_that("with Biomass_genericDisturbance's removal step first, a cut pixel holds only what SCANFI shows", {
  skip_if_not_installed("disturbanceMgmt")
  skip_if_not_installed("timberSupply")
  dir <- file.path(tempdir(), "canlad4")
  toyCanLaD(dir, 2001, list("12" = list(2, 1), "18" = list(2, 1)))
  sim <- toySimInit(2001, params = list(canladUrl = NA_character_))
  file.copy(file.path(dir, sprintf(toyDefault("canladFile"), 2001)), SpaDES.core::inputPath(sim), overwrite = TRUE)
  sim <- suppressMessages(SpaDES.core::spades(sim))
  expect_identical(which(as.vector(terra::values(sim$rstCurrentHarvest)) == 1), c(12L, 18L))
  ## the removal call as in the module's Rmd, with the arguments Biomass_genericDisturbance finds by name
  removed <- disturbanceMgmt::processPartialDist(
    cohortData = sim$cohortData, partialDistTable = sim$partialDistTable,
    pixelGroupMap = sim$pixelGroupMap, replanting = FALSE, currentTime = time(sim),
    partialDistLoc = sim$partialDistLoc
  )
  sim$cohortData <- removed$cohortData
  sim$pixelGroupMap <- removed$pixelGroupMap
  pgm <- as.vector(terra::values(sim$pixelGroupMap))
  expect_false(any(sim$cohortData$age[sim$cohortData$pixelGroup %in% pgm[c(12, 18)]] == 80))  # all cohorts of every age gone
  sim <- timberSupply::plantHistoricHarvest(sim, speciesLayers = toySpeciesLayers())
  pgm <- as.vector(terra::values(sim$pixelGroupMap))
  cd <- sim$cohortData
  expect_true(all(cd$age[cd$pixelGroup %in% pgm[c(12, 18)]] == 2))
  expect_identical(sort(as.character(cd$speciesCode[cd$pixelGroup == pgm[12]])), c("Betu_pap", "Pice_mar"))
  ## an uncut pixel still has its 80-year-old cohorts
  expect_true(all(cd$age[cd$pixelGroup == pgm[11]] == 80))
})
