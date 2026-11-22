test_that("scanfiSnapshotYear is the first 5-year snapshot at or after the year, 2025 at most", {
  expect_identical(scanfiSnapshotYear(c(1985, 1986, 1990, 1991, 2020, 2021, 2024, 2025, 2030)),
                   c(1985L, 1990L, 1990L, 1995L, 2020L, 2025L, 2025L, 2025L, 2025L))
})

test_that("planting table: species above the cover threshold, with maxB in the pixel's ecoregion", {
  ## cell 12 is in ecoregion 1_11 (col 2), cell 18 in 1_12 (col 8)
  tab <- historicHarvestPlantingTable(
    c(12L, 18L), toySpeciesLayers(), toyEcoregionMap(), toySpeciesEcoregion(), currentTime = 2001,
    minCoverThreshold = 5, speciesCodeLevels = toySpp
  )
  data.table::setorderv(tab, c("pixelIndex", "speciesCode"))
  ## Pinu_ban (4%) is below the threshold; Betu_pap has no maxB in 1_12
  expect_identical(as.character(tab$speciesCode), c("Pice_mar", "Betu_pap", "Pice_mar"))
  expect_identical(tab$pixelIndex, c(12L, 12L, 18L))
  expect_identical(as.character(tab$ecoregionGroup), c("1_11", "1_11", "1_12"))
  expect_identical(levels(tab$speciesCode), toySpp)
  expect_true(all(tab$B == 0L))
  ## cover equal to the threshold is not above it
  tab2 <- historicHarvestPlantingTable(12L, toySpeciesLayers(list("12" = c(5, 5, 5))), toyEcoregionMap(),
                                       toySpeciesEcoregion(), 2001, 5, toySpp)
  expect_identical(NROW(tab2), 0L)
})

test_that("plantHistoricHarvest replants the cut pixels with SCANFI's species and leaves the rest alone", {
  sim <- toySimInit(2001)
  cut <- toyRTM(); cut[!is.na(cut)] <- 0; cut[c(12, 13, 18)] <- 1   # 12, 13: pg 1 (1_11); 18: pg 2 (1_12)
  sim$rstCurrentHarvest <- cut
  layers <- toySpeciesLayers(list("13" = c(0, 60, 0)))             # cell 13 is Pinu_ban only
  before <- data.table::copy(sim$cohortData)
  pgmBefore <- as.vector(terra::values(sim$pixelGroupMap))

  sim <- plantHistoricHarvest(sim, speciesLayers = layers)

  pgm <- as.vector(terra::values(sim$pixelGroupMap))
  cd <- sim$cohortData
  ## uncut pixels keep their pixelGroup; the cohorts of those groups are untouched
  uncut <- setdiff(which(!is.na(pgmBefore)), c(12, 13, 18))
  expect_identical(pgm[uncut], pgmBefore[uncut])
  expect_equal(cd[pixelGroup %in% 1:2, names(before), with = FALSE][order(pixelGroup, speciesCode)],
               before[order(pixelGroup, speciesCode)])
  ## planted cohorts: age 2, B = initialB, tracked
  planted <- cd[planted %in% TRUE]
  expect_true(all(planted$age == 2) && all(planted$B == 10))
  species <- function(px) sort(as.character(cd[pixelGroup == pgm[px] & planted %in% TRUE, speciesCode]))
  expect_identical(species(12), c("Betu_pap", "Pice_mar"))      # 1_11: Pice_mar 50, Betu_pap 30
  expect_identical(species(13), "Pinu_ban")
  expect_identical(species(18), "Pice_mar")                     # 1_12: Betu_pap has no maxB
  ## no cohort of the old groups in the cut pixels, and the map and cohortData agree
  expect_false(any(pgm[c(12, 13, 18)] %in% 1:2))
  expect_true(all(pgm[!is.na(pgm) & pgm > 0] %in% cd$pixelGroup))
  expect_true(all(cd$pixelGroup %in% pgm))
})

test_that("a cut pixel with nothing to plant is left without cohorts", {
  sim <- toySimInit(2001)
  cut <- toyRTM(); cut[!is.na(cut)] <- 0; cut[12] <- 1
  sim$rstCurrentHarvest <- cut
  sim <- plantHistoricHarvest(sim, speciesLayers = toySpeciesLayers(list("12" = c(0, 0, 0))))
  expect_equal(as.vector(terra::values(sim$pixelGroupMap))[12], 0)
  expect_false(any(sim$cohortData$planted %in% TRUE))
})

test_that("with Biomass_genericDisturbance's removal step first, a cut pixel holds only what SCANFI shows", {
  skip_if_not_installed("disturbanceMgmt")
  dir <- file.path(tempdir(), "canlad4")
  toyCanLaD(dir, 2001, list("12" = list(2, 1), "18" = list(2, 1)))
  sim <- toySimInit(2001, params = list(canladUrl = NA_character_))
  file.copy(file.path(dir, "canlad_annual_2001_v1_1_20260508.tif"), SpaDES.core::inputPath(sim), overwrite = TRUE)
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
  expect_false(any(sim$cohortData[pixelGroup %in% pgm[c(12, 18)]]$age == 80))  # all cohorts of every age gone
  sim <- plantHistoricHarvest(sim, speciesLayers = toySpeciesLayers())
  pgm <- as.vector(terra::values(sim$pixelGroupMap))
  cd <- sim$cohortData
  expect_true(all(cd[pixelGroup %in% pgm[c(12, 18)]]$age == 2))
  expect_identical(sort(as.character(cd[pixelGroup == pgm[12], speciesCode])), c("Betu_pap", "Pice_mar"))
  ## an uncut pixel still has its 80-year-old cohorts
  expect_true(all(cd[pixelGroup == pgm[11]]$age == 80))
})
