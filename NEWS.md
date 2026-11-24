# historicHarvest (development version)

* First version. Each year the module writes `sim$rstCurrentHarvest` (1 = cut) from CanLaD's annual
  disturbance rasters (harvest class), with the objects `disturbanceMgmt::processPartialDist()` needs to
  remove every cohort in those pixels (`partialDistLoc`, `partialDistTable`). `timberSupply::plantHistoricHarvest()`
  (moved to the timberSupply package, so a `Biomass_genericDisturbance` specification can name it as a string)
  replants the cut pixels with the species SCANFI shows there at the first 5-year snapshot after the
  harvest, through `LandR::updateCohortDataPostHarvest(trackPlanting = TRUE)`.
