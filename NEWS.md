# historicHarvest (development version)

* First version. Each year the module writes `sim$rstCurrentHarvest` (1 = cut) from CanLaD's annual
  disturbance rasters (harvest class), with the objects `disturbanceMgmt::processPartialDist()` needs to
  remove every cohort in those pixels (`partialDistLoc`, `partialDistTable`). `plantHistoricHarvest()`
  replants the cut pixels with the species SCANFI shows there at the first 5-year snapshot after the
  harvest, through `LandR::updateCohortDataPostHarvest(trackPlanting = TRUE)`.
