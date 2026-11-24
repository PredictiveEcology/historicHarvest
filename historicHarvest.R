## Everything in this file and any files in the R directory are sourced during `simInit()`;
## all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used inside any function that was sourced in this module;
## they are namespaced to the module, just like functions in R packages.
## If exact location is required, functions will be: `sim$.mods$<moduleName>$FunctionName`.

defineModule(sim, list(
  name = "historicHarvest",
  description = paste(
    "Replays the forest harvest that actually happened (Canada Landsat Disturbance, CanLaD,",
    "harvest class, 1985 onwards) in a LandR (Biomass_core) simulation, and replants each cut",
    "pixel with the species SCANFI shows growing there at the next 5-year snapshot.",
    "The module only says where harvest happened each year (`rstCurrentHarvest`) and provides",
    "`plantHistoricHarvest()`; the removal of cohorts is done by `Biomass_genericDisturbance`",
    "calling `disturbanceMgmt::processPartialDist()`."
  ),
  keywords = c("harvest", "CanLaD", "SCANFI", "LandR", "disturbance"),
  authors = c(
    person("Eliot", "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca",
           role = c("aut", "cre"))
  ),
  childModules = character(0),
  version = list(historicHarvest = "0.0.0.9002"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "historicHarvest.Rmd"),
  reqdPkgs = list("PredictiveEcology/SpaDES.core@development (>= 3.2.1.9053)",
                  "data.table", "terra",
                  "PredictiveEcology/reproducible@development",
                  "PredictiveEcology/LandR@modsForFireSense (>= 1.2.0.9066)",
                  "PredictiveEcology/pemisc@development",
                  ## the removal step the Rmd tells users to pair this module with
                  "camillegiuliano/disturbanceMgmt@main (>= 0.0.1.9001)"),
  parameters = bindrows(
    defineParameter("canladUrl", "character",
                    paste0("https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                           "canlad_including_insect_defoliation/v1.1/TS_FinalFilter_20260508"),
                    NA, NA,
                    paste("Folder holding CanLaD's annual disturbance rasters (Canada Landsat Disturbance",
                          "including forest insect pest, v1.1, 1985-2025;",
                          "<https://open.canada.ca/data/en/dataset/902801fd-4d9d-4df4-9e95-319e429545cc>).",
                          "`NA` means no download: the file must already be in `inputPath(sim)`.")),
    defineParameter("canladFile", "character", "canlad_annual_%d_v1_1_20260508.tif", NA, NA,
                    "Name of one year's raster in `canladUrl`, with `%d` where the year goes."),
    defineParameter("harvestClasses", "integer", c(2L, 5L), 1L, 8L,
                    paste("CanLaD disturbance classes treated as harvest: 2 (harvest) and 5 (defoliation",
                          "followed by harvest, e.g. salvage logging); both remove the stand.")),
    defineParameter("harvestMinFraction", "numeric", 0.5, 0, 1,
                    paste("CanLaD is 30 m; the simulation's pixels are bigger. A simulation pixel is cut",
                          "when at least this fraction of its CanLaD pixels is harvest.")),
    defineParameter("minCoverThreshold", "numeric", 5, 0, 100,
                    paste("A species is planted in a cut pixel only if its SCANFI cover there is above",
                          "this percentage. Same meaning and default as in `Biomass_borealDataPrep`.")),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    "The column in `sim$sppEquiv` that names species; as in `Biomass_borealDataPrep`."),
    defineParameter("initialB", "numeric", 10, 0, NA,
                    "Biomass (g/m2) of a newly planted cohort. LandR's default is 10."),
    defineParameter(".runInitialTime", "numeric", start(sim), NA, NA,
                    paste("Time of the first `harvest` event. Default: `start(sim)`. A run whose start-year vegetation",
                          "already shows that year's cuts (SCANFI) sets the year after.")),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Used in the names of the cached SCANFI files."),
    defineParameter(".useCache", "character", ".inputObjects", NA, NA,
                    "Which events or functions to Cache.")
  ),
  inputObjects = bindrows(
    expectsInput("rasterToMatch", "SpatRaster",
                 "The simulation's raster. `rstCurrentHarvest` is on its grid.", sourceURL = NA),
    expectsInput("sppEquiv", "data.table",
                 "Species equivalence table, as used by `Biomass_borealDataPrep`.", sourceURL = NA),
    expectsInput("cohortData", "data.table",
                 "Cohorts (from `Biomass_core`). Read for its species and, when planting, updated.",
                 sourceURL = NA),
    expectsInput("pixelGroupMap", "SpatRaster",
                 "Pixel groups (from `Biomass_core`); updated when planting.", sourceURL = NA),
    expectsInput("ecoregionMap", "SpatRaster",
                 "Categorical ecoregion groups (from `Biomass_core`); gives planted cohorts their `ecoregionGroup`.",
                 sourceURL = NA),
    expectsInput("speciesEcoregion", "data.table",
                 "Species by ecoregion traits (from `Biomass_core`); a species is planted only where it has `maxB`.",
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput("rstCurrentHarvest", "SpatRaster",
                  paste("1 where harvest happened this year, 0 elsewhere, `NA` outside the simulation.",
                        "`NULL` in a year without harvest, so a `triggerObj` of it skips that year.")),
    createsOutput("partialDistLoc", "SpatRaster",
                  paste("The same raster, under the name `disturbanceMgmt::processPartialDist()` reads",
                        "for 'where'.")),
    createsOutput("partialDistTable", "data.table",
                  paste("Every species in `cohortData`, any age, with this year as `distYear`: what",
                        "`processPartialDist()` reads to remove all cohorts in the cut pixels."))
  )
))

doEvent.historicHarvest <- function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      ## Before Biomass_regeneration (3), the fire modules (5.13) and growth, so a year's cut is
      ## in place before that year's fire and growth.
      sim <- scheduleEvent(sim, P(sim)$.runInitialTime, "historicHarvest", "harvest", eventPriority = 2)
    },
    harvest = {
      sim <- Harvest(sim)
      if (time(sim) < end(sim)) {
        sim <- scheduleEvent(sim, time(sim) + 1, "historicHarvest", "harvest", eventPriority = 2)
      }
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

Init <- function(sim) {
  ## the same parameter as in Biomass_borealDataPrep must mean the same thing here
  paramCheckOtherMods(sim, "sppEquivCol", ifSetButDifferent = "error")
  paramCheckOtherMods(sim, "minCoverThreshold", ifSetButDifferent = "warning")
  invisible(sim)
}

## One year's harvest: where it happened, and what `processPartialDist()` needs to remove it
Harvest <- function(sim) {
  yr <- as.integer(time(sim))
  cut <- if (!yr %in% canladYears) NULL else getCanLaDHarvest(
    year = yr, rasterToMatch = sim$rasterToMatch,
    url = if (is.na(P(sim)$canladUrl)) NULL else P(sim)$canladUrl,
    file = P(sim)$canladFile, destinationPath = inputPath(sim),
    harvestClasses = P(sim)$harvestClasses, minFraction = P(sim)$harvestMinFraction
  ) |>
    Cache(.functionName = "getCanLaDHarvest",
          userTags = c(currentModule(sim), "canladHarvest", yr, P(sim)$.studyAreaName),
          omitArgs = c("userTags", "destinationPath"))
  if (is.null(cut)) message("historicHarvest: CanLaD has no data for ", yr, " (", min(canladYears), "-", max(canladYears), ")")

  nCut <- if (is.null(cut)) 0L else sum(as.vector(terra::values(cut)) == 1, na.rm = TRUE)
  if (nCut == 0L) {
    ## NULL, not an empty raster: Biomass_genericDisturbance skips a disturbance whose
    ## triggerObj is NULL, and processPartialDist() treats an all-NA raster as "everything"
    sim$rstCurrentHarvest <- sim$partialDistLoc <- sim$partialDistTable <- NULL
    message("historicHarvest: no harvest in ", yr)
    return(invisible(sim))
  }

  sim$rstCurrentHarvest <- sim$partialDistLoc <- cut
  sim$partialDistTable <- data.table::data.table(
    speciesCode = unique(sim$cohortData$speciesCode), distYear = yr,
    ageMin = 0L, ageMax = .Machine$integer.max
  )
  message("historicHarvest: ", nCut, " pixels harvested in ", yr)
  invisible(sim)
}

.inputObjects <- function(sim) {
  ## Everything this module reads comes from the LandR modules; nothing to make here.
  return(invisible(sim))
}
