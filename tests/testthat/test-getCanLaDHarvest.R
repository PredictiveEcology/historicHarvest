test_that("a pixel is cut when at least minFraction of its CanLaD pixels are harvest", {
  dir <- file.path(tempdir(), "canlad1")
  ## cells 12: all harvest; 13: 3/8 harvest; 14: 5/8 harvest; 15: all fire; 16: all class 5
  toyCanLaD(dir, 2001, list("12" = list(2, 1), "13" = list(2, 3/8), "14" = list(2, 5/8),
                            "15" = list(1, 1), "16" = list(5, 1)))
  cut <- toyGetHarvest(2001, dir)
  expect_s4_class(cut, "SpatRaster")
  expect_true(terra::compareGeom(cut, toyRTM(), stopOnError = FALSE))
  expect_identical(which(as.vector(terra::values(cut)) == 1), c(12L, 14L))
  expect_true(is.na(cut[100][[1]]))
})

test_that("harvestClasses decides what counts as harvest, minFraction how much of a pixel", {
  dir <- file.path(tempdir(), "canlad2")
  toyCanLaD(dir, 2002, list("12" = list(2, 1), "13" = list(2, 3/8), "16" = list(5, 1)))
  cut <- toyGetHarvest(2002, dir, harvestClasses = c(2L, 5L), minFraction = 0.3)
  expect_identical(which(as.vector(terra::values(cut)) == 1), c(12L, 13L, 16L))
})

test_that("a year with no harvest gives an all-zero raster, NA outside the simulation", {
  dir <- file.path(tempdir(), "canlad3")
  toyCanLaD(dir, 2003, list())
  cut <- toyGetHarvest(2003, dir)
  expect_equal(unique(as.vector(terra::values(cut))), c(0, NA))
})
