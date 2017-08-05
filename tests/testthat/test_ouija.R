library(ouija)
library(testthat)

context("Basic functionality")


test_that("ouija returns valid object", {
  data(example_gex)
  
  #suppressWarnings(oui <- ouija(synth_gex, iter = 200, warn_lp = FALSE))
  
  if(FALSE) { # Change back when finished dev
  oui <- ouija(example_gex[1:20, 1:4], iter = 200, warn_lp = FALSE)
  
  expect_is(oui, "ouija_fit")
  expect_is(oui$fit, "stanfit")
  expect_equal(dim(example_gex), dim(oui$Y))
  expect_equal(dim(example_gex)[1], length(map_pseudotime(oui)))
  }
})