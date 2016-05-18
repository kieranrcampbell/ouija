
library(ouija)

context("Basic functionality")


test_that("ouija returns valid object", {
  data(synth_gex)
  
  oui <- ouija(synth_gex, iter = 100)
  
  expect_is(oui, "ouija_fit")
  expect_is(oui$fit, "stanfit")
  expect_equal(synth_gex, oui$Y)
  expect_equal(dim(synth_gex)[1], length(map_pseudotime(oui)))
})