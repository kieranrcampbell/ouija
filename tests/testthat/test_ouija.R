library(ouija)
library(testthat)

context("Basic functionality")


test_that("ouija returns valid object", {
  data(example_gex)
  G <- 4
  N <- 50
  small_example_gex <- example_gex[seq_len(N), seq_len(G)] # reduce to small expression matrix

  oui <- ouija(small_example_gex, iter = 500)
  
  expect_is(oui, "ouija_fit")
  expect_is(oui$fit, "stanfit")
  expect_equal(dim(small_example_gex), dim(oui$Y))
  expect_equal(N, length(map_pseudotime(oui)))
  expect_equal(G, length(switch_strengths(oui)))
  expect_equal(G, length(switch_times(oui)))
  
  reg_df <- gene_regulation(oui)
  expect_is(reg_df, "data.frame")
  
  cmat <- consistency_matrix(oui)
  expect_is(cmat, "matrix")
  expect_equal(dim(cmat), c(N, N))
})