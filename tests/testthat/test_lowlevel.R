library(assocr)
context("Low-level fucntions that score functions depend on.")

test_that("nonzero_min returns the right value", {
  expect_equal(nonzero_min(c(1,3,5,2,0)), 1)
  expect_equal(nonzero_min(c(1,1,2,3,5)), 1)
  expect_equal(nonzero_min(c(-1,1,2,3,-15)), 1)
  expect_equal(nonzero_min(c(-1.4,1.23,2.99,3.04,-15)), 1.23)
})

test_that("ceiling_dec properly rounds", {
  # expect_equal(sapply(seq(0,10,0.5), ceiling_dec, 0), c(0,rep(1:10, each=2)))
  expect_equal(ceiling_dec(0.4, 0), 1)
  expect_equal(ceiling_dec(0.2, 1), 0.2)
  expect_equal(ceiling_dec(0.01, 2), 0.02)
  expect_equal(ceiling_dec(0.01, 3), 0.01)
  expect_equal(ceiling_dec(0.02, 2), 0.02)
  expect_equal(ceiling_dec(0.02, 3), 0.02)
})

test_that("count_zeros is accurate", {
  expect_equal(count_zeros(1), -1)
  expect_equal(count_zeros(.1), -1)
  expect_equal(count_zeros(.01), 1)
  expect_equal(count_zeros(.001), 2)
  expect_equal(count_zeros(.0001), 3)

})

test_that("proper rounding up for bounds", {
  expect_equal(ceiling_dec(0.55, count_zeros(0.55) + 1), 1)
  expect_equal(ceiling_dec(0.4, count_zeros(0.4) + 1), 1)
  expect_equal(ceiling_dec(0.055, count_zeros(0.055) + 1), 0.06)
  expect_equal(ceiling_dec(0.06, count_zeros(0.06) + 1), 0.06)
  expect_equal(ceiling_dec(0.05, count_zeros(0.05) + 1), 0.06)
})
