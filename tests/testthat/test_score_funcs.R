library(assocr)
context("Test calculation of measures of association/score functions.")

test_that("attraction case", {
  a <- c(1, 2, 3, 4, 5)
  b <- a + 0.1
  dat <- data.frame(t = c(a,b), m = rep(1:2, each=5))
  expect_equal(calc_score_funcs(dat)$iet.mn, 0.1)
  expect_equal(calc_score_funcs(dat)$iet.md, 0.1)
  expect_equal(calc_score_funcs(dat)$s, -1)
  expect_equal(calc_score_funcs(dat)$m, 1)
})


test_that("repulsion case", {
  a <- seq(0,3,0.1)
  b <- seq(3.5,6.5,0.1)
  dat <- data.frame(t = c(a,b), m = rep(1:2, each=length(a)) )
  expect_equal(calc_score_funcs(dat)$iet.mn, 2)
  expect_equal(calc_score_funcs(dat)$iet.md, 2)
  expect_equal(calc_score_funcs(dat)$s, 1)
  expect_equal(calc_score_funcs(dat)$m, 0)
})


test_that("in-between case", {
  a <- c(0, 1, 1.01, 2)
  b <- c(0.1, 1.9, 3, 3.01)
  dat <- data.frame(t = c(a,b), m = rep(1:2, each=length(a)) )
  expect_equal(calc_score_funcs(dat)$iet.mn, mean(c(.1,.1,1,1.01)))
  expect_equal(calc_score_funcs(dat)$iet.md, median(c(.1,.1,1,1.01)))
  expect_equal(calc_score_funcs(dat)$s, 0)
  expect_equal(calc_score_funcs(dat)$m, 0.5)
})


test_that("W doesn't affect results", {
  a <- seq(0,3,0.1)
  b <- seq(3.5,6.5,0.1)
  dat <- data.frame(t = c(a,b), m = rep(1:2, each=length(a)) )
  expect_equal(calc_score_funcs(dat, W=c(-10,10))$iet.mn, 2)
  expect_equal(calc_score_funcs(dat, W=c(-10,10))$iet.md, 2)
  expect_equal(calc_score_funcs(dat, W=c(-10,10))$s, 1)
  expect_equal(calc_score_funcs(dat, W=c(-10,10))$m, 0)
})

