test_that("reff works", {
  expect_equal(reff(5,1,5), c(1,1,1,1,1))
  expect_equal(reff(5,1,1), c(1,0,0,0,0))
  expect_equal(reff(5,5,1), c(5,0,0,0,0))
})
