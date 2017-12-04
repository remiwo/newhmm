source('R/hmm.R')

test_that("Proper output", {

  output <- H_average(current_scale,H1)

  expect_that( output, is_a("numeric") )
  expect_that( length(output) == 1, is_true() )
  expect_error( H_average(current_scale,'WUIIIIIIIIIDBBBBBB'))

})
