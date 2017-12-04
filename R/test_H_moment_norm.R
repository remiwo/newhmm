source('R/hmm.R')

test_that("Proper output", {

  output <- H_moment_norm(current_scale,H1,100*pi/180)

  expect_that( output, is_a("numeric") )
  expect_that( length(output) < nchar(H1), is_true() )
  expect_error( H_moment_norm(current_scale,'WUIIIIIIIIIDBBBBBB',100*pi/180))
  expect_error( H_moment_norm(current_scale,H1,1000*pi/180))
})
