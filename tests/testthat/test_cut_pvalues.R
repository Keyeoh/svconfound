library(svconfound)

context('Cut p-values')

test_that('cut_pvalues() fails on empty input', {
  expect_error(cut_pvalues())
})

test_that('cut_pvalues() fails on non-numeric input', {
  expect_error(cut_pvalues(NA))
  expect_error(cut_pvalues(LETTERS[1:10]))
  expect_error(cut_pvalues(c(1, 2, 3, 4, 'a', 'b')))
})

test_that('cut_pvalues() output is a factor', {
  expect_is(cut_pvalues(c(0, 0.001, 0.01, 0.1, 1)), 'factor')
})

cut_levels = c('< 1e-10', '< 1e-5', '< 0.001', '< 0.01', '< 0.05', 'NS')

test_that('cut_pvalues() output has all levels', {
  expect_equal(levels(cut_pvalues(c(0, 0.001, 0.01, 0.1, 1))),
               cut_levels)
  expect_equal(levels(cut_pvalues(c(0))),
               cut_levels)
})

test_that('cut_pvalues() works with values outside [0, 1]', {
  expect_equal(cut_pvalues(c(-1, 0.001, 0.01, 0.1, 1)),
               factor(c('< 1e-10', '< 0.001', '< 0.01', 'NS', 'NS'),
                      levels = cut_levels))
  expect_equal(cut_pvalues(c(0, 0.001, 0.01, 3, 9999)),
               factor(c('< 1e-10', '< 0.001', '< 0.01', 'NS', 'NS'),
                      levels = cut_levels))
  expect_equal(cut_pvalues(c(-Inf, 0.001, 0.01, 0.1, 1)),
               factor(c('< 1e-10', '< 0.001', '< 0.01', 'NS', 'NS'),
                      levels = cut_levels))
  expect_equal(cut_pvalues(c(0, 0.001, 0.01, Inf, 1)),
               factor(c('< 1e-10', '< 0.001', '< 0.01', 'NS', 'NS'),
                      levels = cut_levels))
})

test_that('cut_pvalues() works as expected on correct inputs', {
  expect_equal(cut_pvalues(c(0, 0.001, 0.01, 0.1, 1)),
               factor(c('< 1e-10', '< 0.001', '< 0.01', 'NS', 'NS'),
                      levels = cut_levels))

})
