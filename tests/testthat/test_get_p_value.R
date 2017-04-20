library(svconfound)

context('get_p_value()')

test_that('get_p_value() fails on empty input', {
  expect_error(get_p_value())
})

test_that('get_p_value() warns on non-numeric input and returns 1', {
  expect_warning(foo <- get_p_value(NA))
  expect_equal(foo, 1)
  expect_warning(foo <- get_p_value(LETTERS[1:7]))
  expect_equal(foo, 1)
})

test_that('get_p_value() fails on invalid input length', {
  expect_error(get_p_value(c(value = 1)))
  expect_error(get_p_value(c(value = 1, numdf = 2)))
})

test_that('get_p_value() fails on invalid input names', {
  expect_error(get_p_value(c(foo = 1, numdf = 2, dendf = 3)))
  expect_error(get_p_value(c(value = 1, foo = 2, dendf = 3)))
  expect_error(get_p_value(c(value = 1, numdf = 2, foo = 3)))
})

test_that('get_p_Value() works as expected on correct inputs', {
  foo = c(value = 1, numdf = 2, dendf = 3)
  bar = c(value = Inf, numdf = 1, dendf = 1)
  foobar =  c(value = 0, numdf = 77, dendf = 42)

  expect_equal(
    foo_result <- get_p_value(foo),
    1 - pf(foo[['value']], foo[['numdf']], foo[['dendf']])
  )
  expect_gte(foo_result, 0)
  expect_lte(foo_result, 1)
  expect_equal(
    bar_result <- get_p_value(bar),
    1 - pf(bar[['value']], bar[['numdf']], bar[['dendf']])
  )
  expect_gte(bar_result, 0)
  expect_lte(bar_result, 1)
  expect_equal(
    foobar_result <- get_p_value(foobar),
    1 - pf(foobar[['value']], foobar[['numdf']], foobar[['dendf']])
  )
  expect_gte(foobar_result, 0)
  expect_lte(foobar_result, 1)
})
