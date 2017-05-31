library(svconfound)

context('get_var_names()')

test_that('get_var_names() fails on empty input', {
  expect_error(get_var_names())
})

test_that('get_var_names() fails when input is not a data.frame', {
  expect_error(get_var_names(NA))
  expect_error(get_var_names(LETTERS[1:7]))
  expect_error(get_var_names(matrix(1:10, ncol = 2)))
  expect_error(get_var_names(rep(c(TRUE, FALSE), 10)))
  expect_error(get_var_names(1:10))
})

mock_df = data.frame(
  a = 1:10,
  b = LETTERS[1:10],
  c = c(letters[1:5], rep('x', 5)),
  d = 'foo',
  e = rep(c(TRUE, FALSE), 5)
)

mock_output = c(a = 'a', c = 'c', e = 'e')

test_that('get_var_names() works as expected on correct inputs', {
  expect_is(get_var_names(mock_df), 'character')
  expect_equal(get_var_names(mock_df), mock_output)
})

no_valid_df = data.frame(
  b = LETTERS[1:10],
  d = 'foo'
)

test_that('get_var_names() fails when there are no valid var names', {
  expect_error(get_var_names(no_valid_df))
})
