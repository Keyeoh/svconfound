# ==============================================================================
library(minfi)
library(svconfound)

context('get_control_variables()')

# Negative tests ===============================================================
test_that('get_control_variables() fails on empty input', {
  expect_error(get_control_variables())
})

test_that('get_control_variables() fails on wrong input type', {
  expect_error(get_control_variables(matrix(42, nrow = 3, ncol = 2)))
  expect_error(get_control_variables(c(TRUE, FALSE, TRUE)))
  expect_error(get_control_variables(NA))
  expect_error(get_control_variables(NULL))
})

# Positive tests ===============================================================
mock_dimnames = list(
  c('22711390', '22795447', '56682500', '46651360', '24637490', '33665449',
    '54705438', '49720470', '26725400', '57693375', '15700381', '33635504',
    '43720395', '70664314', '71718498', '30724412', '63642461', '47640365',
    '74666473', '31698466', '13643320', '42790394', '28684356', '26772442',
    '21771417'),
  paste0('S', 1:5)
)
mock_green = matrix(2 ^ (1:25), nrow = 25, ncol = 5,
                    dimnames = mock_dimnames)
mock_red = matrix(2 ^ (25:1), nrow = 25, ncol = 5,
                  dimnames = mock_dimnames)

mock_rgset = RGChannelSet(
  Green = mock_green,
  Red = mock_red,
  annotation = 'IlluminaHumanMethylation450k'
)

mock_output = structure(
  c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 19, 19, 19, 19, 19, 18, 18, 18,
    18, 18, 17, 17, 17, 17, 17, 13, 13, 13, 13, 13, 12, 12, 12, 12, 12, 11, 11,
    11, 11, 11, 10, 10, 10, 10, 10, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 23,
    23, 23, 23, 23, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 9, 9, 9, 9, 9, 8, 8,
    8, 8, 8, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20
  ),
  .Dim = c(5L,
           19L),
  .Dimnames = list(
    c("S1", "S2", "S3", "S4", "S5"),
    c(
      "BSC-I C1 Grn",
      "BSC-I C2 Grn",
      "BSC-I C3 Grn",
      "BSC-I C4 Red",
      "BSC-I C5 Red",
      "BSC-I C6 Red",
      "BSC-II C1 Red",
      "BSC-II C2 Red",
      "BSC-II C3 Red",
      "BSC-II C4 Red",
      "Target Removal 1 Grn",
      "Target Removal 2 Grn",
      "Hyb (Low) Grn",
      "Hyb (Medium) Grn",
      "Hyb (High) Grn",
      "Extension (A) Red",
      "Extension (T) Red",
      "Extension (C) Grn",
      "Extension (G) Grn"
    )
  )
)

test_that('get_control_variables() works as expected on correct inputs', {
  expect_equal(get_control_variables(mock_rgset), mock_output)
})
