# ==============================================================================
library(svconfound)

context('svd_analysis()')

# Positive tests ===============================================================
mock_points = structure(
  c(
    3,
    2.94092093904347,
    1.505752286013,
    -0.633974596215561,
    -2.47705771903206,
    -3.16109800568339,
    -2.36602540378444,
    -0.463863220011418,
    1.65534571967039,
    3,
    3,
    2.29813332935693,
    0.520944533000791,
    -1.5,
    -2.81907786235772,
    -2.81907786235773,
    -1.5,
    0.52094453300079,
    2.29813332935693,
    3,
    0,
    0.642787609686539,
    0.984807753012208,
    0.866025403784439,
    0.342020143325669,
    -0.342020143325669,
    -0.866025403784438,
    -0.984807753012208,
    -0.64278760968654,
    -2.44929359829471e-16
  ),
  .Dim = c(10L, 3L)
)

mock_pdata = data.frame(
  g = c(rep('a', 5), rep('b', 5)),
  h = c(rep('c', 2), rep('d', 8)),
  z = 1:10
)

test_that('svd_analysis() works on a basic call (scaling)', {
  mock_result = svd_analysis(
    values = mock_points,
    pdata = mock_pdata,
    center = TRUE,
    scale = TRUE
  )
  expect_equal(mock_result[['center']], TRUE)
  expect_equal(mock_result[['scale']], TRUE)
  expect_equal(dim(mock_result[['variance_explained']]), c(3, 2))
  expect_equal(mock_result[['variance_explained']][['Var']], c(0.66, 0.33, 0),
               tolerance = 0.01)
  expect_is(mock_result, 'SVDAnalysis')
  expect_equal(mock_result[['limit_significant_PC']], 1)
  expect_null(mock_result[['significance_control']])
  expect_equal(mock_result[['significance']][['P_Value']],
               c(0.49, 0, 0.44, 0.07, 0.83, 0.87, 0.6, 0.05, 0.9),
               tolerance = 0.01)
  expect_equal(
    mock_result[['significance']][['Sig']],
    factor(
      c('NS', '< 0.01', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS'),
      levels = c('< 1e-10', '< 1e-5', '< 0.001', '< 0.01', '< 0.05', 'NS')
    )
  )
})

mock_dimnames = list(
  c('22711390', '22795447', '56682500', '46651360', '24637490', '33665449',
    '54705438', '49720470', '26725400', '57693375', '15700381', '33635504',
    '43720395', '70664314', '71718498', '30724412', '63642461', '47640365',
    '74666473', '31698466', '13643320', '42790394', '28684356', '26772442',
    '21771417'),
  paste0('S', 1:10)
)
mock_green = matrix(1:250, nrow = 25, ncol = 10,
                    dimnames = mock_dimnames)
mock_red = matrix(250:1, nrow = 25, ncol = 10,
                  dimnames = mock_dimnames)

mock_rgset = RGChannelSet(
  Green = mock_green,
  Red = mock_red,
  annotation = c(
    array = 'IlluminaHumanMethylation450k',
    annotation = 'ilmn12.hg19'
  )
)

test_that('svd_analysis() works on a basic call (with RGset)', {
  mock_result = svd_analysis(
    values = mock_points,
    pdata = mock_pdata,
    center = TRUE,
    scale = TRUE,
    rgset = mock_rgset
  )
  expect_equal(mock_result[['center']], TRUE)
  expect_equal(mock_result[['scale']], TRUE)
  expect_equal(dim(mock_result[['variance_explained']]), c(3, 2))
  expect_equal(mock_result[['variance_explained']][['Var']], c(0.66, 0.33, 0),
               tolerance = 0.01)
  expect_is(mock_result, 'SVDAnalysis')
  expect_equal(mock_result[['limit_significant_PC']], 1)
  expect_equal(mock_result[['significance_control']][['P_Value']],
               c(0.2, 0.6, 0.96, 0.2, 0.52, 0.98, 0.2, 0.47, 0.99, 0.76, 0.09,
                 0.9, 0.75, 0.09, 0.9, 0.74, 0.09, 0.9, 0.7, 0.1, 0.9, 0.68,
                 0.1, 0.9, 0.67, 0.11, 0.91, 0.66, 0.11, 0.91, 0.24, 0.23, 0.95,
                 0.25, 0.22, 0.94, 0.25, 0.22, 0.94, 0.25, 0.21, 0.94, 0.25,
                 0.21, 0.94, 0.65, 0.12, 0.91, 0.63, 0.12, 0.91, 0.24, 0.24,
                 0.95, 0.24, 0.23, 0.95),
               tolerance = 0.01)
  expect_equal(mock_result[['significance']][['P_Value']][1:9],
               c(0.49, 0, 0.44, 0.07, 0.83, 0.87, 0.6, 0.05, 0.9),
               tolerance = 0.01)
  expect_equal(mock_result[['significance']][['P_Value']][10:66],
               c(0.2, 0.6, 0.96, 0.2, 0.52, 0.98, 0.2, 0.47, 0.99, 0.76, 0.09,
                 0.9, 0.75, 0.09, 0.9, 0.74, 0.09, 0.9, 0.7, 0.1, 0.9, 0.68,
                 0.1, 0.9, 0.67, 0.11, 0.91, 0.66, 0.11, 0.91, 0.24, 0.23, 0.95,
                 0.25, 0.22, 0.94, 0.25, 0.22, 0.94, 0.25, 0.21, 0.94, 0.25,
                 0.21, 0.94, 0.65, 0.12, 0.91, 0.63, 0.12, 0.91, 0.24, 0.24,
                 0.95, 0.24, 0.23, 0.95),
               tolerance = 0.01)
  expect_equal(
    mock_result[['significance']][['Sig']],
    factor(
      c('NS', '< 0.01', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS',
        rep('NS', 57)),
      levels = c('< 1e-10', '< 1e-5', '< 0.001', '< 0.01', '< 0.05', 'NS')
    )
  )
})

# Negative tests ===============================================================
test_that('svd_analysis() fails on empty inputs', {
  expect_error(
    svd_analysis(
      pdata = mock_pdata,
      center = TRUE,
      scale = TRUE
    )
  )
  expect_error(
    svd_analysis(
      values = mock_points,
      center = TRUE,
      scale = TRUE
    )
  )
})

test_that('svd_analysis() fails with wrong pheno data dimensions', {
  expect_error(
    svd_analysis(
      values = mock_points,
      pdata = rbind(mock_pdata, mock_pdata),
      center = TRUE,
      scale = TRUE
    )
  )
})
