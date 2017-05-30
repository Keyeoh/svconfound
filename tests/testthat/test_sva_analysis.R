# ==============================================================================
library(svconfound)

context('sva_analysis()')

# Seed for mock data generation ================================================
set.seed(2032)

# Positive tests ===============================================================
mock_values = matrix(rnorm(1200, mean = 0, sd = 1), nrow = 100, ncol = 12)
mock_values[1:50, 7:12] = mock_values[1:50, 7:12] + 5
mock_values[, 4:6] = mock_values[, 4:6] + 3
mock_values[, 10:12] = mock_values[, 10:12] + 3

mock_pdata = data.frame(
  group = c(rep('a', 6), rep('b', 6)),
  batch = rep(c(rep('x', 3), rep('y', 3)), 2)
)

test_that('sva_analysis() works as expected on a basic call', {
  mock_result = sva_analysis(
    values = mock_values,
    pdata = mock_pdata,
    main_formula = ~ group,
    null_formula = ~ 1
  )
  expect_equal(mock_result[['num_sv']], 1)
  expect_equivalent(mock_result[['mod0']][, 1], rep(1, 12))
  expect_equivalent(mock_result[['mod']][, 1], rep(1, 12))
  expect_equivalent(mock_result[['mod']][, 2], c(rep(0, 6), rep(1, 6)))
  expect_equal(
    mock_result[['significance']][['PC']],
    factor(c('SV-1', 'SV-1'), levels = 'SV-1')
  )
  expect_equal(mock_result[['significance']][['Variable']], c('group', 'batch'))
  expect_equal(
    mock_result[['significance']][['P_Value']],
    c(9.408619e-01, 7.818857e-12),
    tolerance = 0.1
  )
  expect_equal(
    mock_result[['significance']][['Sig']],
    factor(
      c('NS', '< 1e-10'),
      levels = c('< 1e-10', '< 1e-5', '< 0.001', '< 0.01', '< 0.05', 'NS')
    )
  )
})

test_that('sva_analysis() works as expected on a basic call (with vfilter)', {
  mock_result = sva_analysis(
    values = mock_values,
    pdata = mock_pdata,
    main_formula = ~ group,
    null_formula = ~ 1,
    vfilter = 100
  )
  expect_equal(mock_result[['num_sv']], 1)
  expect_equivalent(mock_result[['mod0']][, 1], rep(1, 12))
  expect_equivalent(mock_result[['mod']][, 1], rep(1, 12))
  expect_equivalent(mock_result[['mod']][, 2], c(rep(0, 6), rep(1, 6)))
  expect_equal(
    mock_result[['significance']][['PC']],
    factor(c('SV-1', 'SV-1'), levels = 'SV-1')
  )
  expect_equal(mock_result[['significance']][['Variable']], c('group', 'batch'))
  expect_equal(
    mock_result[['significance']][['P_Value']],
    c(9.408619e-01, 7.818857e-12),
    tolerance = 0.1
  )
  expect_equal(
    mock_result[['significance']][['Sig']],
    factor(
      c('NS', '< 1e-10'),
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
  paste0('S', 1:12)
)
mock_green = matrix(1:300, nrow = 25, ncol = 12,
                    dimnames = mock_dimnames)
mock_red = matrix(300:1, nrow = 25, ncol = 12,
                  dimnames = mock_dimnames)

mock_rgset = RGChannelSet(
  Green = mock_green,
  Red = mock_red,
  annotation = c(
    array = 'IlluminaHumanMethylation450k',
    annotation = 'ilmn12.hg19'
  )
)

test_that('sva_analysis() works as expected on call with RGset', {
  mock_result = sva_analysis(
    values = mock_values,
    pdata = mock_pdata,
    main_formula = ~ group,
    null_formula = ~ 1,
    rgset = mock_rgset
  )
  expect_equal(mock_result[['num_sv']], 1)
  expect_equivalent(mock_result[['mod0']][, 1], rep(1, 12))
  expect_equivalent(mock_result[['mod']][, 1], rep(1, 12))
  expect_equivalent(mock_result[['mod']][, 2], c(rep(0, 6), rep(1, 6)))
  expect_equal(
    mock_result[['significance']][['PC']],
    factor(rep('SV-1', 21), levels = 'SV-1')
  )
  expect_equal(
    mock_result[['significance']][['Variable']],
    factor(
      c('group', 'batch', 'BSC-I C1 Grn', 'BSC-I C2 Grn', 'BSC-I C3 Grn',
        'BSC-I C4 Red', 'BSC-I C5 Red', 'BSC-I C6 Red', 'BSC-II C1 Red',
        'BSC-II C2 Red', 'BSC-II C3 Red', 'BSC-II C4 Red',
        'Target Removal 1 Grn', 'Target Removal 2 Grn', 'Hyb (Low) Grn',
        'Hyb (Medium) Grn', 'Hyb (High) Grn', 'Extension (A) Red',
        'Extension (T) Red', 'Extension (C) Grn', 'Extension (G) Grn'),
      levels = c('group', 'batch', 'BSC-I C1 Grn', 'BSC-I C2 Grn',
                 'BSC-I C3 Grn', 'BSC-I C4 Red', 'BSC-I C5 Red', 'BSC-I C6 Red',
                 'BSC-II C1 Red', 'BSC-II C2 Red', 'BSC-II C3 Red',
                 'BSC-II C4 Red', 'Target Removal 1 Grn',
                 'Target Removal 2 Grn', 'Hyb (Low) Grn', 'Hyb (Medium) Grn',
                 'Hyb (High) Grn', 'Extension (A) Red', 'Extension (T) Red',
                 'Extension (C) Grn', 'Extension (G) Grn')
    )
  )
  expect_equal(
    mock_result[['significance']][['P_Value']],
    c(0.9408619, 7.818857e-12, 0.1842481, 0.1750519, 0.1700909, 0.1381782,
      0.1378948, 0.1376268, 0.1368018, 0.136692, 0.1366418, 0.1366675,
      0.1593007, 0.1594228, 0.1595569, 0.1597011, 0.1598536, 0.1367912,
      0.1370431, 0.1591027, 0.1591931),
    tolerance = 0.1
  )
  expect_equal(
    mock_result[['significance']][['Sig']],
    factor(
      c('NS', '< 1e-10', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS',
        'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS'),
      levels = c('< 1e-10', '< 1e-5', '< 0.001', '< 0.01', '< 0.05', 'NS')
    )
  )
  expect_is(mock_result[['significance_control']], 'data.frame')
})

# Negative tests ===============================================================
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

test_that('sva_analysis() fails when num_sv equals 0', {
  expect_error(
    sva_analysis(
      values = t(mock_points),
      pdata = mock_pdata,
      main_formula = ~ g,
      null_formula = ~ 1
    ),
    regexp = 'Number of estimated surrogate variables is zero.'
  )
})
