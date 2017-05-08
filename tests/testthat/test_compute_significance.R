library(svconfound)

context('compute_significance_data()')

# Describe a mock correct case =================================================
mock_values = structure(
  c(
    -0.420946767931582,
    -0.370776216840104,
    -0.147115353170495,
    0.145382419252616,
    0.369854141961816,
    0.421267001176157,
    0.275564348678966,
    0.00092207487828811,
    -0.274151648005661,
    -0.420946767931582,
    -0.0679852953669871,
    0.247056987459718,
    0.446498560121451,
    0.437018494223607,
    0.223052617958983,
    -0.0952820572023673,
    -0.36903319885662,
    -0.470109605418701,
    -0.351216502919083,
    -0.0679852953669871,
    0.469883046877175,
    0.484702944046869,
    -0.146440939696909,
    0.0981431523827946,
    0.296804972722947,
    0.402430992347165,
    0.3583104554779,
    0.0715412605006976,
    -0.255990052650256,
    0.253212586893879
  ),
  .Dim = c(10L, 3L)
)

mock_pdata = data.frame(
  g = c(rep('a', 5), rep('b', 5)),
  h = c(rep('c', 2), rep('d', 8)),
  z = 1:10
)

mock_c_names = paste0('PC', 1:3)
mock_v_names = c(g = 'g', h = 'h', z = 'z')

mock_true_result_lm = structure(
  list(
    PC = c("PC1", "PC2", "PC3", "PC1", "PC2", "PC3",
           "PC1", "PC2", "PC3"),
    Variable = c("g", "g", "g", "h", "h", "h",
                 "z", "z", "z"),
    P_Value = c(
      0.707932266902524,
      0.0026557834897577,
      0.670562930529479,
      0.0892102803006438,
      0.674355921054121,
      0.0883675042375978,
      0.776582423095358,
      0.0491918582349268,
      0.291630649724286
    ),
    Sig = structure(
      c(6L,
        4L, 6L, 6L, 6L, 6L, 6L, 5L, 6L),
      .Label = c("< 1e-10", "< 1e-5",
                 "< 0.001", "< 0.01", "< 0.05", "NS"),
      class = "factor"
    )
  ),
  class = "data.frame",
  .Names = c("PC",
             "Variable", "P_Value", "Sig"),
  row.names = c(NA, -9L)
)

mock_true_result_kruskal = structure(
  list(
    PC = c("PC1", "PC2", "PC3", "PC1", "PC2", "PC3"),
    Variable = c("g", "g", "g", "h", "h", "h"),
    P_Value = c(
      0.67517361492712,
      0.0119252335930176,
      0.464702099940466,
      0.0886759430447011,
      0.512610438054442,
      0.036713856362704
    ),
    Sig = structure(
      c(6L, 5L, 6L, 6L, 6L, 5L),
      .Label = c("< 1e-10",
                 "< 1e-5", "< 0.001", "< 0.01", "< 0.05", "NS"),
      class = "factor"
    )
  ),
  class = "data.frame",
  .Names = c("PC",
             "Variable", "P_Value", "Sig"),
  row.names = c(NA, -6L)
)

# Positive test ================================================================
test_that('compute_significance_data_var_names() works as expected', {
  expect_is(
    compute_significance_data_var_names(
      values = mock_values,
      pdata = mock_pdata,
      component_names = mock_c_names,
      var_names = mock_v_names,
      method = 'lm'
    ),
    'data.frame'
  )
  expect_equal(
    compute_significance_data_var_names(
      values = mock_values,
      pdata = mock_pdata,
      component_names = mock_c_names,
      var_names = mock_v_names,
      method = 'lm'
    ),
    mock_true_result_lm
  )
  expect_is(
    compute_significance_data_var_names(
      values = mock_values,
      pdata = mock_pdata[, 1:2],
      component_names = mock_c_names,
      var_names = mock_v_names[1:2],
      method = 'kruskal'
    ),
    'data.frame'
  )
  expect_equal(
    compute_significance_data_var_names(
      values = mock_values,
      pdata = mock_pdata[, 1:2],
      component_names = mock_c_names,
      var_names = mock_v_names[1:2],
      method = 'kruskal'
    ),
    mock_true_result_kruskal
  )
})

test_that(
  'compute_significance_data_var_names() works as expected, even with a
  var_names without names', {
  expect_is(
    compute_significance_data_var_names(
      values = mock_values,
      pdata = mock_pdata,
      component_names = mock_c_names,
      var_names = unname(mock_v_names),
      method = 'lm'
    ),
    'data.frame'
  )
  expect_equal(
    compute_significance_data_var_names(
      values = mock_values,
      pdata = mock_pdata,
      component_names = mock_c_names,
      var_names = unname(mock_v_names),
      method = 'lm'
    ),
    mock_true_result_lm
  )
})

test_that('compute_significance_data() works as expected', {
  expect_is(
    compute_significance_data(
      values = mock_values,
      pdata = mock_pdata,
      component_names = mock_c_names,
      method = 'lm'
    ),
    'data.frame'
  )
  expect_equal(
    compute_significance_data(
      values = mock_values,
      pdata = mock_pdata,
      component_names = mock_c_names,
      method = 'lm'
    ),
    mock_true_result_lm
  )
  expect_is(
    compute_significance_data(
      values = mock_values,
      pdata = mock_pdata[, 1:2],
      component_names = mock_c_names,
      method = 'kruskal'
    ),
    'data.frame'
  )
  expect_equal(
    compute_significance_data(
      values = mock_values,
      pdata = mock_pdata[, 1:2],
      component_names = mock_c_names,
      method = 'kruskal'
    ),
    mock_true_result_kruskal
  )
})

# Negative tests ===============================================================
test_that('compute_significance_data_var_names() fails on empty inputs', {
  expect_error(
    compute_significance_data_var_names(
      values = matrix(nrow = 0, ncol = 0),
      pdata = mock_pdata,
      component_names = mock_c_names,
      var_names = mock_v_names,
      method = 'lm'
    )
  )
  expect_error(
    compute_significance_data_var_names(
      values = mock_values,
      pdata = data.frame(),
      component_names = mock_c_names,
      var_names = mock_v_names,
      method = 'lm'
    )
  )
  expect_error(
    compute_significance_data_var_names(
      values = mock_values,
      pdata = mock_pdata,
      component_names = character(0),
      var_names = mock_v_names,
      method = 'lm'
    )
  )
  expect_error(
    compute_significance_data_var_names(
      values = mock_values,
      pdata = mock_pdata,
      component_names = mock_c_names,
      var_names = character(0),
      method = 'lm'
    )
  )
  expect_error(
    compute_significance_data_var_names(
      values = mock_values,
      pdata = mock_pdata,
      component_names = mock_c_names,
      var_names = mock_v_names,
      method = character(0)
    )
  )
})

test_that(
  'compute_significance_data_var_names() fails on Kruskal-Wallis and numeric
  phenotype',
  {
    expect_error(
      compute_significance_data_var_names(
        values = mock_values,
        pdata = mock_pdata,
        component_names = mock_c_names,
        var_names = mock_v_names,
        method = 'kruskal'
      )
    )
    expect_error(
      compute_significance_data(
        values = mock_values,
        pdata = mock_pdata,
        component_names = mock_c_names,
        method = 'kruskal'
      )
    )
  }
)

