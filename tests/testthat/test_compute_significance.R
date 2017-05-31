library(svconfound)

context('compute_significance_data()')

# Describe a mock correct case =================================================
mock_values = structure(
  c(
    -0.420947069874835,
    -0.370772777053382,
    -0.147115870554044,
    0.145384379514911,
    0.369859006417622,
    0.421268916705144,
    0.275562690359924,
    0.000927802204756354,
    -0.274146197132857,
    -0.420947069874835,
    -0.0679865273565225,
    0.247064210147407,
    0.446513024945398,
    0.437009504146823,
    0.223045778371529,
    -0.095271898857681,
    -0.3690229767903,
    -0.470107722301358,
    -0.351224102492757,
    -0.0679865273565224,
    0.0357705857899453,
    0.186793471724375,
    -0.882600158160668,
    0.196852543551029,
    0.0511904683828425,
    -0.118417707273885,
    -0.232623129317155,
    -0.237985132478951,
    -0.131990549646354,
    0.0357705857670688
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
      0.707931400803344,
      0.00265595941039343,
      0.800158781299549,
      0.0892100848401578,
      0.674352123817491,
      0.288994899251423,
      0.77657891016915,
      0.0491895179536767,
      0.940155427188014
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
      0.0162936036210286,
      0.117185087198138,
      0.0886759430447011,
      0.601508134440589,
      0.191694602051888
    ),
    Sig = structure(
      c(6L, 5L, 6L, 6L, 6L, 6L),
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

