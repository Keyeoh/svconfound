# ==============================================================================
library(svconfound)

context('Plot helper functions')

# Seed for mock data generation ================================================
set.seed(2032)

# Positive tests ===============================================================
mock_svd_result =
  structure(
    list(
      variance_explained = structure(
        list(
          PC = structure(
            1:3,
            .Label = c("PC-1",
                       "PC-2", "PC-3"),
            class = "factor"
          ),
          Var = c(0.666666666666667,
                  0.333333333333333, 4.96000104911751e-31)
        ),
        .Names = c("PC", "Var"),
        row.names = c(NA, -3L),
        class = "data.frame"
      ),
      center = TRUE,
      scale = TRUE,
      method = "lm",
      significance = structure(
        list(
          PC = structure(
            c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L),
            .Label = c("PC-1",
                       "PC-2", "PC-3"),
            class = "factor"
          ),
          Variable = c("g",
                       "g", "g", "h", "h", "h", "z", "z", "z"),
          P_Value = c(
            0.493070429486578,
            0.00460859119584767,
            0.447346007945623,
            0.0791687639436954,
            0.83883730021756,
            0.87922897736897,
            0.605341385184567,
            0.0590357113745286,
            0.904273182414704
          ),
          Sig = structure(
            c(6L,
              4L, 6L, 6L, 6L, 6L, 6L, 6L, 6L),
            .Label = c("< 1e-10",
                       "< 1e-5", "< 0.001", "< 0.01", "< 0.05", "NS"),
            class = "factor"
          )
        ),
        class = "data.frame",
        .Names = c("PC",
                   "Variable", "P_Value", "Sig"),
        row.names = c(NA, -9L)
      ),
      significance_control = NULL,
      limit_significant_PC = 1L
    ),
    .Names = c(
      "variance_explained",
      "center",
      "scale",
      "method",
      "significance",
      "significance_control",
      "limit_significant_PC"
    ),
    class = c("list", "SVDAnalysis")
  )

mock_sva_result =
  structure(
    list(
      num_sv = 1,
      mod0 = structure(
        c(1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1),
        .Dim = c(12L, 1L),
        .Dimnames = list(
          c("1",
            "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"),
          "(Intercept)"
        ),
        assign = 0L
      ),
      mod = structure(
        c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
          0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),
        .Dim = c(12L, 2L),
        .Dimnames = list(
          c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
            "11", "12"),
          c("(Intercept)", "groupb")
        ),
        assign = 0:1,
        contrasts = structure(list(group = "contr.treatment"), .Names = "group")
      ),
      significance = structure(
        list(
          PC = structure(c(1L, 1L), .Label = "SV-1", class = "factor"),
          Variable = c("group", "batch"),
          P_Value = c(0.940861896594521,
                      7.81885667322513e-12),
          Sig = structure(
            c(6L, 1L),
            .Label = c("< 1e-10",
                       "< 1e-5", "< 0.001", "< 0.01", "< 0.05", "NS"),
            class = "factor"
          )
        ),
        class = "data.frame",
        .Names = c("PC",
                   "Variable", "P_Value", "Sig"),
        row.names = c(NA, -2L)
      ),
      significance_control = NULL,
      surrogates = structure(
        c(
          0.285859954342453,
          0.246249290001253,
          0.273140070481041,
          -0.282245560554353,
          -0.243137955517904,
          -0.321520334549073,
          0.311789687202384,
          0.304251859733228,
          0.30386525241921,
          -0.304400071056924,
          -0.259932942704041,
          -0.313919249797276
        ),
        .Dim = c(12L, 1L),
        .Dimnames = list(NULL, "SV_1")
      )
    ),
    .Names = c(
      "num_sv",
      "mod0",
      "mod",
      "significance",
      "significance_control",
      "surrogates"
    ),
    class = c("list",
              "SVAAnalysis")
  )

test_that('screeplot() works on simple SVD analysis', {
  mock_scree = screeplot(mock_svd_result)
  expect_is(mock_scree, 'gg')
  expect_is(mock_scree, 'ggplot')
  expect_is(mock_scree[['data']], 'data.frame')
  expect_equal(dim(mock_scree[['data']]), c(3, 2))
  expect_equal(
    names(mock_scree),
    c(
      "data",
      "layers",
      "scales",
      "mapping",
      "theme",
      "coordinates",
      "facet",
      "plot_env",
      "labels"
    )
  )
  expect_equal(length(mock_scree[['layers']]), 3)
  expect_is(mock_scree[['layers']][[1]][['geom']], 'GeomLine')
  expect_is(mock_scree[['layers']][[2]][['geom']], 'GeomVline')
  expect_is(mock_scree[['layers']][[3]][['geom']], 'GeomPoint')
  expect_equal(
    mock_scree[['labels']],
    list(
      x = 'Component',
      y = 'Variance Explained',
      group = '1',
      xintercept = 'xintercept'
    )
  )
})

test_that('plot() works on simple SVD analysis', {
  mock_plot = plot(mock_svd_result)
  expect_is(mock_plot, 'gg')
  expect_is(mock_plot, 'ggplot')
  expect_is(mock_plot[['data']], 'data.frame')
  expect_equal(dim(mock_plot[['data']]), c(9, 4))
  expect_equal(
    mock_plot[['data']],
    mock_svd_result[['significance']]
  )
  expect_equal(
    names(mock_plot),
    c(
      "data",
      "layers",
      "scales",
      "mapping",
      "theme",
      "coordinates",
      "facet",
      "plot_env",
      "labels"
    )
  )
  expect_equal(length(mock_plot[['layers']]), 2)
  expect_is(mock_plot[['layers']][[1]][['geom']], 'GeomVline')
  expect_is(mock_plot[['layers']][[2]][['geom']], 'GeomTile')
  expect_equal(
    mock_plot[['labels']],
    list(
      x = 'Component',
      y = 'Variable',
      xintercept = 'xintercept',
      fill = 'Sig',
      colour = 'Sig'
    )
  )
})

test_that('plot() works on simple SVA analysis', {
  mock_plot = plot(mock_sva_result)
  expect_is(mock_plot, 'gg')
  expect_is(mock_plot, 'ggplot')
  expect_is(mock_plot[['data']], 'data.frame')
  expect_equal(dim(mock_plot[['data']]), c(2, 4))
  expect_equal(
    mock_plot[['data']],
    mock_sva_result[['significance']]
  )
  expect_equal(
    names(mock_plot),
    c(
      "data",
      "layers",
      "scales",
      "mapping",
      "theme",
      "coordinates",
      "facet",
      "plot_env",
      "labels"
    )
  )
  expect_equal(length(mock_plot[['layers']]), 2)
  expect_is(mock_plot[['layers']][[1]][['geom']], 'GeomVline')
  expect_is(mock_plot[['layers']][[2]][['geom']], 'GeomTile')
  expect_equal(
    mock_plot[['labels']],
    list(
      x = 'Component',
      y = 'Variable',
      xintercept = 'xintercept',
      fill = 'Sig',
      colour = 'Sig'
    )
  )
})

# Negative tests ===============================================================
test_that('screeplot() does not work on a SVA analysis', {
  expect_error(
    screeplot(mock_sva_result),
    regexp = 'cannot be called on a SVAAnalysis object'
  )
})
