# svconfound

### A simple SVD/SVA confounding analysis framework

[![Build Status](https://travis-ci.org/Keyeoh/svconfound.svg?branch=master)](https://travis-ci.org/keyeoh/svconfound)
[![codecov](https://codecov.io/gh/Keyeoh/svconfound/branch/master/graph/badge.svg?token=EY9U9yrR0S)](https://codecov.io/gh/Keyeoh/svconfound)

This package implements a simple strategy for analyzing confounding variation on
datasets. The strategy is mostly based on Andrew Teschendorff's book chapter 3
in "Computational and Statistical Epigenomics". It provides two functions for
SVD and SVA analysis of the possible confounding variables among the phenotype
data and control variables, and implementations of the screeplot and plot
methods for them.
