# svconfound

### A simple SVD/SVA confounding analysis framework

[![Build Status](https://travis-ci.org/keyeoh/svconfound.svg?branch=master)](https://travis-ci.org/keyeoh/svconfound)

This package implements a simple strategy for analyzing confounding variation on
datasets. The strategy is mostly based on Andrew Teschendorff's book chapter 3
in "Computational and Statistical Epigenomics". It provides two functions for
SVD and SVA analysis of the possible confounding variables among the phenotype
data and control variables, and implementations of the screeplot and plot
methods for them.
