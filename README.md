# MR_simulations

This repository contains code from singly and multiply robust estimators for making causal survival analysis (see "Multiply robust estimators of causal effects for survival outcomes" by Wen, Hernan and Robins). It contain codes to generate data for the dynamic treatment intervention described in the main paper. In this folders, the files to reproduce the results found in the main paper (and in the supplementary materials) include

datagen.R: contains the function to generate data sets
deterministic_datagen_wide_true.R: code to produce the true parameter estimates
program.R: codes to produce the parameter estimates and standard errors from each of the estimators described in the paper
