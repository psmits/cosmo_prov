#!/bin/bash

cd R

nohup nice R CMD BATCH --vanilla ../R/cosmo_analysis.r
nohup nice R CMD BATCH --vanilla ../R/trait_analysis.r
nohup nice R CMD BATCH --vanilla ../R/diversity.r
