#!/bin/bash

cd R

nohup nice R CMD BATCH --vanilla ../R/cosmo_analysis.r
