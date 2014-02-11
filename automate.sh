#!/bin/bash

cd ./R

nohup nice R CMD BATCH --vanilla ../R/time_plots.r
nohup nice R CMD BATCH --vanilla ../R/surv_analysis.r

cd ..

echo 'analysis complete' | mail -s 'mam' psmits@uchicago.edu
