#!/bin/bash
# all taxa models
FILES=../data/data_dump/surv*.data.R
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ./weibull_phy_surv sample num_samples=15000 num_warmup=15000 thin=30 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/wei_surv_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    ./exp_phy_surv sample num_samples=15000 num_warmup=15000 thin=30 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/exp_surv_${i}.csv &
  done
  wait
done
