#!/bin/bash
FILES=../data/meta_dump/*
for f in $FILES;
do
  n=${f//[^0-9]/};
  for i in `seq 1 4`;
  do
    # this is the poisson model
    ./degree_phy_model sample random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/meta_dump/meta_sample/samples_${n}_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    # this is the negative binomial model
    ./deg_phy_over sample random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/meta_dump/meta_sample/samples_${n}_${i}.csv &
  done
  wait
done
