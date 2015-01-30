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
      output file=../data/meta_dump/meta_sample/pois_samp_${n}_${i}.csv &
  done
  wait
  grep lp__ ../data/meta_dump/meta_sample/pois_samp_${n}_1.csv > \
    ../data/meta_dump/meta_sample/combined_pois_${n}.csv
  sed '/^[#1]/d' ../data/meta_dump/meta_sample/pois_samp_{$n}_*.csv >> \
    ../data/meta_dump/meta_sample/combined_pois_${n}.csv
  for i in `seq 1 4`;
  do
    # this is the negative binomial model
    ./deg_phy_over sample num_samples=5000 num_warmup=5000 thin=5 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/meta_dump/meta_sample/nb_samp_${n}_${i}.csv &
  done
  wait
  grep lp__ ../data/meta_dump/meta_sample/nb_samp_${n}_1.csv > \
    ../data/meta_dump/meta_sample/combined_nb_${n}.csv
  sed '/^[#1]/d' ../data/meta_dump/meta_sample/nb_samp_{$n}_*.csv >> \
    ../data/meta_dump/meta_sample/combined_nb_${n}.csv
done
