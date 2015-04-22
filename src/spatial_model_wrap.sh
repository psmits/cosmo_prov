#!/bin/bash
FILES=../data/data_dump/spatial*
for f in $FILES;
do
  n=${f//[^0-9]/};
  for i in `seq 1 4`;
  do
    # this is the poisson model
    ./spatial_model sample num_samples=5000 num_warmup=5000 random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/meta_dump/meta_sample/spatial_samp_${n}_${i}.csv &
  done
  wait
#  for i in `seq 1 4`;
#  do
#    # this is the poisson model
#    ./spatial_zip sample num_samples=5000 num_warmup=5000 random seed=420 \
#      id=$i \
#      data file=$f \
#      output file=../data/meta_dump/meta_sample/zip_samp_${n}_${i}.csv &
#  done
#  wait
#  for i in `seq 1 4`;
#  do
#    # this is the poisson model
#    ./spatial_hurdle sample num_samples=5000 num_warmup=5000 random seed=420 \
#      id=$i \
#      data file=$f \
#      output file=../data/meta_dump/meta_sample/hurdle_samp_${n}_${i}.csv &
#  done
#  wait
done
