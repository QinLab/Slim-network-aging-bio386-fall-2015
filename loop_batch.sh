#!/bin/bash
for p in 1.0 0.9 0.8
do
   for LOne in 0.001 0.0015
   do
      echo "$p $LOne"
      R --vanilla --slave -f 20151101-net-sim-ginppi.R --args $Lone 0.0002 5 $p 5
   done
done
