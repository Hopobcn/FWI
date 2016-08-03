#!/bin/bash

# @ job_name      = fwi.openmp
# @ initialdir    = .
# @ output        = Results/%j.out
# @ error         = Results/%j.err
# @ total_tasks   = 1
# @ cpus_per_task = 16
# # @ gpus_per_node = 4
# @ wall_clock_limit = 06:00:00
# @ features = k80
ulimit -c unlimited

source environment_mt.sh
export OMP_NUM_THREADS=16
#bind threads to real cores
export MP_BIND=yes
export MP_BLIST=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15

i=16

./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies_small.txt  2> fwi.$i.small.out

mv fwi.log fwi.$i.small.log
cat fwi.$i.small.log | grep STATS | grep Maingrid | grep GLOBAL   > fwi.$i.small.global.log
cat fwi.$i.small.log | grep STATS | grep Maingrid | grep STRESS   > fwi.$i.small.stress.log
cat fwi.$i.small.log | grep STATS | grep Maingrid | grep VELOCITY > fwi.$i.small.velocity.log

count=0;
total=0;
for j in $( awk '{ print $10; }' fwi.$i.small.global.log );
do
    total=$(echo $total+$j | bc )
    ((count++))
done
echo "scale=2; $total / $count" | bc > fwi.$i.small.stats

count=0;
total=0;
for j in $( awk '{ print $10; }' fwi.$i.small.stress.log );
do
    total=$(echo $total+$j | bc )
    ((count++))
done
echo "scale=2; $total / $count" | bc >> fwi.$i.small.stats

count=0;
total=0;
for j in $( awk '{ print $10; }' fwi.$i.small.velocity.log );
do
    total=$(echo $total+$j | bc )
    ((count++))
done
echo "scale=2; $total / $count" | bc >> fwi.$i.small.stats

#############################################################################################################
./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies_medium.txt 2> fwi.$i.medium.out

mv fwi.log fwi.$i.medium.log
cat fwi.$i.medium.log | grep STATS | grep Maingrid | grep GLOBAL   > fwi.$i.medium.global.log
cat fwi.$i.medium.log | grep STATS | grep Maingrid | grep STRESS   > fwi.$i.medium.stress.log
cat fwi.$i.medium.log | grep STATS | grep Maingrid | grep VELOCITY > fwi.$i.medium.velocity.log

count=0;
total=0;
for j in $( awk '{ print $10; }' fwi.$i.medium.global.log );
do
    total=$(echo $total+$j | bc )
    ((count++))
done
echo "scale=2; $total / $count" | bc > fwi.$i.medium.stats

count=0;
total=0;
for j in $( awk '{ print $10; }' fwi.$i.medium.stress.log );
do
    total=$(echo $total+$j | bc )
    ((count++))
done
echo "scale=2; $total / $count" | bc >> fwi.$i.medium.stats

count=0;
total=0;
for j in $( awk '{ print $10; }' fwi.$i.medium.velocity.log );
do
    total=$(echo $total+$j | bc )
    ((count++))
done
echo "scale=2; $total / $count" | bc >> fwi.$i.medium.stats


./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies_large.txt  2> fwi.$i.large.out 

mv fwi.log fwi.$i.large.log
cat fwi.$i.large.log | grep STATS | grep Maingrid | grep GLOBAL   > fwi.$i.large.global.log
cat fwi.$i.large.log | grep STATS | grep Maingrid | grep STRESS   > fwi.$i.large.stress.log
cat fwi.$i.large.log | grep STATS | grep Maingrid | grep VELOCITY > fwi.$i.large.velocity.log

count=0;
total=0;
for j in $( awk '{ print $10; }' fwi.$i.large.global.log );
do
    total=$(echo $total+$j | bc )
    ((count++))
done
echo "scale=2; $total / $count" | bc > fwi.$i.large.stats

count=0;
total=0;
for j in $( awk '{ print $10; }' fwi.$i.large.stress.log );
do
    total=$(echo $total+$j | bc )
    ((count++))
done
echo "scale=2; $total / $count" | bc >> fwi.$i.large.stats

count=0;
total=0;
for j in $( awk '{ print $10; }' fwi.$i.large.velocity.log );
do
    total=$(echo $total+$j | bc )
    ((count++))
done
echo "scale=2; $total / $count" | bc >> fwi.$i.large.stats


