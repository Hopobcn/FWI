#!/bin/bash

# @ job_name      = fwi.mpi
# @ initialdir    = .
# @ output        = Results/%j.out
# @ error         = Results/%j.err
# @ total_tasks   = 16
# @ cpus_per_task = 4
# @ gpus_per_node = 4
# @ wall_clock_limit = 06:00:00
# @ features = k80
ulimit -c unlimited

NGPUS=16
source environment_mt.sh
export OMP_NUM_THREADS=4

for i in $NGPUS
do
    ###########################################################################################################
    mpirun -np $i ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies_small.txt  2> fwi.$i.small.out
  
    count_global=0;
    total_global=0;
    count_stress=0;
    total_stress=0;
    count_veloci=0;
    total_veloci=0;
    for k in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
    do
        mv fwi.$k.log fwi.$i.$k.small.log
        cat fwi.$i.$k.small.log | grep STATS | grep Maingrid | grep GLOBAL   > fwi.$i.$k.small.global.log
        cat fwi.$i.$k.small.log | grep STATS | grep Maingrid | grep STRESS   > fwi.$i.$k.small.stress.log
        cat fwi.$i.$k.small.log | grep STATS | grep Maingrid | grep VELOCITY > fwi.$i.$k.small.velocity.log

        
        for j in $( awk '{ print $10; }' fwi.$i.$k.small.global.log );
        do
            total_global=$(echo $total_global+$j | bc )
     
            ((count_global++))
        done
                
        for j in $( awk '{ print $10; }' fwi.$i.$k.small.stress.log );
        do
            total_stress=$(echo $total_stress+$j | bc )
            ((count_stress++))
        done
            
        for j in $( awk '{ print $10; }' fwi.$i.$k.small.velocity.log );
        do
            total_veloci=$(echo $total_veloci+$j | bc )
            ((count_veloci++))
        done

    done
    echo "scale=2; (16 * $total_global) / $count_global" | bc > fwi.$i.small.stats
    echo "scale=2; (16 * $total_stress) / $count_stress" | bc >> fwi.$i.small.stats
    echo "scale=2; (16 * $total_veloci) / $count_veloci" | bc >> fwi.$i.small.stats

    #############################################################################################################
    mpirun -np $i ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies_medium.txt  2> fwi.$i.medium.out
    
    count_global=0;
    total_global=0;
    count_stress=0;
    total_stress=0;
    count_veloci=0;
    total_veloci=0;
    for k in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
    do
        mv fwi.$k.log fwi.$i.$k.medium.log
        cat fwi.$i.$k.medium.log | grep STATS | grep Maingrid | grep GLOBAL   > fwi.$i.$k.medium.global.log
        cat fwi.$i.$k.medium.log | grep STATS | grep Maingrid | grep STRESS   > fwi.$i.$k.medium.stress.log
        cat fwi.$i.$k.medium.log | grep STATS | grep Maingrid | grep VELOCITY > fwi.$i.$k.medium.velocity.log

        
        for j in $( awk '{ print $10; }' fwi.$i.$k.medium.global.log );
        do
            total_global=$(echo $total_global+$j | bc )
     
            ((count_global++))
        done
                
        for j in $( awk '{ print $10; }' fwi.$i.$k.medium.stress.log );
        do
            total_stress=$(echo $total_stress+$j | bc )
            ((count_stress++))
        done
            
        for j in $( awk '{ print $10; }' fwi.$i.$k.medium.velocity.log );
        do
            total_veloci=$(echo $total_veloci+$j | bc )
            ((count_veloci++))
        done

    done
    echo "scale=2; (16 * $total_global) / $count_global" | bc > fwi.$i.medium.stats
    echo "scale=2; (16 * $total_stress) / $count_stress" | bc >> fwi.$i.medium.stats
    echo "scale=2; (16 * $total_veloci) / $count_veloci" | bc >> fwi.$i.medium.stats

    
    mpirun -np $i ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies_large.txt  2> fwi.$i.large.out 

    count_global=0;
    total_global=0;
    count_stress=0;
    total_stress=0;
    count_veloci=0;
    total_veloci=0;
    for k in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15
    do
        mv fwi.$k.log fwi.$i.$k.large.log
        cat fwi.$i.$k.large.log | grep STATS | grep Maingrid | grep GLOBAL   > fwi.$i.$k.large.global.log
        cat fwi.$i.$k.large.log | grep STATS | grep Maingrid | grep STRESS   > fwi.$i.$k.large.stress.log
        cat fwi.$i.$k.large.log | grep STATS | grep Maingrid | grep VELOCITY > fwi.$i.$k.large.velocity.log

        
        for j in $( awk '{ print $10; }' fwi.$i.$k.large.global.log );
        do
            total_global=$(echo $total_global+$j | bc )
     
            ((count_global++))
        done
                
        for j in $( awk '{ print $10; }' fwi.$i.$k.large.stress.log );
        do
            total_stress=$(echo $total_stress+$j | bc )
            ((count_stress++))
        done
            
        for j in $( awk '{ print $10; }' fwi.$i.$k.large.velocity.log );
        do
            total_veloci=$(echo $total_veloci+$j | bc )
            ((count_veloci++))
        done

    done
    echo "scale=2; (16 * $total_global) / $count_global" | bc >  fwi.$i.large.stats
    echo "scale=2; (16 * $total_stress) / $count_stress" | bc >> fwi.$i.large.stats
    echo "scale=2; (16 * $total_veloci) / $count_veloci" | bc >> fwi.$i.large.stats


done



