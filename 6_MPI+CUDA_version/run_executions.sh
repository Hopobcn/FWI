#!/bin/bash

source environment_hulk.sh
export OMP_NUM_THREADS=4
#bind threads to real cores
export MP_BIND=no

for i in 1
do
    ###########################################################################################################
    mpirun -np $i ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies_small.txt  2> fwi.$i.small.out
  
    count_global=0;
    total_global=0;
    count_stress=0;
    total_stress=0;
    count_veloci=0;
    total_veloci=0;
    for k in 0  
    do
        mv fwi.0$k.log fwi.$i.0$k.small.log
        cat fwi.$i.0$k.small.log | grep STATS | grep Maingrid | grep GLOBAL   > fwi.$i.0$k.small.global.log
        cat fwi.$i.0$k.small.log | grep STATS | grep Maingrid | grep STRESS   > fwi.$i.0$k.small.stress.log
        cat fwi.$i.0$k.small.log | grep STATS | grep Maingrid | grep VELOCITY > fwi.$i.0$k.small.velocity.log

        
        for j in $( awk '{ print $10; }' fwi.$i.0$k.small.global.log );
        do
            total_global=$(echo $total_global+$j | bc )
     
            ((count_global++))
        done
                
        for j in $( awk '{ print $10; }' fwi.$i.0$k.small.stress.log );
        do
            total_stress=$(echo $total_stress+$j | bc )
            ((count_stress++))
        done
            
        for j in $( awk '{ print $10; }' fwi.$i.0$k.small.velocity.log );
        do
            total_veloci=$(echo $total_veloci+$j | bc )
            ((count_veloci++))
        done

    done
    echo "scale=2; $total_global / $count_global" | bc > fwi.$i.small.stats
    echo "scale=2; $total_stress / $count_stress" | bc >> fwi.$i.small.stats
    echo "scale=2; $total_veloci / $count_veloci" | bc >> fwi.$i.small.stats

    #############################################################################################################
    mpirun -np $i ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies_medium.txt  2> fwi.$i.medium.out
    
    count_global=0;
    total_global=0;
    count_stress=0;
    total_stress=0;
    count_veloci=0;
    total_veloci=0;
    for k in 0  
    do
        mv fwi.0$k.log fwi.$i.0$k.medium.log
        cat fwi.$i.0$k.medium.log | grep STATS | grep Maingrid | grep GLOBAL   > fwi.$i.0$k.medium.global.log
        cat fwi.$i.0$k.medium.log | grep STATS | grep Maingrid | grep STRESS   > fwi.$i.0$k.medium.stress.log
        cat fwi.$i.0$k.medium.log | grep STATS | grep Maingrid | grep VELOCITY > fwi.$i.0$k.medium.velocity.log

        
        for j in $( awk '{ print $10; }' fwi.$i.0$k.medium.global.log );
        do
            total_global=$(echo $total_global+$j | bc )
     
            ((count_global++))
        done
                
        for j in $( awk '{ print $10; }' fwi.$i.0$k.medium.stress.log );
        do
            total_stress=$(echo $total_stress+$j | bc )
            ((count_stress++))
        done
            
        for j in $( awk '{ print $10; }' fwi.$i.0$k.medium.velocity.log );
        do
            total_veloci=$(echo $total_veloci+$j | bc )
            ((count_veloci++))
        done

    done
    echo "scale=2; $total_global / $count_global" | bc > fwi.$i.medium.stats
    echo "scale=2; $total_stress / $count_stress" | bc >> fwi.$i.medium.stats
    echo "scale=2; $total_veloci / $count_veloci" | bc >> fwi.$i.medium.stats

    
    mpirun -np $i ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies_large.txt  2> fwi.$i.large.out 

    count_global=0;
    total_global=0;
    count_stress=0;
    total_stress=0;
    count_veloci=0;
    total_veloci=0;
    for k in 0 
    do
        mv fwi.0$k.log fwi.$i.0$k.large.log
        cat fwi.$i.0$k.large.log | grep STATS | grep Maingrid | grep GLOBAL   > fwi.$i.0$k.large.global.log
        cat fwi.$i.0$k.large.log | grep STATS | grep Maingrid | grep STRESS   > fwi.$i.0$k.large.stress.log
        cat fwi.$i.0$k.large.log | grep STATS | grep Maingrid | grep VELOCITY > fwi.$i.0$k.large.velocity.log

        
        for j in $( awk '{ print $10; }' fwi.$i.0$k.large.global.log );
        do
            total_global=$(echo $total_global+$j | bc )
     
            ((count_global++))
        done
                
        for j in $( awk '{ print $10; }' fwi.$i.0$k.large.stress.log );
        do
            total_stress=$(echo $total_stress+$j | bc )
            ((count_stress++))
        done
            
        for j in $( awk '{ print $10; }' fwi.$i.0$k.large.velocity.log );
        do
            total_veloci=$(echo $total_veloci+$j | bc )
            ((count_veloci++))
        done

    done
    echo "scale=2; $total_global / $count_global" | bc >  fwi.$i.large.stats
    echo "scale=2; $total_stress / $count_stress" | bc >> fwi.$i.large.stats
    echo "scale=2; $total_veloci / $count_veloci" | bc >> fwi.$i.large.stats


done

