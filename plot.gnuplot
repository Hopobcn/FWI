#!/usr/bin/gnuplot

set title 'FWI Strong Scalability'

set key inside center top horizontal box

set yrange [0:*]
set xrange [2.5:22.5]

set xlabel 'freq [hz]'
set ylabel 'Mcells/s'

set xtics (2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5)
set ytics 50 nomirror

set datafile separator ","


#plot "0_OpenMP_version/fwi.openmp.csv" using 1:2 title "OpenMP" w linesp lt 1, "2_MPI_version/fwi.mpi.csv" using 1:2 title "MPI" w linesp lt 3, "4_OpenACC_version/fwi.openacc.csv" using 1:2 title "OpenACC" w linesp lt 4, "5_CUDA_version/fwi.cuda.csv" using 1:2 title "CUDA" w linesp lt 2, "6_MPI+CUDA_version/fwi.mpicuda.csv" using 1:2 title "MPI+CUDA" w linesp lt 5

plot "fwi.all.csv" using 1:2 title "OpenMP" w linesp lt 1, "fwi.all.csv" using 1:3 title "MPI" w linesp lt 3, "fwi.all.csv" using 1:4 title "OpenACC" w linesp lt 4, "fwi.all.csv" using 1:5 title "CUDA" w linesp lt 2, "fwi.all.csv" using 1:6 title "MPI+CUDA" w linesp lt 5

pause -1

set title 'FWI Strong Scalability - Speedup'

set key outside center top horizontal box

set yrange [0:30]
set xrange [2.5:22.5]

set xlabel 'freq [hz]'
set ylabel 'Speedup'

set xtics (2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5)
set ytics 1 nomirror


plot "fwi.all.csv" using 1:7 title "OpenMP" w linesp lt 1, "fwi.all.csv" using 1:8 title "MPI" w linesp lt 3, "fwi.all.csv" using 1:9 title "OpenACC" w linesp lt 4, "fwi.all.csv" using 1:10 title "CUDA" w linesp lt 2, "fwi.all.csv" using 1:11 title "MPI+CUDA" w linesp lt 5


pause -1
