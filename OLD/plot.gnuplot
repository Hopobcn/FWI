#!/usr/bin/gnuplot

set title 'FWI Strong Scalability'

set key outside center right vertical box

set yrange [0:2000]
set xrange [2.5:80]

set xlabel 'freq [hz]'
set ylabel 'Mcells/s'

set xtics (2.5, 5.0, 10.0, 15.0, 20.0, 30, 40, 50, 60, 70, 80)
set ytics 50 nomirror

set datafile separator ","


plot "fwi.all.csv" using 1:2 title "OpenMP(16)"       w linesp lt 1, \
     "fwi.all.csv" using 1:3 title "MPI(16)"          w linesp lt 2, \
     "fwi.all.csv" using 1:4 title "MPI(32)"          w linesp lt 3, \
     "fwi.all.csv" using 1:5 title "MPI(64)"          w linesp lt 4, \
     "fwi.all.csv" using 1:7 title "OpenACC(4)"       w linesp lt 6, \
     "fwi.all.csv" using 1:8 title "CUDA(4)"          w linesp lt 7, \
     "fwi.all.csv" using 1:9  title "MPI+CUDA(4)"     w linesp lt 8, \
     "fwi.all.csv" using 1:10 title "MPI+CUDA(8)"     w linesp lt 9, \
     "fwi.all.csv" using 1:11 title "MPI+CUDA(16)"    w linesp lt 10,\
     "fwi.all.csv" using 1:12 title "MPI+CUDA(32)"    w linesp lt 11

pause -1

set title 'FWI Strong Scalability - Speedup'

set key outside center top horizontal box

set yrange [0:30]
set xrange [2.5:80]

set xlabel 'freq [hz]'
set ylabel 'Speedup'

set xtics (2.5, 5.0, 10.0, 15.0, 20.0, 30, 40, 50, 60, 70, 80)
set ytics 1 nomirror


plot "fwi.all.csv" using 1:13 title "OpenMP(16)"   w linesp lt 1, \
     "fwi.all.csv" using 1:14 title "MPI(16)"      w linesp lt 2, \
     "fwi.all.csv" using 1:15 title "MPI(32)"      w linesp lt 3, \
     "fwi.all.csv" using 1:16 title "MPI(64)"      w linesp lt 4, \
     "fwi.all.csv" using 1:18 title "OpenACC(4)"   w linesp lt 6, \
     "fwi.all.csv" using 1:19 title "CUDA(4)"      w linesp lt 7, \
     "fwi.all.csv" using 1:20 title "MPI+CUDA(4)"  w linesp lt 8, \
     "fwi.all.csv" using 1:21 title "MPI+CUDA(8)"  w linesp lt 9, \
     "fwi.all.csv" using 1:22 title "MPI+CUDA(16)" w linesp lt 10, \
     "fwi.all.csv" using 1:23 title "MPI+CUDA(32)" w linesp lt 11


pause -1
