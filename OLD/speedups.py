#!/usr/bin/python
import sys
import csv
import argparse
import numpy as np


def count_nlines(fname):
    freqs = []
    with open(fname) as f:
        for i, l in enumerate(f):
            freqs.append(l)
            pass
    return (i+1, freqs)

def parse_args():
    parser = argparse.ArgumentParser(description='Gather time statistics and produce plots.')
    parser.add_argument('path',   metavar='path', type=str, help='FWI freq params path', default='./SetupParams/fwi_frequencies.txt', nargs='?')
    args = parser.parse_args()
    
    return args.path
   
def read_csv_metrics(fname):
    r = []
    with open(fname, 'rb') as csvfile:
        data = csv.reader(csvfile)
        for row in data:
            r.append( row[1] )
    return r

def max_array_lenght( *args):
    max_len = 0
    for arg in args:
        length = len(arg)
        if length > max_len:
            max_len = length
    return max_len

def calc_speedup(rbase, ropti, fid):
    metric_base = rbase[fid] if fid < len(rbase) else ""
    metric_opti = ropti[fid] if fid < len(ropti) else ""
    
    speedup = float(metric_opti)/float(metric_base) if metric_opti != "" and metric_base != "" else ""

    return (metric_opti, speedup)

def write_csv_metrics_and_speedups(fname, r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, nfreqs, freqs):
    with open(fname, 'w') as csvfile:
        writer = csv.writer(csvfile) 
        
        print nfreqs
        print len(r0)

        for fid in range(0, nfreqs):
            (m0, s0) = calc_speedup(r1, r0, fid)
            (m1, s1) = calc_speedup(r1, r1, fid)
            (m2, s2) = calc_speedup(r1, r2, fid)
            (m3, s3) = calc_speedup(r1, r3, fid)
            (m4, s4) = calc_speedup(r1, r4, fid)
            (m5, s5) = calc_speedup(r1, r5, fid)
            (m6, s6) = calc_speedup(r1, r6, fid)
            (m7, s7) = calc_speedup(r1, r7, fid)
            (m8, s8) = calc_speedup(r1, r8, fid)
            (m9, s9) = calc_speedup(r1, r9, fid)
            (m10, s10) = calc_speedup(r1, r10, fid)

            writer.writerow( (str(freqs[fid]).rstrip(), m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10) )

def main():
    freqfile = parse_args()

    (nfreqs, freqs) = count_nlines(freqfile)

    r0 = read_csv_metrics('0_OpenMP_version/fwi.openmp.csv')

    r1 = read_csv_metrics('2_MPI_version/fwi.mpi.16.csv')
    r2 = read_csv_metrics('2_MPI_version/fwi.mpi.32.csv')
    r3 = read_csv_metrics('2_MPI_version/fwi.mpi.64.csv')
    r4 = [] # r4 = read_csv_metrics('2_MPI_version/fwi.mpi.128.csv')

    r5 = read_csv_metrics('4_OpenACC_version/fwi.openacc.csv')
    r6 = read_csv_metrics('5_CUDA_version/fwi.cuda.csv')

    r7 = read_csv_metrics('6_MPI+CUDA_version/fwi.mpicuda.4.csv')
    r8 = read_csv_metrics('6_MPI+CUDA_version/fwi.mpicuda.8.csv') 
    r9 = read_csv_metrics('6_MPI+CUDA_version/fwi.mpicuda.16.csv')
    r10 = read_csv_metrics('6_MPI+CUDA_version/fwi.mpicuda.32.csv')
 
    write_csv_metrics_and_speedups('fwi.all.csv', r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, nfreqs, freqs)


if __name__ == "__main__":
    main()
