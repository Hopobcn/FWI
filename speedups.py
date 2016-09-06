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

def write_csv_metrics_and_speedups(fname, r0, r2, r4, r5, r6, nfreqs, freqs):
    with open(fname, 'w') as csvfile:
        writer = csv.writer(csvfile) 

        for fid in range(0, nfreqs):
            m0 = r0[fid]
            m2 = r2[fid]
            m4 = r4[fid]
            m5 = r5[fid]
            m6 = r6[fid]
            
            s0 = 1.0
            s2 = float(m2)/float(m0)
            s4 = float(m4)/float(m0)
            s5 = float(m5)/float(m0)
            s6 = float(m6)/float(m0)

            writer.writerow( (str(freqs[fid]).rstrip(), m0, m2, m4, m5, m6, s0, s2, s4, s5, s6) )

def main():
    freqfile = parse_args()

    (nfreqs, freqs) = count_nlines(freqfile)

    r0 = read_csv_metrics('0_OpenMP_version/fwi.openmp.csv')
    r2 = read_csv_metrics('2_MPI_version/fwi.mpi.csv')
    r4 = read_csv_metrics('4_OpenACC_version/fwi.openacc.csv')
    r5 = read_csv_metrics('5_CUDA_version/fwi.cuda.csv')
    r6 = read_csv_metrics('6_MPI+CUDA_version/fwi.mpicuda.csv')

    write_csv_metrics_and_speedups('fwi.all.csv', r0, r2, r4, r5, r6, nfreqs, freqs)


if __name__ == "__main__":
    main()
