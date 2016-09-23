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
    parser.add_argument('nranks', metavar='n',    type=int, help='# MPI ranks launched (eq. to #GPUs)')
    parser.add_argument('path',   metavar='path', type=str, help='FWI freq params path', default='../SetupParams/fwi_frequencies.txt', nargs='?')
    args = parser.parse_args()
    
    return (args.nranks, args.path)
    
def main():
    (numprocs, freqfile) = parse_args()

    (nfreqs, freqs) = count_nlines(freqfile)

    print 'num freqs: '+str(nfreqs)
    print 'p  f  type   time [s] Mcells/s'

    nshots = 1
    ngrads = 1

    nsteps = ngrads * nshots * 2 # RTM algorithm executes a FORWARD + BACKWARD propagation
    ntests = 0 #nshots
    ntotal = nsteps + ntests
    print 'num total: ', ntotal

    rt = [[] for y in range(nfreqs)] # temporal array for storing parsed data
    rm = [[0 for x in range(ntotal)] for y in range(nfreqs)]
    for i in range(0, numprocs):
        procid = '0'+str(i) if i < 10 else str(i)
        fname  = 'fwi.'+procid+'.log'



        
        j = 0
        f = 0    
        for line in open(fname, 'r'):
            if ('STATS' in line) and ('GLOBAL' in line):
                if j == ntotal: 
                    j = 0
                    f += 1
                
                istr   = '0'+str(i) if i < 10 else str(i)
                jstr   = '0'+str(j) if j < 10 else str(j)
                istest = '----' if j < nsteps else 'TEST'
                sline  = line.split();
                
                rt[f].append(float(sline[6]))
                rm[f][j] += float(sline[9])

                sys.stdout.write(istr+' '+jstr+' '+istest+' '+sline[6]+' '+sline[9]+'\n')
                
                j += 1
    

    with open('fwi.mpicuda.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)

        for fid in range(0, nfreqs):
            time   = np.array(rt[fid]).astype(np.float)
            metric = np.array(rm[fid]).astype(np.float)

            writer.writerow( (str(freqs[fid]).rstrip(), str(np.mean(metric))) )

            #print >> sys.stderr, 'freq: ', freq, 'TIME     - Mean:', '  {0:.8f}'.format(np.mean(time)), 'Median:', '  {0:.8f}'.format(np.median(time)), 'Max:', '  {0:.8f}'.format(np.amax(time)), 'Min:', '  {0:.8f}'.format(np.amin(time)), 'STD:', np.std(time)
            print >> sys.stderr, 'Freq: ', freqs[fid].rstrip(), 'Mcells/s - Mean:', '{0:.8f}'.format(np.mean(metric)), 'Median:', '{0:.8f}'.format(np.median(metric)), 'Max:', '{0:.8f}'.format(np.amax(metric)), 'Min:', '{0:.8f}'.format(np.amin(metric)), 'STD:', np.std(metric)



if __name__ == "__main__":
    main()
