
import argparse
import hifive
import sys
import glob
import h5py
import numpy

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import rcParams

def load_data(infile, chroms, resolutions):
    starts = infile['starts'][...]
    chromosomes = infile['chromosomes'][...]
    data = {}
    for res in resolutions:
        data[res] = {}
        for i, chrom in enumerate(chromosomes):
            if chrom not in chroms:
                continue
            start = (starts[i] / res) * res
            dist = infile['dist.%s.%i' % (chrom, res)][...]
            valid_rows = infile['valid.%s.%i' % (chrom, res)][...]
            corr = infile['corr.%s.%i' % (chrom, res)][...]
            valid = numpy.zeros(corr.shape, dtype=numpy.bool)
            N, M = corr.shape
            valid = numpy.zeros((N, M), dtype=numpy.int32)
            for i in range(min(N - 1, M)):
                P = N - i - 1
                valid[:P, i] = valid_rows[(i + 1):] * valid_rows[:P]
            temp = corr * dist
            valid[numpy.where(numpy.abs(temp) == numpy.inf)] = False
            data[res][chrom] = [start, temp, valid]
    return data

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--transform')
    parser.add_argument('--out')
    args = parser.parse_args()
    
    infile1 = h5py.File(args.transform, 'r')
    resolutions = infile1['resolutions'][...]
    chroms = infile1['chromosomes'][...]
    data1 = load_data(infile1, chroms, resolutions)
    infile1.close()

    '''
    #for now, don't plot this
    for resolution in data1.keys():
        for chromo in chroms:
            N = data1[resolution][chromo][1].shape[0]
            full=numpy.empty((N,N))
            #full=full/0
            for i in range(100):
                temp1 = numpy.arange(N - i - 1)
                temp2 = numpy.arange(i+1, N)
                full[temp1, temp2] = data1[resolution][chromo][1][temp1, i]
                full[temp2, temp1] = full[temp1, temp2]
            x=0.8
            plt.matshow(full,cmap='seismic',vmin=-x,vmax=x)
            plt.colorbar()
            plt.show()
            plt.savefig(args.out+'.res'+str(resolution)+'.chr'+chromo+'.pdf')    
   '''
main()
