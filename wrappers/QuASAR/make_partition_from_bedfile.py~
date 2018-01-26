import argparse
import sys
import os

import hifive
import h5py
import numpy
import gzip
import re
import subprocess as subp

class Encode_Fend(hifive.Fend):

    def _load_binned_from_bed(self, fname):
        chromosomes = []
        data = {}
        feature_names = []
        input = open(fname, 'r')
        for line in input:
            temp = line.strip('\n').split('\t')
            start = int(temp[1])
            chrom = temp[0]
            if chrom not in data:
                data[chrom] = []
            stop = int(temp[2])
            idx = int(temp[3])
            data[chrom].append((start, stop, idx))
        input.close()
        chromosomes = data.keys()
        chromosomes.sort()
        chromosomes = numpy.array(chromosomes)
        dtypes = [('chr', numpy.int32), ('start', numpy.int32), ('stop', numpy.int32), ('mid', numpy.int32),
                  ('idx', numpy.int32)]
        dtypes2 = [('start', numpy.int32), ('stop', numpy.int32), ('idx', numpy.int32)]
        count = 0
        for i, chrom in enumerate(chromosomes):
            data[chrom].sort()
            data[chrom] = numpy.array(data[chrom], dtype=numpy.dtype(dtypes2))
            count += data[chrom].shape[0]
        bins = numpy.empty(count, dtype=numpy.dtype(dtypes))
        pos = 0
        for i, chrom in enumerate(chromosomes):
            data_len = data[chrom].shape[0]
            bins['chr'][pos:(pos + data_len)] = i
            bins['start'][pos:(pos + data_len)] = data[chrom]['start'][:]
            bins['stop'][pos:(pos + data_len)] = data[chrom]['stop'][:]
            bins['mid'][pos:(pos + data_len)] = (data[chrom]['start'][:] + data[chrom]['stop'][:]) / 2
            bins['idx'][pos:(pos + data_len)] = data[chrom]['idx'][:]
            pos += data_len
        return [bins, chromosomes]

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--nodes')
    parser.add_argument('--partition')
    parser.add_argument('--resolution',type=int)
    parser.add_argument('--subset_chromosomes',default='NA')
    args = parser.parse_args()

    nodes_modified_file=args.partition+'.tmp'
    subp.check_output(['bash','-c','zcat -f '+args.nodes+' | sort -k1,1 -k2,2n | gzip > '+nodes_modified_file+'.sorted'])
    nodes_modified=open(nodes_modified_file,'w')
    for line in gzip.open(nodes_modified_file+'.sorted','r'):
        items=line.strip().split('\t')
        chromo,start,end,name=items[0],items[1],items[2],items[3]
        if args.subset_chromosomes!='NA':
            if chromo not in args.subset_chromosomes.split(','):
                continue
        nodes_modified.write(re.sub('chr','',chromo)+'\t'+start+'\t'+end+'\t'+name+'\n')
    nodes_modified.close()

    #uniform bins
    myfends=Encode_Fend(args.partition, mode='w', binned=int(args.resolution)) 
    myfends.load_bins(nodes_modified_file, format='bed') 
    myfends.save() 

    os.remove(nodes_modified_file)
    os.remove(nodes_modified_file+'.sorted')

main()
