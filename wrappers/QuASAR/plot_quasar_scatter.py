#!/usr/bin/env python

import sys
import glob
import argparse

import h5py
import numpy

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import rcParams


def main():
    parser = generate_parser()
    args = parser.parse_args()
    infile1 = h5py.File(args.input1, 'r')
    infile2 = h5py.File(args.input2, 'r')
    resolutions = numpy.intersect1d(infile1['resolutions'][...], infile2['resolutions'][...])
    chroms = numpy.intersect1d(infile2['chromosomes'][...], infile2['chromosomes'][...])
    results = {}
    data1 = load_data(infile1, chroms, resolutions)
    data2 = load_data(infile2, chroms, resolutions)
    infile1.close()
    infile2.close()
    results = {}
    results[(args.input1.split('/')[-1].strip('.quasar'), args.input2.split('/')[-1].strip('.quasar'))] = correlate_samples(data1, data2)
    for resolution in data1.keys():
        for chromo in chroms:
            plt.scatter(data1[resolution][chromo][1].flatten(),data2[resolution][chromo][1].flatten(),alpha=0.1,color='red')
            plt.show()
            plt.savefig(args.output+'.res'+str(resolution)+'.chr'+chromo+'.pdf')

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

def correlate_samples(data1, data2):
    resolutions = data1.keys()
    resolutions.sort()
    chroms = data1[resolutions[0]].keys()
    chroms.sort()
    results = numpy.zeros((len(resolutions), len(chroms) + 1), dtype=numpy.float64)
    for i, res in enumerate(resolutions):
        temp = numpy.zeros(6, dtype=numpy.float64)
        for j, chrom in enumerate(chroms):
            start1 = data1[res][chrom][0]
            start2 = data2[res][chrom][0]
            if start2 > start1:
                start1 = (start2 - start1) / res
                start2 = 0
            else:
                start2 = (start1 - start2) / res
                start1 = 0
            stop1 = data1[res][chrom][1].shape[0] - start1
            stop2 = data2[res][chrom][1].shape[0] - start2
            if stop2 > stop1:
                stop2 = stop1 + start2
                stop1 = start1 + stop1
            else:
                stop1 = stop2 + start1
                stop2 = start2 + stop2
            trans1 = data1[res][chrom][1][start1:stop1, :]
            valid1 = data1[res][chrom][2][start1:stop1, :]
            trans2 = data2[res][chrom][1][start2:stop2, :]
            valid2 = data2[res][chrom][2][start2:stop2, :]
            valid = numpy.where(valid1 & valid2)
            if valid[0].shape[0] == 0:
                continue
            X = numpy.sum(trans1[valid])
            Y = numpy.sum(trans2[valid])
            X2 = numpy.sum(trans1[valid] ** 2.0)
            Y2 = numpy.sum(trans2[valid] ** 2.0)
            XY = numpy.sum(trans1[valid] * trans2[valid])
            N = valid[0].shape[0]
            temp += [X, Y, X2, Y2, XY, N]
            Xmu = X / N
            Ymu = Y / N
            X2mu = X2 / N
            Y2mu = Y2 / N
            XYmu = XY / N
            if Xmu ** 2.0 > X2mu or Ymu ** 2.0 > Y2mu:
                continue
            Xstd = (X2mu - Xmu ** 2.0) ** 0.5
            Ystd = (Y2mu - Ymu ** 2.0) ** 0.5
            results[i, j] = (XYmu - Xmu * Ymu) / (Xstd * Ystd)
        Xmu = temp[0] / temp[5]
        Ymu = temp[1] / temp[5]
        X2mu = temp[2] / temp[5]
        Y2mu = temp[3] / temp[5]
        XYmu = temp[4] / temp[5]
        Xstd = (X2mu - Xmu ** 2.0) ** 0.5
        Ystd = (Y2mu - Ymu ** 2.0) ** 0.5
        results[i, -1] = (XYmu - Xmu * Ymu) / (Xstd * Ystd)
    return results

def write_results(results, resolutions, chroms, args):
    resolutions.sort()
    chroms.sort()
    output = open(args.output, 'w')
    temp = "Sample1\tSample2\tResolution\tAll"
    for chrom in chroms:
        temp += "\t%s" % chrom
    print >> output, temp
    keys = results.keys()
    keys.sort()
    for r, res in enumerate(resolutions):
        for key in keys:
            temp = [key[0], key[1], str(res), str(results[key][r, -1])]
            for i in range(results[key].shape[1] - 1):
                temp.append(str(results[key][r, i]))
            print >> output, '\t'.join(temp)
    output.close()

def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Find a quality score for a HiC dataset"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(dest="input1", type=str, action='store', help="First quasar file name")
    parser.add_argument(dest="input2", type=str, action='store', help="Second quasar file name")
    parser.add_argument(dest="output", type=str, action='store', help="Results file name")
    return parser

if  __name__ == "__main__":
    main()
