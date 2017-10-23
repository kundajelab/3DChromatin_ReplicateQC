#!/usr/bin/env python

import sys
import os

import hifive
import h5py
import numpy
import gzip

class Encode_Data(hifive.HiCData):

    def load_data_from_int(self, fendfilename, filelist):
        self.history += "HiCData.load_data_from_int(fendfilename='%s', filelist=%s, ) - " % (fendfilename, str(filelist))
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, \
                ("The fend file %s was not found. No data was loaded.\n") % (fendfilename),
            self.history += "Error: '%s' not found\n" % fendfilename
            return None
        self.fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                       os.path.dirname(self.file)), os.path.basename(fendfilename))
        self.fends = h5py.File(fendfilename, 'r')
        idx2index = {}
        if 'binned' in self.fends['/'].attrs and self.fends['/'].attrs['binned'] is not None:
            self.binned = True
            f = self.fends['bins'][...]
            for i in range(f.shape[0]):
                idx2index[f['idx'][i]] = (f['chr'][i], i)
            chr_indices = self.fends['bin_indices'][...]
        else:
            self.binned = False
            f = self.fends['fends'][...]
            for i in range(f.shape[0]):
                idx2index[f['idx'][i]] = (f['chr'][i], i)
            chr_indices = self.fends['chr_indices'][...]
        if 'fends' in self.fends and self.fends['fends'] is not None:
            self.re = True
        else:
            self.re = False
        self.history = self.fends['/'].attrs['history'] + self.history
        self.chr2int = {}
        chroms = self.fends['chromosomes'][...]
        for i, j in enumerate(chroms):
            self.chr2int[j] = i
        # load data from all files, skipping if chromosome not in the fend file.
        if isinstance(filelist, str):
            filelist = [filelist]
        int_filelist = []
        for filename in filelist:
            int_filelist.append("%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(filename)),
                                           os.path.dirname(self.file)), os.path.basename(filename)))
        self.int_filelist = ",".join(int_filelist)
        fend_pairs = []
        for i in range(len(chroms)):
            fend_pairs.append([])
            for j in range(i + 1):
                fend_pairs[i].append({})
        total_reads = 0
        data = []
        for i in range(chroms.shape[0]):
            data.append([])
            for j in range(i + 1):
                data[-1].append({})
        for fname in filelist:
            if not os.path.exists(fname):
                if not self.silent:
                    print >> sys.stderr, ("The file %s was not found...skipped.\n") % (fname.split('/')[-1]),
                self.history += "'%s' not found, " % fname
                continue
            if not self.silent:
                print >> sys.stderr, ("Loading data from %s...") % (fname.split('/')[-1]),
            input = gzip.open(fname, 'r', 1)
            new_reads = 0
            for line in input:
                temp = line.strip('\n').split('\t')
                try:
                    chrint1, bin1 = idx2index[int(temp[1])]
                    chrint2, bin2 = idx2index[int(temp[3])]
                    count = int(temp[4])
                except:
                    continue
                if chrint2 < chrint1:
                    chrint1, chrint2, bin1, bin2 = chrint2, chrint1, bin2, bin1
                elif chrint1 == chrint2 and bin2 < bin1:
                    bin1, bin2 = bin2, bin1
                key = (bin1 - chr_indices[chrint1], bin2 - chr_indices[chrint2])
                if key not in data[chrint2][chrint1]:
                    data[chrint2][chrint1][key] = 0
                data[chrint2][chrint1][key] += count
                new_reads += count
            input.close()
            total_reads += new_reads
            if not self.silent:
                print >> sys.stderr, ("\r%s\r%i validly-mapped read pairs loaded.\n") % (' ' * 50, new_reads),
        total_bin_pairs = 0
        for i in range(len(data)):
            for j in range(len(data[i])):
                total_bin_pairs += len(data[i][j])
        if total_bin_pairs == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            self.history += "Error: no valid data loaded\n"
            return None
        if not self.silent:
            print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i valid bin pairs\n") %\
                             (total_reads, total_bin_pairs),
        # write fend pairs to h5dict
        self._parse_fend_pairs(data)
        self.history += 'Success\n'
        return None

def main():
    int_fname, fend_fname, out_fname = sys.argv[1:4]
    data = Encode_Data(out_fname, 'w')
    data.load_data_from_int(fend_fname, int_fname)
    data.save()


if __name__ == "__main__":
    main()
