#!/usr/bin/env python

import sys
import glob
import argparse
import re
import os

def main():
    scorefile = open( sys.argv[1], 'r')
    sample1, sample2 = sys.argv[2:4]
    lines=scorefile.readlines()
    header = None
    for i in range(len(lines)):
        if lines[i] == 'Replicate Score Results\n':
            header = i + 2
            break
    header_line = lines[header].rstrip('\n').split('\t')
    data_line = lines[header + 1].rstrip('\n').split('\t')
    d={}
    for i in range(3, len(header_line)):
        d[header_line[i]]=data_line[i]    
    outfile=open(os.path.dirname(sys.argv[1])+'/'+re.sub('QuASAR-QC','',re.sub('QuASAR-Rep.scores.','',os.path.basename(sys.argv[1]))),'w')
    for chromo, score in d.iteritems():
        outfile.write(re.sub('.quasar_transform','',sample1)+'\t'+re.sub('.quasar_transform','',sample2)+'\t'+'chr'+re.sub('chr','',str(chromo))+'\t'+score+'\n')
    outfile.close()

if  __name__ == "__main__":
    main()
