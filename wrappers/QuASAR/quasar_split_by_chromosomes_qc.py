#!/usr/bin/env python

import sys
import glob
import argparse
import re
import os

def main():
    scorefile = open( sys.argv[1], 'r')
    samplename=sys.argv[2]
    lines=scorefile.readlines()
    header = None
    for i in range(len(lines)):
        if lines[i] == 'Quality Score Results\n':
            header = i + 2
            break
    header_line = lines[header].rstrip('\n').split('\t')
    data_line = lines[header + 1].rstrip('\n').split('\t')
    d={}
    for i in range(2, len(header_line)):
        d[header_line[i]]=data_line[i]
    for chromo, score in d.iteritems():
        outfile=open(os.path.dirname(sys.argv[1])+'/'+'chr'+chromo+'.'+re.sub('.QuASAR-QC','',re.sub('.QuASAR-Rep.','',os.path.basename(sys.argv[1]))),'w')
        outfile.write(samplename+'\t'+score+'\n')
        outfile.close()

if  __name__ == "__main__":
    main()
