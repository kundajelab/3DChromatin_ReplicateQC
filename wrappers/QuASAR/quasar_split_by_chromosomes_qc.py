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
    titles,scores=lines[0],lines[1]
    scorelist=scores.strip().split('\t')
    d={}
    l=titles.strip().split('\t')
    for i in range(len(l)):
        if l[i] in ['All','Sample1','Resolution']:
            continue
        d[l[i]]=i    
    print d
    print scorelist
    for chromo in d.keys():
        print chromo
        outfile=open(os.path.dirname(sys.argv[1])+'/'+'chr'+chromo+'.'+re.sub('.QuASAR-QC','',re.sub('.QuASAR-Rep.','',os.path.basename(sys.argv[1]))),'w')
        outfile.write(samplename+'\t'+scorelist[d[chromo]]+'\n')
        outfile.close()

if  __name__ == "__main__":
    main()
