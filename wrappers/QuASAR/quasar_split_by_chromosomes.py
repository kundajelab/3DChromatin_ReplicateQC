#!/usr/bin/env python

import sys
import glob
import argparse
import re
import os

def main():
    scorefile = open( sys.argv[1], 'r')
    lines=scorefile.readlines()
    titles,scores=lines[0],lines[1]
    scorelist=scores.strip().split('\t')
    d={}
    l=titles.strip().split('\t')
    for i in range(len(l)):
        if l[i] in ['All','Sample1','Sample2','Resolution']:
            continue
        d[l[i]]=i    
    for chromo in d.keys():
        outfile=open(os.path.dirname(sys.argv[1])+'/'+'chr'+chromo+'.'+re.sub('QuASAR-QC','',re.sub('QuASAR-Rep.','',os.path.basename(sys.argv[1]))),'w')
        outfile.write(re.sub('.quasar_transform','',scorelist[0])+'\t'+re.sub('.quasar_transform','',scorelist[1])+'\t'+scorelist[d[chromo]]+'\n')
        outfile.close()

if  __name__ == "__main__":
    main()
