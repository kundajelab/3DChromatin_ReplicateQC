
from __future__ import print_function
import argparse
import os
import sys
#from software.genomedisco.genomedisco import concordance_utils

def main():
    command_methods = {'preprocess': concordance_utils.preprocess,
                         'concordance': concordance_utils.concordance,
                         'summary': concordance_utils.summary,
                         'cleanup':concordance_utils.clean_up,
                       'run_all': concordance_utils.run_all}
    
    command, args = concordance_utils.parse_args_replicateqc()
    command_methods[command](**args)


if __name__ == "__main__":
    main()
