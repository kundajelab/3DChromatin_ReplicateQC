#!/bin/bash

d="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#make metadata for the examples
metadata_samples=${d}/metadata.samples
printf "HIC001\t${d}/HIC001.res50000.gz\n" > ${metadata_samples}
printf "HIC002\t${d}/HIC002.res50000.gz\n" >> ${metadata_samples}

metadata_pairs=${d}/metadata.pairs
printf "HIC001\tHIC002\n" > ${metadata_pairs}
