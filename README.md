3DChromatin_ReplicateQC
===
Welcome! This repository will allow you to measure the quality and reproducibility of 3D genome data.
It computes the following:

**Quality scores per sample** using
- QuASAR-QC (part of the hifive suite at http://github.com/bxlab/hifive)

**Reproducibility for pairs of samples** using 
- HiCRep (http://github.com/qunhualilab/hicrep)
- GenomeDISCO (http://github.com/kundajelab/genomedisco)
- HiC-Spector (http://github.com/gersteinlab/HiC-spector)
- QuASAR-Rep (part of the hifive suite at http://github.com/bxlab/hifive)

![ScreenShot](https://github.com/kundajelab/3DChromatin_ReplicateQC/blob/master/Fig1_outline_2017_08_23.png)

Installation
=====

1. Install [Anaconda2](https://www.continuum.io/downloads), which contains python as well as a set of commonly used packages. 3DChromatin_ReplicateQC is compatible with Python 2.

2. Run the installation script as below:

```
git clone http://github.com/kundajelab/3DChromatin_ReplicateQC
cd 3DChromatin_ReplicateQC
install_scripts/install_3DChromatin_ReplicateQC.sh
```
**Note if you are installing these locally**: There are a few parameters you can provide to the installation script, to point it to your desired python installation (where you installed anaconda, e.g. `/home/anaconda/bin/python`), R installation, R library, modules and bedtools installation. Thus, you can run the above script as follows:

```
install_scripts/install_3DChromatin_ReplicateQC.sh --pathtopython /path/to/your/python --pathtor /path/to/your/R --rlib /path/to/your/Rlibrary --modules modulename --pathtobedtools path/to/your/bedtools
```

Quick start
====

Say you want to compare 2 contact maps. For this example, we will use a subset of datasets from Rao et al., 2014. 
First, configure the files used in the example:

```
examples/configure_example.sh
```

Then run all methods (both QC and reproducibility as follows):

```
python 3DChromatin_ReplicateQC.py run_all --metadata_samples examples/metadata.samples --metadata_pairs examples/metadata.pairs --bins examples/Nodes.w40000.bed.gz --outdir examples/output 
```

Output
======
The scores are summarized in the output directory under `results/summary/`.

Inputs
=============
Before running 3DChromatin_ReplicateQC, make sure to have the following files:

- **contact map** For each of your samples, you need a file containing the counts assigned to each pair of bins in your contact map, and should have the format `chr1 bin1 chr2 bin2 value`. Note: 3DChromatin_ReplicateQC assumes that this file contains the contacts for all chromosomes, and will split it into individual files for each chromosome.

- **bins** This file contains the full set of genomic regions associated with your contact maps, in the format `chr start end name` where name is the name of the bin as used in the contact map files above. 3DChromatin_ReplicateQC supports both fixed-size bins and variable-sized bins (e.g. obtained by partitioning the genome into restriction fragments). 

3DChromatin_ReplicateQC takes the following inputs:

- `--metadata_samples` Information about the samples being compared. Tab-delimited file, with columns "samplename", "samplefile". Note: each samplename should be unique. Each samplefile listed here should follow the format "chr1 bin1 chr2 bin2 value

- `--metadata_pairs` Each row is a pair of sample names to be compared, in the format "samplename1 samplename2". Important: sample names used here need to correspond to the first column of the --metadata_samples file.

- `--bins` A (gzipped) bed file of the all bins used in the analysis. It should have 4 columns: "chr start end name", where the name of the bin corresponds to the bins used in the contact maps.

- `--re_fragments` Add this flag if the bins are not uniform bins in the genome (e.g. if they are restriction-fragment-based).By default, the code assumes the bins are of uniform length.

- `--methods` Which method to use for measuring concordance or QC. Comma-delimited list. Possible methods: "GenomeDISCO", "HiCRep", "HiC-Spector", "QuASAR-Rep", "QuASAR-QC". By default all methods are run

- `--parameters_file` File with parameters for reproducibility and QC analysis. For details see ["Parameters file"](#parameters-file)

- `--outdir` Name of output directory. DEFAULT: replicateQC

- `--running_mode` The mode in which to run the analysis. This allows you to choose whether the analysis will be run as is, or submitted as a job through sge or slurm. Available options are: "NA" (default, no jobs are submitted). Coming soon: "sge", "slurm"

- `--concise_analysis` Set this flag to obtain a concise analysis, which means replicateQC is measured but plots that might be more time/memory consuming are not created. This is useful for quick testing or running large-scale analyses on hundreds of comparisons.

- `--subset_chromosomes` Comma-delimited list of chromosomes for which you want to run the analysis. By default the analysis runs on all chromosomes for which there are data. This is useful for quick testing

Analyzing multiple dataset pairs
======
To analyze multiple pairs of contact maps (or multiple contact maps if just computing QC), all you need to do is add any additional datasets you want to analyze to the `--metadata_samples` file and any additional pairs of datasets you want to compare to the `--metadata_pairs` files. 

Parameters file
======

The parameters file specifies the parameters to be used with 3DChromatin_ReplicateQC. The format of the file is: `method_name parameter_name parameter_value`. The default parameters file used by 3DChromatin_ReplicateQC is:

```
GenomeDISCO	subsampling	lowest
GenomeDISCO	tmin		3
GenomeDISCO	tmax	3
GenomeDISCO	norm	sqrtvc
HiCRep		h	5
HiCRep		maxdist	5000000
HiC-Spector	n	20
```
Note: all of the above parameters need to be specified in the parameters file.

More questions about this repository?
====
Contact Oana Ursu

oursu@stanford.edu

Thanks to Michael Sauria for providing wrapper scripts around the QuASAR method, Tao Yang for his assistance in integrating HiCRep into this repository, and Koon-Kiu Yan for his assistance in integrating HiC-Spector into this repository.