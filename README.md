3DChromatin_ReplicateQC
1;95;0c===
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
1. Make sure you have R>=3.4.0 (needed for HiCRep).

1. Install [Anaconda2](https://www.continuum.io/downloads), which contains python as well as a set of commonly used packages. 3DChromatin_ReplicateQC is compatible with Python 2.

2. First clone the repository:

```
git clone http://github.com/kundajelab/3DChromatin_ReplicateQC
```

Then, run the following installation command:

```
3DChromatin_ReplicateQC/install_scripts/install_3DChromatin_ReplicateQC.sh
```
**Note if you are installing locally**: There are a few parameters you can provide to the installation script, to point it to your desired python installation (where you installed anaconda, e.g. `/home/anaconda/bin/python`), R installation, R library and modules. Thus, you can run the above as described in the following example.

Assume the following:
- path to your anaconda installation of python is `/home/my_anaconda2/bin/python`
- path to your R installation is `/home/my_R/bin/R`
- path to your R libraries is `/home/my_R_libraries`
- you wish to load any modules that could be loaded as `module load <module name>`

Then, your installation commands would look as below:
```
3DChromatin_ReplicateQC/install_scripts/install_3DChromatin_ReplicateQC.sh --pathtopython /home/my_anaconda2/bin/python --pathtor /home/my_R/bin/R --rlib /home/my_R_libraries --modules module1,module2` 
```

Quick start
====

Say you want to compare 2 contact maps. For this example, we will use a subset of datasets from Rao et al., 2014. 
First, configure the files used in the example:

```
3DChromatin_ReplicateQC/examples/configure_example.sh
```

Then run all methods (both QC and reproducibility as follows):

```
cd 3DChromatin_ReplicateQC
3DChromatin_ReplicateQC run_all --metadata_samples examples/metadata.samples --metadata_pairs examples/metadata.pairs --bins examples/Bins.w40000.bed.gz --outdir examples/output 
```

Output
======
The scores are summarized in the output directory under `scores/`.

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
GenomeDISCO|subsampling	lowest
GenomeDISCO|tmin		3
GenomeDISCO|tmax	3
GenomeDISCO|norm	sqrtvc
GenomeDISCO|scoresByStep	no
GenomeDISCO|removeDiag		yes
GenomeDISCO|transition		yes
HiCRep|h	5
HiCRep|maxdist	5000000
HiC-Spector|n	20
QuASAR|rebinning	resolution
SGE|text	"-l h_vmem=3G"
slurm|text	"--mem 3G"
```
Note: all of the above parameters need to be specified in the parameters file.

Here are details about setting these parameters:

**GenomeDISCO parameters**
- `GenomeDISCO|subsampling` This allows subsampling the datasets to a specific desired sequencing depth. Possible values are: `lowest` (subsample to the depth of the sample with the lower sequencing depth from the pair being compared), `<samplename>` where <samplename> is the name of the sample that is used to determine the sequencing depth to subsample from. 

- `GenomeDISCO|tmin` The minimum number of steps of random walk to perform. Integer, > 0.

- `GenomeDISCO|tmax` The max number of steps of random walk to perform. Integer, > tmin.
 
- `GenomeDISCO|norm` The normalization to use on the data when running GenomeDISCO. Possible values include: `uniform` (no normalization), `sqrtvc`.

- `GenomeDISCO|scoresByStep` Whether to report the score at each t. By default (GenomeDISCO|scoresByStep no), only the final reproducibility score is returned.

- `GenomeDISCO|removeDiag` Whether to set the diagonal to entries in the contact map to 0. By default (GenomeDISCO|removeDiag yes), the diagonal entries are set to 0.

- `GenomeDISCO|transition` Whether to convert the normalized contact map to an appropriate transition matrix before running the random walks. By default (GenomeDISCO|transition yes) the normalized contact map is converted to a proper transition matrix, such that all rows sum to 1 exactly.

**HiCRep parameters**
- `HiCRep|h` The h parameter in HiCRep that determines the extent of 2D smoothing. See the HiCRep paper (http://genome.cshlp.org/content/early/2017/08/30/gr.220640.117) for details. Integer, >=0.

- `HiCRep|maxdist` The maximum distance to consider when computing the HiCRep score. Integer, should be a amultiple of the resolution of the data.

**HiC-Spector parameters**
- `HiC-Spector|n` The number of eigenvectors to use for HiC-Spector. Integer, > 0.

**QuASAR parameters**
- `QuASAR|rebinning` The rebinning distance. See the QuASAR paper (https://www.biorxiv.org/content/early/2017/10/17/204438) for details. Integer.

**Job submission parameters**
- `SGE|text` Text to append to the job submission for SGE. The default is "-l h_vmem=3G".

- `slurm|text` Text to append to the job submission for slurm. The default is "--mem 3G". 

**Note about normalization**: At the moment, the different methods operate on different types of normalizations. For GenomeDISCO, the user can specify the desired normalization. For HiCRep and HiC-Spector the scores are computed on the provided data, without normalization.
Thus, if you have normalized data, then you can provide that as an input, and set `GenomeDISCO|norm` to uniform. If you have raw data, then your HiCRep and HiC-Spector scores will be run on the raw data, and GenomeDISCO will be run on the normalization you specify with `GenomeDISCO|norm`.

Running 3DChromatin_ReplicateQC step by step
============================================
3DChromatin_ReplicateQC consists of multiple steps, which are run in sequence by default. However, the user may decide to run the steps individually, which can be useful for instance when running 3DChromatin_ReplicateQC with job submission engines that runs the comparisons in parallel as separate jobs.

**3DChromatin_ReplicateQC steps**

**preprocess**

Preprocesses all datasets provided in `--metadata_samples`.

Example command: 
```
3DChromatin_ReplicateQC preprocess --metadata_samples examples/metadata.samples --bins examples/Bins.w40000.bed.gz --outdir examples/output --parameters_file examples/example_parameters.txt
```

**qc**

Runs QC methods on all samples provided in `--metadata_samples`. Note that this step is only compatible with `--methods QuASAR-QC`, since this is the only method that provides sample-level quality scores. If you provide any other methods to this step in addition to QuASAR-QC, they will not be run.

Example command: 
```
3DChromatin_ReplicateQC qc --metadata_samples examples/metadata.samples --outdir examples/output --methods QuASAR-QC
```

**concordance**

Runs reproducibility methods on all samples pairs provided in `--metadata_pairs`. Note that the only methods that can be run with this step are: GenomeDISCO, HiCRep, HiC-Spector and QuASAR-Rep. If QuASAR-QC is provided here, it will not be run. 

Example command: 
```
3DChromatin_ReplicateQC concordance --metadata_pairs examples/metadata.pairs --outdir examples/output --methods GenomeDISCO,HiCRep,HiC-Spector,QuASAR-Rep 
```

**summary**

Summarizes scores across all comparisons.

Example command: 
```
3DChromatin_ReplicateQC summary --metadata_samples examples/metadata.samples --metadata_pairs examples/metadata.pairs --bins examples/Bins.w40000.bed.gz --outdir examples/output --methods GenomeDISCO,HiCRep,HiC-Spector,QuASAR-Rep,QuASAR-QC
```

**cleanup**

Clean up superfluous files, leaving only the scores.

Example command: 
```
3DChromatin_ReplicateQC cleanup --outdir examples/output
```

Running this code with job submission engines
============

It is possible to run 3DChromatin_ReplicateQC with job submission engines, specifically either SGE or slurm.
To do so, modify the parameters `SGE|text` or `slurm|text` respectively, to add any additional parameters to the job run.

Then, run the steps sequentially (that is, wait for all jobs of a given step to complete before launching the next step), while specifying `--running_mode` to either `sge` or `slurm`.

For instance, an example analysis workflow for SGE would be:
```
3DChromatin_ReplicateQC preprocess --running_mode sge --metadata_samples examples/metadata.samples --bins examples/Bins.w40000.bed.gz --outdir examples/output --parameters_file examples/example_parameters.txt
3DChromatin_ReplicateQC qc --running_mode sge --metadata_samples examples/metadata.samples --outdir examples/output --methods QuASAR-QC
3DChromatin_ReplicateQC concordance --running_mode sge --metadata_pairs examples/metadata.pairs --outdir examples/output --methods GenomeDISCO,HiCRep,HiC-Spector,QuASAR-Rep
3DChromatin_ReplicateQC summary --running_mode sge --metadata_samples examples/metadata.samples --metadata_pairs examples/metadata.pairs --bins examples/Bins.w40000.bed.gz --outdir examples/output --methods GenomeDISCO,HiCRep,HiC-Spector,QuASAR-Rep,QuASAR-QC
3DChromatin_ReplicateQC cleanup --running_mode sge --outdir examples/output
```

Similarly, for slurm, change sge to slurm for the `--running_mode`.

More questions about this repository?
====
Contact Oana Ursu (oursu@stanford.edu)

Or submit an issue for this repository.

Thanks
===
**Code**

This repository was put together by Oana Ursu. Thanks to Michael Sauria for providing wrapper scripts around the QuASAR method, Tao Yang for his assistance in integrating HiCRep into this repository, and Koon-Kiu Yan for his assistance in integrating HiC-Spector into this repository.

**Testing**

Thanks to the Noble lab (Gurkan Yardminci, Jie Liu, Charles Grant), as well as Michael Sauria for testing the code out and for suggestions for improvement.

**docker** (coming soon)

Thanks to Anna Shcherbina for help making this code dockerized.