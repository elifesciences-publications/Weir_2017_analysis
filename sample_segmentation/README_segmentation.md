# README for segmentation sample scripts
### Author: Nicholas Weir, Harvard University


The scripts included here are intended to achieve two goals:
1. Provide the user with scripts demonstrating how to perform your own segmentation using our pyto_segmenter package (batch_mito_seg.py, batch_pex_seg.py, and mito_and_pex_seg_nooverlap.py)
2. Provide one example of a complete segmentation and YFP signal measurement pipeline for a figure (fig_1_analysis.py). This pipeline produces identical data to fig_1C_S1D_data.csv (see figure_1 folder).

__WARNING__: For our actual analyses in the paper, segmentation and YFP signal detection was parallelized using the Harvard University Research Computing Odyssey cluster. As we assume most readers won't have access to this cluster, we are providing sample scripts for performing these analyses without parallelization. These scripts may therefore take a _very_ long time (days) to run. Parallelized analysis scripts using the slurm distributed computing framework are available upon request.

__Description of files in this directory__:
- batch_mito_seg.py: Script for segmenting mitochondria from a set of mitochondrial marker images.
- batch_pex_seg.py: Script for segmenting peroxisomes from a set of peroxisome marker images.
- mito_and_pex_seg_nooverlap.py: Script for segmenting both peroxisomes and mitochondria from a set of images containing both marker channels, which additionally removes regions of mitochondria that overlap with segmented peroxisomes. This was critical for the production of Figure 1.
- norm_camera.py: This script was used to normalize both within and across images to account for uneven illumination and fluctuations in camera noise that occurred during time-lapse images. _All images provided in the Dryad repository for this paper have already had this normalization performed._
- mito_and_pex_seg_nooverlap.sh: A shell script to execute segmentation using mito_and_pex_seg_nooverlap.py. Provided as an example of how to launch segmentation on a linux/osx system for a directory of images.

__To perform sample segmentation and analysis of Figure 1 data using the scripts provided:__

1. Download the source data from the Dryad repository for this paper. Unzip the figure_1c_s1d data file.
2. Open the file mito_and_pex_seg_nooverlap.sh in a text editor and change the path to mito_and_pex_seg_nooverlap.py to point to the script file.
3. Make sure that the pyto_segmenter package is in your PYTHONPATH.
4. In a terminal or bash shell, run the following commands, changing the paths provided below to point to your clone of the data and scripts:

/path/to/mito_and_pex_seg_nooverlap.sh -d /path/to/fig_1c_s1d/ -ht 750 -lt 375

python3 /path/to/yfp_analysis.py -d /path/to/fig_1c_s1d/

__Completing the above steps will produce the csv file provided for Figure 1C and S1D. See the associated R script for further analysis and plotting.__

_last updated 5.29.2017_
