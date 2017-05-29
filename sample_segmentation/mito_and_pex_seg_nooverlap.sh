#!/bin/bash
img_dir=$1
high_threshold=$2
low_threshold=$3
cd $img_dir
# change the following path to point to your clone of the segmentation scripts
python3 /path/to/Weir_2017_analysis/sample_segmentation/mito_and_pex_seg_nooverlap.py \
    -d $img_dir -ht $high_threshold -lt $low_threshold
