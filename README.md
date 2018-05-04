
This code is associated with the paper from Weir et al., "The AAA protein Msp1 mediates clearance of excess tail-anchored proteins from the peroxisomal membrane". eLife, 2017. http://dx.doi.org/10.7554/eLife.28507


# README for the Weir et al. 2017 paper data

This repository contains the following:
- Sample scripts for segmentation of peroxisomes and mitochondria followed by fluorescence intensity measurement using [pyto_segmenter](https://github.com/deniclab/pyto_segmenter) (see sample_segmentation folder)
- Post-segmentation and fluorescence intensity measurement data, in csv format, for each quantitative microscopy experiment
- R scripts for generation of quantitative microscopy

Two other repositories are critical for complete re-production of the analyses in this paper:
- pyto_segmenter package (above)
- Dryad repository of raw microscopy data

__Before trying to re-capitulate segmentation using the scripts in sample_segmentation, be sure to read the readme within that folder.__

The R scripts contained in this repository require the following packages to be installed:
- readr
- plyr
- dplyr
- ggplot2
- gridExtra
- Cairo
- minpack.lm

See [pyto_segmenter](https://github.com/deniclab/pyto_segmenter) for the list of python dependencies.

_Last updated 9.5.17_
