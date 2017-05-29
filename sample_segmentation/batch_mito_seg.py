# import dependencies
import os
import sys
import argparse
# the following line assumes that the pyto_segmenter package is already
# installed and is in the user's PYTHONPATH. See
# https://github.com/deniclab/pyto_segmenter to acquire this segmentation
# package.
from pyto_segmenter import MitoSegment



parser = argparse.ArgumentParser(description = 'Segment mitochondria from \
                                 images and return pickled objects.')
parser.add_argument('-d', '--img_dir', required = True, 
                    help = 'directory containing images to segment.')
parser.add_argument('images', nargs = '*',
                    help = 'filenames for images, can be full path or just \
                    the image filename.')
# process arguments provided to this script.
args = parser.parse_args()
print(args)
img_dir = args.img_dir
images = [img[2:] for img in args.images]

# perform segmentation on peroxisome marker image names provided as an
# argument to this script.
for img in images:
    os.chdir(img_dir)
    print('SEGMENTING ' + img)
    mito_segmenter = MitoSegment.MitoSegmenter(img, seg_method = 'canny',
                                               high_threshold = 250,
                                               low_threshold = 125,
                                               min_cutoff = 2300)       
    mito_obj = mito_segmenter.segment()
    mito_obj.rm_border_objs()
    mito_obj.pickle(output_dir = img_dir + '/pickles')
    del mito_obj
