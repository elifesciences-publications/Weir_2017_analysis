# import dependencies
import os
import sys
import argparse
# the following line assumes that the pyto_segmenter package is already
# installed and is in the user's PYTHONPATH. See
# https://github.com/deniclab/pyto_segmenter to acquire this segmentation
# package.
from pyto_segmenter import PexSegment, MitoSegment



parser = argparse.ArgumentParser(description = 'Segment peroxisomes from \
                                 images and return pickled objects.')
parser.add_argument('-d', '--img_dir', required = True, 
                    help = 'directory containing images to segment.')
parser.add_argument('-ht', '--high_threshold', required = True,
                    help = 'high threshold for canny segmentation.')
parser.add_argument('-lt', '--low_threshold', required = True,
                    help = 'low threshold for canny segmentation.')

args = parser.parse_args()
print(args)
img_dir = args.img_dir
high_threshold = int(args.high_threshold)
low_threshold = int(args.low_threshold)

# perform segmentation on both peroxisome and mitochondrial marker images,
# removing overlapping peroxisomal regions from mitochondrial objects.
os.chdir(img_dir)
flist = os.listdir()
imgs = [f for f in flist if '.tif' in f.lower()]
pex_imgs = [im for im in imgs if '594' in im]
mito_imgs = [im for im in imgs if '447' in im]
pex_imgs.sort()
mito_imgs.sort()
if len(pex_imgs) != len(mito_imgs):
    raise ValueError('Length of pex and mito img sets do not match.')
mito_list = mito_imgs 
pex_list = pex_imgs
for i in range(0,len(pex_list)):
    os.chdir(img_dir)
    # segment the relevant pex marker image 
    print('SEGMENTING ' + pex_list[i] + ' and ' + mito_list[i])
    pex_segmenter = PexSegment.PexSegmenter(pex_list[i], seg_method = 'canny',
                                           high_threshold = high_threshold,
                                           low_threshold = low_threshold)
    pex_obj = pex_segmenter.segment()
    pex_obj.rm_border_objs()
    os.chdir(img_dir)
    # segment the mitochondrial marker from the same stage position
    mito_segmenter = MitoSegment.MitoSegmenter(mito_list[i], seg_method = 'canny',
                                               high_threshold = 250,
                                               low_threshold = 125,
                                               min_cutoff = 2300)
    mito_obj = mito_segmenter.segment()
    mito_obj.rm_border_objs()
    # remove mitochondrial object regions that overlap with peroxisomes
    mito_obj.mitochondria[pex_obj.threshold_img == 1] = 0 # remove overlap
    pex_obj.peroxisomes[mito_obj.threshold_img == 1] = 0 # remove overlap
    # save mito objects for use in measuring YFP signal later
    mito_obj.pickle(output_dir = img_dir + '/pickles', filename = 
                    mito_obj.filename[0:mito_obj.filename.index('.tif')] + '_mito.pickle' )
    # save pex objects for use in measuring YFP signal later
    pex_obj.pickle(output_dir = img_dir + '/pickles', filename =
                   pex_obj.filename[0:pex_obj.filename.index('.tif')] + '_pex.pickle' )
    del pex_obj
    del mito_obj
