# import dependencies
import os
import sys
import argparse
# the following line assumes that the pyto_segmenter package is already
# installed and is in the user's PYTHONPATH. See
# https://github.com/deniclab/pyto_segmenter to acquire this segmentation
# package.
from pyto_segmenter import PexSegment, MitoSegment
import numpy as np
import pandas as pd
from skimage import io
import pickle
import re
# parse arguments from command line
parser = argparse.ArgumentParser(description = 'Measure fluorescence \
                                 intensity at pexs and mitochondria.')
parser.add_argument('-d', '--img_dir', required = True, 
                    help = 'directory containing images to segment.')

args = parser.parse_args()
print(args)
img_dir = args.img_dir

os.chdir(img_dir)
# get the filenames for all images in the directory
files = [f for f in os.listdir() if '.tif' in f.lower()]
# extract the YFP images based on the laser channel identifier
yfp_imgs = [y for y in files if '515' in y]
if expt_type == 'yfp':
    output_frame = pd.DataFrame({'img': [],
                                 'obj_channel': [],
                                 'obj_number': [],
                                 'volume': [],
                                 'yfp': [],
                                })
pickles = get_pickle_set(img_dir)
yfp_ids = get_img_ids(yfp_imgs)
yfp_ids = {v: k for k, v in yfp_ids.items()}
# begin measuring YFP signal for a given segmented object file
for p in pickles:
    os.chdir(img_dir + '/pickles')
    cfile = open(p,'rb')
    cpickle = pickle.load(cfile)
    print('current segmented object image: ' + cpickle.filename)
    # get the stage position identifier and the laser channel using get_img_ids
    (pickle_id, pickle_channel) = get_img_ids([cpickle.filename],
                                              return_channel = True)
    pickle_id = pickle_id[cpickle.filename]
    # use the attribute type to determine if its mitochondria or pexs
    pickle_channel = pickle_channel[cpickle.filename]
    if hasattr(cpickle, 'mitochondria'):
        obj_type = 'mito'
    else:
        obj_type = 'pex'
    print('current segmented object identifier: ' + pickle_id)
    print('current segmented object channel: ' + pickle_channel)
    os.chdir(img_dir)
    print('current YFP image file: ' + yfp_ids[pickle_id])
    # read in the YFP fluorescence image
    yfp_img = io.imread(yfp_ids[pickle_id])
    # initialize variables to hold values measured next
    yfp_vals = {}
    volumes_v2 = {}
    for obj in cpickle.obj_nums:
        print('     current obj number: ' + str(obj))
        if obj_type == 'mito':
            # measure the total YFP signal that overlaps with the object
            yfp_vals[obj] = sum(yfp_img[cpickle.mitochondria == obj])
            # get the object volume
            volumes_v2[obj] = len(np.flatnonzero(cpickle.mitochondria == obj))
        else:
            yfp_vals[obj] = sum(yfp_img[cpickle.peroxisomes == obj])
            volumes_v2[obj] = len(np.flatnonzero(cpickle.peroxisomes == obj))
        print('     yfp intensity: ' + str(yfp_vals[obj]))
        print('')
    print('Appending data to output...')
    # use a pandas dataframe to generate the table for later output
    currimg_data = pd.DataFrame({'img': pd.Series(data =
                                                  [cpickle.filename]*len(cpickle.obj_nums),
                                                  index = cpickle.obj_nums),
                                 'obj_channel': pd.Series(data =
                                                          [pickle_channel]*len(cpickle.obj_nums),
                                                      index = cpickle.obj_nums),
                                 'obj_number': pd.Series(data = cpickle.obj_nums,
                                                         index = cpickle.obj_nums),
                                 'volume': volumes_v2,
                                 'yfp': yfp_vals,
                                })
    output_frame = pd.concat([output_frame, currimg_data])
print('')
print('-----------------------------------------------------------------')
print('-----------------------------------------------------------------')
print('')
print('saving data...')
if not os.path.isdir(img_dir + '/analysis_output'):
    os.mkdir(img_dir + '/analysis_output')
output_frame.to_csv(img_dir + '/analysis_output/' + str(array_n) +
                    '_analysis_output.csv')

## helper functions ##
def get_pickle_set(img_dir):
    '''Get the subset of pickles for analysis by a given instance.'''
    os.chdir(img_dir + '/pickles')
    pickle_list = [p for p in os.listdir() if '.pickle' in p]
    pickle_list.sort()
    return(pickle_list)

def get_img_ids(img_files, return_channel = False):
    '''Extracts image filenames lacking wavelength identifiers.'''
    channel_re = re.compile('^w\d+[A-Za-z]*[ .]')
    img_ids = []
    channels = []
    for img in img_files:
        print('_______________________________________________________')
        print('     generating image identifier for ' + img)
        split_im = img.split('_')
        rm_channel = '_'.join([i for i in split_im if re.search(channel_re, i)
                               == None])
        channel = [i for i in split_im if re.search(channel_re, i)]
        if len(channel) > 1:
            print('WARNING: more than one match to channel ID string!')
        channel = channel[0].split('.')[0]
        channel = channel[-3:]
        print('     image identifier: ' + rm_channel)
        print('     channel: ' + channel)
        img_ids.append(rm_channel)
        channels.append(channel)
    fname_id_dict = dict(zip(img_files, img_ids))
    channel_dict = dict(zip(img_files, channels))
    print('')
    print('done extracting image identifiers.')
    if not return_channel:
        return(fname_id_dict)
    if return_channel:
        return((fname_id_dict, channel_dict))
