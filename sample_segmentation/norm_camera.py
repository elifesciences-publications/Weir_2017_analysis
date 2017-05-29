import os
import sys
import argparse
import numpy as np
import re
import skimage.io as io

parser = argparse.ArgumentParser(description = 'Normalize images to camera \
                                 background signal image.')
parser.add_argument('-d', '--img_dir', required = True, 
                    help = 'directory containing images to segment.')

args = parser.parse_args()
img_dir = args.img_dir

def get_img_ids(img_files):
    '''Extracts image filenames lacking wavelength identifiers.'''
    channel_re = re.compile('^w\d+[A-Za-z]*[ .]')
    img_ids = []
    for img in img_files:
        print('_______________________________________________________')
        print('     generating image identifier for ' + img)
        split_im = img.split('_')
        rm_channel = '_'.join([i for i in split_im if re.search(channel_re, i)
                               == None])
        print('     image identifier: ' + rm_channel)
        img_ids.append(rm_channel)
    fname_id_dict = dict(zip(img_files, img_ids))
    print('')
    print('done extracting image identifiers.')
    return(fname_id_dict)

def main(img_dir):
    os.chdir(img_dir)
    if not os.path.isdir(img_dir+'/normalized'):
        os.mkdir(img_dir+'/normalized')
    if not os.path.isdir(img_dir+'/normalized/bg_ims'):
        os.mkdir(img_dir+'/normalized/bg_ims')
    print('Running normalization on ' + img_dir)
    img_files = [i for i in os.listdir(img_dir) if '.tif' in i.lower()]
    print('')
    print('List of image files in the directory:')
    print(img_files)
    print('')
    print('retrieving image identifiers from filenames...')
    fname_id_dict = get_img_ids(img_files)
    print('filename IDs:')
    print(fname_id_dict.values())
    bg_imgs = [img for img in img_files if 'noise' in img.lower()]
    print('')
    print('list of background image files:')
    print(bg_imgs)
    bg_stack = np.empty(shape = (len(bg_imgs),512,512))
    for i in range(0,len(bg_imgs)):
        bg_stack[i,:,:] = io.imread(bg_imgs[i])
    io.imsave(img_dir+'/normalized/bg_ims/bg_stack.tif',bg_stack.astype('uint16'))
    bg_stack_mean = np.mean(bg_stack, dtype = np.float64)
    bg_stack_mean = int(bg_stack_mean)
    print('mean of all background images: ' + str(bg_stack_mean))
    bg_img_means = np.mean(bg_stack, axis = (1,2), dtype = np.float64)
    bg_img_means = bg_img_means.astype(int)
    print('mean of each background image:')
    print(bg_img_means)
    bg_means_diff = bg_img_means - bg_stack_mean
    bg_stack_norm = np.copy(bg_stack)
    print('')
    print('normalizing background images...')
    for i in range(0,len(bg_imgs)):
        bg_stack_norm[i,:,:] = bg_stack_norm[i,:,:] - bg_means_diff[i]
    print('saving normalized camera background image stack...')
    io.imsave(img_dir+'/normalized/bg_ims/bg_stack_norm.tif',bg_stack_norm.astype('uint16'))
    print('calculating mean 2d background image...')
    mean_2d_img = np.mean(bg_stack_norm, axis = 0, dtype = np.float64)
    print('saving mean 2d background image...')
    io.imsave(img_dir+'/normalized/bg_ims/mean_bg_img.tif',mean_2d_img.astype('uint16'))
    print('generating 2d difference-from-mean image...')
    diff_from_2d_mean = mean_2d_img - np.mean(mean_2d_img,
                                              dtype = np.float64).astype(int)
    print('saving 2d difference-from-mean image...')
    io.imsave(img_dir+'/normalized/bg_ims/diff_from_mean_bg_img.tif',
              diff_from_2d_mean.astype('int16'))
    bg_ids = []
    print('')
    print('generating image id:background mean difference dict...')
    for img in bg_imgs:
        bg_ids.append(fname_id_dict[img])
    bg_mean_dict = dict(zip(bg_ids, bg_means_diff.tolist()))
    print('image id:background mean difference dict:')
    print(bg_mean_dict)
    print('')
    print('normalizing images...')
    for f in img_files:
        raw_img = io.imread(f)
        print('     normalizing ' + f)
        print('     background normalizing factor: ' +
              str(bg_mean_dict[fname_id_dict[f]]))
        norm_img = raw_img - bg_mean_dict[fname_id_dict[f]] - diff_from_2d_mean
        io.imsave('normalized/'+f,norm_img.astype('uint16'))

if __name__ == '__main__':
    main(img_dir)
