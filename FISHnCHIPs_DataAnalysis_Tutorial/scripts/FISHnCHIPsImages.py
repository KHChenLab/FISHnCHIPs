#!/usr/bin/env python
# coding: utf-8
"""
Module for FISHnCHIPs image processing
This script is for FISHnCHIPs image processing. It contains classes and functions for image segmentation, registration, and mask creation. 
The script can segment DAPI images and register all other images to DAPI. 
It can also segment bleach-subtracted images and create masks for each cell. 
The masks contain information such as cell type, image name, FOV, area, centroid, and intensity. 


@updated: 23Mar2023 - Xinrui
Update the processing:
    1. For each DAPI FOV, segment on contrasted JPG 
        Default: vmax at 99% percentile & manual adjust specific FOVs
    2. For all hyb images, register to DAPI, and keep records of all shifts 
        Remove the outliers and get the median / mean of the shifts 
        Apply THE SHIFT to all FOVs of hyb images and bleach images 
    3. Overlay hyb images on DAPI to check the registration is reasonable 
    4. Same strategy to apply masks, and parse masks locations and intensity

"""
 
import os
import csv
import pickle
import skimage.io
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from scripts.colorizeImages import make_composite

from skimage import measure
from tifffile import imwrite                
from scripts.segmentationFunctions import watershed_segmentation, cellpose_segmentation
from scripts.segmentationFunctions import dilate_mask_func, show_dilate_mark, dilate_mask_func_watershed
from scripts.registerFunction import register_slice, register_slice_with_shift

matplotlib.use('Qt5Agg')

class yml_class:

    def __init__(self, user_dict, output_path, **kwargs):
        self.output_path = output_path
        self.num_fov = user_dict['fovs']
        self.mainpath = user_dict['mainpath']
        self.schema = user_dict['schema']
        self.dapiname = user_dict['dapiname']
        self.overlap_in_px = user_dict['overlap_in_px']
        self.segmentation_type = user_dict['segmentation']
        self.segment_DAPI_format = user_dict['segment_DAPI_format']
        
        self.imagename_dict = self._getschemadict()
        self._createoutputpath_full()
        self.bleach_namedict = self._get_bleach_imagename_dict()
        for keys, values in kwargs.items():
            setattr(self, keys, values)

    def _getschemadict(self):
        self.imagename_dict = {}
        with open(self.schema , mode='r', encoding='utf-8-sig') as file:
            reader = csv.reader(file)
            ### Format of the schema: dye_hyb, moduleName ---------------
            # ----------------- Example ---------------------------------
            # Cy7_00, Module1
            # Cy7_01, Module2
            # Alexa594_00, Module3
            # -----------------------------------------------------------
            self.imagename_dict = {rows[0]:rows[1] for rows in reader} 
        return self.imagename_dict

    def _createoutputpath_full(self):
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)
            
    def _get_bleach_imagename_dict(self):
        self.bleach_namedict = {}
        for dye_hyb in self.imagename_dict:
            self.bleach_namedict[dye_hyb] = {}
            
            channel_list = dye_hyb.split('_')
            dye = channel_list[0]
            hyb = channel_list[1]
            self.bleach_namedict[dye_hyb]= '%s_Bleach_%s' % (dye, hyb)

        return self.bleach_namedict
    
    def get_dapi_image_by_fov(self, fov):
        # Example of DAPI image filename: DAPI_0_99.tif
        self.dapiname_fov = os.path.join(self.mainpath, self.dapiname + fov + self.segment_DAPI_format)
        print(self.dapiname_fov)
        return skimage.io.imread(self.dapiname_fov) #loads image from filename
    
    def get_image_by_fov(self, fov, dye_hyb):
        # Example of image filename: Cy5_00_99.tif
        imagename_fov = os.path.join(self.mainpath, '%s_%s.tif' % (dye_hyb , fov))
        print(imagename_fov)
        return skimage.io.imread(imagename_fov)
    
    def get_bleach_image_by_fov(self, fov, dye_hyb):
        # Exampel of bleach image filename: Cy5_Bleach_00_99.tif
        channel_list = dye_hyb.split('_')
        bleach_image_fn = os.path.join(self.mainpath, '%s_Bleach_%s_%s.tif' % (channel_list[0], channel_list[1], fov))
        
        if os.path.exists(bleach_image_fn):
            return skimage.io.imread(bleach_image_fn)
        else:
            print('Bleach file not exist. ')
            return None
    

class image_class(yml_class):
    def __init__(self, image_arr, fov, user_dict, outputpathold, **kwargs):
        super().__init__(user_dict, outputpathold)
        self.outputpathold = outputpathold
        self.image_fov = image_arr
        self.fov = fov
        self.user_dict = user_dict
        self.output_path = self.output_path
    
    def remove_n_px(self, n_px, verbose=True):
        new_image = self.image_fov[n_px: (2048-n_px), n_px: (2048-n_px)]
        if verbose:
            print("Cropping the edge of image. Final image shape is", new_image.shape)
        return image_class(new_image, self.fov, self.user_dict, self.outputpathold)
    
    def register_to_dapi(self, dapi_img):
        self.image_fov, shifts = register_slice(dapi_img, self.image_fov)
        return self.image_fov, shifts 
    
    def apply_shifts(self, shifts):
        shifted_img = register_slice_with_shift(self.image_fov, shifts)
        return image_class(shifted_img, self.fov, self.user_dict, self.outputpathold)
    
    
    def register_to_dapi_and_remove_background(self, dapi_img, dye_hyb, n_px):
        subtracted_fn = self.output_path + dye_hyb + '_subtracted_' + self.fov + '.tif'
        
        if os.path.exists(subtracted_fn):
            removebg_img = skimage.io.imread(subtracted_fn)
        else:
            self.image_fov, shifts = register_slice(dapi_img, self.image_fov)
            
            bleach_image_fn = os.path.join(self.mainpath, '%s_%s.tif' % (self.bleach_namedict[dye_hyb], self.fov))
            if os.path.exists(bleach_image_fn):
                bleach_image = skimage.io.imread(self.mainpath + self.bleach_namedict[dye_hyb] + '_' + self.fov + '.tif')
                bleach_image = register_slice_with_shift(bleach_image[n_px: (2048-n_px), n_px: (2048-n_px)], shifts)
                subtracted = self.image_fov.astype(np.int32) - bleach_image.astype(np.int32)
            else:
                # if there is no bleach image, subtract nothing.
                subtracted = self.image_fov.astype(np.int32)
            removebg_img = np.clip(subtracted, 0, subtracted.max()).astype(np.uint16)
            imwrite(subtracted_fn, removebg_img)
        
        return image_class(removebg_img, self.fov, self.user_dict, self.outputpathold)

      
    def segmentation(self, max_val = None, min_val = None):
        if self.segmentation_type == 'watershed':
            self.kernel_size = self.user_dict['kernel_size']
            self.cutoff_threshold = self.user_dict['cutoff_threshold']
            self.opening_threshold = self.user_dict['opening_threshold']
            self.filtermask = self.user_dict['filtermask']
            self.filtermask_size = self.user_dict['filtermask_size']
            self.markers = watershed_segmentation(self.image_fov, self.cutoff_threshold, self.opening_threshold, max_val, min_val, 
                                                  kernel_size = self.kernel_size,  filtermask=self.filtermask, filtermask_size = self.filtermask_size)
            
        elif self.segmentation_type == 'cellpose':
            self.cellsize = self.user_dict['cellsize']
            self.mode = self.user_dict['segment_mode']
            self.filtermask = self.user_dict['filtermask']
            self.filtermask_size = self.user_dict['filtermask_size']
            self.markers =  cellpose_segmentation(self.image_fov, self.cellsize, self.mode,  
                                                  filtermask = self.filtermask, 
                                                  filtermask_size = self.filtermask_size)

        else:
            print('Did you type anything wrong ?')
        return self.markers


    def create_props_mask(self, markers, celltype):
        cellmask_list = []
        number_of_mark = np.unique(markers).shape[0]
        print("Number of markers:" + str(number_of_mark))

        if number_of_mark == None:
            self.inten_image_in_mask_dict = None
        props_mask = measure.regionprops(markers,  intensity_image=self.image_fov)

        for region in props_mask:
            mean_intensity = region['mean_intensity']
            sum_intensity = region['image_intensity'].sum()
            area = region['area']
            max_intensity = region['image_intensity'].max()
            label = region['label']
            xycentroid = region['centroid'] # y, x = centroid
            inten_per_99_9 = np.percentile(region['image_intensity'],99.9)
            inten_per_99_5 = np.percentile(region['image_intensity'],99.5)
            inten_per_99_95 = np.percentile(region['image_intensity'],99.95)
            median_val = np.median(region['image_intensity'])

            img_name = self.imagename_dict[celltype]
            cellmask_list.append(Cell_mask(celltype, img_name, self.fov, 
                                           max_intensity,mean_intensity,sum_intensity,
                                           area,label,xycentroid, 
                                           inten_per_99_9 = inten_per_99_9,
                                           inten_per_99_95 = inten_per_99_95, 
                                           inten_per_99_5 = inten_per_99_5, 
                                           median_val = median_val ))
        return cellmask_list

class Cell_mask:
    def __init__(self, celltype, img_name, fov, max_inten, mean_inten, sum_inten, area, label, xycentroid, **kwargs):
        self.celltype = celltype
        self.img_name = img_name
        self.fov = fov
        self.max_inten = max_inten
        self.mean_inten = mean_inten
        self.sum_inten = sum_inten
        self.area = area
        self.label = label
        self.xycentroid = xycentroid
        for keys, values in kwargs.items():
            setattr(self, keys, values)

    def __repr__(self):
        return "Mask"
    

class Cell_masks:
    def __init__(self, list_of_cellmask, **kwargs):
        self.num_masks = len(list_of_cellmask)
        self.all_masks = list_of_cellmask
        self.all_fov = self._get_number_fov()
        self.all_celltype = self._get_celltype_name()

    def __repr__(self):
        return "Masks"

    def get_by_celltype(self, celltype):
        mask_by_celltype = []
        for mask in self.all_masks:
            if mask.celltype == celltype:
                mask_by_celltype.append(mask)
        return Cell_masks(mask_by_celltype)

    def get_by_fov(self, fov):
        mask_by_fov = []
        for mask in self.all_masks:
            if mask.fov == fov:
                mask_by_fov.append(mask)
        return Cell_masks(mask_by_fov)

    def _get_number_fov(self):
        number_fov_list = []
        for mask in self.all_masks:
            number_fov_list.append(mask.fov)
        number_fov = np.unique(number_fov_list)
        return list(number_fov)

    def _get_celltype_name(self):
        celltype_list = []
        for mask in self.all_masks:
            celltype_list.append(mask.celltype)
        celltypes = np.unique(celltype_list)
        return list(celltypes)


### ------------------------------------------------------------------------------------------------------------
def register_one_slice(user_dict, output_path, dye_hyb, fov, save_registered_img = True):
    print('Registering --- %s \t %s' % (dye_hyb, fov))
    image_obj_0 = yml_class(user_dict, output_path)
    overlap_in_px = user_dict['overlap_in_px']
    dapi_fn = os.path.join(output_path, 'DAPI_rm%i_%s.tif' % (overlap_in_px, str(fov)))
    if os.path.exists(dapi_fn):
        img_dapi_fov = skimage.io.imread(dapi_fn)
    else:
        img_dapi = image_class(image_obj_0.get_dapi_image_by_fov(fov), fov, user_dict, output_path)
        img_dapi = img_dapi.remove_n_px(n_px = overlap_in_px)
        img_dapi_fov = img_dapi.image_fov
        imwrite(dapi_fn, img_dapi_fov)
    
    img_hyb = image_class(image_obj_0.get_image_by_fov(fov, dye_hyb), fov, user_dict, output_path)
    img_hyb = img_hyb.remove_n_px(n_px = overlap_in_px)
            
    registered, shifts = img_hyb.register_to_dapi(img_dapi_fov)
    registered = np.clip(registered, 0, registered.max()).astype(np.uint16)
    
    if save_registered_img:
        registered_fn = os.path.join(output_path, '%s_%s_reg.tif' % (dye_hyb, fov)) 
        imwrite(registered_fn, registered)
    
    shifts_dict = {'FOV': [fov], 'shift_x': [shifts[0]], 'shift_y': [shifts[1]]}
    shifts_df = pd.DataFrame(shifts_dict)

    return shifts_df


def register_one_module(user_dict, output_path, dye_hyb, save_registered_img=True):
    ### Register one module all fovs to DAPI, calculate average shifts and apply to all 
    print(dye_hyb)
    overlap_in_px = user_dict['overlap_in_px']
    
    # Initialization
    image_obj_0 = yml_class(user_dict, output_path)
    ndigits = len(str(user_dict['fovs'])) # 2, 3 etc 
    fov_list = [str(x).zfill(ndigits) for x in range(user_dict['fovs'])]
    
    shifts_dict = {}
    for fov in fov_list:
        
        dapi_fn = os.path.join(output_path, 'DAPI_rm%i_%s.tif' % (overlap_in_px, str(fov)))
        if os.path.exists(dapi_fn):
            img_dapi_fov = skimage.io.imread(dapi_fn)
        else:
            img_dapi = image_class(image_obj_0.get_dapi_image_by_fov(fov), fov, user_dict, output_path)
            img_dapi = img_dapi.remove_n_px(n_px = overlap_in_px)
            img_dapi_fov = img_dapi.image_fov
            imwrite(dapi_fn, img_dapi_fov)
        
        img_hyb = image_class(image_obj_0.get_image_by_fov(fov, dye_hyb), fov, user_dict, output_path)
        img_hyb = img_hyb.remove_n_px(n_px = overlap_in_px)
                
        registered, shifts = img_hyb.register_to_dapi(img_dapi_fov)
        registered = np.clip(registered, 0, registered.max()).astype(np.uint16)
        
        shifts_dict[fov] = shifts
        
        if save_registered_img:
            registered_fn = os.path.join(output_path, '%s_%s_reg.tif' % (dye_hyb, fov)) 
            imwrite(registered_fn, registered)
    
    shifts_dict_fn = os.path.join(output_path, '%s_reg_shifts.csv' % dye_hyb)
    with open(shifts_dict_fn, 'w') as f:
        f.write('FOV,shift_x,shift_y\n')
        for fov in shifts_dict.keys():
            f.write('%s,%f,%f\n' % (str(fov), shifts_dict[fov][0], shifts_dict[fov][1]))
        
            
    return shifts_dict

def register_one_module_TOOSlow(user_dict, output_path, dye_hyb, save_registered_img=True):
    ### Register one module all fovs to DAPI, calculate average shifts and apply to all 
    print(dye_hyb)
    overlap_in_px = user_dict['overlap_in_px']
    
    # Initialization
    image_obj_0 = yml_class(user_dict, output_path)
    ndigits = len(str(user_dict['fovs'])) # 2, 3 etc 
    fov_list = [str(x).zfill(ndigits) for x in range(user_dict['fovs'])]
    
    shifts_dict = {}
    for fov in fov_list:
        
        dapi_fn = os.path.join(output_path, 'DAPI_rm%i_%s.tif' % (overlap_in_px, str(fov)))
        if os.path.exists(dapi_fn):
            img_dapi_fov = skimage.io.imread(dapi_fn)
        else:
            img_dapi = image_class(image_obj_0.get_dapi_image_by_fov(fov), fov, user_dict, output_path)
            img_dapi = img_dapi.remove_n_px(n_px = overlap_in_px)
            img_dapi_fov = img_dapi.image_fov
            imwrite(dapi_fn, img_dapi_fov)
        
        img_hyb = image_class(image_obj_0.get_image_by_fov(fov, dye_hyb), fov, user_dict, output_path)
        img_hyb = img_hyb.remove_n_px(n_px = overlap_in_px)
                
        registered, shifts = img_hyb.register_to_dapi(img_dapi_fov)
        registered = np.clip(registered, 0, registered.max()).astype(np.uint16)
        
        shifts_dict[fov] = shifts
        
        if save_registered_img:
            registered_fn = os.path.join(output_path, '%s_%s_reg.tif' % (dye_hyb, fov)) 
            imwrite(registered_fn, registered)
    
    shifts_dict_fn = os.path.join(output_path, '%s_reg_shifts.csv' % dye_hyb)
    with open(shifts_dict_fn, 'w') as f:
        f.write('FOV,shift_x,shift_y\n')
        for fov in shifts_dict.keys():
            f.write('%s,%f,%f\n' % (str(fov), shifts_dict[fov][0], shifts_dict[fov][1]))
        
            
    return shifts_dict

def plot_shifts_distribution(shifts_x, shifts_y, output_file = './distribution_shifts.jpg'):
        
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 3))
    ax1.hist(shifts_x, label='shifts_x', bins=20)
    ax1.set_title('shifts_x')
    ax2.hist(shifts_y, label='shifts_y', bins=20)
    ax2.set_title('shifts_y')
    plt.savefig(output_file, bbox_inches='tight', dpi=120)
    plt.close()


def get_shifts_in_common(shifts_x, shifts_y, max_shift = 50, plot_distribution=True):
    # shifts_x = []
    # shifts_y = []
    # for item in shifts_dict:
    #     shifts_x.append(shifts_dict[item][0])
    #     shifts_y.append(shifts_dict[item][1])
            
    shifts_x = [x for x in shifts_x if abs(x)<=max_shift]
    shifts_y = [y for y in shifts_y if abs(y)<=max_shift]
    
    return np.mean(shifts_x), np.mean(shifts_y)

def process_one_slice(user_dict, img_path, dapi_path, dye_hyb, fov, shifts, 
                      subtract_bleach = True, save_overlay = True):
    # Initialization
    overlap_in_px = user_dict['overlap_in_px']
    image_obj_0 = yml_class(user_dict, img_path)
    
    
    # remove edges and shift 
    img_hyb = image_class(image_obj_0.get_image_by_fov(fov, dye_hyb), fov, user_dict, img_path)
    img_hyb = img_hyb.remove_n_px(n_px = overlap_in_px)
    shifted_img_hyb = register_slice_with_shift(img_hyb.image_fov, shifts)
    
    if subtract_bleach:
        # bleach image
        bleach_img_fn = os.path.join(image_obj_0.mainpath, '%s_%s.tif' % (image_obj_0.bleach_namedict[dye_hyb], fov))
        bleach_img = skimage.io.imread(bleach_img_fn)
        shifted_bleach_img = register_slice_with_shift(bleach_img[overlap_in_px: (2048-overlap_in_px), overlap_in_px: (2048-overlap_in_px)], shifts)
        subtracted = shifted_img_hyb.astype(np.int32) - shifted_bleach_img.astype(np.int32)
        subtracted = np.clip(subtracted, 0, subtracted.max()).astype(np.uint16)
    else:
        subtracted = shifted_img_hyb
        
    processed_fn = os.path.join(img_path, '%s_processed_%s.tif' % (dye_hyb, fov))
    imwrite(processed_fn, subtracted)
    
    if save_overlay:
        dapi_fn = os.path.join(dapi_path, 'DAPI_rm%i_%s.tif' % (overlap_in_px, fov))
        dapi_img = skimage.io.imread(dapi_fn)
        colors = [(0, 0, 1), (1, 1, 1)] # blue and white
        composite = make_composite([dapi_img, subtracted], colors)
        composite_fn = os.path.join(img_path, 'overlay_DAPI_%s_%s.jpg' % (dye_hyb, fov))
        plt.imshow(composite)
        plt.axis(False)
        plt.gca().set_aspect('equal')
        plt.tight_layout()
        plt.savefig(composite_fn, bbox_inches='tight', dpi = 300)
    
    return image_class(subtracted, fov, user_dict, img_path)
    

### ------------------------------------------------------------------------------------------------------------


def segment_dapi_one_fov(user_dict, output_path, fov, min_image=None, max_image=None):
    overlap_in_px = user_dict['overlap_in_px']
    image_obj_0 = yml_class(user_dict, output_path)
    
    img_dapi = image_class(image_obj_0.get_dapi_image_by_fov(fov), fov, user_dict, output_path)
    img_dapi = img_dapi.remove_n_px(n_px = overlap_in_px)
    dapi_fn = output_path + 'DAPI_rm%i_fov_%s.tif' % (overlap_in_px, str(fov))
    
    if not os.path.exists(dapi_fn):
        imwrite(dapi_fn, img_dapi.image_fov)
    
    markers = img_dapi.segmentation(max_val = max_image, min_val = min_image)
    if user_dict['segmentation'] == 'cellpose':
        dilate_mask = dilate_mask_func(markers, k=user_dict['mask_dilate_factor'])
    elif user_dict['segmentation'] == 'watershed':
        dilate_mask = dilate_mask_func_watershed(markers, k=user_dict['mask_dilate_factor'])
    
    mask_fn = output_path + 'dilateMask_fov_' + fov + '.tif'
    imwrite(mask_fn, dilate_mask)
    
    mask_jpg = output_path + '/dapi' + '_overlay_mask_fov_' + fov + '.jpg'
    show_dilate_mark(img_dapi.image_fov, dilate_mask, overlay = True, filename = mask_jpg)
    
    
    return True


def segment_one_fov(user_dict, output_path, fov):
    overlap_in_px = user_dict['overlap_in_px']
    list_all_cellmask = []
    image_obj_0 = yml_class(user_dict, output_path)
    celltype_to_analyse = image_obj_0.imagename_dict.keys()
    print(str(len(celltype_to_analyse)) + " Modules to analyse")
    
    img_dapi = image_class(image_obj_0.get_dapi_image_by_fov(fov), fov, user_dict, output_path)    
    img_dapi = img_dapi.remove_n_px(n_px = overlap_in_px)
    dapi_rm_fn = output_path + 'DAPI_rm%i_fov_%s.tif' % (overlap_in_px, str(fov))
    
    if not os.path.exists(dapi_rm_fn):
        imwrite(dapi_rm_fn, img_dapi.image_fov)
        
    filename_dilatemask = output_path + 'dilateMask_fov_' + str(fov) + '.tif'
    if not os.path.exists(filename_dilatemask):
        markers = img_dapi.segmentation()
        print("segmentation... Done")
        #after segmentation, it might shrink the thing that we are segmenting in size --> dilate
        
        if user_dict['segmentation'] == 'cellpose':
            dilate_mask = dilate_mask_func(markers, k=user_dict['mask_dilate_factor'])
        elif user_dict['segmentation'] == 'watershed':
            dilate_mask = dilate_mask_func_watershed(markers, k=user_dict['mask_dilate_factor'])
        imwrite(filename_dilatemask, dilate_mask)
        show_dilate_mark(img_dapi.image_fov, dilate_mask, overlay = True, 
                         filename = output_path + '/dapi' + '_overlay_mask_fov_' + fov + '.jpg')
    else:
        dilate_mask = skimage.io.imread(filename_dilatemask)
        print("Loading segmentation.. Done")
        
    for celltype in celltype_to_analyse:
        cell_image = image_class(image_obj_0.get_image_by_fov(fov, celltype), fov, user_dict, output_path).remove_n_px(n_px = overlap_in_px)
        cell_image = cell_image.register_to_dapi_and_remove_background(img_dapi.image_fov, celltype, overlap_in_px)
        plt.close('all')
        list_cell_mask = cell_image.create_props_mask(dilate_mask, celltype)
        list_all_cellmask = list_all_cellmask  + list_cell_mask
        plt.close('all')
        
    return list_all_cellmask


def subtract_one_image(user_dict, output_path, celltype, fov):
    overlap_in_px = user_dict['overlap_in_px']
    image_obj_0 = yml_class(user_dict, output_path)
    
    img_dapi_fn = output_path + 'DAPI_rm%i_fov_%s.tif' % (overlap_in_px, str(fov))
    if os.path.exists(img_dapi_fn):
        image_dapi = skimage.io.imread(img_dapi_fn)
    else:
        img_dapi = image_class(image_obj_0.get_dapi_image_by_fov(fov), fov, user_dict, output_path)
        img_dapi = img_dapi.remove_n_px(n_px = overlap_in_px)
        image_dapi = img_dapi.image_fov
    celltype_img = image_class(image_obj_0.get_image_by_fov(fov, celltype), fov, user_dict, output_path).remove_n_px(n_px = overlap_in_px)
    celltype_img = celltype_img.register_to_dapi_and_remove_background(image_dapi, celltype, overlap_in_px)
    
    return celltype_img
    
def load_masks(image_path, mask_path, user_dict, dye_hyb, fov):
    fn_mask = os.path.join(mask_path, 'dilateMask_fov_%s.tif' % fov)
    dilate_mask = skimage.io.imread(fn_mask)
    print('loading masks for %s' % (fov))
    
    fn_image = os.path.join(image_path, '%s_adjusted_%s.tif' % (dye_hyb, fov))
    image = skimage.io.imread(fn_image)
    
    cell_image = image_class(image, fov, user_dict, image_path)
    list_cell_mask = cell_image.create_props_mask(dilate_mask, dye_hyb)
    
    return list_cell_mask






### Need to update 
def segment_img_by_fov(user_dict, fov_list, output_path, max_image = None, min_image = None, filename_all_cellmask = None):
    ### Currently only used for hCRC data analyis
    
    # Initialization 
    image_obj_0 = yml_class(user_dict, output_path)
    print('Image to anchor: ' + user_dict['anchor_name'])
    print('celltype: ' + user_dict['anchor_celltype'])
    overlap_in_px = user_dict['overlap_in_px']
    if(filename_all_cellmask is None):
        list_all_cellmask = []
        filename_all_cellmask = output_path + 'list_all_cellmask.pkl'
    else:
        list_all_cellmask = pickle.load(filename_all_cellmask)
        
    # Find minimum and maximum intensity of stitched image
    fusedpath = user_dict['fusedpath']
    max_image = np.max(skimage.io.imread(fusedpath))
    min_image = np.min(skimage.io.imread(fusedpath))

    
    for fov in fov_list:
        # Segment bleach subtracted images 
        print("processing ..", fov)
        img_anchor= image_class(image_obj_0.get_anchor_image_by_fov(fov, anchor_name=user_dict['anchor_name']), fov, user_dict, output_path)
        img_anchor = img_anchor.remove_n_px(n_px = overlap_in_px)        
        img_bleach = image_class(image_obj_0.get_bleach_image_by_fov(fov, anchor_name=user_dict['anchor_name']), fov, user_dict, output_path)
        
        if img_bleach is not None:
            img_bleach = img_bleach.remove_n_px(n_px = overlap_in_px)
            subtracted = img_anchor.image_fov.astype(np.int32) - img_bleach.image_fov.astype(np.int32)
        else:
            subtracted = img_anchor.image_fov.astype(np.int32)
            
        removebg_img = np.clip(subtracted, 0, subtracted.max()).astype(np.uint16)
        
        img_subtracted = image_class(removebg_img, fov, user_dict, output_path)
        markers = img_subtracted.segmentation(max_val = max_image, min_val = min_image)
        print("segmentation.. Done")
        if(user_dict['segmentation'] == 'watershed'):
            print(user_dict['segmentation'] + " ---- no force dilation ---------")
            dilate_mask = markers.astype(np.int16)
        else:
            dilate_mask = dilate_mask_func(markers, k=user_dict['mask_dilate_factor'])
        
        show_dilate_mark(img_anchor.image_fov, dilate_mask, overlay = True, 
                         filename = output_path + user_dict['anchor_celltype'] + '_overlay_mask_fov_' + fov + '.jpg')
        
        
        list_all_cellmask = list_all_cellmask + img_anchor.create_props_mask(dilate_mask, user_dict['anchor_celltype'])
        plt.close('all')
        

    return list_all_cellmask


