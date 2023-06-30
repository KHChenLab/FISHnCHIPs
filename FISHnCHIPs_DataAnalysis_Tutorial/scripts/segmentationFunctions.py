import os
import cv2
import skimage.io

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skimage.segmentation import watershed
from cellpose import models
from cv2 import dilate
from skimage import measure


def norm_image_func(image, max, min):
    """" Params: image (2 or 3D arr), maximum intensity and mininum intensity
         Return: normalise image arr by min = 0 and max = 1 (2 or 3D arr)
    """
    image_max = max
    image_min = min
    return (image - image_min) / (image_max - image_min)

def watershed_segmentation(image_arr, cutoff_threshold, opening_threshold, max_val, min_val, 
                           kernel_size=5, filtermask=True, filtermask_size = 3000):
    norm_img = 255 * norm_image_func(image_arr, max_val, min_val)
    img = norm_img.astype(np.uint8)

    _, thresh_img = cv2.threshold(img, cutoff_threshold, 255, cv2.THRESH_BINARY)

    kernel = np.ones((kernel_size, kernel_size), np.uint8)
    opening = cv2.morphologyEx(thresh_img, cv2.MORPH_OPEN, kernel, iterations = opening_threshold)  

    dist_transform = cv2.distanceTransform(opening, cv2.DIST_L2, 3)
    _, sure_fg = cv2.threshold(dist_transform, 0.01*dist_transform.max(),255,0) 
    sure_fg = np.uint8(sure_fg)

    # Marker labelling
    ret, markers = cv2.connectedComponents(sure_fg)
 
    # Add one to all labels so that sure background is not 0, but 1
    markers = markers+1
    markers = watershed(img, markers)

    if filtermask==True:
        print("filtering mask...")
        num, freq = np.unique(markers, return_counts =True)
        for number, frequency in zip(num, freq):
            if number ==1:
                print("removing", frequency)
                index = np.where(markers.flat == number)
                markers.flat[index] = 0
            if int(frequency) < filtermask_size:
                index = np.where(markers.flat == number)
                markers.flat[index] = 0
    print(type(markers))
    return markers

def cellpose_segmentation(imgs,  cellsize, mode, filtermask=True, filtermask_size = 3000):
    model_choice = mode #@param ["Cytoplasm","Cytoplasm2", "Cytoplasm2_Omnipose", "Bacteria_Omnipose", "Nuclei"]

    if model_choice=="Cytoplasm":
        model_type="cyto"

    elif model_choice=="Cytoplasm2":
        model_type="cyto2"

    elif model_choice=="Cytoplasm2_Omnipose":
        model_type="cyto2_omni"

    elif model_choice=="Bacteria_Omnipose":
        model_type="bact_omni"
        cellsize = 0
    elif model_choice=="Nuclei":
        model_type="nuclei"

    else:
        model_type="cyto"

    model = models.Cellpose(gpu=False, model_type=model_type)
    masks, flows, styles, diams = model.eval(imgs, diameter=cellsize, channels=[0,0], do_3D=False)
    masks  = np.array(masks).squeeze()
    if filtermask==True:
        print("filtering mask..")
        num, freq = np.unique(masks, return_counts =True)
        for number, frequency in zip(num, freq):
            if frequency < filtermask_size:
                index = np.where(masks.flat == number)
                masks.flat[index] = 0
    print(type(masks))
    return masks

def dilate_mask_func_watershed(markers, k=15):
    kernel = np.ones((k,k)).astype('uint8')
    dilate_masks = dilate(markers.astype('uint8'), kernel,1)
    return dilate_masks

def dilate_mask_func(markers, k=15):
    kernel = np.ones((k,k)).astype('uint8')
    dilate_masks = dilate(markers, kernel,1)
    return dilate_masks

def show_dilate_mark(dapi_image, dilate_masks, overlay = True, filename = "diate_mask.jpg"):
    #function from skimage to produce data on properties of each nucleus identified in the cell (e.g. centroid, orientation etc.)
    props = measure.regionprops(dilate_masks) 
    # centroid_list = []
    label_list = []
    for region in props:
        label = region['label']
        # centroid = region['centroid']  # y, x
        # centroid_list.append(centroid)
        label_list.append(label)
    # centroid_list = [(t[1], t[0]) for t in centroid_list] # x, y
    
    if overlay:
        # Overlay
        fig = plt.figure()
        ax = plt.gca()
        ax.imshow(dapi_image, cmap='gray', interpolation='none') 
        ax.imshow(dilate_masks, cmap = 'jet', alpha=0.4*(dilate_masks>0), interpolation='none') 
        ax.set_title(str(len(label_list)) + " masks")
        ax.axis('off')
    
    else:
        fig, axes = plt.subplots(1,2)
        ax = axes.ravel()
        ax[0].imshow(dapi_image, vmax= np.percentile(dapi_image, 99.5), cmap = 'gray')
        ax[0].set_title("Image")
        ax[0].axis('off')
        ax[1].imshow(dilate_masks, cmap = 'jet')
        ax[1].set_title(str(len(label_list)) + " masks")
        ax[1].axis('off')
    fig.tight_layout()
    fig.savefig(filename, dpi = 500)
    plt.close(fig)
    

def get_mask_intensity_matrix(All_masks):
    cell_bit_mean_intendict = {}
    cell_bit_sum_intendict = {}
    cell_bit_max_intendict = {}
    cell_bit_median_intendict = {}
    
    for i, celltype in enumerate(All_masks.all_celltype):
        mean_inten_list = []
        sum_inten_list = []
        max_inten_list = []
        inten_median_list = []
        mask_celltype = All_masks.get_by_celltype(celltype).all_masks
        for mask in mask_celltype:
            sum_inten = mask.sum_inten
            max_inten = mask.max_inten
            mean_inten = mask.mean_inten
            medians = mask.median_val
    
            if sum_inten != 0:
                sum_inten_list.append(sum_inten)
                max_inten_list.append(max_inten)
                mean_inten_list.append(mean_inten)
            else:
                sum_inten_list.append(0)
                max_inten_list.append(0)
                mean_inten_list.append(0)
            inten_median_list.append(medians)
    
        cell_bit_sum_intendict[celltype] = sum_inten_list
        cell_bit_max_intendict[celltype] = max_inten_list
        cell_bit_mean_intendict[celltype] = mean_inten_list
        cell_bit_median_intendict[celltype] = inten_median_list
    
    df_mean = pd.DataFrame(cell_bit_mean_intendict)
    df_sum = pd.DataFrame(cell_bit_sum_intendict)
    df_max = pd.DataFrame(cell_bit_max_intendict)
    df_median= pd.DataFrame(cell_bit_median_intendict)
    
    return df_mean, df_sum, df_max, df_median



def get_mask_info(All_masks):
    # save the masks info (x,y, centroid , area, inten, fov, etc.)
    df_fov_list = []
    df_celltype_list = []
    df_imgname_list = []
    df_mean_list = []
    df_max_list = []
    df_sum_list = []
    df_area_list = []
    df_centroid = []
    
    
    for mask in All_masks.all_masks:
        fov = mask.fov
        area = mask.area
        celltype = mask.celltype
        img_name = mask.img_name
    
        if mask.mean_inten != 0:
            mean_inten = mask.mean_inten
        else:
            mean_inten = 0.1
    
        if mask.max_inten != 0:
            max_inten = mask.max_inten
        else:
            max_inten = 0.1
    
        if mask.sum_inten != 0:
            sum_inten = mask.sum_inten
        else:
            sum_inten = 0.1
    
        xycentroid = mask.xycentroid
    
        df_fov_list.append(fov)
        df_celltype_list.append(celltype)
        df_imgname_list.append(img_name)
        df_mean_list.append(mean_inten)
        df_max_list.append(max_inten)
        df_sum_list.append(sum_inten)
        df_area_list.append(area)
        df_centroid.append(xycentroid)
    
    df = {}
    df['Fov'] = df_fov_list
    df['Celltype']= df_celltype_list
    df['Img_name'] = df_imgname_list
    df['Mean_inten']=df_mean_list
    df['Sum_inten']=df_sum_list
    df['Max_inten']=df_max_list
    df['Area']=df_area_list
    df['Centroid']=df_centroid # y, x 
    df_frame = pd.DataFrame(df)
    
    return df_frame


def get_centroids(image_path, mask_prefix, dapi_prefix, fov=None):
    # Input: segmented masks 
    # Ouput: list of centroid positions: x, y 
    # 
    # Added more properties for QC on masks
    # Updated on 26 Mar 2023 by Xinrui 
    mask_img_fn = os.path.join(image_path, ('%s_%s.tif' % (mask_prefix, fov)))
    dapi_img_fn = os.path.join(image_path, ('%s_%s.tif' % (dapi_prefix, fov)))
    mask_image = skimage.io.imread(mask_img_fn)
    dapi_image = skimage.io.imread(dapi_img_fn)
    props = measure.regionprops(mask_image, intensity_image=dapi_image)
    
    centroids = {}
    centroid_list = []
    label_list = []
    area_list =[]
    mean_inten_list = []
    max_inten_list = []
    median_inten_list = []
    sum_inten_list = []
    inten_per_99_list = []
    
    
    for region in props:
        label = region['label']
        centroid = region['centroid'] # y, x = centroid
        area = region['area']
        mean_intensity = region['mean_intensity']
        sum_intensity = region['image_intensity'].sum()
        max_intensity = region['image_intensity'].max()
        median_val = np.median(region['image_intensity'])
        inten_per_99 = np.percentile(region['image_intensity'],99)
        
        
        centroid_list.append(centroid)
        label_list.append(label)
        area_list.append(area)
        
        mean_inten_list.append(mean_intensity)
        max_inten_list.append(max_intensity)
        median_inten_list.append(median_val)
        sum_inten_list.append(sum_intensity)
        inten_per_99_list.append(inten_per_99)
        
    centroid_list = [(t[1], t[0]) for t in centroid_list] # x, y 
    
    if fov is not None:
        fov_list = [fov]*len(centroid_list)
        
    centroids['Centroid'] = centroid_list
    centroids['label'] = label_list
    centroids['FOV'] = fov_list
    centroids['Area'] = area_list
    
    centroids['DAPI_Mean_inten'] = mean_inten_list
    centroids['DAPI_Max_inten'] = max_inten_list
    centroids['DAPI_Median_inten'] = median_inten_list
    centroids['DAPI_sum_inten'] = sum_inten_list
    centroids['DAPI_99pct_inten'] = inten_per_99_list
    
    df_centroids = pd.DataFrame(centroids)
    
    return df_centroids


def get_mask_positions_inFOV(mask_info):
    data_one_celltype = mask_info.loc[mask_info['Celltype'] == mask_info['Celltype'][0]]
    centroid_list =data_one_celltype['Centroid'].to_list()
    list_fov = data_one_celltype['Fov'].to_list()
    list_x = [t[1] for t in centroid_list]
    list_y = [t[0] for t in centroid_list]
    
    return list_fov, list_x, list_y

def get_mask_positions_inFuse(list_fov, list_x, list_y, fov_x, fov_y, image_size = 1648):
    # create position file for whole stitched fovs
    # update for rectangulars @Xinrui 14 March 2023
        
    if fov_y%2 == 0:
        start = fov_x * (fov_y-1)
        size = fov_x
        x_scale = {}
        y_scale = {}
    
        for i in range(fov_y):
            if i%2 ==0:
                for j in range(fov_x):
                    fov = str(start-size*int(i)+j)
                    x_scale[fov] = j
                    y_scale[fov] = i
                    print('row: %i col: %i fov: %s' % (i, j, fov) )
            if i%2 ==1:
                for j in range(1,size+1):
                    fov = str(start-size*(int(i-1))-j)
                    x_scale[fov] = (j-1)
                    y_scale[fov] = i
                    print( '%i %i %s' % (i, j, fov) )
                    
    
    if fov_y%2 == 1:
        start = fov_x * (fov_y-1)
        size = fov_x 
        x_scale = {}
        y_scale = {}
        
        for i in range(fov_y):
            if i%2 ==0:
                for j in range(1,size+1):
                    fov = str(start-size*int(i-1)-j)
                    x_scale[fov] = j
                    y_scale[fov] = i
                    print('row: %i col: %i fov: %s' % (i, j, fov) )
            if i%2 ==1:
                for j in range(size):
                    fov = str(start-size*(int(i))+j)
                    x_scale[fov] = j+1
                    y_scale[fov] = i
                    print( '%i %i %s' % (i, j, fov) )
    
    new_x_list = []
    new_y_list = []
    for fov, x, y in zip(list_fov, list_x, list_y):
        fov_s = str(int(fov))
        new_y = y_scale[fov_s] * image_size + float(y)
        new_x = x_scale[fov_s] * image_size + float(x)
        new_x_list.append(new_x)
        new_y_list.append(new_y)
    
    dict_pos = {}
    dict_pos['x'] = new_x_list
    dict_pos['y'] = new_y_list
    df_pos = pd.DataFrame(dict_pos)
    # df_pos.to_csv(output_path + 'xy_position_' + 'allfovs.csv')
    
    return df_pos
