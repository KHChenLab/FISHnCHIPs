# [FISHnCHIPs: Fluorescence In Situ Hybridization of Cellular HeterogeneIty and Gene Expression Programs](https://www.biorxiv.org/content/10.1101/2023.04.11.536345v1)

## About
We present FISHnCHIPs for highly sensitive in situ profiling of cell types and gene expression programs. FISHnCHIPs achieves this by simultaneously imaging ∼2-35 co-expressed genes that are spatially co-localized in tissues, resulting in similar spatial information as single-gene FISH, but at ∼2-20-fold higher sensitivity. See (https://www.biorxiv.org/content/10.1101/2023.04.11.536345v1). This software guides users to design and evaluate their own gene panel using input from scRNA-seq. We also provide the FISHnCHIPs data analysis pipeline that we used to process the raw FISHnCHIPs data and cluster the module-cell matrix to define cell types.

![Workflow of panel design](/assets/images/FISHnCHIPs-Workflow-PanelDesign_v3.jpg)

## Contents
* [Prerequisites](https://github.com/ChengDHow/FISHnCHIPs#prerequisites)
    * [System requirements](https://github.com/ChengDHow/FISHnCHIPs#system-requirements)
    * [Software requirements](https://github.com/ChengDHow/FISHnCHIPs#software-requirements)
    * [Tested dependencies](https://github.com/ChengDHow/FISHnCHIPs#tested-dependencies)
* [Getting started](https://github.com/ChengDHow/FISHnCHIPs#getting-started)
    * [Installation](https://github.com/ChengDHow/FISHnCHIPs#installation)
    * [Gene Panel Design Demo](https://github.com/ChengDHow/FISHnCHIPs#gene-panel-design-demo)
        * [Generation of panel design](https://github.com/ChengDHow/FISHnCHIPs#generation-of-panel-design)
        * [Panel Evaluation](https://github.com/ChengDHow/FISHnCHIPs#panel-evaluation)
    * [Data Analysis Demo](https://github.com/ChengDHow/FISHnCHIPs#data-analysis-demo)
        * [Generation of cell mask images](https://github.com/ChengDHow/FISHnCHIPs#generation-of-cell-mask-images)
        * [Calculate cell positions & brightness intensity](https://github.com/ChengDHow/FISHnCHIPs#calculate-cell-positions--brightness-intensity)
    * [FISHnCHIPs manuscript figure](https://github.com/ChengDHow/FISHnCHIPs#fishnchips-manuscript-figure)
* [Parameter settings](https://github.com/ChengDHow/FISHnCHIPs#parameter-settings)
* [Citations and code dependencies](https://github.com/ChengDHow/FISHnCHIPs#citation-and-code-dependencies)
* [Authors](https://github.com/ChengDHow/FISHnCHIPs#authors)
* [License](https://github.com/ChengDHow/FISHnCHIPs#license)
* [Acknowledgements](https://github.com/ChengDHow/FISHnCHIPs#license) 

## Prerequisites
### System requirements:
A computer that can run Python and/or R, with at least 16 GB of RAM. No non-standard hardware is required.

### Software requirements:
* Python 3.9.16
* pandas 1.5.2
* numpy 1.23.5
* matplotlib 3.6.2
* numba 0.56.4
* seaborn 0.12.2
* scipy 1.10.0

***For Gene Panel Design:***
* FISHnCHIPs 0.1.0 gene panel design package
* igraph 0.10.4
* leidenalg 0.9.1

***For Image Processing:***
* FISHnCHIPs 0.1.0 image processing package
* cellpose 2.2
* scikit-image 0.19.3

### Tested dependencies
***Packages common to both gene panel design and data processing:***
* Python 3.9.16
* pandas 1.5.2
* numpy 1.23.5
* matplotlib 3.6.2
* numba 0.56.4
* seaborn 0.12.2
* scipy 1.10.0
* scanpy 1.9.2

***For Gene Panel Design:***
* igraph 0.10.4
* leidenalg 0.9.1

***For Image Processing:***
* cellpose 2.2
* scikit-image 0.19.3

**For manuscript figures (R packages):**
* Seurat 4.3.0
* tidyverse 2.0.0
* dittoSeq 1.10.0
* ComplexHeatmap 2.14.0
* ggpubr 0.6.0
* ggrepel 0.9.3



## Getting started
### Installation
1. Download and install [Anaconda](https://www.anaconda.com/products/distribution#download-section).
2. Spyder and the common packages will be installed with Anaconda. Installation of **Anaconda** typically takes less than 30 minutes.
3. For image processing, install [cellpose](https://www.cellpose.org/) and [scikit-image](https://scikit-image.org/) using either pip or conda install:
    **Cellpose**
    * pip install cellpose
    * conda install -c conda-forge cellpose
    **scikit-image**
    * pip install scikit-image
    * conda install -c anaconda scikit-image
4. For gene panel design, install [igraph](https://python.igraph.org/en/stable/) and [leidenalg](https://leidenalg.readthedocs.io/en/stable/install.html) using either pip or conda install:
    **igraph**
    * pip install igraph
    * conda install -c conda-forge python-igraph
    **leidenalg**
    * pip install leidenalg
    * conda install -c conda-forge leidenalg
5. Code has been tested to run on Spyder and Jupyter Notebook but also can also be run from any other Python 3 compatible IDE.

### <ins>Gene Panel Design Demo</ins>
The tutorial for the workflow of gene panel design and evaluation are provided as a jupyter notebook file **`Gene panel design tutorial.ipynb`** in the **FISHnCHIPs_GenePanelDesign_Tutorial** folder.
In the tutorial, various functions were called from the **`FISHnCHIPS_0.1.0.py`** package in the **scripts folder**.

### Input files provided (csv format):
* Correlation matrix
* Gene-Cell expression matrix
* Reference markers
* scRNA celltype

### Functions:
#### <ins>Generation of panel design</ins>
1. **`get_panel`** (used for both cell-centric and gene-centric panel design)   
    * Inputs: 
        * correlation matrix
        * reference markers
    * Hyperparameters: 
        * correlation threshold
        * minimum number of genes required from each marker
    <p> Function returns a DataFrame containing the genes selected for the panel based on the hyperparameters, the correlation of each gene with the reference marker gene and the cluster that it belongs to.</p>

2. **`get_filtered_genes`** (only for gene-centric panel design)
    * Inputs: 
        * Gene-cell expression matrix
    * Hyperparameters: 
        * Minimum gene expression level
        * Minimum number of cells that gene is present in
        * Maximum number of cells that gene is present in 
        * Name pattern of genes to exclude from panel
    <p> Function for prefiltering the genes with low gene expression level, specific naming patterns, minimum and maximum number of cells that gene is present in to ensure that they are adequately expressed genes.</p>
    <p>Returns the list of genes that passes the filtering requirements.</p>

3. **`leiden_corr_matrix`** (only for gene-centric panel design)
    * Inputs: 
        * Post-filtered correlation matrix
    * Hyperparameters: 
        * Edge correlation threshold
        * Partition type (Default: ModularityVertexPartition)
    <p>As the reference marker file is not available for gene-centric panel design, the clustering of genes will have to be conducted using algorithms based on the gene-gene correlation.</p>
    <p>Function uses the leiden algorithm to cluster the genes and removes genes that have correlation lower than the threshold with all other genes, returning a tuple containing a dataframe of the cluster that each selected gene belongs to and the cluster network graph of the selected genes.</p>

#### <ins>Panel Evaluation</ins>

4. **`get_cumulative_signals`** (only for cell-centric panel evaluation)
    * Inputs: 
        * Gene-cell expression matrix
        * Panel design information
        * scRNA celltype
    * Hyperparameters: 
        * Usage of multiple probes
    <p>Function calculates the signal gain, conserved signal gain and Signal Specificity Ratio (SSR) of each gene using the gene-cell expression matrix with genes selected during the panel design and returning the results as a DataFrame.</p>

5. **`get_cell_bit_matrix`** (only for gene-centric panel evaluation)
    * Inputs: 
        * Gene-cell expression matrix
        * Panel design information
    * Hyperparameters: 
        * Name of column indicating gene cluster
        * Usage of multiple probes
    <p>Function returns a DataFrame that contains the expression level of each cluster of genes in each cell (cell-bit matrix).</p>

6. **`evaluate_gene_centric_panel`** (only for gene-centric panel evaluation)
    * Inputs: 
        * Gene-cell expression matrix
        * Panel design information
        * Cell-bit matrix
    * Hyperparameters: 
        * Signal threshold
    <p>Function returns a DataFrame containing the signal gain for each gene and the cumulative signal gain for genes in the same cluster. The signal threshold hyperparameter was used for differentiating signal cells from background noise.</p>


### <ins>Data Analysis Demo</ins>
The tutorial for the analysis of FISHnCHIPs image data is provided as a jupyter notebook file **`FISHnCHIPs data analysis tutorial.ipynb`** in the **FISHnCHIPs_DataAnalysis_Tutorial** folder.
In the tutorial, various functions were called from the **`FISHnCHIPS_DataAnalysis.py`** package <mark>(To be packaged, containing FISHnCHIPsImages.py, registerFunction.py and segmentationFunction.py)</mark> in the **scripts folder**.

### Input files provided:
* DAPI Images (TIF format)
* Hyb Images (TIF format)
* Bleached Hyb Images (TIF format)
* Schema file (CSV file)

### Functions:
#### <ins>Generation of cell mask images **`FISHnCHIPsImages.py`**</ins> 
1. **`segment_dapi_one_fov`**   
    * Inputs: 
        * DAPI Image (TIF format)
    * Hyperparameters:
        * Parameters from yaml file
        * Output path 
        * List of DAPI FOVs to segment
    <p> Function saves a DAPI image with overlapping border removed (TIF format), dilated cell masks image (TIF format) and an overlay of cell masks onto the DAPI image (JPG format) in the output path folder.</p>

2. **`subtract_one_image`**   
    * Inputs: 
        * Hyb Images (TIF format)
        * Bleached Hyb Images (TIF format)

    <p> Function returns the bleach-subtracted Hyb images in TIF format.</p>

3. **`segment_one_fov`**
    * Inputs: 
        * DAPI Image (TIF format)
        * Hyb Images (TIF format)
        * Bleached Hyb Images (TIF format)

    An all-in-one function that performs the functions of **`segment_dapi_one_fov`** and **`subtract_one_image`** when all DAPI, Hyb and bleached Hyb images are ready. It also returns list of all cell mask for spatial and mask intensity analysis.

#### <ins>Calculate cell positions & brightness intensity **`segmentationFunctions.py`**</ins>

4. **`get_centroids`** 
    * Inputs: 
        * Dilated cell mask image (TIF format)
        * DAPI Image after border removal (TIF format)
    * Hyperparameters:
        * Output path 
        * Prefix of dilated cell mask image
        * Prefix of DAPI image after border removal
        * List of FOVs to segment
    <p> Function returns the coordinates of the cell masks of each FOV as a DataFrame.</p>

5. **`get_mask_positions_inFuse`**  
    * Inputs: 
        * List of FOVs that cell masks belongs to
        * List of x-coordinates of corresponding cell masks
        * List of y-coordinates of corresponding cell masks
        * Number of FOVs along x-axis
        * Number of FOVs along y-axis

    <p> Function returns the coordinates of the cell masks in the context where all FOVs are fused as a DataFrame.</p>

6. **`get_mask_intensity_matrix`**
    * Inputs: 
        * List of all cell masks

    Function takes in the list of all cell masks from the **`segment_one_fov`** function as input and returns 4 separate DataFrames with the mean, median, maximum and summation of mask intensity of each cell for each Hyb.

7. **`get_mask_info`** 
    * Inputs: 
        * List of all cell masks

    Function takes in the list of all cell masks from the **`segment_one_fov`** function as input and returns various information of the cell masks including list of FOVs analysed, cell types identified, mask intensity, area of masks and mask spatial positions as a DataFrame.

8. **`get_mask_positions_inFOV`**
    * Inputs: 
        * Cell mask information

    Function takes in the cell mask information from the **`get_mask_info`** function as input and returns a tuple of fov, x-coordinate and y-coordinate of each cell position in their respective FOV which can be fed into the **`get_mask_positions_inFuse`** function to obtain the spatial position of cell masks in the context where all FOVs are fused.



### <ins>FISHnCHIPs manuscript figure</ins>
Figures 3 and 4 are produced using R data visualizaion packages while Figure 5 and 6 are produced with Python visualization tools. For Figure 5 and 6, the packages used are the same as the tutorial, hence simply run the jupyter notebook to reproduce the figures. For Figure 3 and 4, functions from the standalone package **`capFISHImage`** were used. Please install the package provided in the **package** folder using the following command and run the script provided accordingly:

    install.packages('./package/capFISHImage_0.1.0.zip', repos=NULL, type='source')



## Parameter settings
An explanation of each of the parameters used in **`Gene panel design tutorial.ipynb`** and **`FISHnCHIPs data analysis tutorial.ipynb`**.

<details><summary><i>Gene panel design tutorial</i></summary>

* **min_corr**: Minimum correlation for genes to be considered as highly correlated
* **min_ngenes**: Minimum number of genes to include in the panel design from each reference marker celltype
* **min_expression**: Expression threshold of genes in cells to binarize gene expression level
* **min_cells**: Minimum number of cells that the genes are expressed in
* **max_cells**: Maximum number of cells that the genes are expressed in
* **filt_name_pattern**: Naming pattern of mitochondrial or pseudo genes to remove

</details>

<details><summary><i>FISHnCHIPs data analysis tutorial</i></summary>

<details><summary><i>Parameters in yaml file:</i></summary>

* **mainpath**: The file path where DAPI and hyb image TIF file are stored
* **schema**: The file path where the schema csv file is located. Schema file should contain the cell type to analyse and its corresponding dye
* **dapiname**: Prefix name of dapi file
* **fovs**: Number of field of views (FOVs) to process
* **fov_x**: Number of FOVs along x-axis
* **fov_y**: Number of FOVs along y-axis
* **overlap_in_px**: Number of pixel that overlaps between each FOV
* **background**: Default "bleach"; Type of background image used to offset background noise
* **remove_background**: Default 'subtract'; How background noise is removed
* **show_dapiseg**: Boolean; Whether to run segmentation function
* **segmentation**: Select segmentation method; 'cellpose' or 'watershed'
* **segment_mode**: 'Cytoplasm', 'Cytoplasm2', 'Cytoplasm2_Omnipose', 'Nuclei' for cellpose
* **cellsize**: Cell size in μm
* **flow_threshold**: Default 0.4; Increase threshold if cellpose is not returning as many ROIs as you’d expect.
* **cellprob_threshold**: Default 0.0; Decrease this threshold if cellpose is not returning as many ROIs as you’d expect.
* **mask_dilate_factor**: Amount of dilation applied to cell mask
* **filtermask**: Boolean; Whether to remove small masks
* **filtermask_size**: Minimum mask size
* **anchor_name**: Prefix of hyb image to segment e.g.'Cy7_3_'
* **anchor_celltype**: Celltype label of the hyb image e.g.'CAF_1'

*Configuration for watershed segmentation*

* **fusedpath**: File path of fused image to determine the minimum and maximum intensity at a certain percentile
* **cutoff_threshold**: Threshold value which is used to classify the pixel values into foreground and background classes, creating a binary image
* **opening_threshold**: Number of iterations of erosion and dilation that the image goes through
* **kernel_size**: Size of kernel that slides through the image which will determine how much the image is eroded or dilated


</details>

<details><summary><i>Parameters in tutorial:</i></summary>

* **segmentDAPI**: Boolean input. True if segmentation step is desired.
* **dapi_prefix**: DAPI filename prefix
* **prefix**: Post-segmentation image file prefix
* **yml_file**: yaml filename

</details>

</details>

<p></p>

## Citation and code dependencies
* The raw scRNA-seq count matrix data was preprocessed using the [Seurat](https://satijalab.org/seurat/) pipeline to produce the cell-scaled gene-cell matrix used in this tutorial. The clustermap shown in Figure 3c-j were clustered using functions from the Seurat package as well.

    *Hao and Hao et al. Integrated analysis of multimodal single-cell data. Cell (2021) [Seurat V4]*


**The function provided in the tutorials uses the following packages that are included in the standard Anaconda distribution:**

<details><summary><ins>Python</ins></summary>

* [Pandas](https://pandas.pydata.org/)
* [Numpy](https://numpy.org/)
* [Matplotlib](https://matplotlib.org/)
* [Numba](https://numba.pydata.org/)
* [Seaborn](https://seaborn.pydata.org/)
* [Scipy](https://scipy.org/)
* [Scanpy](https://scanpy.readthedocs.io/en/stable/)
* [igraph](https://igraph.org/) 
* [Leidenalg](https://github.com/vtraag/leidenalg)
* [cellpose](https://www.cellpose.org/)
* [scikit-image](https://scikit-image.org/)

</details>

<details><summary><ins>R</ins></summary>

* [Tidyverse](https://www.tidyverse.org/)
* [Dittoseq](https://github.com/dtm2451/dittoSeq) 
* [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) 

</details>


## Authors
* ZHOU Xinrui *-Main code/script author, pipeline design and code troubleshooting-* (https://github.com/Xinrui0523)
* CHENG Teh How *-GitHub setup, tutorial design and code troubleshooting-* (https://github.com/ChengDHow)

## License
See full license [here](https://github.com/ChengDHow/FISHnCHIPs/blob/main/LICENSE).

## Acknowledgements
* Shyam Prabhakar
* SEOW Wan Yi
* Norbert HA How Ong
* Jeeranan Boonruangkan

## Questions
* CHEN Kok Hao, chenkh@gis.a-star.edu.sg
* CHOU Shijie Nigel, Nigel_Chou@gis.a-star.edu.sg
* ZHOU Xinrui, Zhou_Xinrui@gis.a-star.edu.sg
