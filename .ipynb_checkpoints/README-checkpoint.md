# FISHnCHIPs Gene Panel Design & Image Processing
* insert description here

## Contents
* Prerequisites
    * system requirements
    * software requirements
    * tested dependencies
* Getting started
    * Installation
    * Gene Panel Design Demo
        * Time taken to run
        * Generation of panel design
        * Clustering of candidate genes
    * Image Processing Demo
        * Time taken to run
        * Generation of cell mask images
        * Calculate cell positions & brightness intensity
        * Stitching multiple FOVs of cell mask images
* Parameter settings
    * Gene Panel Design
    * Image Processing
* Citations and code dependencies
* Authors
* License
* Acknowledgements 

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
***Common packages used in both gene panel design and image processing:***
* Python 3.9.16
* pandas 1.5.2
* numpy 1.23.5
* matplotlib 3.6.2
* numba 0.56.4
* seaborn 0.12.2
* scipy 1.10.0

***For Gene Panel Design:***
* igraph 0.10.4
* leidenalg 0.9.1

***For Image Processing:***
* cellpose 2.2
* scikit-image 0.19.3

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
5. Code has been tested to run on Spyder but also can also be run from command line (run MainScript.py) or any other Python 3 compatible IDE.

### Gene Panel Design Demo
Insert example codes and description on their funtionality

### Image Processing Demo
Insert example codes and description on their funtionality

## Parameter settings
Explain parameters and file naming conventions used in demo codes

## Citation and code dependencies

## Authors

## License
See full license [here]

## Acknowledgements

## Questions
