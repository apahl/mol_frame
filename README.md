# MolFrames

Working with Dask and Pandas DataFrames that can handle chemical structures.  

**[Update 29-Jan-2018]**: The module has been rewritten to now handle **Dask** and Pandas DataFrames.
Most functions and methods (except for the SDF reader) now return a MolFrame where the `data` object is a **Dask** DataFrame to allow handling of arbitrarily large files and datasets. To explicitly get a Pandas DataFrame, call `.compute()` on a MolFrame instance.  
A `.get_data_type()`  method was added to show the data type of the underlying `data` obejct.

## Features:
Output as sortable and optionally selectable HTML Tables:
![HTML Tables](res/molframe.png)  
(this is an extension of bluenote10's NimData [html browser](https://github.com/bluenote10/NimData/blob/master/src/nimdata/html.nim) 
which uses JQuery Datatables.)

Interactive HoloViews plots with structure tooltips:
![Scatter Plot](res/scatter_test.png)

Operations on molecules. For long operations, this displays a progress bar in the Notebook.

## Requirements
The recommended way to install the dependencies is via [conda](https://www.anaconda.com/download/).
* Python 3 (come on, you can do the switch to 3!)
* Jupyter Notebook
* [RDKit](http://rdkit.org/)
* [Dask](https://dask.pydata.org/)
* [HoloViews](http://holoviews.org/) (optionally, used for plotting)

The code is written and tested on Ubuntu 17.04 / 17.10 and Python 3.6 and is intended for use in the Jupyter Notebook.

The main class is the MolFrame, which is a wrapper around a Pandas dataframe, exposing all DataFrame methods and extending it with some chemical functionality from the RDKit.  
The underlying DataFrame is contained in the MolFrame.data object and can always be accessed directly, if necessary.

See the accompanying [Tutorial](tutorials/tutorial1.ipynb) notebook for further examples.

## Installation
After installing the requirements, clone this repo, then the module can be used by including the project's base directory (mol_frame) in Python's import path (I actually prefer this to using setuptools, because a simple git pull will get you the newest version).
This can be achieved as follows:

Put a file with the extension `.pth`, e.g. `my_packages.pth`, into one of the `site-packages` directories of your Python installation (e.g. `~/anaconda3/envs/chem/lib/python3.6/site-packages/`) and put the path to the base directory of this project (mol_frame) into it.
(I have the path to a dedicated folder on my machine included in such a .pth file and link all my development projects to that folder. This way, I need to create / modify the .pth file only once.)

## Disclaimer
This is work in progress, compatibility-breaking changes will happen.
