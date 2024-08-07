"""
Expected folder structure is:

path/
├── sample1/
|   sample1.h5ad
│   ├── raw_feature_bc_matrix/
│   ├── spatial/
|----metadata.csv
....
├── samplek/
|   samplek.h5ad
│   ├── raw_feature_bc_matrix/
│   ├── spatial/     
"""
import scanpy as sc
import anndata
import pandas as pd
import imageio as im
import gzip
import json
from io import BytesIO

def dirpaste(*args):
    """
    A function to concatenate strings to form a directory path.
    """
    return '/'.join(args)

def gex_anndata(path, sample):
    """
    A function to load 10x Genomics single-cell gene expression data into anndata object.
    
    Expects a folder of the form:
    
    path/
    ├── sample1/
    |   sample1.h5ad
    │   ├── raw_feature_bc_matrix/
    ....
    ├── samplek/
    |   samplek.h5ad
    │   ├── raw_feature_bc_matrix/
    
    Parameters
    ----------
    path: str
        Path to the directory containing the sample data.
    sample: str
        Name of the sample.
        
    Returns
    -------
    adata: anndata.AnnData
        An anndata object containing the loaded data.
    """
    f = sample + '.h5ad'
    adata = sc.read_h5ad(dirpaste(path, sample, f))

def gex_spatial(path, sample,
                spatial_key='spatial',
                feature_folder='raw_feature_bc_matrix',
                spatial_folder='spatial',
                scalefactors='scalefactors_json.json.gz',
                tissue_image='tissue_hires_image.png.gz',
                tissue_positions='tissue_positions_list.csv.gz',
                gzip_compression=True):
    """
    A function to load 10x Visium spatial data into anndata object.
    
    Expects a folder of the form:
    
    path/
    ├── sample1/
    |   sample1.h5ad
    │   ├── raw_feature_bc_matrix/
    │   ├── spatial/
    |----metadata.csv
    ....
    ├── samplek/
    |   samplek.h5ad
    │   ├── raw_feature_bc_matrix/
    │   ├── spatial/ 
    
    Parameters
    ----------
    path: str
        Path to the directory containing the sample data.
    sample: str
        Name of the sample.
    spatial_key: str
        Key to store spatial data in anndata object.
    feature_folder: str
        Name of the folder containing the feature barcode matrix.
    spatial_folder: str
        Name of the folder containing the spatial data.
    scalefactors: str
        Name of the file containing scalefactors.
    tissue_image: str
        Name of the file containing the tissue image.
    tissue_positions: str
        Name of the file containing the tissue positions.
    gzip_compression: bool
        Whether the files are compressed
        
    Returns
    -------
    adata: anndata.AnnData
        An anndata object containing the loaded data.
    """
    if gzip_compression == False:
        scalefactors.replace('.gz', '')
        tissue_image.replace('.gz', '')
        tissue_positions.replace('.gz', '')
        
    adata = sc.read_10x_mtx(dirpaste(path, sample, feature_folder))
    
    spot_factors = dirpaste(path, sample, spatial_folder, scalefactors)
    imf = dirpaste(path, sample, spatial_folder, tissue_image)
    tissue_positions = dirpaste(path, sample, spatial_folder, tissue_positions)
     
    # handle compression specific commands:
    if gzip_compression:
        with gzip.open(spot_factors, 'r') as file:
            spots = json.load(file)
            
        with gzip.open(imf, 'rb') as f:
            image_file = BytesIO(f.read())
    else:
        with open(spot_factors, 'r') as file:
            spots = json.load(file)
        image_file = imf
            
    # load spatial data into adata object:        
    library_id = sample
    adata.uns[spatial_key] = {library_id: {}}
    adata.uns[spatial_key][library_id]["images"] = {}
    adata.uns[spatial_key][library_id]["metadata"] = {}
    adata.uns[spatial_key][library_id]["metadata"]['source_image_path'] = imf
    adata.uns[spatial_key][library_id]["images"] = {"hires": im.v2.imread(image_file)}
    adata.uns[spatial_key][library_id]["scalefactors"] = {
        "tissue_hires_scalef": spots["tissue_hires_scalef"],
        "spot_diameter_fullres": spots['spot_diameter_fullres'],
    }

    # load tissue positions metadata:
    coords = pd.read_csv(tissue_positions, header=None, index_col=0, compression='gzip')
    coords.columns = ["in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]
    coords.set_index(coords.index.astype(adata.obs.index.dtype), inplace=True)
    
    # merge tissue positions with adata:
    adata.obs = pd.merge(adata.obs, coords, how="left", left_index=True, right_index=True)
    adata.obs['in_tissue'] = adata.obs['in_tissue'].astype('category')
    adata.obsm[spatial_key] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values
    adata.obs.drop(columns=["pxl_row_in_fullres", "pxl_col_in_fullres"], inplace=True)
    return adata
