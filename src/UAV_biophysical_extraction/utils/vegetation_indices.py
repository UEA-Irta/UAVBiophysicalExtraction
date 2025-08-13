
import numpy as np
import rasterio
from PIL import Image



def calculate_VI(bands, variable, crs=None, raster_transform=None, save=False, save_format='tif'):
    # Ensure required bands are present
    required_bands = {
        'NDVI': ['NIR', 'Red'],
        'NDWI': ['NIR', 'Green'],
        'CIg': ['NIR', 'Green'],
        'SAVI': ['NIR', 'Red'],
        'MTVI2': ['NIR', 'Red', 'Green'],
        'NDRE': ['NIR', 'RedEdge'],
        'VARI': ['Green', 'Red', 'Blue'],
        'IronOxide': ['Red', 'Blue'],
        'NGRDI': ['Green', 'Red'],
    }

    if variable in required_bands:
        for band in required_bands[variable]:
            if band not in bands or bands[band] is None:
                raise ValueError(f"Band {band} is required for {variable} but is missing.")

    # Compute indices
    if variable == 'B':
        index = bands['Blue']
    elif variable == 'G':
        index = bands['Green']
    elif variable == 'R':
        index = bands['Red']
    elif variable == 'RE':
        index = bands['RedEdge']
    elif variable == 'NIR':
        index = bands['NIR']
    elif variable == 'VARI':
        index = (bands['Green'] - bands['Red']) / (bands['Green'] + bands['Red'] - bands['Blue'])
    elif variable == 'IronOxide':
        index = bands['Red'] / bands['Blue']
    elif variable == 'NGRDI':
        index = (bands['Green'] - bands['Red']) / (bands['Green'] + bands['Red'])
    elif variable == 'NDVI':
        index = (bands['NIR'] - bands['Red']) / (bands['NIR'] + bands['Red'])
    elif variable == 'NDWI':
        index = (bands['NIR'] - bands['Green']) / (bands['NIR'] + bands['Green'])
    elif variable == 'CIg':
        index = (bands['NIR'] / bands['Green']) - 1
    elif variable == 'SAVI':
        L = 0.5
        index = ((bands['NIR'] - bands['Red']) / (bands['NIR'] + bands['Red'] + L)) * (1 + L)
    elif variable == 'MTVI2':
        index = (1.5 * (1.2 * (bands['NIR'] - bands['Green']) - 2.5 * (bands['Red'] - bands['Green']))) / \
                ((2 * bands['NIR'] + 1) ** 2 - (6 * bands['NIR'] - 5 * bands['Red'] ** 0.5) - 0.5) ** 0.5
    elif variable == 'NDRE':
        index = (bands['NIR'] - bands['RedEdge']) / (bands['NIR'] + bands['RedEdge'])
    else:
        raise ValueError(f"Variable {variable} not recognized")

    # Save the result
    if save:
        if save and save_format == 'jpg':
            normalized_index = (255*(index - np.nanmin(index)) / (np.nanmax(index) - np.nanmin(index))).astype(
                np.uint8)
            img = Image.fromarray(normalized_index, mode='L')
            img.save(f'{variable}', format="JPEG")

        if save and save_format == 'tif':
            if crs is None or raster_transform is None:
                raise ValueError("CRS and raster transform are required to save the file")

            with rasterio.open(
                    f'{variable}.tif',
                    'w',
                    driver='GTiff',
                    height=index.shape[0],
                    width=index.shape[1],
                    count=1,
                    dtype=index.dtype,
                    crs=crs,
                    transform=raster_transform
            ) as dst:
                dst.write(index, 1)

    return index