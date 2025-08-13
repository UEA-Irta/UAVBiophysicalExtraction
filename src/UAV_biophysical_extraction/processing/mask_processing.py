
from utils.vegetation_indices import calculate_VI
from utils.thresholds import threshold_apply
import numpy as np
from geopandas import GeoDataFrame
import rasterio
from rasterio.mask import mask
from rasterio import features
from rasterio.io import MemoryFile
from PIL import Image
import fiona
from shapely.geometry import shape


def mask_processing(general_params,  height_params=None, spectral_params=None):
    """
    Generates a mask to isolate relevant pixels based on specified thresholds. It can be used as a vegetation or ground mask.

    Parameters
    ----------
    general_params : dict
        Parameters applicable to all mask types:
            - grid : str, optional
                Path to a shapefile for masking.
            - apply_mask : bool
                Whether to save as a shapefile (True) or raster (False). Default is True.
            - output : str
                Output file name (without extension). Default is 'mask'.

    height_params : dict, optional
        Parameters specific to height data:
            - height_image_file : str or .tif
                Path to the spectral image file or image.tif.
            - lower_threshold : float or bool
                Lower threshold for the height mask. Default is True (apply threshold using otsu method).
            - upper_threshold : float or bool
                Upper threshold for the height mask. Default is False.

    spectral_params : dict, optional
        Parameters specific to spectral data:
            - spectral_image_file : str or .tif
                Path to the spectral image file or image.tif.
            - spectral_data : str
                Spectral type data: RGB or multispectral data.
            - spectral_order : dict
                Bands order in .tif or .jpg imagery.
            - lower_threshold : float or bool
                Lower threshold for the mask. Default is True (apply threshold using otsu method).
            - upper_threshold : float or bool
                Upper threshold for the mask. Default is False.
            - vi: str
                Vegetation index to calculate. Options include:
                - 'NGRDI' (Normalized Green-Red Difference Index) for **RGB** images.
                          Threshold: 0 (above 0 for vegetation, below 0 for ground).
                - 'SAVI' (Soil Adjusted Vegetation Index) for **multispectral** images with red and near-infrared bands.
                          Threshold: 0.2 (values above 0.2 indicate vegetation, values below 0.2 indicate ground). (Quintanilla-Albornoz, M., et al. 2023. https://doi.org/10.1007/s00271-023-00888-1)

    Returns
    -------
    None or GeoDataFrame
        Saves the output as a raster or shapefile based on `apply_mask`.
    """

    # general part

    grid = general_params.get('grid')
    apply_mask = general_params.get('apply_mask', True)
    output = general_params.get('output', 'mask')

    height_mask = None
    spectral_mask = None


    # height part

    if height_params:
        height_image = height_params.get('height_image_file')
        lower_threshold = height_params.get('lower_threshold', True)
        upper_threshold = height_params.get('upper_threshold', False)

        if height_image.endswith('.tif'):
            with rasterio.open(height_image, 'r') as height_raster:
                height_crs = height_raster.crs
                if grid:
                    with fiona.open(grid, 'r') as shapefile:
                        grid_data = [feature['geometry'] for feature in shapefile]
                    height_data, height_transform = mask(height_raster, grid_data, crop=True)
                else:
                    height_data = height_raster.read()
                    height_transform = height_raster.transform

        else:
            raise ValueError(f"Unsupported image format: {height_image}")

        mask_threshold, up_limit, low_limit = threshold_apply(height_data, upper_threshold, lower_threshold, 'height')

        if apply_mask:
            mask_image_height = mask_threshold
        else:
            mask_image_height = mask_threshold * height_data

        print(f'height lower threshold: {low_limit}, height upper threshold: {up_limit}')
        mask_image_height = np.squeeze(mask_image_height)


    # spectral part

    if spectral_params:
        spectral_image = spectral_params.get('spectral_image_file')
        spectral_data = spectral_params.get('spectral_data')
        bands_ordered = spectral_params.get('bands_ordered')
        lower_threshold = spectral_params.get('lower_threshold', True)
        upper_threshold = spectral_params.get('upper_threshold', False)
        blue_order = bands_ordered.get('Blue', 1)
        red_order = bands_ordered.get('Red', 2)
        green_order = bands_ordered.get('Green', 3)
        red_edge_order = bands_ordered.get('RedEdge', 4)
        nir_order = bands_ordered.get('NIR', 5)

        if spectral_data == 'RGB':
            vi = spectral_params.get('vi', 'NGRDI')

        elif spectral_data == 'multispectral':
            vi = spectral_params.get('vi', 'SAVI')

        else:
            f'error: {spectral_data} is not a valid spectral data type'


        if spectral_image.endswith('.tif'):

            with rasterio.open(spectral_image, 'r') as spectral_raster:
                spectral_crs = spectral_raster.crs
                if grid:
                    with fiona.open(grid, 'r') as shapefile:
                        grid_data = [feature['geometry'] for feature in shapefile]
                    spectral_data, spectral_transform = mask(spectral_raster, grid_data, crop=True)
                else:
                    spectral_data = spectral_raster.read()
                    spectral_transform = spectral_raster.transform

                if blue_order is not None:
                    B = spectral_data[blue_order - 1, :, :]
                if red_order is not None:
                    R = spectral_data[red_order - 1, :, :]
                if green_order is not None:
                    G = spectral_data[green_order - 1, :, :]
                if red_edge_order is not None:
                    RedEdge = spectral_data[red_edge_order - 1, :, :]
                if nir_order is not None:
                    NIR = spectral_data[nir_order - 1, :, :]

                bands = {
                    'Blue': B,
                    'Green': G,
                    'Red': R,
                    'RedEdge': RedEdge,
                    'NIR': NIR
                }

                VI_raster = calculate_VI(bands, vi, crs=spectral_crs, raster_transform=spectral_transform, save=True)
                mask_threshold, up_limit, low_limit = threshold_apply(VI_raster, upper_threshold, lower_threshold,'spectral')

                if apply_mask:
                    mask_image_spectral = mask_threshold
                else:
                    mask_image_spectral = mask_threshold * VI_raster

                print(f'{vi} lower threshold: {low_limit}, {vi} upper threshold: {up_limit}')
                mask_image_spectral = np.squeeze(mask_image_spectral)

        elif spectral_image.endswith('.jpg') or spectral_image.endswith('.jpeg'):
            img = Image.open(spectral_image)
            if img.mode == 'L':
                spectral_data = np.array(img, dtype=np.float32) / 255.0
                R = G = B = spectral_data

            else:
                img = img.convert('RGB')
                spectral_data = np.array(img, dtype=np.float32) / 255.0
                R = spectral_data[:, :, 0]
                G = spectral_data[:, :, 1]
                B = spectral_data[:, :, 2]

                bands = {
                    'Blue': B,
                    'Green': G,
                    'Red': R
                }
                VI_raster = calculate_VI(bands, vi, crs=None, raster_transform=None, save=True, save_format='jpg')
                VI_raster = (255 * (VI_raster - np.nanmin(VI_raster)) / (np.nanmax(VI_raster) - np.nanmin(VI_raster))).astype(np.uint8)
                mask_threshold, up_limit, low_limit = threshold_apply(VI_raster, upper_threshold, lower_threshold,'spectral')

                if apply_mask:
                    mask_image_spectral = mask_threshold
                    mask_image_spectral = np.nan_to_num(mask_image_spectral, nan=0)
                    mask_image_spectral[mask_image_spectral == 1] = 255
                    if np.unique(mask_image_spectral).size == 2:
                        mask_image_spectral = (255 * (mask_image_spectral - np.nanmin(mask_image_spectral)) / (np.nanmax(mask_image_spectral) - np.nanmin(mask_image_spectral))).astype(np.uint8)

                else:
                    mask_image_spectral = mask_threshold * VI_raster
                    if mask_image_spectral.dtype != np.uint8:
                        mask_image_spectral = (255 * (mask_image_spectral - np.nanmin(mask_image_spectral)) / (np.nanmax(mask_image_spectral) - np.nanmin(mask_image_spectral))).astype(np.uint8)

                img = Image.fromarray(mask_image_spectral, mode='L')
                img.save(f"{output}.jpg", format="JPEG")

        else:
            raise ValueError(f"Unsupported image format: {spectral_image}")


    # save

    if apply_mask:
        if height_params is not None and spectral_params is None:
            final_mask = GeoDataFrame({'geometry': [shape(geom) for geom, val in features.shapes(mask_image_height.astype('float32'),transform=height_transform) if val == 1]}, crs=height_crs)

        elif spectral_params is not None and height_params is None:
            final_mask = GeoDataFrame({'geometry': [shape(geom) for geom, val in features.shapes(mask_image_spectral.astype('float32'),transform=spectral_transform) if val == 1]}, crs=spectral_crs)

        else:
            gdf_height = GeoDataFrame({'geometry': [shape(geom) for geom, val in features.shapes(mask_image_height.astype('float32'),transform=height_transform) if val == 1]}, crs=height_crs)

            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=mask_image_spectral.shape[0], width=mask_image_spectral.shape[1], count=1, dtype=mask_image_spectral.dtype, crs=spectral_crs, transform=spectral_transform) as final:
                    final.write(mask_image_spectral, 1)

                    cut, cut_transform = mask(final, gdf_height.geometry, crop=True)

            final_mask = GeoDataFrame({'geometry': [shape(geom) for geom, val in features.shapes(cut.astype('float32'),transform=cut_transform) if val == 1]}, crs=spectral_crs)

        final_mask.to_file(f'{output}.shp')

    else:
        if height_params is not None and spectral_params is None:
            with rasterio.open(f'{output}.tif', 'w', driver='GTiff', count=1, dtype='float32', crs=height_crs, transform=height_transform, width=mask_image_height.shape[1], height=mask_image_height.shape[0]) as dst:
                dst.write(mask_image_height, 1)

        elif spectral_params is not None and height_params is None:
            with rasterio.open(f'{output}.tif', 'w', driver='GTiff', count=1, dtype='float32', crs=spectral_crs, transform=spectral_transform, width=mask_image_spectral.shape[1], height=mask_image_spectral.shape[0]) as dst:
                dst.write(mask_image_spectral, 1)
        else:
            final_mask = np.stack([mask_image_height, mask_image_spectral], axis=-1)
            with rasterio.open(f'{output}.tif', 'w', driver='GTiff', count=2, dtype='float32', crs=height_crs, transform=height_transform , width=final_mask.shape[1], height=final_mask.shape[0]) as dst:
                dst.write(final_mask[:, :, 0], 1)
                dst.write(final_mask[:, :, 1], 2)

                dst.update_tags(1, description="height")
                dst.update_tags(2, description=f"{vi}")
