
from utils.vegetation_indices import calculate_VI
import os
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.mask import mask
import rasterio as rio
from rasterstats import zonal_stats
import fiona
import warnings



def zonalstats_processing(general_params, structure_params, VI_params, LST_params):
    """
    Generate a dataset with biophysical traits based on the following parameters.

    Parameters
    ----------
    general_params : dict
        Parameters applicable to all mask types:
            - grid : str
                Path to a shapefile for masking.
            - vegetation_mask : str
                Path to the vegetation mask shapefile.
            - statistics : list of str
                List of statistical operations to apply on the data. Options include:
                - "mean", "max", "min", "std", "percentile_25", "percentile_50", "percentile_75".
            - output : str
                Output file name (without extension). Default is "UAV_biophysical_stats".

    structure_params : dict, optional
        Parameters specific to vegetation type and height data:
            - height_image_file : str or .tif
                Path to the height data image file.
            - vegetation_type : dict
                A dictionary specifying the vegetation types and their respective parameters:
                    - herbaceous : dict
                        - enabled : bool
                            Whether herbaceous vegetation is included. Default is True.
                    - woody : dict
                        - enabled : bool
                            Whether woody vegetation is included. Default is False.
                        - distance_trees : float
                            Distance between trees for woody vegetation (in meters). Default is 1.5.
                        - distance_rows : float
                            Distance between rows for woody vegetation (in meters). Default is 3.

    VI_params : dict, optional
        Parameters specific to spectral vegetation index (VI) data:
            - spectral_image_file : str or .tif
                Path to the spectral image file.
            - spectral_data : str
                Type of spectral data: "multispectral" or "RGB". Default is "multispectral".
            - bands_ordered : dict
                A dictionary specifying the band order for multispectral images:
                    - Example: {"Blue": 1, "Green": 2, "Red": 3, "RedEdge": 4, "NIR": 5}.
            - vi : list of str
                List of vegetation indices to calculate from the bands. Options include:
                    - "NDVI" (Normalized Difference Vegetation Index)
                    - "NDWI" (Normalized Difference Water Index)
                    - "CIg" (Green Chlorophyll Index)
                    - "SAVI" (Soil Adjusted Vegetation Index)
                    - "MTVI2" (Modified Triangular Vegetation Index)
                    - "NDRE" (Normalized Difference Red Edge)
                    - "VARI" (Visible Atmospherically Resistant Index)
                    - "IronOxide", and others.

    LST_params : dict, optional
        Parameters specific to Land Surface Temperature (LST) data:
            - temperature_image_file : str or .tif
                Path to the temperature image file.
            - soil_mask : str or None
                Path to the soil mask shapefile (optional). Default is None.
            - vegetation_cover_mask : str or None
                Path to the vegetation cover mask shapefile (optional). Default is None.

    Returns
    -------
    None
        This function processes the input parameters and generates a dataset of biophysical traits.

    structural_traits :
            - Fvc : float
                Fractional green cover, expressed as a percentage (%/100).
            - Hc : float
                Vegetation height (m). Includes all calculated statistics (e.g., mean, max, min, std).
            - Ac : float
                Vegetation zenithal area (m²).
            - Vc : float
                Vegetation volume (m³).
            - Wc : float, optional
                Vegetation width (m). Includes all calculated statistics. Extracted only for woody crops.
            - Hc/Wc : float, optional
                The ratio between maximum height (hc) and width (wc) (m/m). Extracted only for woody crops.

    VI_traits :
            - vi : list of floats
                Calculated vegetation indices based on the specified options.

    LST_traits :
            - Tc : float
                Temperature of vegetation (K).
            - Trad : float
                Radiometric temperature (K).
            - Ts : float, optional
                Temperature of soil (K). Included only if a soil mask is provided.
            - Tvc : float, optional
                Temperature of the vegetation cover (K). Included only if a vegetation cover mask is provided.


    """

    # general part

    grid = general_params.get('grid')
    vegetation_mask = general_params.get('vegetation_mask')
    statistics = general_params.get('statistics', ['mean', 'max'])
    output = general_params.get('output', 'UAV_biophysical_stats')

    with fiona.open(grid) as grid_shp:
        grid_crs = grid_shp.crs
        grid_data = gpd.read_file(grid, crs=grid_crs)
        grid_dataset = None

    with fiona.open(vegetation_mask) as veg_mask_shp:
        veg_mask_crs = veg_mask_shp.crs
        veg_mask_data = gpd.read_file(vegetation_mask, crs=veg_mask_crs)
        veg_mask_shapes = [feature['geometry'] for feature in veg_mask_shp]


    # structure part

    if structure_params:

        height_file = structure_params.get("height_image_file")
        veg_params = structure_params["vegetation_type"]

        if veg_params["woody"]["enabled"] and veg_params["herbaceous"]["enabled"]:
            raise ValueError("Vegetation type cannot be both woody and herbaceous.")

        with rio.open(height_file, 'r') as raster_data:
            height_data, transform_height = mask(raster_data, veg_mask_shapes, crop=True, nodata=np.nan)

        if height_data.ndim == 3:
            height_data = height_data[0]

        with rasterio.open('height.tif', 'w', driver='GTiff', height=height_data.shape[0], width=height_data.shape[1],
                           count=1, dtype=height_data.dtype, crs=raster_data.crs, transform=transform_height) as dst:
            dst.write(height_data, 1)

        grid_dataset = grid_data.copy()

        for statistic in statistics:
            values_index = zonal_stats(grid, 'height.tif', nodata=np.nan, stats=statistic)
            stats_index = [zone[statistic] for zone in values_index]
            grid_dataset[f'Hc_{statistic}'] = stats_index

        os.remove('height.tif')

        fc_list = np.zeros(len(grid_data))
        ac_list = np.zeros(len(grid_data))
        is_herbaceous = veg_params["herbaceous"]["enabled"]
        is_woody = veg_params["woody"]["enabled"]

        if is_woody:
            distance_trees = veg_params["woody"]["distance_trees"]
            distance_rows = veg_params["woody"]["distance_rows"]

        for i in range(len(grid_data)):
            single = grid_data.iloc[i]
            area_poly = single["geometry"].area
            if veg_mask_data.crs != grid_crs:  # Ensure they match
                veg_mask_data = veg_mask_data.to_crs(grid_crs)
            poly = gpd.overlay(veg_mask_data, gpd.GeoDataFrame(geometry=[single["geometry"]], crs=grid_crs))

            ac = poly.area.sum()
            ac_list[i] = ac

            if is_herbaceous:
                fc = ac / area_poly
            elif is_woody:
                fc = ac / (distance_trees * distance_rows)
            else:
                fc = 1

            fc_list[i] = fc
            ac_list[i] = ac

        grid_dataset['Fvc'] = fc_list
        grid_dataset['Ac'] = ac_list



        if veg_params["woody"]["enabled"]:
            if veg_params["woody"]["row_crop"]:
                grid_dataset['Wc'] = distance_rows * grid_dataset['Fvc']
            else:
                grid_dataset['Wc'] = 2 * np.sqrt((grid_dataset['Fvc'] * distance_trees * distance_rows) / 3.1416)

        elif veg_params["herbaceous"]["enabled"]:
            grid_dataset['Wc'] = grid_dataset['Ac'] / grid_dataset['Fvc']

        else:
            print("The training system is not valid.")

        if 'Hc_max' in grid_dataset.columns or 'Hc_mean' in grid_dataset.columns:
            hc_column = 'Hc_max' if 'Hc_max' in grid_dataset.columns else 'Hc_mean'
            grid_dataset['Vc'] = grid_dataset[hc_column] * grid_dataset['Ac']
            hc_column = 'Hc_mean' if 'Hc_mean' in grid_dataset.columns else 'Hc_max'
            grid_dataset['Hc_Wc'] = grid_dataset[hc_column] / grid_dataset['Wc']
        else:
            raise ValueError("statistics must contain either 'Hc_max' or 'Hc_mean'.")

        cols_to_move = ['Fvc', 'Ac', 'Vc', 'Wc', 'Hc_Wc']
        other_cols = [col for col in grid_dataset.columns if col not in cols_to_move]
        origin_cols = [col for col in grid_data.columns if col not in cols_to_move]
        cols_dif = list(set(other_cols) - set(origin_cols))

        grid_dataset = grid_dataset[origin_cols + cols_to_move + cols_dif]
        print("The structural features have now been successfully calculated!")



    if VI_params:

        image_file = VI_params.get("spectral_image_file")
        bands_ordered = VI_params.get("bands_ordered")
        vi = VI_params.get("vi")
        blue_order = bands_ordered.get('Blue', 1)
        red_order = bands_ordered.get('Red', 2)
        green_order = bands_ordered.get('Green', 3)
        red_edge_order = bands_ordered.get('RedEdge', 4)
        nir_order = bands_ordered.get('NIR', 5)

        if grid_dataset is not None:
            warnings.filterwarnings("ignore", message="Column names longer than 10 characters will be truncated")
            grid_dataset.to_file('dataset.shp')
            grid_extract = 'dataset.shp'

        else:
            grid_dataset = grid_data.copy()
            grid_extract = grid

        with rasterio.open(image_file, 'r') as raster_data:
            crs = raster_data.crs
            spectral_data, spectral_transform = mask(raster_data, veg_mask_shapes, crop=True)
            spectral_data[spectral_data == 0] = np.nan

            B = None
            R = None
            G = None
            RedEdge = None
            NIR = None

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

            for index in vi:
                VI_raster = calculate_VI(bands, index, crs=crs, raster_transform=spectral_transform, save=True)

                for statistic in statistics:
                    values_index = zonal_stats(grid_extract, f'{index}.tif', nodata=np.nan, stats=statistic)
                    stats_index = [zone[statistic] for zone in values_index]
                    grid_dataset[f'{index}_{statistic}'] = stats_index

                os.remove(f'{index}.tif')

        if grid_extract == 'dataset.shp':
            os.remove('dataset.shp')

        print("The VI features have now been successfully calculated!")



    if LST_params:

        temp_file = LST_params.get("temperature_image_file")
        soil_mask_file = LST_params.get("soil_mask_file", None)
        cc_mask_file = LST_params.get("vegetation_cover_mask_file", None)

        if grid_dataset is not None:
            warnings.filterwarnings("ignore", message="Column names longer than 10 characters will be truncated")
            grid_dataset.to_file('dataset.shp')
            grid_extract = 'dataset.shp'
        else:
            grid_dataset = grid_data.copy()
            grid_extract = grid


        if soil_mask_file is not None:
            with fiona.open(soil_mask_file) as soil_mask_shp:
                soil_mask_shapes = [feature['geometry'] for feature in soil_mask_shp]


        if cc_mask_file is not None:
            with fiona.open(cc_mask_file) as cc_mask_shp:
                cc_mask_shapes = [feature['geometry'] for feature in cc_mask_shp]

        with rasterio.open(temp_file, 'r') as raster_data:
            temp_crs = raster_data.crs
            veg_mask, veg_mask_transform = mask(raster_data, veg_mask_shapes, crop=True, nodata=np.nan)
            veg_mask = np.squeeze(veg_mask)

            with rasterio.open('veg_mask.tif', 'w', driver='GTiff', height=veg_mask.shape[0], width=veg_mask.shape[1], count=1, dtype=veg_mask.dtype, crs=temp_crs, transform=veg_mask_transform) as veg_raster:
                veg_raster.write(veg_mask, 1)

            for statistic in statistics:
                values_index = zonal_stats(grid_extract, 'veg_mask.tif', nodata=np.nan, stats=statistic)
                stats_index = [zone[statistic] for zone in values_index]
                stats_index_with_offset = [(value + 273.15) if value is not None else None for value in stats_index]
                grid_dataset[f'Tc_{statistic}'] = stats_index_with_offset

            if soil_mask_file is not None and cc_mask_file is not None:
                warnings.filterwarnings("ignore", message="Column names longer than 10 characters will be truncated")
                grid_dataset.to_file('Tc.shp')

            if soil_mask_file is not None:
                soil_mask, soil_mask_transform = mask(raster_data, soil_mask_shapes, crop=True, nodata=np.nan)
                soil_mask = np.squeeze(soil_mask)
                with rasterio.open('soil_mask.tif', 'w', driver='GTiff', height=soil_mask.shape[0], width=soil_mask.shape[1], count=1, dtype=soil_mask.dtype, crs=temp_crs, transform=soil_mask_transform) as soil_raster:
                    soil_raster.write(soil_mask, 1)

                for statistic in statistics:
                    values_index = zonal_stats('Tc.shp', 'soil_mask.tif', nodata=np.nan, stats=statistic)
                    stats_index = [zone[statistic] for zone in values_index]
                    stats_index_with_offset = [(value + 273.15) if value is not None else None for value in stats_index]
                    grid_dataset[f'Ts_{statistic}'] = stats_index_with_offset

            if cc_mask_file is not None:
                warnings.filterwarnings("ignore", message="Column names longer than 10 characters will be truncated")
                grid_dataset.to_file('Tc_Ts.shp')
                if soil_mask_file is not None:
                    grid_extract = 'Tc_Ts.shp'
                else:
                    grid_extract = 'Tc.shp'

            if cc_mask_file is not None:
                cc_mask, cc_mask_transform = mask(raster_data, cc_mask_shapes, crop=True, nodata=np.nan)
                cc_mask = np.squeeze(cc_mask)
                with rasterio.open('cc_mask.tif', 'w', driver='GTiff', height=cc_mask.shape[0], width=cc_mask.shape[1], count=1, dtype=cc_mask.dtype, crs=temp_crs, transform=cc_mask_transform) as soil_raster:
                    soil_raster.write(cc_mask, 1)

                for statistic in statistics:
                    values_index = zonal_stats(grid_extract, 'cc_mask.tif', nodata=np.nan, stats=statistic)
                    stats_index = [zone[statistic] for zone in values_index]
                    stats_index_with_offset = [(value + 273.15) if value is not None else None for value in stats_index]
                    grid_dataset[f'Tcc_{statistic}'] = stats_index_with_offset

        if grid_extract == 'dataset.shp':
            os.remove('dataset.shp')

        os.remove('veg_mask.tif')
        if soil_mask_file is not None and cc_mask_file is not None:
            os.remove('Tc.shp')
            os.remove('soil_mask.tif')
        if cc_mask_file is not None:
            os.remove('Tc_Ts.shp')
            os.remove('cc_mask.tif')

        Trad = zonal_stats(grid, temp_file, stats='mean')
        Trad_list = [zone['mean'] for zone in Trad]
        stats_index_with_offset = [(value + 273.15) if value is not None else None for value in Trad_list]
        grid_dataset['T_R'] = stats_index_with_offset
        print("The LST features have now been successfully calculated!")

    grid_dataset.to_file(output + '.shp')
    grid_dataset.to_excel(output + '.xlsx')


