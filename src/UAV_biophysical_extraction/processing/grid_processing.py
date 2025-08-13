
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
from shapely.affinity import translate, rotate


def grid_processing(grid_params):
    """
    Generates a grid based on the parameters specified in the grid_params dictionary.

    Parameters
    ----------
    grid_params : dict
        A dictionary containing the following keys:
        - x_origin, y_origin : float
            Coordinates of the center of the polygon at the lower-left origin/upper-right origin.
        - width, height : float
            Width and height of the polygon. Defaults to 1.
        - epsg : int
            EPSG code representing the reference system. Defaults to 32631 (UTM Zone 31N, WGS84).
        - horizontal_distance, vertical_distance : float
            Horizontal and vertical distance between polygons. Defaults to 1.
        - num_columns, num_rows : int
            Number of columns and rows. Defaults to 1.
        - azimuth : float
            Azimuth angle in degrees (measured clockwise from North). Defaults to 0.
        - buffer : float
            Buffer distance for polygons. Defaults to 0.
        - grid_file : str
            File name without extension for the output file. Default is None (no file saved).

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame containing the grid geometries.
    """

    x_origin = grid_params.get('x_origin', 0)
    y_origin = grid_params.get('y_origin', 0)
    width = grid_params.get('width', 1)
    height = grid_params.get('height', 1)
    epsg = grid_params.get('epsg', 32631)
    horizontal_distance = grid_params.get('horizontal_distance', 1)
    vertical_distance = grid_params.get('vertical_distance', 1)
    num_columns = grid_params.get('num_columns', 1)
    num_rows = grid_params.get('num_rows', 1)
    azimuth = grid_params.get('azimuth', 0)
    buffer = grid_params.get('buffer', 0)
    grid_file = grid_params.get('grid_file', None)

    base_polygon = Polygon([(0, 0), (width, 0), (width, height), (0, height), (0, 0)])

    wkt_list = []
    for i in range(num_rows):
        for j in range(num_columns):
            x_offset = horizontal_distance * j - width / 2
            y_offset = vertical_distance * i - height / 2

            current_polygon = translate(base_polygon, xoff=x_offset, yoff=y_offset)
            current_polygon = current_polygon.buffer(buffer)
            current_polygon = rotate(current_polygon, azimuth, origin=(0, 0))
            current_polygon = translate(current_polygon, xoff=x_origin, yoff=y_origin)

            wkt = current_polygon.wkt
            wkt_list.append(wkt)

    df = pd.DataFrame({'geometry': wkt_list})
    gdf = gpd.GeoDataFrame(df, geometry=gpd.GeoSeries.from_wkt(wkt_list), crs=f"EPSG:{epsg}")

    if grid_file:
        gdf.to_file(grid_file + '.shp')

    return gdf