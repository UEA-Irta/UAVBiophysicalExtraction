from processing.grid_processing import grid_processing
from processing.mask_processing import mask_processing
from processing.zonalstats_processing import zonalstats_processing
from config.grid_params import grid_params
from config.mask_params import mask_params, height_params, spectral_params
from config.zonalstats_params import stats_params, structure_params, VI_params, LST_params

def main(grid=True, mask=True, zonalstats=True):

    # GRID
    if grid:
        grid_processing(grid_params)

    # MASK
    if mask:
        mask_processing(mask_params, height_params, spectral_params)

    # ZONAL STATS
    if zonalstats:
        zonalstats_processing(stats_params, structure_params, VI_params, LST_params)



if __name__ == "__main__":
    main(grid=True, mask=True, zonalstats=True)
