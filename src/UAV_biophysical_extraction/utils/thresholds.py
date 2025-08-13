
import numpy as np
from skimage import filters



def threshold_apply(data, upper_threshold, lower_threshold, data_type):

    def calculate_threshold(threshold, image_data, type):
        if threshold is True:
            if type == 'height':
                valid_values = (image_data - np.nanmin(image_data)) / (np.nanmax(image_data) - np.nanmin(image_data))
                return filters.threshold_otsu(valid_values)
            elif type == 'spectral':
                valid_values = image_data[~np.isnan(image_data)]
                return filters.threshold_otsu(valid_values)
        elif threshold is False:
            return None
        elif isinstance(threshold, (int, float)):
            return threshold
        else:
            raise ValueError("Threshold must be int, float, True, or False.")

    up_limit = calculate_threshold(upper_threshold, data, data_type)
    low_limit = calculate_threshold(lower_threshold, data, data_type)

    if up_limit is not None and low_limit is None:
        mask_threshold = np.where(data <= up_limit, 1, np.nan)
    elif low_limit is not None and up_limit is None:
        mask_threshold = np.where(data >= low_limit, 1, np.nan)
    elif up_limit is not None and low_limit is not None:
        mask_threshold = np.where((data >= low_limit) & (data <= up_limit), 1, np.nan)
    else:
        mask_threshold = np.full(data.shape, np.nan)

    return mask_threshold, up_limit, low_limit