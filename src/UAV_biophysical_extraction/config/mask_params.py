# mask_processing

mask_params =   {"grid":         "grid.shp",
                  "apply_mask":   True, # or False
                  "output":       "mask"}


height_params = {"height_image_file": "height.tif",
                 "lower_threshold":   True, # or False or threshold value
                 "upper_threshold":   False} # or True or threshold value


spectral_params ={"spectral_image_file":       "multispectral.tif",
                  "spectral_data":             "multispectral", # or RGB
                  "bands_ordered":             {"Blue": 1,
                                                "Green": 2,
                                                "Red": 3,
                                                "RedEdge": 4,
                                                "NIR": 5},
                  "vi":                         "SAVI", # or NGRDI among others
                  "lower_threshold":            True, # or False or threshold value
                  "upper_threshold":            False} # or True or threshold value