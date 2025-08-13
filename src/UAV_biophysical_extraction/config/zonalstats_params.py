

# zonalstats_processing


stats_params = {"grid":               "grid.shp",
                  "vegetation_mask":    "mask_crop.shp",
                  "statistics":         ["mean", "max", "min", "std", "percentile_25", "percentile_50", "percentile_75"],
                  "output":             "UAV_biophysical_traits"}


structure_params = {"height_image_file": "height.tif",
                    "vegetation_type": {
                        "herbaceous": {
                            "enabled": False,
                        },
                        "woody": {
                            "enabled": True,
                            "distance_trees": 5,
                            "distance_rows": 6,
                            "row_crop": False
                        }}}


VI_params = {"spectral_image_file":       "multispectral.tif",
             "spectral_data":             "multispectral", # or RGB
             "bands_ordered":             {"Blue": 1,
                                           "Green": 2,
                                           "Red": 3,
                                           "RedEdge": 4,
                                           "NIR": 5},
             "vi":                         ["NDVI", "SAVI"]  #["B", "G", "R", "RE", "NIR", "NDVI", "NDWI", "CIg", "SAVI", "MTVI2", "NDRE", "VARI", "IronOxide"]
             }


LST_params = {"temperature_image_file":        "temperature.tif",
              "soil_mask_file":                "mask_soil.shp",
              "vegetation_cover_mask_file":    "mask_cc.shp"
             }
