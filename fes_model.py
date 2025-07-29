# fes_model.py (Final Version)
import numpy as np
import pyfes
import yaml
import os
import traceback

def get_fes_tide(lons, lats, times, fes_data_path, template_yaml_path):
    """
    Calculates FES tides using the modern YAML configuration method,
    including the necessary longitude conversion.
    """
    try:
        
        os.environ['DATASET_DIR'] = fes_data_path
        
        
        handlers = pyfes.load_config(template_yaml_path)

        
        lons_360 = np.mod(lons, 360)

        
        tide, lp, _ = pyfes.evaluate_tide(handlers['tide'],
                                          times,
                                          lons_360,
                                          lats)
        
        
        pure_tide_cm = tide + lp
        pure_tide_m = pure_tide_cm / 100.0
        
        return pure_tide_m

    except Exception as e:
        print(f"An error occurred in the FES model: {e}")
        traceback.print_exc()
        # Return an array of NaNs that matches the input size in case of failure
        return np.full_like(lons, np.nan)
