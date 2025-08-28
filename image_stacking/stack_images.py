import os
import re
import sys
import skimage
import pandas as pd
import warnings
import argparse
import numpy as np
from natsort import natsorted
import tifffile
import logging
from tifffile import TiffWriter

warnings.filterwarnings("ignore")

# Read in parameters
parser = argparse.ArgumentParser(description="Stack images pipeline")
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument("-i", "--input", action="store", help="Input folder containing CA-videos in tiff format", required=True)
requiredNamed.add_argument("-o", "--output", action="store", help="Output folder for output files")

# Parse arguments
args = parser.parse_args()

# Assign the parsed arguments to variables
param_input_folder = args.input
param_output_folder = args.output if args.output else param_input_folder  # Default output is the input folder if not provided


# Try block to ensure the script finishes with proper logging, even if an error occurs
try:
    # Get folders from input directory
    exp_folders = [folder for folder in os.listdir(param_input_folder) if os.path.isdir(os.path.join(param_input_folder, folder))]
    print(exp_folders)
    if os.path.isdir(param_output_folder):
        if 'image_stacks' in os.listdir(param_output_folder):
            logging.error(f'\nERROR: Output folder exists. Image stacks already generated for {param_input_folder} folder.\n')
            sys.exit()

    os.makedirs(param_output_folder, exist_ok=True)

    # Set up logging to log to both terminal and file
    log_file_path = os.path.join(param_output_folder, "log.txt")
    #os.makedirs(os.path.join(param_input_folder, "image_stacks"), exist_ok=True)

    logging.basicConfig(level=logging.INFO, 
                        format='%(asctime)s - %(levelname)s - %(message)s', 
                        handlers=[logging.StreamHandler(), logging.FileHandler(log_file_path)])
    

    stack_l = []

    # Loop over input folders
    for exp in exp_folders:
        print(exp)
        
        # Get the path to the subfolder (e.g., cv8000) inside each exp folder
        exp_path = os.path.join(param_input_folder, exp)
        subfolder = [f for f in os.listdir(exp_path) if os.path.isdir(os.path.join(exp_path, f))][0]

        subfolder_path = os.path.join(exp_path, subfolder)

        # Get all .tif files from the selected subfolder
        image_paths = natsorted([
            os.path.join(subfolder_path, file)
            for file in os.listdir(subfolder_path)
            if file.lower().endswith('.tif') and subfolder in file
        ])

        well_list = [os.path.basename(path).rsplit('_', 2)[-2] for path in image_paths]

        # Generate df containing image path and well id
        df = pd.DataFrame({'tiff_paths': image_paths, 'well': well_list, 'plate_id': exp})

        # Generate output folder
        output_folder = os.path.join(param_output_folder, 'image_stacks', exp)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Loop unique well ids
        for value in df.well.unique():
            # Get file path of all tiffs associated with the well id and generate a stack out of it
            image_stack = np.array([skimage.io.imread(path, plugin='tifffile') for path in df[df.well == value].tiff_paths.values])

            logging.info(f'Well: {value} | Stack dtype: {image_stack.dtype} | Stack shape: {image_stack.shape}')

            skimage.io.imsave(os.path.join(output_folder, f'{value}.tiff'), image_stack)
            stack_l.append({'plate_id': exp, 'filename': value, 'filepath': os.path.join(output_folder, f'{value}.tiff')})
            

        
        logging.info('')
        logging.info('')


        
    # Save the dataframe
    #df.to_csv(os.path.join(param_output_folder, 'image_stacking.csv'))
    stack_df = pd.DataFrame(stack_l)
    stack_df.to_csv(os.path.join(param_output_folder, 'image_stack_paths.csv'), sep=';', index=None)


    logging.info('Script finished without errors')

except Exception as e:
    logging.error(f"ERROR: {e}")
    sys.exit(1)