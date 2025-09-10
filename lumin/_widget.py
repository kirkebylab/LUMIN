import warnings
#warnings.filterwarnings("ignore")
#warnings.filterwarnings("ignore", category=SyntaxWarning)

import napari
from magicgui import magicgui
from lumin import utils
import numpy as np
from qtpy.QtWidgets import QComboBox, QMessageBox
import pandas as pd
from cellpose import models
from stardist.models import StarDist2D
from lumin.segmentation import run_cellpose, run_stardist, run_manual_selection
from lumin import plot
from lumin import activity
import os
import shutil
import random
from skimage import io
from skimage import exposure
from lumin import preprocess
from qtpy.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import seaborn as sns
import traceback
import pickle
import gc
from datetime import datetime
from time import perf_counter

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"


# List of pre-trained cellpose models
CP_models = ['-- Select --', "cyto3", "cyto2", "cyto", "nuclei", "tissuenet_cp3", 
             "livecell_cp3", "yeast_PhC_cp3", "yeast_BF_cp3", "bact_phase_cp3", 
             "bact_fluor_cp3", "deepbacs_cp3", "cyto2_cp3"]

# Takes widget name as input and it removes
def remove_widget_if_exists(viewer, widget_name):
    for dock_widget in viewer.window._dock_widgets.values():
        widget = dock_widget.widget() 
        if widget.objectName() == widget_name:
            viewer.window.remove_dock_widget(dock_widget)
            break


# Disable or enable widget value 
def disable_enable_value(widget, mode, value_to_disable):
    combo_box = widget.native
    if isinstance(combo_box, QComboBox):
        index = combo_box.findText(value_to_disable)
        if index != -1:
            item = combo_box.model().item(index)
            if item is not None and mode == 'disable':
                item.setEnabled(False)
            if item is not None and mode == 'enable':
                item.setEnabled(True)


# Remove layers from napari window before adding new
def remove_layers(viewer):
    for layer in list(viewer.layers):
        viewer.layers.remove(layer)


# Disable placeholder of widget
def disable_placeholder(widget):
    combo_box = widget.native
    if isinstance(combo_box, QComboBox):
        combo_box.model().item(0).setEnabled(False)



def segmentation_widget():

    # Removes quantification widgets 
    viewer = napari.current_viewer()
    remove_widget_if_exists(viewer, 'Trace quantification')
    remove_widget_if_exists(viewer, 'Baseline quantification')
    remove_widget_if_exists(viewer, 'Baseline normalization')
    remove_widget_if_exists(viewer, 'Spontaneous activity')

    # Widget layout
    @magicgui(
        layout='vertical',
        input_file=dict(widget_type='FileEdit', label='Input file:',value='' , tooltip='Specify image stack metadata file.'),
        project_dir=dict(widget_type='FileEdit',value='', label='Project directory:', mode='d', tooltip='Specify project directory for pipeline output.'),
        seg_label=dict(widget_type='Label', label='<div style="text-align: center; display: block; width: 100%;"><b>———Segmentation settings———</b></div>'),
        selection_mode = dict(widget_type='ComboBox',name = 'selection_mode', label='ROI segmentation mode', value='-- Select --', choices=['-- Select --','Automated', 'Manual selection'],  tooltip='Specify segmentation method.'),

        # nuclear channel first channel of image
        nuclear_stain = dict(widget_type='ComboBox',name='nuclear_stain',  label='Nuclear stain', value='-- Select --', choices=['-- Select --', 'First frame', 'None'],  tooltip='Specify whether the stack contains nuclear stain as first frame.'),
        stain_to_segment = dict(widget_type='ComboBox',name = 'stain_to_segment', label='Stain to segment',value='-- Select --', choices=['-- Select --', 'Nuclear (StarDist)', 'Cytoplasmic (Cellpose)', 'Nuclear (StarDist) and cytoplasmic (Cellpose)'],  tooltip='Segmentation approach to create ROI mask.'),

        # Stardist settings
        label_sd=dict(widget_type='Label', label='<div style="text-align: center; display: block; width: 100%;">StarDist parameters</div>'),
        prob_thresh_sd = dict(widget_type="FloatSpinBox",label="Probability/Score Threshold",min=0.0,max=1.0,step=0.05,value=0.48, tooltip='Higher values lead to fewer segmented objects.'),
        overlap_thresh_sd = dict(widget_type="FloatSpinBox",label="Overlap threshold",min=0.0,max=1.0,step=0.05,value=0.3, tooltip='Higher values allow segmented objects to overlap substantially.'),

        # Cellpose settings
        label_cp=dict(widget_type='Label', label='<div style="text-align: center; display: block; width: 100%;">Cellpose parameters</div>'),
        model_cp = dict(widget_type='ComboBox',name='model_cp', label='Model',value = '-- Select --',  choices=[*CP_models],  tooltip='Specify Cellpose model for segmentation.'),
        diameter_cp = dict(widget_type='SpinBox', name='diameter_cp',label='Diameter', value=30, min=0, max=100,step=5, tooltip='Approximate diameter of cells in pixels to be segmented.'),
        cellprob_threshold_cp = dict(widget_type='FloatSpinBox', name='cellprob_threshold_cp',label='Cell probability threshold', value=0.0, min=-8.0, max=8.0, step=0.5, tooltip='Cell probability threshold (set lower to get more cells and larger cells).'),
        flow_threshold_cp = dict(widget_type='FloatSpinBox', name='flow_threshold_cp',label='Flow threshold', value=0.4, min=0.0, max=5.0, step=0.2, tooltip='Threshold on maximum allowed error (set higher to get more cells, or to zero to turn off).'),

        # Manual segmentation settings
        co_stain = dict(widget_type='ComboBox',name='co_stain',  label='Co-stain', value='-- Select --', choices=['-- Select --', True, False],  tooltip='Specify whether annotating immunolabelled cells.'),
        marker_name = dict(widget_type='LineEdit',name = 'marker_name', label='Marker', value='',  tooltip='Cell population name to be annotated.'),
        annotation_point_size = dict(widget_type='IntSlider', name='annotation_point_size',label='Point size', value=20, min=0, max=70, step=1, tooltip='Point size for manual annotation.'),
        optimize_button  = dict(widget_type='PushButton', text='Test settings on random image', tooltip='Samples a random recording from input to test segmentation parameters.', enabled=False),
        
        # Post-filtering settings
        nuclear_overlap = dict(widget_type='FloatSlider', name='nuclear_overlap',label='Nuclear overlap', value = 0, min=0, max=1, step=0.01, tooltip='Minimum fraction of nuclear mask overlapping with cytoplasmic mask.'),
        nuclear_area = dict(widget_type='RangeSlider', name='nuclear_area',label='Nuclear area', value = (0,1), min=0, max=1, step=1, tooltip='Range of nuclear mask size in pixels.'),
        cell_area = dict(widget_type='RangeSlider', name='cell_area',label='Cell area', value = (0,1), min=0, max=1, step=1, tooltip='Range of cytoplasmic mask size in pixels.'),
        ca_intensity = dict(widget_type='RangeSlider', name='ca_intensity',label='Fluorescence intensity', value = (0,1), min=0, max=1, step=1, tooltip='Average raw pixel intensity across the ROI.'),
        filters_checkbox = dict(widget_type='CheckBox', name='filters_checkbox',label='Apply filters', value = True,  tooltip='Test impact with and without the filters (the checkbox is only for visualization purposes. Values from sliders are applied when running the pipeline).'),

    )

    def widget(input_file, project_dir,  seg_label, selection_mode, nuclear_stain, co_stain, marker_name, stain_to_segment,label_sd,  prob_thresh_sd, overlap_thresh_sd ,label_cp, model_cp,diameter_cp,  cellprob_threshold_cp, flow_threshold_cp, optimize_button, nuclear_overlap,nuclear_area,cell_area,  ca_intensity, filters_checkbox, annotation_point_size):
        pass

    widget.call_button.enabled = False

    widget.native.setObjectName('Cell segmentation')

    class WidgetState:
        def __init__(self):
            self.original_mask = None
            self.original_mask_nuclear = None
            self.original_df = None
            self.skip_update = False 
            self.cp_model_dict = {}
            self.model_sd = None 

    widget_state = WidgetState()

    def _get_stardist_model():
        if widget_state.model_sd is None:
            try:
                widget_state.model_sd = StarDist2D.from_pretrained('2D_versatile_fluo')
            except Exception as e:
                QMessageBox.critical(None, "StarDist Model Error", f"Failed to load StarDist model: {e}")
                return None
        return widget_state.model_sd

    def _get_cellpose_model(model_name_cp):
        if model_name_cp == '-- Select --':
            QMessageBox.warning(None, "Cellpose Model", "Please select a Cellpose model.")
            return None
        if model_name_cp not in widget_state.cp_model_dict:
            try:
                widget_state.cp_model_dict[model_name_cp] = models.CellposeModel(model_type=model_name_cp, gpu=True)
            except Exception as e:
                QMessageBox.critical(None, "Cellpose Model Error", f"Failed to load Cellpose model {model_name_cp}: {e}")
                del widget_state.cp_model_dict[model_name_cp] # Ensure partially loaded model is removed
                return None
        return widget_state.cp_model_dict[model_name_cp]
    



    # Function that loops the input folders and feeds the images to the pipeline functions
    @widget.call_button.clicked.connect
    def _pipeline():
        
        # Pop up box to confirm pipeline execution
        msg_box = QMessageBox()
        msg_box.setWindowTitle("Confirm Action")
        msg_box.setText("Are you sure you want to proceed?")
        msg_box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        msg_box.setDefaultButton(QMessageBox.No)
        
        result = msg_box.exec_()

        try: 
            if result == QMessageBox.Yes:
                start_time = perf_counter()
                print('Running Segmentation and signal extraction pipeline...')

                widget.call_button.enabled = False
                viewer = napari.current_viewer()
                remove_layers(viewer)

                # Extract values from widgets
                input_file = widget.input_file.value
                selection_mode = widget.selection_mode.value
                co_stain = widget.co_stain.value
                if co_stain == '-- Select --': co_stain = False
                image_df, annotated_image_df = utils.parse_input_output(input_file = input_file, project_dir = widget.project_dir.value, selection_mode = selection_mode, co_stain = co_stain)
                project_dir = widget.project_dir.value
                nuclear_stain = widget.nuclear_stain.value
                stain_to_segment = widget.stain_to_segment.value
                prob_thresh_sd = widget.prob_thresh_sd.value
                overlap_thresh_sd = widget.overlap_thresh_sd.value
                model_cp = widget.model_cp.value
                diameter_cp = widget.diameter_cp.value
                cellprob_threshold_cp = widget.cellprob_threshold_cp.value
                flow_threshold_cp = widget.flow_threshold_cp.value
                marker_name = widget.marker_name.value
                annotation_point_size = widget.annotation_point_size.value
                cell_area_min, cell_area_max = widget.cell_area.value
                nuclear_area_min, nuclear_area_max = widget.nuclear_area.value
                intensity_min, intensity_max = widget.ca_intensity.value
                nuclear_overlap = widget.nuclear_overlap.value

                # Generate parameter list to be written to file
                parameter_list = [f'Input file: {input_file}', f'Project directory: {project_dir}', f'ROI selection mode: {selection_mode}',
                                  f'Segmentation - Nuclear stain: {nuclear_stain}']

                # Determine if first frame of stack is nuclear channel or part of the CA video
                if nuclear_stain == 'None':
                    first_frame = 0
                else: 
                    first_frame = 1

                for index, row in image_df.iterrows():

                    # Try except to catch if required columns not there
                    try:
                        image_id = row['image_id']
                        biological_replicate = row['biological_replicate']
                        stimulation = row['stimulation']
                        filename = row['filename']
                        plate_id = row['plate_id']
                        filepath = row['filepath']
                    
                    except Exception as e: 
                        traceback.print_exc()

                    plots_output_dir = os.path.join(project_dir, f'Segmentation/Plots/{plate_id}/{filename}_{stimulation}_{biological_replicate}')
                    mask_output_dir = os.path.join(project_dir, f'Segmentation/Masks/{plate_id}')
                    mask_name = f'{filename}_{stimulation}_{biological_replicate}_final_mask.tiff'



                    if not os.path.exists(plots_output_dir):
                        os.makedirs(plots_output_dir)
                    
                    if not os.path.exists(mask_output_dir):
                        os.makedirs(mask_output_dir)

                    
                    print(f'Processing image {index+1}/{len(image_df)}: Reading from {filepath} and projecting to its max intensity...')

                    # Get image
                    image_projected, image_stack = utils.read_and_project_image(filepath = filepath, first_frame=first_frame)

                    if selection_mode == 'Manual selection':
                        print(f'Running manual segmentation...')

                        if index == 0:
                            parameter_list.extend([
                                f"Segmentation - Co-stain: {co_stain}",
                                *( [f"Segmentation - Marker: {marker_name}"] if co_stain else [] ),
                                f"Segmentation - Point size: {annotation_point_size}"
                            ])
                        if co_stain == True:
                            markers = annotated_image_df.loc[annotated_image_df['image_id'] == image_id, 'marker_name'].tolist()
                        else: markers = []

                        if (image_id not in annotated_image_df['image_id'].to_list() or (co_stain == True and marker_name not in markers)):
                            if image_id  in annotated_image_df['image_id'].to_list():
                                max_label = annotated_image_df[annotated_image_df.image_id == image_id].max_label.values[0]
                            else: max_label = 0

                            image_projected, mask = run_manual_selection(image = image_projected, point_size = annotation_point_size, first_label=max_label) # Segment
                            cell_properties_df = utils.get_cell_properties(mask = mask, image = image_projected) # Get properties of projected image
                            cell_properties_df = utils.extract_raw_traces(image_stack = image_stack[first_frame:], mask = mask, cell_properties_df = cell_properties_df) # Get raw traces
                            if co_stain: cell_properties_df['marker'] = marker_name 
                            
                            # Save results
                            if co_stain: annotated_image_df.loc[len(annotated_image_df.index)] = [image_id,filename, biological_replicate,stimulation,plate_id, marker_name, cell_properties_df.label.max(), filepath,  os.path.join(mask_output_dir,  mask_name)]
                            else: annotated_image_df.loc[len(annotated_image_df.index)] = [image_id,filename, biological_replicate,stimulation,plate_id, cell_properties_df.label.max(), filepath,  os.path.join(mask_output_dir,  mask_name)]
                            annotated_image_df.to_csv(f'{project_dir}/annotated_images.csv', sep=';')

                            io.imsave(os.path.join(mask_output_dir,  mask_name), mask)
                            plot.segmentation(image = image_projected, mask = mask,  title_image = 'Calcium image', title_mask = 'Mask outline', output_path = plots_output_dir, file_name = 'calcium_mask')
                            plot.overlay_labels(image = image_projected, mask = mask, cell_properties_df = cell_properties_df, output_path = plots_output_dir, file_name = 'labelled_mask')

                        else: # if image is annotated, skip it and move to the next one
                            print(f'\nOmitting {image_id}: Image is already annotated.\n')
                            continue

                        del image_stack, image_projected, mask # Clean memory
                        
    
                    elif selection_mode == 'Automated' and stain_to_segment == 'Nuclear (StarDist)' and nuclear_stain == 'First frame':

                        if index == 0:
                            parameter_list.extend([f'Segmentation - Stain to segment: {stain_to_segment}', f'Stardist - Probability/Score Threshold: {prob_thresh_sd}', 
                                                f'Stardist - Overlap threshold: {overlap_thresh_sd}', f'Filtering - Nuclear area: {nuclear_area_min} - {nuclear_area_max}',
                                                f'Filtering - Fluorescence intensity: {intensity_min} - {intensity_max}'])

                        first_frame_image = image_stack[:1][0]

                        print(f'Running StarDist segmentation...')


                        # Run segmentation and get label properties
                        mask = run_stardist(image = first_frame_image, model_sd = _get_stardist_model(), prob_thresh_sd = prob_thresh_sd, overlap_thresh_sd = overlap_thresh_sd)

                        cell_properties_df = utils.get_cell_properties(mask = mask, image = image_projected)
                        cell_properties_df = cell_properties_df.rename(columns={'area': 'nuclear_area'})

                        # Filtering
                        filtered_mask, cell_properties_df = utils.filter_labels(mask = mask, cell_properties_df = cell_properties_df, nuclear_area_min = nuclear_area_min, nuclear_area_max = nuclear_area_max, intensity_min = intensity_min, intensity_max = intensity_max)

                        cell_properties_df = utils.extract_raw_traces(image_stack = image_stack[first_frame:], mask = filtered_mask, cell_properties_df = cell_properties_df)

                        io.imsave(os.path.join(mask_output_dir,  mask_name), filtered_mask)

                        # Visualization 
                        if exposure.is_low_contrast(first_frame_image): 
                            first_frame_image = exposure.equalize_adapthist(first_frame_image, clip_limit=0.03)

                        if exposure.is_low_contrast(image_projected): 
                            image_projected = exposure.equalize_adapthist(image_projected, clip_limit=0.03)

                        # Visualization 
                        plot.segmentation(image = first_frame_image, mask = mask,  title_image = 'Nuclear image', title_mask = 'Mask outline', output_path = plots_output_dir, file_name = 'nuclear_mask')

                        plot.segmentation(image = image_projected, mask = filtered_mask,  title_image = 'Calcium image', title_mask = 'Filtered mask outline', output_path = plots_output_dir, file_name = 'calcium_mask_filtered')

                        plot.overlay_labels(image = image_projected, mask = filtered_mask, cell_properties_df = cell_properties_df, output_path = plots_output_dir, file_name = 'labelled_mask_filtered')

                        del image_stack, image_projected, mask, filtered_mask, first_frame_image

                    elif selection_mode == 'Automated' and stain_to_segment == 'Cytoplasmic (Cellpose)' and nuclear_stain == 'None':
                        if index == 0:
                            parameter_list.extend([f'Segmentation - Stain to segment: {stain_to_segment}', f'Cellpose - Model: {model_cp}', f'Cellpose - Diameter: {diameter_cp}', 
                                                f'Cellpose - Cell probability threshold: {cellprob_threshold_cp}', f'Cellpose - Flow threshold: {flow_threshold_cp}', 
                                                f'Filtering - Cell area: {cell_area_min} - {cell_area_max}',f'Filtering - Fluorescence intensity: {intensity_min} - {intensity_max}'])

                        # Run segmentation and get label properties
                        cp_model_obj = _get_cellpose_model(model_cp)

                        print(f'Running Cellpose segmentation...')

                        mask = run_cellpose(image = image_projected, model_cp = cp_model_obj, diameter_cp = diameter_cp, cellprob_threshold_cp = cellprob_threshold_cp, flow_threshold_cp = flow_threshold_cp)
                        cell_properties_df = utils.get_cell_properties(mask = mask, image = image_projected)

                        cell_properties_df = cell_properties_df.rename(columns={'area': 'cell_area'})

                        # Filtering
                        filtered_mask, cell_properties_df = utils.filter_labels(mask = mask, cell_properties_df=cell_properties_df, cell_area_min = cell_area_min, cell_area_max=cell_area_max, intensity_min = intensity_min, intensity_max = intensity_max)

                        cell_properties_df = utils.extract_raw_traces(image_stack = image_stack, mask = filtered_mask, cell_properties_df = cell_properties_df)

                    
                        io.imsave(os.path.join(mask_output_dir,  mask_name), filtered_mask)

                        if exposure.is_low_contrast(image_projected): 
                            image_projected = exposure.equalize_adapthist(image_projected, clip_limit=0.03)

                        plot.segmentation(image = image_projected, mask = mask,  title_image = 'Calcium image', title_mask = 'Mask outline', output_path = plots_output_dir, file_name = 'calcium_mask')
                        plot.segmentation(image = image_projected, mask = filtered_mask,  title_image = 'Calcium image', title_mask = 'Filtered mask outline', output_path = plots_output_dir, file_name = 'calcium_mask_filtered')
                        plot.overlay_labels(image = image_projected, mask = filtered_mask, cell_properties_df = cell_properties_df, output_path = plots_output_dir, file_name = 'labelled_mask_filtered')

                        del image_stack, image_projected, mask, filtered_mask
                        

                    elif selection_mode == 'Automated' and stain_to_segment == 'Nuclear (StarDist) and cytoplasmic (Cellpose)':
                        if index == 0:
                            parameter_list.extend([f'Segmentation - Stain to segment: {stain_to_segment}',f'Stardist - Probability/Score Threshold: {prob_thresh_sd}', 
                                f'Stardist - Overlap threshold: {overlap_thresh_sd}', f'Cellpose - Model: {model_cp}', f'Cellpose - Diameter: {diameter_cp}', 
                                f'Cellpose - Cell probability threshold: {cellprob_threshold_cp}', f'Cellpose - Flow threshold: {flow_threshold_cp}', 
                                f'Filtering - Nuclear overlap: {nuclear_overlap}', f'Filtering - Nuclear area: {nuclear_area_min} - {nuclear_area_max}',
                                f'Filtering - Cell area: {cell_area_min} - {cell_area_max}',f'Filtering - Fluorescence intensity: {intensity_min} - {intensity_max}'])

                        # Segmentation
                        cp_model_obj = _get_cellpose_model(model_cp)
                        first_frame_image = image_stack[:1][0]
                        print(f'Running StarDist and Cellpose segmentation...')

                        mask_nuclear = run_stardist(image = first_frame_image, model_sd = _get_stardist_model(), prob_thresh_sd = prob_thresh_sd, overlap_thresh_sd = overlap_thresh_sd)
                        mask_cyto = run_cellpose(image = image_projected,  model_cp = cp_model_obj, diameter_cp = diameter_cp, cellprob_threshold_cp = cellprob_threshold_cp, flow_threshold_cp = flow_threshold_cp)

                        cell_properties_nuclear_df = utils.get_cell_properties(mask = mask_nuclear, image = first_frame_image)
                        cell_properties_cyto_df = utils.get_cell_properties(mask = mask_cyto, image = image_projected)

                        overlap_df, filtered_mask_cyto, filtered_mask_nuclear = utils.nuclei_cell_intersection(mask_nuclear = mask_nuclear, df_nuclear = cell_properties_nuclear_df, mask_cyto = mask_cyto, df_cyto = cell_properties_cyto_df)

                        cell_properties_df = utils.get_cell_properties(mask = filtered_mask_cyto, image = image_projected)

                        cell_properties_df = cell_properties_df.merge(overlap_df[['label', 'overlap_fraction_nuclear', 'nuclear_area', 'cell_area', 'nuclear_id']],on='label', how='left') 
                        cell_properties_df = cell_properties_df.drop("area", axis='columns')

                        filtered_mask_cyto, filtered_mask_nuclear, cell_properties_df = utils.filter_labels(mask = filtered_mask_cyto, cell_properties_df=cell_properties_df, mask_nuclear = filtered_mask_nuclear , nuclear_area_min = nuclear_area_min, nuclear_area_max=nuclear_area_max,  cell_area_min = cell_area_min, cell_area_max=cell_area_max, nuclear_overlap=nuclear_overlap,  intensity_min = intensity_min, intensity_max = intensity_max)
                        
                        cell_properties_df = utils.extract_raw_traces(image_stack = image_stack[first_frame:], mask = filtered_mask_cyto, cell_properties_df = cell_properties_df)

                        io.imsave(os.path.join(mask_output_dir,  mask_name), filtered_mask_cyto)

                        # Visualization 
                        if exposure.is_low_contrast(first_frame_image): 
                            first_frame_image = exposure.equalize_adapthist(first_frame_image, clip_limit=0.02)

                        if exposure.is_low_contrast(image_projected): 
                            image_projected = exposure.equalize_adapthist(image_projected, clip_limit=0.008)

                        plot.segmentation(image = first_frame_image, mask = mask_nuclear,  title_image = 'Nuclear image', title_mask = 'Mask outline', output_path = plots_output_dir, file_name = 'nuclear_mask')
                        plot.segmentation(image = image_projected, mask = mask_cyto,  title_image = 'Calcium image', title_mask = 'Mask outline', output_path = plots_output_dir, file_name = 'calcium_mask')
                        plot.segmentation(image = image_projected, mask = filtered_mask_cyto,  title_image = 'Calcium image', title_mask = 'Filtered mask outline', output_path = plots_output_dir, file_name = 'calcium_mask_filtered')
                        plot.overlay_labels(image = image_projected, mask = mask_cyto, mask_nuclear = mask_nuclear, cell_properties_df = cell_properties_df, output_path = plots_output_dir, file_name = 'nuclear_calcium')
                        plot.overlay_labels(image = image_projected, mask = filtered_mask_cyto, mask_nuclear = filtered_mask_nuclear, cell_properties_df = cell_properties_df, output_path =plots_output_dir, file_name = 'nuclear_calcium_filtered')
                        plot.overlay_labels(image = image_projected, mask = filtered_mask_cyto, cell_properties_df = cell_properties_df, output_path = plots_output_dir, file_name = 'labelled_mask_filtered')
                        
                        del image_stack, image_projected, mask_nuclear, mask_cyto, filtered_mask_cyto, filtered_mask_nuclear
                        del cell_properties_nuclear_df, cell_properties_cyto_df, overlap_df

                    print(f'Saving output...\n')

                    # Add metadata
                    cell_properties_df['image_id'] = image_id
                    cell_properties_df['biological_replicate'] = biological_replicate
                    cell_properties_df['stimulation'] = stimulation
                    cell_properties_df['filename'] = filename
                    cell_properties_df['plate_id'] = plate_id
                    cell_properties_df['mask_path'] = os.path.join(mask_output_dir,  mask_name)
                    
                    columns = set(row.keys()).difference(cell_properties_df.columns)

                    for column in columns:
                        cell_properties_df[column] = row[column]

                    if len(cell_properties_df) > 0:

                        # Plot traces
                        #os.makedirs(os.path.join(output_dir,  'Raw_traces'))
                        plot.overlaid_traces(cell_properties_df = cell_properties_df, trace='raw')
                        plt.savefig(os.path.join(plots_output_dir,  'raw_traces.pdf'), bbox_inches='tight')

                        #for img in cell_properties_df.image_id.unique():
                        #    plot.cellwise_traces(cell_properties_df = cell_properties_df[cell_properties_df.image_id == img].copy(), trace='raw', baseline=False, spikes = True, spikes_mode = 'all', output_path = os.path.join(output_dir, 'Raw_traces'))

                        table_dir = f'{project_dir}/Segmentation/Tables'
                        if not os.path.exists(table_dir): os.makedirs(table_dir)

                        with open(os.path.join(table_dir, 'cell_properties_signal_extraction.pkl'), 'ab') as f:
                            pickle.dump(cell_properties_df, f)

                        del cell_properties_df

                        gc.collect()
                    
                    

                df_list = []
                with open(os.path.join(table_dir, 'cell_properties_signal_extraction.pkl'), 'rb') as f:
                    while True:
                        try:
                            df_list.append(pickle.load(f))
                        except EOFError:
                            break

                cell_properties_df = pd.concat(df_list, ignore_index=True)

                cell_properties_df['plate_id'] = cell_properties_df['plate_id'].astype('category')
                cell_properties_df['stimulation'] = cell_properties_df['stimulation'].astype('category')
                cell_properties_df['biological_replicate'] = cell_properties_df['biological_replicate'].astype('category')
                if 'marker' in cell_properties_df.columns:
                    cell_properties_df['marker'] = cell_properties_df['marker'].astype('category')

                cell_properties_df['plate_id_biological_replicate'] = cell_properties_df['plate_id'].astype(str) + '_' + cell_properties_df['biological_replicate'].astype(str)


                if 'Unnamed: 0' in cell_properties_df.columns:
                    cell_properties_df = cell_properties_df.drop(columns='Unnamed: 0')

                cell_properties_df.to_pickle(os.path.join(table_dir, 'cell_properties_signal_extraction.pkl'))
                cell_properties_df.to_csv(os.path.join(table_dir, 'cell_properties_signal_extraction.csv'), sep=';')


                time = datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                end_time = perf_counter()
                runtime_seconds = end_time - start_time

                parameter_list.extend([f'',f'Date and time: {time}', f'Pipeline runtime (s): {runtime_seconds}' ])
                
                with open(os.path.join(project_dir, f'Segmentation', 'segmentation_signal_extraction_config.txt'), 'w') as f:
                    for line in parameter_list:
                        f.write(f'{line}\n')

                print('All images annotated... Stopping the pipeline...\n')
            

            else:
                print("Cancelled.")

        except Exception as e: 
            print("\nAn error occurred:", e)
            traceback.print_exc()
        
        widget.call_button.enabled = True
        # Reset settings
        widget.input_file.value = ''
        widget.project_dir.value = ''
        widget.selection_mode.value = '-- Select --'
        widget.nuclear_stain.value = '-- Select --'
        widget.stain_to_segment.value = '-- Select --'
        widget.prob_thresh_sd.value = 0.48
        widget.overlap_thresh_sd.value = 0.3
        widget.model_cp.value = '-- Select --'
        widget.diameter_cp.value = 30
        widget.cellprob_threshold_cp.value = 0.0
        widget.flow_threshold_cp.value = 0.4
        widget.co_stain.value = '-- Select --'
        widget.marker_name.value = ''
        widget.annotation_point_size.value = 20
        widget.nuclear_overlap.value = 0
        widget.nuclear_area.value = (0,1)
        widget.cell_area.value = (0,1)
        widget.ca_intensity.value = (0,1)
        widget.filters_checkbox.value = True



        

        

    ######################### Input directory and files ##################################
            
    # Hiding widgets
    widget.stain_to_segment.visible = False
    widget.label_sd.visible = False
    widget.prob_thresh_sd.visible = False
    widget.overlap_thresh_sd.visible = False
    widget.label_cp.visible = False
    widget.model_cp.visible = False
    widget.diameter_cp.visible = False
    widget.cellprob_threshold_cp.visible = False
    widget.flow_threshold_cp.visible = False
    widget.optimize_button.visible = False
    widget.nuclear_overlap.visible = False
    widget.cell_area.visible = False
    widget.nuclear_area.visible = False
    widget.ca_intensity.visible = False
    widget.filters_checkbox.visible = False
    widget.annotation_point_size.visible = False
    widget.co_stain.visible = False
    widget.marker_name.visible = False


    '''    # Sets defg
    @widget.selection_mode.changed.connect
    def _reset_to_defaults_mode():
        widget.nuclear_stain.value = '-- Select --'
        widget.stain_to_segment.value = '-- Select --'
        widget.annotation_point_size.value = 15'''

    # Function reseting to default segmentation settings
    def _reset_to_defaults_segment(boolean_reset_stain_to_segment):
        widget.prob_thresh_sd.value = 0.48
        widget.overlap_thresh_sd.value = 0.3
        widget.model_cp.value = '-- Select --'
        widget.diameter_cp.value = 30
        widget.cellprob_threshold_cp.value = 0
        widget.flow_threshold_cp.value = 0.4
        if boolean_reset_stain_to_segment == True:
            widget.stain_to_segment.value = '-- Select --'




    # Display widgets depending of whether automated or manual selection is selected
    @widget.selection_mode.changed.connect
    def _selection_mode_change():
        widget.nuclear_stain.value = '-- Select --'
        widget.stain_to_segment.value = '-- Select --'
        widget.annotation_point_size.value = 20

        disable_placeholder(widget.selection_mode)
        if widget.selection_mode.value == 'Automated' :
            widget.stain_to_segment.visible = True
            widget.annotation_point_size.visible = False
            widget.co_stain.visible = False
            widget.marker_name.visible = False
            widget.call_button.enabled = False


        elif widget.selection_mode.value == 'Manual selection' :
            widget.stain_to_segment.visible = False
            widget.label_sd.visible = False
            widget.prob_thresh_sd.visible = False
            widget.overlap_thresh_sd.visible = False
            widget.label_cp.visible = False
            widget.model_cp.visible = False
            widget.diameter_cp.visible = False
            widget.cellprob_threshold_cp.visible = False
            widget.flow_threshold_cp.visible = False
            widget.optimize_button.visible = False
            widget.nuclear_overlap.visible = False
            widget.cell_area.visible = False
            widget.nuclear_area.visible = False
            #widget.mask_area.visible = False
            widget.ca_intensity.visible = False
            widget.filters_checkbox.visible = False
            widget.annotation_point_size.visible = True
            widget.nuclear_stain.visible = True
            widget.co_stain.visible = True
            widget.call_button.enabled = True

    @widget.co_stain.changed.connect
    def _co_stain_change():
        disable_placeholder(widget.co_stain)
        widget.marker_name.value = ''
        widget.marker_name.visible = widget.co_stain.value == True

        if widget.co_stain.value == True:
            widget.marker_name.visible = True
        else:
            widget.marker_name.visible = False


    # Display settings depeding if nuclear stain in stack or not
    @widget.nuclear_stain.changed.connect
    def _nuclear_stain_change():
        disable_placeholder(widget.nuclear_stain)
        _reset_to_defaults_segment(True)
        
        if widget.nuclear_stain.value == 'None':
            disable_enable_value(widget.stain_to_segment, 'disable', 'Nuclear (StarDist) and cytoplasmic (Cellpose)')
            disable_enable_value(widget.stain_to_segment, 'disable', 'Nuclear (StarDist)')
            disable_enable_value(widget.stain_to_segment, 'enable', 'Cytoplasmic (Cellpose)')

        
        else:
            disable_enable_value(widget.stain_to_segment, 'enable', 'Nuclear (StarDist) and cytoplasmic (Cellpose)')
            disable_enable_value(widget.stain_to_segment, 'enable', 'Nuclear (StarDist)')
            disable_enable_value(widget.stain_to_segment, 'disable', 'Cytoplasmic (Cellpose)')


        widget.label_sd.visible = False
        widget.prob_thresh_sd.visible = False
        widget.overlap_thresh_sd.visible = False
        widget.label_cp.visible = False
        widget.model_cp.visible = False
        widget.diameter_cp.visible = False
        widget.cellprob_threshold_cp.visible = False
        widget.flow_threshold_cp.visible = False
        widget.optimize_button.visible = False
        widget.nuclear_overlap.visible = False
        widget.cell_area.visible = False
        widget.nuclear_area.visible = False
        #widget.mask_area.visible = False
        widget.ca_intensity.visible = False
        widget.filters_checkbox.visible = False

    # Once cellpose model is selected enable to run optimization
    @widget.model_cp.changed.connect
    def _model_cp_changed():
        disable_placeholder(widget.model_cp)
        widget.optimize_button.enabled = True
        widget.call_button.enabled = True
        

    # Controls image frame to segment
    # Resets post processing settings if changed
    @widget.stain_to_segment.changed.connect
    def _stain_to_segment_change():
        disable_placeholder(widget.stain_to_segment)
        _reset_to_defaults_segment(False)
        widget.optimize_button.enabled = False
        widget.call_button.enabled = False


        widget_state.skip_update = True
        widget['nuclear_overlap'].value = 0
        widget['nuclear_area'].value = (0, 1)
        widget['cell_area'].value = (0, 1)
        widget['ca_intensity'].value = (0, 1)
        widget_state.skip_update = False

        # Show stardist parameters
        if widget.stain_to_segment.value == 'Nuclear (StarDist)'  :
            widget.optimize_button.enabled = True
            widget.call_button.enabled = True
            widget.stain_to_segment.visible = True
            widget.nuclear_stain.visible = True
            widget.label_sd.visible = True
            widget.label_cp.visible = False
            widget.prob_thresh_sd.visible = True
            widget.overlap_thresh_sd.visible = True
            widget.model_cp.visible = False
            widget.diameter_cp.visible = False
            widget.cellprob_threshold_cp.visible = False
            widget.flow_threshold_cp.visible = False
            widget.optimize_button.visible = True
            widget.nuclear_overlap.visible = False
            widget.cell_area.visible = False
            widget.nuclear_area.visible = False
            widget.ca_intensity.visible = False
            widget.filters_checkbox.visible = False

        # Show cellpose parameters 
        elif widget.stain_to_segment.value == 'Cytoplasmic (Cellpose)':
            widget.stain_to_segment.visible = True
            widget.label_sd.visible = False
            widget.prob_thresh_sd.visible = False
            widget.overlap_thresh_sd.visible = False
            widget.label_cp.visible = True
            widget.model_cp.visible = True
            widget.diameter_cp.visible = True
            widget.cellprob_threshold_cp.visible = True
            widget.flow_threshold_cp.visible = True
            widget.optimize_button.visible = True
            widget.nuclear_overlap.visible = False
            widget.cell_area.visible = False
            widget.nuclear_area.visible = False
            widget.ca_intensity.visible = False
            widget.filters_checkbox.visible = False

        # Show parameters for dual segmentation
        elif widget.stain_to_segment.value == 'Nuclear (StarDist) and cytoplasmic (Cellpose)':
            widget.stain_to_segment.visible = True
            widget.label_sd.visible = True
            widget.prob_thresh_sd.visible = True
            widget.overlap_thresh_sd.visible = True
            widget.label_cp.visible = True
            widget.model_cp.visible = True
            widget.diameter_cp.visible = True
            widget.cellprob_threshold_cp.visible = True
            widget.flow_threshold_cp.visible = True
            widget.optimize_button.visible = True
            widget.nuclear_overlap.visible = False
            widget.cell_area.visible = False
            widget.nuclear_area.visible = False
            #widget.mask_area.visible = False
            widget.ca_intensity.visible = False
            widget.filters_checkbox.visible = False


    # Function allowing user to optimize segmentation
    @widget.optimize_button.clicked.connect
    def _optimize_segmentation():

        print('Testing segmentation settings on random image...')

        try:
            # Function setting the ranges for post-processing settings
            def set_slider_range(column, scaler):

                if widget[column].max < int(widget_state.original_df[column].max() * scaler):

                    if int(widget_state.original_df[column].max() * scaler) > widget[column].value[1]:
                        widget[column].max = int(widget_state.original_df[column].max() * scaler)

                    else: 
                        widget[column].max = int(widget_state.original_df[column].max() * scaler)
                        widget[column].value = (0, int(widget_state.original_df[column].max() * scaler))

                if widget[column].value[0] == 0 and widget[column].value[1] == 1:
                    widget[column].value = (0, widget[column].max)

            # Extract values from widgets
            selection_mode = widget.selection_mode.value

            image_df, _ = utils.parse_input_output(input_file = widget.input_file.value, project_dir = widget.project_dir.value, selection_mode = selection_mode)
            nuclear_stain = widget.nuclear_stain.value
            stain_to_segment = widget.stain_to_segment.value
            prob_thresh_sd = widget.prob_thresh_sd.value
            overlap_thresh_sd = widget.overlap_thresh_sd.value
            model_cp = widget.model_cp.value
            diameter_cp = widget.diameter_cp.value
            cellprob_threshold_cp = widget.cellprob_threshold_cp.value
            flow_threshold_cp = widget.flow_threshold_cp.value

            # Reads the scalers from the settings file
            config_df = pd.read_csv('configs/settings.csv', sep=None)
            nuclear_area_scaler = float(config_df[config_df.parameter == 'nuclear_area_scaler'].value.values[0])
            cell_area_scaler = float(config_df[config_df.parameter == 'cell_area_scaler'].value.values[0])
            intensity_scaler = float(config_df[config_df.parameter == 'intensity_scaler'].value.values[0])

            viewer = napari.current_viewer()

            remove_layers(viewer)

            # Reading the input image
            if nuclear_stain == 'None':
                first_frame = 0
            else: 
                first_frame = 1

            # Draw random image
            i = random.randint(0, len(image_df)-1)
            filepath = image_df.iloc[i].filepath
            print(f'Reading image from {filepath} and projecting to its max intensity...')

            image_projected, image_stack = utils.read_and_project_image(filepath = filepath, first_frame = first_frame)

            # Calling segmentation functions depending of settings specified by the user
            if selection_mode == 'Automated' and stain_to_segment == 'Nuclear (StarDist)':

                widget.nuclear_overlap.visible = False
                widget.cell_area.visible = False
                widget.nuclear_area.visible = True
                widget.ca_intensity.visible = True
                widget.filters_checkbox.visible = True

                first_frame_image = image_stack[:1][0]
                print(f'Running StarDist segmentation...')

                mask = run_stardist(image = first_frame_image, model_sd = _get_stardist_model(), prob_thresh_sd = prob_thresh_sd, overlap_thresh_sd = overlap_thresh_sd)

                viewer.add_image(image_stack, name='CA video')
                viewer.add_image(first_frame_image, name='Nuclear')
                viewer.add_image(image_projected, name='Cyto')
                viewer.add_labels(mask, name='Mask final')

                widget_state.original_mask = mask
                widget_state.original_df = utils.get_cell_properties(mask = widget_state.original_mask, image = image_projected)
                widget_state.original_df = widget_state.original_df.rename(columns={'area': 'nuclear_area'})

                set_slider_range('nuclear_area', nuclear_area_scaler)
                set_slider_range('ca_intensity', intensity_scaler )
                

            elif selection_mode == 'Automated' and stain_to_segment == 'Cytoplasmic (Cellpose)':

                widget.nuclear_overlap.visible = False
                widget.cell_area.visible = True
                widget.nuclear_area.visible = False
                widget.filters_checkbox.visible = True
                widget.ca_intensity.visible = True

                cp_model_obj = _get_cellpose_model(model_cp)
                print(f'Running Cellpose segmentation...')
                mask = run_cellpose(image = image_projected, model_cp = cp_model_obj, diameter_cp = diameter_cp, cellprob_threshold_cp = cellprob_threshold_cp, flow_threshold_cp = flow_threshold_cp)
                
                viewer.add_image(image_stack, name='CA video')
                viewer.add_image(image_projected, name='Cyto')
                viewer.add_labels(mask, name='Mask final')

                widget_state.original_mask = mask
                widget_state.original_df = utils.get_cell_properties(mask = widget_state.original_mask, image = image_projected)

                widget_state.original_df = widget_state.original_df.rename(columns={'area': 'cell_area'})

                set_slider_range('cell_area', cell_area_scaler)
                set_slider_range('ca_intensity', intensity_scaler )
                viewer.dims.current_step = (0,)


            
            elif selection_mode == 'Automated' and stain_to_segment == 'Nuclear (StarDist) and cytoplasmic (Cellpose)':

                widget.nuclear_overlap.visible = True
                widget.cell_area.visible = True
                widget.nuclear_area.visible = True
                widget.ca_intensity.visible = True
                widget.filters_checkbox.visible = True

                cp_model_obj = _get_cellpose_model(model_cp)

                first_frame_image = image_stack[:1][0]
                print(f'Running Stardist and Cellpose segmentation...')

                mask_nuclear = run_stardist(image = first_frame_image, model_sd = _get_stardist_model(), prob_thresh_sd = prob_thresh_sd, overlap_thresh_sd = overlap_thresh_sd)
                mask_cyto = run_cellpose(image = image_projected,  model_cp = cp_model_obj, diameter_cp = diameter_cp, cellprob_threshold_cp = cellprob_threshold_cp, flow_threshold_cp = flow_threshold_cp)
            
                viewer.add_image(image_stack, name='CA video')
                viewer.add_image(first_frame_image, name='Nuclear')
                viewer.add_image(image_projected, name='Cyto')
                viewer.add_labels(mask_nuclear, name='Mask nuclear', colormap={label: [0.0, 1.0, 1.0, 1.0] for label in np.unique(mask_nuclear) if label != 0})
                viewer.layers['Mask nuclear'].contour = 3
                viewer.layers['Mask nuclear'].opacity = 1

                viewer.add_labels(mask_cyto, name='Mask cyto')

                original_df_nuclear = utils.get_cell_properties(mask = mask_nuclear, image = first_frame_image)
                original_df_cyto = utils.get_cell_properties(mask = mask_cyto, image = image_projected)

                overlap_df, filtered_mask_cyto, filtered_mask_nuclear = utils.nuclei_cell_intersection(mask_nuclear = mask_nuclear, df_nuclear = original_df_nuclear, mask_cyto = mask_cyto, df_cyto = original_df_cyto)
                
                viewer.add_labels(filtered_mask_nuclear, name='Mask nuclear filtered', colormap={label: [1.0, 0.5, 0.0, 1.0] for label in np.unique(filtered_mask_nuclear) if label != 0})
                viewer.layers['Mask nuclear filtered'].contour = 3
                viewer.layers['Mask nuclear filtered'].opacity = 1

                viewer.add_labels(filtered_mask_cyto, name='Mask final')

                widget_state.original_df = utils.get_cell_properties(mask = filtered_mask_cyto, image = image_projected)

                widget_state.original_df = widget_state.original_df.merge(overlap_df[['label', 'overlap_fraction_nuclear', 'nuclear_area', 'cell_area', 'nuclear_id']],on='label', how='left') 

                widget_state.original_df = widget_state.original_df.drop("area", axis='columns')


                widget_state.original_mask = filtered_mask_cyto
                widget_state.original_mask_nuclear = filtered_mask_nuclear

                viewer.layers['Mask nuclear filtered'].visible = False
                viewer.layers['Mask cyto'].visible = False
                viewer.layers['Mask nuclear'].visible = False

                set_slider_range('nuclear_area', nuclear_area_scaler )
                set_slider_range('cell_area', cell_area_scaler )
                set_slider_range('ca_intensity', intensity_scaler )

                viewer.dims.current_step = (0,)

            print(f'Segmentation done, try post-filtering settings to optimize the segmentation... and test the segmentation with another image...\n')

            apply_intensity_area_filters()

        except Exception as e: 
            print("\nAn error occurred:", e)
            traceback.print_exc()




    def apply_intensity_area_filters():
        if widget_state.skip_update:
            return
        
        widget.filters_checkbox.value = True

        viewer = napari.current_viewer()
        cell_properties_df = widget_state.original_df.copy()
        mask = widget_state.original_mask.copy()
        nuclear_area_min, nuclear_area_max = widget.nuclear_area.value
        intensity_min, intensity_max = widget.ca_intensity.value
        cell_area_min, cell_area_max = widget.cell_area.value
        nuclear_overlap = widget.nuclear_overlap.value

        if widget.stain_to_segment.value == 'Nuclear (StarDist)':
            mask, _ = utils.filter_labels( mask = mask, cell_properties_df=cell_properties_df, nuclear_area_min=nuclear_area_min, nuclear_area_max = nuclear_area_max, intensity_min=intensity_min, intensity_max = intensity_max)

        elif widget.stain_to_segment.value == 'Cytoplasmic (Cellpose)':
            mask,  _ = utils.filter_labels( mask = mask, cell_properties_df=cell_properties_df ,cell_area_min = cell_area_min, cell_area_max=cell_area_max , intensity_min=intensity_min, intensity_max = intensity_max)

        # Determine which labels to keep
        elif widget.stain_to_segment.value == 'Nuclear (StarDist) and cytoplasmic (Cellpose)':
            mask, mask_nuclear, _ = utils.filter_labels( mask = mask, cell_properties_df=cell_properties_df, mask_nuclear = widget_state.original_mask_nuclear , nuclear_area_min = nuclear_area_min, nuclear_area_max=nuclear_area_max,  cell_area_min = cell_area_min, cell_area_max=cell_area_max, nuclear_overlap=nuclear_overlap, intensity_min=intensity_min, intensity_max = intensity_max )

            viewer.layers['Mask nuclear filtered'].data = mask_nuclear

        viewer.layers['Mask final'].data = mask

    # Connect both sliders to the same handler
    widget.ca_intensity.changed.connect(apply_intensity_area_filters)
    widget.nuclear_area.changed.connect(apply_intensity_area_filters)
    widget.nuclear_overlap.changed.connect(apply_intensity_area_filters)
    widget.cell_area.changed.connect(apply_intensity_area_filters)

    # Functionality for filters_checkbox
    @widget.filters_checkbox.changed.connect
    def _filters_checkbox_change():
        if widget.filters_checkbox.value == False:
            viewer = napari.current_viewer()

            if widget.stain_to_segment.value == 'Nuclear (StarDist) and cytoplasmic (Cellpose)':

                viewer.layers['Mask nuclear filtered'].data = widget_state.original_mask_nuclear.copy()
                   
            viewer.layers['Mask final'].data =  widget_state.original_mask.copy()

        else:
            apply_intensity_area_filters()

    return widget


class PlotViewer(QWidget):
    def __init__(self, figures):
        super().__init__()
        self.figures = figures
        self.current_index = 0

        # Setup UI layout
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        # Use one persistent canvas attached to a single figure
        self.fig, self.ax = plt.subplots( figsize=(10, 12))
        self.canvas = FigureCanvas(self.fig)
        self.layout.addWidget(self.canvas)

        # Buttons and label
        button_layout = QHBoxLayout()
        self.prev_button = QPushButton("Previous")
        self.next_button = QPushButton("Next")
        self.index_label = QLabel()
        button_layout.addWidget(self.prev_button)
        button_layout.addWidget(self.index_label)
        button_layout.addWidget(self.next_button)
        self.layout.addLayout(button_layout)

        # Connect buttons
        self.prev_button.clicked.connect(self.prev_plot)
        self.next_button.clicked.connect(self.next_plot)

        # Initial plot
        self.update_plot()

    def update_plot(self):
        # Clear the current figure completely (this resets zoom but it's intended)
        self.fig.clf()

        # Copy the entire contents of the target figure onto this one
        source_fig = self.figures[self.current_index]
        num_axes = len(source_fig.axes)

        for i, source_ax in enumerate(source_fig.axes):
            ax = self.fig.add_subplot(num_axes, 1, i + 1)

            for line in source_ax.lines:
                ax.plot(
                    *line.get_data(),
                    color=line.get_color(),
                    linestyle=line.get_linestyle(),
                    linewidth=line.get_linewidth(),
                    marker=line.get_marker(),
                )

            ax.set_xlim(source_ax.get_xlim())
            ax.set_ylim(source_ax.get_ylim())
            ax.tick_params(axis='y', labelsize=6)
            ax.tick_params(axis='x', labelsize=8)

            if i < num_axes - 1:
                ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            else:
                ax.tick_params(axis='x', labelsize=8)

            ax.yaxis.tick_right()
            ax.set_ylabel(source_ax.get_ylabel(), rotation=0, fontsize=6, labelpad=5)
            ax.set_title(source_ax.get_title())
        
        self.fig.tight_layout( h_pad=0.3, w_pad=0.5)
        self.canvas.draw()
        self.index_label.setText(f"{self.current_index + 1} / {len(self.figures)}")

    def next_plot(self):
        if self.current_index < len(self.figures) - 1:
            self.current_index += 1
            self.update_plot()

    def prev_plot(self):
        if self.current_index > 0:
            self.current_index -= 1
            self.update_plot()


class PlotViewerBaseline(QWidget):
    def __init__(self, figure):
        super().__init__()
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.canvas = FigureCanvas(figure)
        self.layout.addWidget(self.canvas)

def get_colors(cell_properties_df, project_dir):
    unique_stimulations = cell_properties_df["stimulation"].cat.categories

    # Read user specified palette
    if os.path.isfile(os.path.join(project_dir, 'colors.csv')):
        colors = pd.read_csv(os.path.join(project_dir, 'colors.csv'), sep=None)
        palette = dict(zip(colors['stimulation'], colors['color']))
        
        # Check all categories has color specified
        for cond in unique_stimulations:
            if cond not in palette.keys():
                raise ValueError(f'Color group {cond} is not present in color table')

    else:
        colors = sns.color_palette(n_colors=len(unique_stimulations))
        hex_colors = ["#{:02X}{:02X}{:02X}".format(int(color[0]*255), int(color[1]*255), int(color[2]*255)) for color in colors]
        palette = dict(zip(unique_stimulations, hex_colors))

    palette = {k: v for k, v in palette.items() if k in cell_properties_df['stimulation'].cat.categories} 
    #palette = {cat:palette[cat] for cat in cell_properties_df.stimulation.cat.categories} # Reorder palette according to the categories

    return palette

def get_cell_properties_df(project_dir):
    if not os.path.isfile(os.path.join(project_dir, 'Segmentation', 'Tables', 'cell_properties_signal_extraction.pkl')):
        raise ValueError(f"\ncell_properties_signal_extraction.pkl file not found from {os.path.join(project_dir, 'Tables')} ")
    else:
        return pd.read_pickle(os.path.join(project_dir, 'Segmentation', 'Tables', 'cell_properties_signal_extraction.pkl'))
    


# Reset settings
def reset_settings_sc_analysis(widget):
    # project_dir analysis_mode activity_type control_condition 
    widget.stimulation_frame.value = 0
    widget.kcl_frame.value = -1
    widget.sliding_window_size.value = 75
    widget.percentile_threshold.value = 15
    widget.spike_prominence_threshold.value = 0
    widget.spike_amplitude_width_ratio.value = 0
    widget.analysis_window_start.value = 0
    widget.analysis_window_end.value = 0
    widget.baseline_std_threshold.value = 0
    widget.imaging_interval.value = 0
    widget.n_clusters.value = 5




# ANALYSIS
def single_cell_widget():

    viewer = napari.current_viewer()
    remove_layers(viewer)
    remove_widget_if_exists(viewer, 'Cell segmentation')

    @magicgui(
        layout='vertical',
        project_dir=dict(widget_type='FileEdit',value='', label='Project directory:', mode='d', tooltip='Specify project directory for pipeline output'),
        analysis_mode = dict(widget_type='ComboBox', name = 'analysis_mode', label='Analysis mode', value='-- Select --', choices=['-- Select --','Compound-evoked activity', 'Spontaneous activity'],  tooltip='Experimental type.'),
        activity_type = dict(widget_type='ComboBox', name = 'activity_type', label='Activity type', value='-- Select --', choices=['-- Select --','Spontaneous', 'Baseline change'],  tooltip='Type of anticipated cellular activity.'),
        control_condition = dict(widget_type='LineEdit',name = 'control_condition', label='Control condition', value='',  tooltip='Name of stimulation to be used as control'),
        norm_label=dict(widget_type='Label', label='<div style="text-align: center; display: block; width: 100%;"><b>———Normalization settings———</b></div>'),
        normalization_mode = dict(widget_type='ComboBox',name = 'normalization_mode', label='Normalization', value='-- Select --', choices=['-- Select --','Sliding window', 'Pre-stimulus window'],  tooltip='Normalization method (use pre-stimulus window only when recording contains a stable baseline signal).'),
        stimulation_frame = dict(widget_type="SpinBox",label="Stimulation frame",value=0, step = 1, tooltip='Image frame for stimulation administration.'),
        kcl_frame = dict(widget_type="SpinBox",label="KCl stimulation frame",value=-1, step=1, tooltip='Image frame for KCl administration. If no KCl was added specify -1.'), # make dynamic
        sliding_window_size = dict(widget_type="SpinBox",label="Sliding window size",value=75, step=1, min=1, tooltip='Sliding window size in frames'), # make dynamic
        percentile_threshold = dict(widget_type='IntSlider', name='percentile_threshold',label='Percentile threshold', value=15, min=1, max=50, step=1, tooltip='Percentile of fluorescence values classified as baseline. Increase the percentile threshold for low activity recordings.'),
        spike_label=dict(widget_type='Label', label='<div style="text-align: center; display: block; width: 100%;"><b>———Trace quantification———</b></div>'),
        spike_prominence_threshold = dict(widget_type='FloatSpinBox', name='spike_prominence_threshold',label='Prominence', value=0, step=0.01, tooltip = 'Threshold to identify peaks based on how much local maximum stand out from the surrounding baseline.'),
        spike_amplitude_width_ratio = dict(widget_type='FloatSpinBox', name='spike_amplitude_width_ratio',label='Amplitude width ratio', value=0, step=0.005, tooltip = 'Parameter to exclude high width low amplitude peaks (If value is 0 nothing is excluded).'),
        analysis_window_start=dict(widget_type='SpinBox', label='Analysis window start', value=0, tooltip = 'Specify analysis window start frame.'),
        analysis_window_end=dict(widget_type='SpinBox', label='Analysis window end', value=0, tooltip = 'Specify analysis window end frame.'),
        baseline_std_threshold = dict(widget_type='SpinBox', name='baseline_std_threshold',label='Standard deviation threshold', value=0, step=0.1, tooltip = 'Specify how many standard deviations AUC of a cell needs to exceed the mean AUC of the control to be classified as repsonding.'),
        imaging_interval=dict(widget_type='FloatSpinBox', label='Imaging interval (s)', value=0, step=0.1, tooltip = 'Imaging interval in seconds.'),
        downstream_label=dict(widget_type='Label', label='<div style="text-align: center; display: block; width: 100%;"><b>———Downstream analysis———</b></div>'),
        n_clusters = dict(widget_type='IntSlider', name='n_clusters',label='Number of clusters', value=5, min=1, max=10, step=1, tooltip='Number of clusters for k-means clustering.'),
        optimize_button  = dict(widget_type='PushButton', text='Test settings on random recording', tooltip='Samples a random recording from input to test quantification parameters.', enabled=True),

    )

    def widget(project_dir, analysis_mode, activity_type, control_condition,  norm_label, normalization_mode, stimulation_frame,  sliding_window_size, percentile_threshold  , spike_label, spike_prominence_threshold, spike_amplitude_width_ratio, analysis_window_start, analysis_window_end, baseline_std_threshold, imaging_interval,kcl_frame, downstream_label, n_clusters, optimize_button):
        pass

    widget.native.setObjectName('Trace quantification')

    class WidgetState:
        def __init__(self):
            self.cell_properties_df = None
            #self.control_condition = None
            #self.treatment_conditions = None

    widget_state = WidgetState()

    @widget.call_button.clicked.connect
    def _run_analysis():
        msg_box = QMessageBox()
        msg_box.setWindowTitle("Confirm Action")
        msg_box.setText("Are you sure you want to proceed?")
        msg_box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        msg_box.setDefaultButton(QMessageBox.No)
        
        result = msg_box.exec_()
        try: 
            if result == QMessageBox.Yes:
                start_time = perf_counter()

                viewer = napari.current_viewer()
                widget.call_button.enabled = False

                # Extract values from widgets
                project_dir = widget.project_dir.value
                analysis_mode = widget.analysis_mode.value
                activity_type = widget.activity_type.value
                normalization_mode = widget.normalization_mode.value    
                control_condition =  widget.control_condition.value
                #analysis_setup = widget.analysis_setup.value
                stimulation_frame = widget.stimulation_frame.value    
                kcl_frame = widget.kcl_frame.value
                sliding_window_size = widget.sliding_window_size.value
                percentile_threshold = widget.percentile_threshold.value
                spike_prominence_threshold = widget.spike_prominence_threshold.value
                spike_amplitude_width_ratio = widget.spike_amplitude_width_ratio.value
                baseline_std_threshold = widget.baseline_std_threshold.value
                analysis_window_start = widget.analysis_window_start.value  
                analysis_window_end = widget.analysis_window_end.value 
                imaging_interval = widget.imaging_interval.value
                n_clusters = widget.n_clusters.value

                plt.close()
                plt.clf()

                widget_state.cell_properties_df = get_cell_properties_df(project_dir)
                cell_properties_df = widget_state.cell_properties_df.copy()

                if control_condition not in cell_properties_df.stimulation.tolist():
                    raise ValueError(f"Condition '{control_condition}' not in input data")
                
                if not imaging_interval > 0:
                    raise ValueError(f"Imaging interval needs to be above 0")

                for col in cell_properties_df.select_dtypes(include='category').columns:
                    cell_properties_df[col] = cell_properties_df[col].cat.remove_unused_categories()

                # Generate parameter list to be written to file
                parameter_list = [f'Project directory: {project_dir}', f'Analysis mode: {analysis_mode}',
                                f'Control condition: {control_condition}', f'Normalization mode: {normalization_mode}']
                
                ##plot_dir = f'{project_dir}/Quantification/Plots'
                ##table_dir = f'{project_dir}/Quantification/Tables'

                if os.path.exists(os.path.join(project_dir, 'Quantification')): shutil.rmtree(os.path.join(project_dir, 'Quantification'), ignore_errors=True)


                if 'marker' not in cell_properties_df.columns or analysis_mode == 'Spontaneous activity' or (analysis_mode == 'Compound-evoked activity' and activity_type == 'Spontaneous'):
                    plot_dir = f'{project_dir}/Quantification/Plots'
                    if not os.path.exists(plot_dir): os.makedirs(plot_dir)

                table_dir = f'{project_dir}/Quantification/Tables'
                if not os.path.exists(table_dir): os.makedirs(table_dir)


                print('Starting signal quantification...')
                
                # Normalization
                if normalization_mode == 'Sliding window':
                    print('Running sliding window normalization...')
                    parameter_list.extend([f'Normalization - Sliding window size: {sliding_window_size}', f'Normalization - Percentile threshold: {percentile_threshold}'])
                    cell_properties_df = preprocess.sliding_window(cell_properties_df = cell_properties_df, sliding_window_size = sliding_window_size, percentile_threshold = percentile_threshold)

                elif normalization_mode == 'Pre-stimulus window':
                    print('Running pre-stimulus window normalization...')
                    parameter_list.extend([f'Normalization - Stimulation frame: {stimulation_frame}'])
                    cell_properties_df = preprocess.pre_stimulation(cell_properties_df = cell_properties_df, stimulation_frame = stimulation_frame)

                palette = get_colors(cell_properties_df, project_dir)

                # Sort the df for plotting
                stimulations = [control_condition] + [c for c in palette.keys() if c != control_condition]
                cell_properties_df['stimulation'] = pd.Categorical(cell_properties_df['stimulation'], categories=stimulations, ordered=True)
                cell_properties_df = cell_properties_df.sort_values('stimulation')

                palette = {cat:palette[cat] for cat in cell_properties_df.stimulation.cat.categories} # Order the palette
                cell_properties_df['colors'] = cell_properties_df['stimulation'].map(palette) # Map colors to stimulation groups

                if 'marker' not in cell_properties_df.columns:

                    # Plot number of cells
                    cell_count_df = cell_properties_df.groupby(['biological_replicate',  'image_id', 'stimulation'], observed=True).size().reset_index(name='count') 
                    ax_barplot = plot.all_conditions_barplot(dataframe = cell_count_df, palette = palette, xcolumn='biological_replicate', ycolumn='count', hue='stimulation' )  
                    plt.ylabel('Cells per image', fontsize='xx-large')
                    plt.xlabel('')
                    plt.title('')
                    plt.savefig(os.path.join(plot_dir, 'barplot_cell_count_per_well.pdf'), bbox_inches='tight')
                    plt.close()              

                    # Quantify and plot KCL response
                    if kcl_frame > 0:
                        cell_properties_df = utils.compute_auc(cell_properties_df = cell_properties_df,  start_frame = kcl_frame, column = 'AUC_kcl')
                        mean_kcl_df = cell_properties_df.groupby(['image_id', 'stimulation', 'biological_replicate'], observed=True)['AUC_kcl'].mean().reset_index() 
                        ax_barplot = plot.all_conditions_barplot(dataframe = mean_kcl_df, palette = palette, xcolumn='biological_replicate', ycolumn='AUC_kcl', hue='stimulation' )  
                        plt.ylabel('Mean KCl AUC', size='xx-large')
                        plt.xlabel('')
                        plt.savefig(os.path.join(plot_dir, 'barplot_kcl_auc_per_well.pdf'), bbox_inches='tight')
                        plt.close()   

                        mean_kcl_df.to_csv(os.path.join(table_dir, 'kcl_auc_per_well.csv'), sep=';')
                        mean_kcl_df.to_pickle(os.path.join(table_dir, 'kcl_auc_per_well.pkl'))


                def spontaneous_activity_output():
                    # Visualize baseline estimation and peak detection
                    #plot_list = plot.cellwise_traces(cell_properties_df = cell_properties_df, trace='raw',baseline=True)
                    #plot_list_spikes = plot.cellwise_traces(cell_properties_df = cell_properties_df, trace='dff',baseline=False, spikes=True, spikes_mode='all')

                    # Generate output folder
                    #output_dir_quantification = os.path.join(project_dir, 'Quantification_output')

                    
                    #category_string = '_'.join(cell_properties_df.stimulation.cat.categories.tolist())
                    #category_string = category_string.replace('_', '_vs_', 1)
                    #output_dir_quantification = os.path.join(project_dir, 'Quantification', f'{category_string}_output')



                    #if os.path.exists(plot_dir): shutil.rmtree(plot_dir, ignore_errors=True)
                    #if os.path.exists(table_dir): shutil.rmtree(table_dir, ignore_errors=True)


                    #output_dir_quantification = os.path.join(project_dir, 'Quantification')


                    #os.makedirs(output_dir_quantification)
                    #os.makedirs(f'{output_dir_quantification}/Tables')

                    # Dicts defining ylabels and titles for swarmplots
                    ylabel_dict = {'amplitude': '$\Delta F/F_0$', 'prominence':'$\Delta F/F_0$', 'frequency':'Spikes/min', 'width':'Time (s)','rise_time':'Time (s)', 'decay_time':'Time (s)'}
                    title_dict = {'amplitude': 'Amplitude','prominence':'Prominence', 'frequency':'Frequency', 'width':'Width','rise_time':'Rise time','decay_time':'Decay time', 'cluster':'Cluster'}

                    print('Plotting and saving output...')
                    for exp_replicate in cell_properties_df.plate_id_biological_replicate.unique():

                        output_dir_exp_replicate = os.path.join(plot_dir,exp_replicate)
                        os.makedirs(output_dir_exp_replicate)

                        subset_df = cell_properties_df[cell_properties_df.plate_id_biological_replicate == exp_replicate].copy()
                        subset_df["stimulation"] = subset_df["stimulation"].cat.remove_unused_categories()

                        # Heatmap
                        if 'marker' not in cell_properties_df.columns:
                            ax_heatmap = plot.heatmap(cell_properties_df=subset_df,   imaging_interval=imaging_interval,  cmap = 'plasma', palette = palette, minmax_bool=False)
                            plt.savefig(os.path.join(output_dir_exp_replicate, f'heatmap.pdf'),  bbox_inches='tight')
                            plt.close()

                        # Plot properties for each plate
                        for property in ['amplitude', 'width', 'rise_time', 'decay_time', 'frequency', 'prominence']:
                            if 'marker' in cell_properties_df.columns:
                                ax_beeswarm = plot.beeswarm(cell_properties_df = subset_df[subset_df.frequency > 0].copy(), control_condition = control_condition, palette=palette, y=property, x='marker', hue='stimulation',hue_separation=0.3, separation_between_plots=2, max_plot_width=0.3)

                            else:
                                ax_beeswarm = plot.beeswarm(cell_properties_df = subset_df[subset_df.frequency > 0].copy(), control_condition = control_condition, palette=palette, y=property, x='stimulation')
                            plt.ylabel(ylabel_dict[property])
                            plt.title(title_dict[property])
                            plt.savefig(os.path.join(output_dir_exp_replicate, f'beeswarm_{property}.pdf'),  bbox_inches='tight')
                            plt.close()              

                    # Count number of active cells, plot barplot and save tables
                    response_perc_well_df, response_perc_rep_df = utils.percentage_responding(cell_properties_df = cell_properties_df, analysis_type='spontaneous')
                    #response_perc_well_df, response_perc_rep_df = utils.percentage_responding_spontaneous(cell_properties_df = cell_properties_df)

                    if 'marker' not in cell_properties_df.columns:
                        ax_barplot_well = plot.all_conditions_barplot(dataframe=response_perc_well_df, palette=palette, ycolumn="proportion_responding", xcolumn='biological_replicate', hue = 'stimulation')
                        ax_barplot_well.set_ylabel(r"% active cells (Ca$^{2+}$)", fontsize='xx-large')
                        plt.xlabel('')
                        plt.title('')
                        plt.savefig(os.path.join(plot_dir, 'active_cells_per_well_barplot.pdf'), bbox_inches='tight')
                        plt.close()
                    
                    if 'marker' in cell_properties_df.columns:
                        ax_barplot_rep = plot.all_conditions_barplot(dataframe=response_perc_rep_df, palette=palette, ycolumn="proportion_responding", xcolumn='marker', hue='stimulation')
                    else: 
                        ax_barplot_rep = plot.all_conditions_barplot(dataframe=response_perc_rep_df, palette=palette, ycolumn="proportion_responding", xcolumn='stimulation')
                    ax_barplot_rep.set_ylabel(r"% active cells (Ca$^{2+}$)", fontsize='xx-large')
                    plt.xlabel('')
                    plt.title('')
                    plt.savefig(os.path.join(plot_dir, 'active_cells_per_replicate_barplot.pdf'), bbox_inches='tight')
                    plt.close()

                    response_perc_well_df.to_csv(os.path.join(table_dir, 'active_cells_per_well.csv'), sep=';')
                    response_perc_well_df.to_pickle(os.path.join(table_dir, 'active_cells_per_well.pkl'))
                    response_perc_rep_df.to_csv(os.path.join(table_dir, 'active_cells_per_biological_replicate.csv'), sep=';')
                    response_perc_rep_df.to_pickle(os.path.join(table_dir, 'active_cells_per_biological_replicate.pkl'))

                
                    # Plot traces one by one
                    '''for img in cell_properties_df.image_id.unique():

                        plate = cell_properties_df[cell_properties_df.image_id == img]['plate_id'].iloc[0]
                        filename = cell_properties_df[cell_properties_df.image_id == img]['filename'].iloc[0]
                        stimulation = cell_properties_df[cell_properties_df.image_id == img]['stimulation'].iloc[0]
                        biological_replicate = cell_properties_df[cell_properties_df.image_id == img]['biological_replicate'].iloc[0]

                        output_norm_dir = os.path.join(project_dir, 'Per_image_output', plate, f'{filename}_{stimulation}_{biological_replicate}', 'Baseline_traces' )
                        output_spikes_dir = os.path.join(project_dir, 'Per_image_output', plate, f'{filename}_{stimulation}_{biological_replicate}', 'Spikes_traces' )
                        
                        if os.path.exists(output_norm_dir): shutil.rmtree(output_norm_dir, ignore_errors=True)
                        os.makedirs(output_norm_dir)

                        if os.path.exists(output_spikes_dir): shutil.rmtree(output_spikes_dir, ignore_errors=True)
                        os.makedirs(output_spikes_dir)

                        #plot.cellwise_traces(cell_properties_df = cell_properties_df[cell_properties_df.image_id == img].copy(), trace='raw', baseline=True, spikes = False, spikes_mode = 'all', output_path = output_norm_dir) 
                        #plot.cellwise_traces(cell_properties_df = cell_properties_df[cell_properties_df.image_id == img].copy(), trace='dff', baseline=False, spikes = True, spikes_mode = 'all', output_path = output_spikes_dir) 
                    '''
                    cell_properties_filtered_df = cell_properties_df[cell_properties_df.frequency > 0].copy()


                    # Scale spike properties for pca, clustering and heatmap, 
                    cell_properties_filtered_df = utils.scale_spike_properties(cell_properties_df = cell_properties_filtered_df)
                    # Clustering
                    array = cell_properties_filtered_df[['frequency_scaled', 'width_scaled', 'rise_time_scaled', 'decay_time_scaled', 'amplitude_scaled']].to_numpy()
                    cell_properties_filtered_df = utils.k_means_clustering(cell_properties_df = cell_properties_filtered_df, array=array, n_clusters = n_clusters)

                    ax_pca, cell_properties_filtered_df = plot.biplot(cell_properties_df = cell_properties_filtered_df, palette = palette)

                    if 'marker' not in cell_properties_df.columns:

                        # PCA
                        plt.savefig(os.path.join(plot_dir, 'spike_properties_pca.pdf'), bbox_inches='tight')
                        plt.close()    
                            
                        # Plot heatmap
                        ax_heatmap = plot.cluster_heatmap(cell_properties_df = cell_properties_filtered_df, cbar=True)
                        plt.savefig(os.path.join(plot_dir, 'cluster_heatmap.pdf'), bbox_inches='tight')
                        plt.close()        
                        
                        valid_combinations = cell_properties_filtered_df[['stimulation', 'biological_replicate']].drop_duplicates().assign(combined=lambda d: d['stimulation'].astype(str) + "_" + d['biological_replicate'].astype(str))['combined'].values

                        cluster_percentages_df = (
                            cell_properties_filtered_df
                            .groupby(['stimulation', 'biological_replicate'])['cluster']
                            .value_counts(normalize=True)   
                            .mul(100)                       
                            .rename('percentage')
                            .reset_index()
                            .assign(combined=lambda d: d['stimulation'].astype(str) + "_" + d['biological_replicate'].astype(str))
                        )
                        cluster_percentages_df = cluster_percentages_df[cluster_percentages_df.combined.isin(valid_combinations)]
                        
                        ax_barplot = plot.all_conditions_barplot(dataframe = cluster_percentages_df, palette = palette, ycolumn='cluster', xcolumn='percentage', hue='stimulation')
                        ax_barplot.set_ylabel(r"% of cells", fontsize='xx-large')
                        plt.xlabel('')
                        plt.title('')
                        plt.savefig(os.path.join(plot_dir, 'cluster_barplot.pdf'), bbox_inches='tight')
                        plt.close()       

                        cluster_percentages_df.to_csv(os.path.join(table_dir, 'cluster_percentages.csv'), sep=';') 
                        cluster_percentages_df.to_pickle(os.path.join(table_dir, 'cluster_percentages.pkl')) 

                    # Plot properties on beeswarm and on pca plot
                    for property in ['cluster','amplitude', 'width', 'rise_time', 'decay_time', 'frequency']:
                        if property != 'cluster':
                            if 'marker' in cell_properties_df.columns:   
                                ax_beeswarm = plot.beeswarm(cell_properties_df = cell_properties_filtered_df, control_condition = control_condition, palette=palette, y=property, x='marker', hue='stimulation',hue_separation=0.3, separation_between_plots=2, max_plot_width=0.3)
                            else:   
                                ax_beeswarm = plot.beeswarm(cell_properties_df = cell_properties_filtered_df, control_condition = control_condition, palette=palette, y=property, x='stimulation')
                            plt.ylabel(ylabel_dict[property])
                            plt.title(title_dict[property])
                            plt.savefig(os.path.join(plot_dir, f'beeswarm_{property}_pooled.pdf'),  bbox_inches='tight')
                            plt.close()    

                        if 'marker' not in cell_properties_df.columns:   

                            ax_pca_property = plot.pca_property(cell_properties_df = cell_properties_filtered_df, color_by = property)
                            plt.title(title_dict[property])
                            plt.savefig(os.path.join(plot_dir, f'pca_{property}.pdf'),  bbox_inches='tight')
                            plt.close()

                    # Save results
                    cell_properties_df.to_pickle(os.path.join(table_dir, 'cell_properties_spontaneous.pkl'))
                    cell_properties_filtered_df.to_pickle(os.path.join(table_dir, 'cell_properties_spontaneous_active_cells.pkl'))

                    cell_properties_df.to_csv(os.path.join(table_dir, 'cell_properties_spontaneous.csv'), sep=';')
                    cell_properties_filtered_df.to_csv(os.path.join(table_dir, 'cell_properties_spontaneous_active_cells.csv'), sep=';')

                    

                #if analysis_mode == 'Baseline change':
                if analysis_mode == 'Compound-evoked activity':
                    if activity_type == 'Baseline change':
                        print('Running baseline change analysis...')
                        parameter_list.extend([f'Quantification - Analysis window start: {analysis_window_start}', f'Quantification - Analysis window end: {analysis_window_end}', f'Quantification - Standard deviation threshold: {baseline_std_threshold}' , f'Quantification - Imaging interval: {imaging_interval}', f'Quantification - KCl stimulation frame: {kcl_frame}'])

                        cell_properties_df = utils.compute_auc(cell_properties_df = cell_properties_df,  start_frame = analysis_window_start, end_frame = analysis_window_end, column='AUC')


                        # Activity
                        cell_properties_df = activity.baseline_change(cell_properties_df=cell_properties_df, control_condition=control_condition, std_threshold=baseline_std_threshold)
                        
                        response_perc_well_df, response_perc_rep_df = utils.percentage_responding(cell_properties_df = cell_properties_df,analysis_type='baseline')

                        #response_perc_well_df, response_perc_rep_df = utils.percentage_responding_baseline(cell_properties_df = cell_properties_df)

                        print('Plotting and saving output...')
                        
                        response_perc_well_df.to_csv(os.path.join(table_dir, 'percentage_responding_per_well.csv'), sep=';')
                        response_perc_well_df.to_pickle(os.path.join(table_dir, 'percentage_responding_per_well.pkl'))

                        response_perc_rep_df.to_csv(os.path.join(table_dir, 'percentage_responding_per_replicate.csv'), sep=';')
                        response_perc_rep_df.to_pickle(os.path.join(table_dir, 'percentage_responding_per_replicate.pkl'))

                        # Clustering
                        array = np.array(cell_properties_df.dff.to_list())[:, :(kcl_frame if kcl_frame >= 0 else None)]
                        cell_properties_df = utils.k_means_clustering(cell_properties_df=cell_properties_df,array=array, n_clusters=n_clusters)

                        # Cluster centroids
                        cluster_dict, colors_dict = utils.cluster_centroids(cell_properties_df=cell_properties_df)
                        # initialize list of lists 

                        # Create the pandas DataFrame 
                        centroid_df = pd.DataFrame({
                            "cluster": list(cluster_dict.keys()),
                            "centroid": [arr.tolist() for arr in cluster_dict.values()]
                        })
                        centroid_df.to_csv(os.path.join(table_dir, 'cluster_centroids.csv'), sep=';')
                        centroid_df.to_pickle(os.path.join(table_dir, 'cluster_centroids.pkl'))
                        
                        if 'marker' not in cell_properties_df.columns:

                            # Plot overview of the data
                            ax_barplot = plot.all_conditions_barplot(dataframe=response_perc_well_df, palette=palette, ycolumn="proportion_positive_cells", xcolumn='biological_replicate', hue = 'stimulation')
                            ax_barplot.set_ylabel(r"% responding cells (Ca$^{2+}$)", fontsize='xx-large')
                            plt.xlabel('')
                            plt.title('')
                            plt.savefig(os.path.join(plot_dir, 'barplot_responding_cells_per_well.pdf'), bbox_inches='tight')
                            plt.close()

                            ax_barplot = plot.all_conditions_barplot(dataframe=response_perc_rep_df, palette=palette, ycolumn="proportion_positive_cells", xcolumn='stimulation')
                            ax_barplot.set_ylabel(r"% responding cells (Ca$^{2+}$)", fontsize='xx-large')
                            plt.xlabel('')
                            plt.title('')
                            plt.savefig(os.path.join(plot_dir, 'barplot_responding_cells_per_replicate.pdf'), bbox_inches='tight')
                            plt.close()

                            # Plot cluster centroids
                            ax_cluster_lineplot = plot.cluster_centroids(cluster_dict=cluster_dict,palette=colors_dict, imaging_interval = imaging_interval)
                            plt.savefig(os.path.join(plot_dir, 'cluster_centroids.pdf'), bbox_inches='tight')
                            plt.close()

                            # Plot cluster-to-treatment barplot
                            cluster_df = cell_properties_df.groupby(['stimulation','biological_replicate'])['cluster'].value_counts(normalize=True).mul(100).reset_index()
                            ax_cluster_barplot = plot.all_conditions_barplot(dataframe=cluster_df, palette=colors_dict, ycolumn='proportion',xcolumn='stimulation',hue='cluster')
                            plt.savefig(os.path.join(plot_dir, 'stimulation_to_cluster.pdf'), bbox_inches='tight')
                            plt.close()
                            cluster_df['cluster_color'] = cluster_df['cluster'].astype(int).map(colors_dict)
                            cluster_df.to_csv(os.path.join(table_dir, 'stimulation_to_cluster.csv'), sep=';')
                            cluster_df.to_pickle(os.path.join(table_dir, 'stimulation_to_cluster.pkl'))

                        
                        # Valid diffs
                        for treatment_condition in cell_properties_df.stimulation.cat.categories.tolist():
                            if treatment_condition != control_condition:
                                valid_comparisons = cell_properties_df[cell_properties_df.stimulation.isin([treatment_condition, control_condition])].groupby(['plate_id_biological_replicate']).stimulation.nunique()
                                valid_comparisons = valid_comparisons[valid_comparisons == 2].index  
                                comparison_df = cell_properties_df[(cell_properties_df['plate_id_biological_replicate'].isin(valid_comparisons)) & (cell_properties_df.stimulation.isin([treatment_condition, control_condition]))]
                                comparison_df["stimulation"] = comparison_df["stimulation"].cat.remove_unused_categories()
                                # Generate output folder
                                if 'marker' in comparison_df.columns:
                                    comparison_table_dir = os.path.join(project_dir, 'Quantification',f'{control_condition}_vs_{treatment_condition}','Tables')
                                    comparison_plot_dir = os.path.join(project_dir, 'Quantification',f'{control_condition}_vs_{treatment_condition}','Plots')

                                else:   
                                    comparison_table_dir = os.path.join(project_dir, 'Quantification',f'{control_condition}_vs_{treatment_condition}','Tables')
                                    comparison_plot_dir = os.path.join(project_dir, 'Quantification',f'{control_condition}_vs_{treatment_condition}','Plots')

                                if not os.path.exists(comparison_table_dir): os.makedirs(comparison_table_dir)
                                if not os.path.exists(comparison_plot_dir): os.makedirs(comparison_plot_dir)

                                if 'marker' not in comparison_df.columns:

                                    for exp_replicate in comparison_df.plate_id_biological_replicate.unique():
                                        subset_df = comparison_df[comparison_df.plate_id_biological_replicate == exp_replicate].copy()
                                        subset_df["stimulation"] = subset_df["stimulation"].cat.remove_unused_categories()

                                        # Beeswarm
                                        with plt.rc_context({"figure.dpi": 350, "figure.figsize": (1.8, 2.2)}):

                                            ax_beeswarm = plot.beeswarm(cell_properties_df=subset_df, control_condition=control_condition,  std_threshold=2,  palette = palette, brace=True, x='stimulation', y='AUC', control_condition_mean=True)
                                            plt.ylabel(r"Ca$^{2+}$ response (AUC)")
                                            plt.savefig(os.path.join(comparison_plot_dir, f'beeswarm_{exp_replicate}.pdf'),  bbox_inches='tight')
                                            plt.close()              

                                        # Heatmap
                                        ax_heatmap = plot.heatmap(cell_properties_df=subset_df,   imaging_interval=imaging_interval,  cmap = 'plasma', palette = palette)
                                        plt.savefig(os.path.join(comparison_plot_dir, f'heatmap_{exp_replicate}.pdf'),  bbox_inches='tight')
                                        plt.close()              
                                        plt.close()              

                                        # Trace
                                        ax_trace = plot.overlaid_traces_two_groups(cell_properties_df = subset_df,  control_condition = control_condition, treatment_condition= treatment_condition, trace='dff', mean=True, start_frame= analysis_window_start, end_frame = analysis_window_end, stimulation_frame = stimulation_frame,  palette=palette, imaging_interval=imaging_interval)
                                        plt.savefig(os.path.join(comparison_plot_dir, f'trace_{exp_replicate}.pdf'),  bbox_inches='tight')
                                        plt.close()              

                                # barplot
                                filtered_response_perc_df = response_perc_well_df[response_perc_well_df.image_id.isin(comparison_df.image_id.unique())].copy()
                                filtered_response_perc_df["stimulation"] = filtered_response_perc_df["stimulation"].cat.remove_unused_categories()
                                
                                if 'marker' in filtered_response_perc_df.columns:
                                    response_perc_mean_df = filtered_response_perc_df.groupby(['stimulation', 'biological_replicate', 'marker'])['proportion_positive_cells'].mean().reset_index()
                                else:
                                    response_perc_mean_df = filtered_response_perc_df.groupby(['stimulation', 'biological_replicate'])['proportion_positive_cells'].mean().reset_index()


                                if 'marker' in response_perc_mean_df.columns:
                                    ax_two_conditions_barplot = plot.two_conditions_barplot(response_perc_mean_df = response_perc_mean_df, palette=palette, x = 'marker', y = 'proportion_positive_cells', hue = 'stimulation')
                                    plt.ylabel(r"% responding (Ca$^{2+}$)", fontsize='large')
                                    #plt.title(marker, fontsize='large')
                                    plt.savefig(os.path.join(comparison_plot_dir, f'barplot_responding_cells_per_replicate.pdf'),  bbox_inches='tight')
                                    plt.close()              

                                    ax_swarmplot = plot.beeswarm(comparison_df, x="marker", y='AUC',hue='stimulation',  palette=palette, ax=None, separation_between_plots=1.8, max_plot_width=0.5)
                                    plt.savefig(os.path.join(comparison_plot_dir, 'beeswarm_plot_auc.pdf'),  bbox_inches='tight')
                                    plt.close()              

                                    

                                else:
                                    with plt.rc_context({"figure.dpi": 350, "figure.figsize": (1.8, 2.2)}):
                                        ax_two_conditions_barplot = plot.two_conditions_barplot(response_perc_mean_df = response_perc_mean_df, palette=palette, x='stimulation', y='proportion_positive_cells')
                                        plt.ylabel(r"% responding (Ca$^{2+}$)", fontsize='large')
                                        plt.savefig(os.path.join(comparison_plot_dir, 'barplot_responding_cells_per_replicate.pdf'),  bbox_inches='tight')
                                        plt.close()              

                                
                                comparison_df.to_csv(os.path.join(comparison_table_dir,f'{control_condition}_vs_{treatment_condition}_cell_properties_baseline_change.csv'), sep=';')
                                comparison_df.to_pickle(os.path.join(comparison_table_dir,f'{control_condition}_vs_{treatment_condition}_cell_properties_baseline_change.pkl'))

                                filtered_response_perc_df.to_csv(os.path.join(comparison_table_dir,f'{control_condition}_vs_{treatment_condition}_percentage_responding_per_well.csv'), sep=';')
                                filtered_response_perc_df.to_pickle(os.path.join(comparison_table_dir,f'{control_condition}_vs_{treatment_condition}_percentage_responding_per_well.pkl'))

                                response_perc_mean_df.to_csv(os.path.join(comparison_table_dir,f'{control_condition}_vs_{treatment_condition}_percentage_responding_per_replicate.csv'), sep=';')
                                response_perc_mean_df.to_pickle(os.path.join(comparison_table_dir,f'{control_condition}_vs_{treatment_condition}_percentage_responding_per_replicate.pkl'))

                                
                        
                        cell_properties_df.to_csv(os.path.join(table_dir,f'cell_properties_baseline_change.csv'), sep=';')
                        cell_properties_df.to_pickle(os.path.join(table_dir,f'cell_properties_baseline_change.pkl'))


                    elif activity_type == 'Spontaneous':
                        parameter_list.extend([f'Quantification - Prominence threshold: {spike_prominence_threshold}', f'Quantification - Amplitude width ratio: {spike_amplitude_width_ratio}',   f'Quantification - Imaging interval: {imaging_interval}', f'Quantification - KCl stimulation frame: {kcl_frame}'])

                        cell_properties_df, start_frame, end_frame = activity.spike(cell_properties_df = cell_properties_df, prominence = spike_prominence_threshold,  amplitude_width_ratio=spike_amplitude_width_ratio, imaging_interval=imaging_interval, start_frame=analysis_window_start, end_frame = analysis_window_end)
                        spontaneous_activity_output()


                elif analysis_mode == 'Spontaneous activity':
                    print('Running spontaneous activity analysis...')
                    parameter_list.extend([f'Quantification - Prominence threshold: {spike_prominence_threshold}', f'Quantification - Amplitude width ratio: {spike_amplitude_width_ratio}',   f'Quantification - Imaging interval: {imaging_interval}', f'Quantification - KCl stimulation frame: {kcl_frame}'])

                    # Detect events
                    cell_properties_df, start_frame, end_frame = activity.spike(cell_properties_df = cell_properties_df, prominence = spike_prominence_threshold,  amplitude_width_ratio=spike_amplitude_width_ratio, imaging_interval=imaging_interval, end_frame = kcl_frame)
                    spontaneous_activity_output()
                    
                time = datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                end_time = perf_counter()
                runtime_seconds = end_time - start_time

                parameter_list.extend([f'Downstream - Number of clusters:{n_clusters}',f'',f'Date and time: {time}', f'Pipeline runtime (s): {runtime_seconds}' ])

                with open(os.path.join(project_dir,f'Quantification', 'normalization_signal_quantification_config.txt'), 'w') as f:
                    for line in parameter_list: f.write(f'{line}\n')

                print('Analysis ready...')

                

        except Exception as e: 
            print("\nAn error occurred:", e)
            traceback.print_exc()

        widget.call_button.enabled = True
        widget.project_dir.value = ''
        widget.analysis_mode.value = '-- Select --'
        widget.activity_type.value = '-- Select --'
        widget.control_condition.value = ''
        widget.normalization_mode.value = '-- Select --'
        reset_settings_sc_analysis(widget)
        # project_dir analysis_mode activity_type control_condition normalization_mode



    # Hiding widgets
    widget.stimulation_frame.visible = False
    widget.kcl_frame.visible = False
    widget.sliding_window_size.visible = False
    widget.spike_label.visible = False
    widget.percentile_threshold.visible = False
    widget.spike_prominence_threshold.visible = False
    widget.spike_amplitude_width_ratio.visible = False
    widget.baseline_std_threshold.visible = False
    widget.analysis_window_start.visible = False
    widget.analysis_window_end.visible = False
    widget.imaging_interval.visible = False
    widget.activity_type.visible = False
    widget.downstream_label.visible = False
    widget.n_clusters.visible = False


    @widget.project_dir.changed.connect
    def _project_dir_changed():

        widget.analysis_mode.value = '-- Select --'
        disable_enable_value(widget.analysis_mode, 'disable', '-- Select --')


    @widget.analysis_mode.changed.connect
    def _analysis_mode_changed():
        widget.spike_label.visible = True
        widget.normalization_mode.value = '-- Select --'
        widget.activity_type.value = '-- Select --'
        reset_settings_sc_analysis(widget)
        disable_enable_value(widget.analysis_mode, 'disable', '-- Select --')
        disable_enable_value(widget.normalization_mode, 'enable', '-- Select --')
        disable_enable_value(widget.normalization_mode, 'enable','Sliding window')
        disable_enable_value(widget.normalization_mode, 'enable', 'Pre-stimulus window')
        disable_enable_value(widget.activity_type, 'enable', '-- Select --')


        widget.sliding_window_size.visible = False
        widget.percentile_threshold.visible = False
        widget.stimulation_frame.visible = False


        if widget.analysis_mode.value == 'Compound-evoked activity':
            widget.activity_type.visible = True

        else:   
            disable_enable_value(widget.normalization_mode, 'enable', 'Sliding window')
            disable_enable_value(widget.normalization_mode, 'disable', 'Pre-stimulus window')
            widget.activity_type.visible = False
            widget.spike_prominence_threshold.visible = True
            widget.spike_amplitude_width_ratio.visible = True
            widget.imaging_interval.visible = True
            widget.kcl_frame.visible = True
            widget.baseline_std_threshold.visible = False
            widget.analysis_window_start.visible = False
            widget.analysis_window_end.visible = False
            widget.downstream_label.visible = True
            widget.n_clusters.visible = True

    @widget.activity_type.changed.connect
    def _activity_mode_connect():
        disable_enable_value(widget.activity_type, 'disable', '-- Select --')
        widget.normalization_mode.value = '-- Select --'
        reset_settings_sc_analysis(widget)
        disable_enable_value(widget.normalization_mode, 'enable', '-- Select --')
        widget.sliding_window_size.visible = False
        widget.percentile_threshold.visible = False
        widget.stimulation_frame.visible = False

        if widget.activity_type.value == 'Baseline change':
            disable_enable_value(widget.normalization_mode, 'disable', 'Sliding window')
            disable_enable_value(widget.normalization_mode, 'enable', 'Pre-stimulus window')
            widget.baseline_std_threshold.visible = True
            widget.analysis_window_start.visible = True
            widget.analysis_window_end.visible = True
            widget.imaging_interval.visible = True
            widget.kcl_frame.visible = True
            widget.spike_prominence_threshold.visible = False
            widget.spike_amplitude_width_ratio.visible = False
            widget.downstream_label.visible = True
            widget.n_clusters.visible = True



        elif widget.activity_type.value == 'Spontaneous':
            disable_enable_value(widget.normalization_mode, 'enable', 'Sliding window')
            disable_enable_value(widget.normalization_mode, 'enable', 'Pre-stimulus window')
            widget.baseline_std_threshold.visible = False
            widget.analysis_window_start.visible = True
            widget.analysis_window_end.visible = True
            widget.imaging_interval.visible = True
            widget.kcl_frame.visible = True
            widget.spike_prominence_threshold.visible = True
            widget.spike_amplitude_width_ratio.visible = True
            widget.downstream_label.visible = True
            widget.n_clusters.visible = True


    @widget.normalization_mode.changed.connect
    def _normalization_mode_changed():

        disable_enable_value(widget.normalization_mode, 'disable', '-- Select --')
        
        if widget.normalization_mode.value == 'Sliding window':
            widget.sliding_window_size.visible = True
            widget.percentile_threshold.visible = True
            widget.stimulation_frame.visible = False

        if widget.normalization_mode.value == 'Pre-stimulus window':
            widget.sliding_window_size.visible = False
            widget.percentile_threshold.visible = False
            widget.stimulation_frame.visible = True

        reset_settings_sc_analysis(widget)


    @widget.optimize_button.clicked.connect
    def _optimize_quantification():

        remove_layers(viewer)

        try:
            project_dir = widget.project_dir.value
            analysis_mode = widget.analysis_mode.value
            normalization_mode = widget.normalization_mode.value    
            control_condition =  widget.control_condition.value
            #analysis_setup = widget.analysis_setup.value
            stimulation_frame = widget.stimulation_frame.value    
            kcl_frame = widget.kcl_frame.value
            sliding_window_size = widget.sliding_window_size.value
            percentile_threshold = widget.percentile_threshold.value
            spike_prominence_threshold = widget.spike_prominence_threshold.value
            spike_amplitude_width_ratio = widget.spike_amplitude_width_ratio.value
            baseline_std_threshold = widget.baseline_std_threshold.value
            analysis_window_start = widget.analysis_window_start.value  
            analysis_window_end = widget.analysis_window_end.value  
            activity_type = widget.activity_type.value
            imaging_interval = widget.imaging_interval.value

            widget_state.cell_properties_df = get_cell_properties_df(project_dir)
            cell_properties_df = widget_state.cell_properties_df.copy()


            if control_condition not in cell_properties_df.stimulation.tolist():
                raise ValueError(f"Condition '{control_condition}' not in input data")
            
            if not imaging_interval > 0:
                raise ValueError(f"Imaging interval needs to be above 0")


            print('Testing quantification settings')

            if (analysis_mode == 'Spontaneous activity') or (analysis_mode == 'Compound-evoked activity' and activity_type == 'Spontaneous'):
                print('RANDOM1')
                samples = cell_properties_df['image_id'].unique().tolist()
                sample = random.choice(samples)
                temp_df = cell_properties_df[cell_properties_df['image_id'] == sample].copy()
                print(f"Sampled {cell_properties_df[cell_properties_df['image_id'] == sample].filename.values[0]} from {cell_properties_df[cell_properties_df['image_id'] == sample].plate_id.values[0]}")

            elif analysis_mode == 'Compound-evoked activity' and activity_type == 'Baseline change':
                print('RANDOM2')

                temp_df = cell_properties_df.copy()
                temp_df['plate_id_biological_replicate'] = temp_df['plate_id'].astype(str) + '_' + temp_df['biological_replicate'].astype(str)
                samples = temp_df['plate_id_biological_replicate'].unique().tolist()
                sample = random.choice(samples)
                temp_df = temp_df[temp_df['plate_id_biological_replicate'] == sample].copy()
                treatment_condition = random.choice([c for c in set(temp_df.stimulation.tolist()) if c != control_condition])
                temp_df = temp_df[temp_df['stimulation'].isin([treatment_condition, control_condition])].copy()

                # Sort the df for plotting
                stimulations = sorted(temp_df['stimulation'].unique())
                stimulations = [control_condition] + [c for c in stimulations if c != control_condition]
                temp_df['stimulation'] = pd.Categorical(temp_df['stimulation'], categories=stimulations, ordered=True)
                temp_df = temp_df.sort_values('stimulation')


            if normalization_mode == 'Sliding window':
                print('Sliding')

                temp_df = preprocess.sliding_window(cell_properties_df = temp_df, sliding_window_size = sliding_window_size, percentile_threshold = percentile_threshold)

            elif normalization_mode == 'Pre-stimulus window':
                print('Pre-stimulus')
                
                temp_df = preprocess.pre_stimulation(cell_properties_df = temp_df, stimulation_frame = stimulation_frame)
                
                #plot_list = plot.cellwise_traces(cell_properties_df = temp_df, trace='raw',baseline=False)
            print(temp_df.columns)
            palette = get_colors(cell_properties_df, project_dir)


            #remove_widget_if_exists(viewer, 'Baseline')
            remove_widget_if_exists(viewer, 'Baseline normalization')
            remove_widget_if_exists(viewer, 'Spontaneous activity')
            remove_widget_if_exists(viewer, 'Baseline quantification') 
            remove_widget_if_exists(viewer, 'Trace quantification') 


            if analysis_mode == 'Compound-evoked activity' and activity_type == 'Baseline change':
                temp_df = utils.compute_auc(cell_properties_df = temp_df,  start_frame = analysis_window_start, end_frame = analysis_window_end, column='AUC')
                temp_df = activity.baseline_change(cell_properties_df=temp_df, control_condition = control_condition, std_threshold = baseline_std_threshold)

                # Plot results

                fig, ax = plt.subplots(2, 1, figsize=(6, 10))
                if kcl_frame > 0:
                    ax[0] = plot.overlaid_traces_two_groups(cell_properties_df = temp_df,  control_condition = control_condition, treatment_condition= treatment_condition, trace='dff', mean=True, start_frame= analysis_window_start, end_frame = analysis_window_end, stimulation_frame = stimulation_frame, ax=ax[0], palette=palette, imaging_interval=1, kcl_frame=kcl_frame)
                else:
                    ax[0] = plot.overlaid_traces_two_groups(cell_properties_df = temp_df,  control_condition = control_condition, treatment_condition= treatment_condition, trace='dff', mean=True, start_frame= analysis_window_start, end_frame = analysis_window_end, stimulation_frame = stimulation_frame, ax=ax[0], palette=palette, imaging_interval=1)
                ax[0].set_xlabel('')  # Adjust labels  
                

                ax = plot.beeswarm(cell_properties_df = temp_df, control_condition = control_condition,  std_threshold=baseline_std_threshold, ax=ax[1], palette=palette, y='AUC',x='stimulation', control_condition_mean=True, brace=True,)
                plt.ylabel(r"Ca$^{2+}$ response (AUC)")
                plt.tight_layout()
            
                widget_baseline = PlotViewerBaseline(fig)
                widget_baseline.setObjectName('Baseline quantification')
                dw_baseline = viewer.window.add_dock_widget(widget_baseline, area='right', name='Baseline quantification')

                # Tabify dw1 and dw2 with trace_quant_widget
                dw_main = viewer.window.add_dock_widget(widget, area='right', name='Trace quantification')
                viewer.window._qt_window.tabifyDockWidget(dw_main, dw_baseline)


            elif (analysis_mode == 'Spontaneous activity') or (analysis_mode == 'Compound-evoked activity' and activity_type == 'Spontaneous'):
                if kcl_frame < 0: kcl_frame = None
                if analysis_mode == 'Spontaneous activity':
                    cell_properties_df, start_frame, end_frame = activity.spike(cell_properties_df = temp_df, prominence = spike_prominence_threshold,  amplitude_width_ratio=spike_amplitude_width_ratio, imaging_interval=imaging_interval, end_frame = kcl_frame)
                elif analysis_mode == 'Compound-evoked activity':
                    cell_properties_df, start_frame, end_frame = activity.spike(cell_properties_df = temp_df, prominence = spike_prominence_threshold,  amplitude_width_ratio=spike_amplitude_width_ratio, imaging_interval=imaging_interval, start_frame=analysis_window_start, end_frame = analysis_window_end)


                plot_list = plot.cellwise_traces(cell_properties_df = temp_df, trace='raw',baseline=True)
                plot_list_spikes = plot.cellwise_traces(cell_properties_df = temp_df, trace='dff',baseline=False, spikes=True, spikes_mode='all')


                stack_events, mask_events = plot.overlay_events(cell_properties_df = temp_df)
                mask = io.imread(temp_df.mask_path.values[0])

                widget_normalization = PlotViewer(plot_list)
                widget_normalization.setObjectName('Baseline normalization')
                dw1 = viewer.window.add_dock_widget(widget_normalization, area='right', name='Baseline normalization')

                widget_activity = PlotViewer(plot_list_spikes)
                widget_activity.setObjectName('Spontaneous activity')
                dw2 = viewer.window.add_dock_widget(widget_activity, area='right', name='Spontaneous activity')

                # Tabify dw1 and dw2 with trace_quant_widget
                dw_main = viewer.window.add_dock_widget(widget, area='right', name='Trace quantification')
                viewer.window._qt_window.tabifyDockWidget(dw_main, dw1)
                viewer.window._qt_window.tabifyDockWidget(dw_main, dw2)
                viewer.add_image(stack_events, name='CA video')

                viewer.add_labels(mask, name='Labels', colormap={label: [0.0, 1.0, 1.0, 1.0] for label in np.unique(mask) if label != 0})
                viewer.layers['Labels'].contour = 3
                viewer.layers['Labels'].visible = False

                points = list(zip(temp_df['centroid-0'], temp_df['centroid-1']))

                text = {
                    'string': '{label_id}',
                    'size': 7,
                    'color': '#FF00FF',
                    'anchor': 'center',
                    'translation': [0, 0],
                    'properties': {'label_id': temp_df['label'].tolist()}
                }

                viewer.add_points(points, name='Label id', size=10, face_color='transparent', border_color='transparent', text=text)
                viewer.add_labels(mask_events, name='Activity', colormap={100:[0.9882, 0.6275, 0.1686, 1.0]})
                

                viewer.dims.current_step = (1,)

            print(f'Done, try different threshold and test quantification with another image...\n')



        except Exception as e: 
            print("\nAn error occurred:", e)
            traceback.print_exc()

    return widget


def napari_experimental_provide_dock_widget():
    return segmentation_widget, single_cell_widget, {"name": "My Pipeline Launcher"}



