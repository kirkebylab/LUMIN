import warnings
warnings.filterwarnings("ignore")
import  sys
from skimage import util, segmentation 
import numpy as np
import napari
from napari_blob_detection import points_to_labels
from scipy import ndimage as ndi
from csbdeep.utils import normalize



# from napari_scripts.utils import read_and_project_image

# Notes: check that each diff has a control and stimulated sample
# concat image name, condition, diff for unique id

# parse_input function takes CA videos input folder and returns dataframe containing image, sample and condition information
#  returns annotated_image_df from previous run if output folder exists, alternatively returns empty annotated_image_df





# set the values to defaults ?????
def run_stardist(image = None,  model_sd = None,  prob_thresh_sd = None, overlap_thresh_sd = None):
    mask, _ = model_sd.predict_instances(normalize(image,  1,99.8,axis=(0,1)), nms_thresh=overlap_thresh_sd, prob_thresh=prob_thresh_sd)

    del model_sd

    return mask


def run_cellpose(image = None, model_cp = None, diameter_cp = None, cellprob_threshold_cp = None, flow_threshold_cp = None):

    segmentation_params = {
        'diameter': int(diameter_cp),
        'cellprob_threshold': cellprob_threshold_cp,
        'flow_threshold': flow_threshold_cp,
        'resample': False,  
        'do_3D': False,  
        'stitch_threshold': 0.0,  
    }

    mask, _, _ = model_cp.eval(image,**segmentation_params)

    del model_cp

    return mask

def run_manual_selection(image: np.ndarray,  point_size: int, first_label: int = 0):

    label_text = {"string": "label",'size': 6,'color': '#eb1717'}
    roi_color = '#66abd9'
    all_labels = set() # reset all_labels variable
    
    # Initialize napari viewer with the image
    viewer, image_layer = napari.imshow(image)
    viewer.window._qt_viewer.layerButtons.hide() # Hide unnecessary buttons
    #state = viewer.camera.get_state()
    #viewer.window.qt_viewer.view.camera.set_state(state)

    features = {'label': np.empty(0, dtype=int)}
    points_layer = viewer.add_points(features = features,name='points', size=point_size, face_color=roi_color, symbol='disc',text=label_text) # Adding points layer to viewer

    @points_layer.events.data.connect
    def update_points_layer():
        current_labels = {label for label in points_layer.properties.get('label', [])}
        # Update the global set with the current labels
        all_labels.update(current_labels)
        # Assign the next label ID
        next_label = max(all_labels, default=first_label) + 1  # Use default=0 for an empty set
        points_layer.feature_defaults['label'] = next_label

    # fire function set_point_labels once at the start
    update_points_layer()

    # Start the event loop and show the viewer
    viewer.show(block=True)
    
    if len(vars(points_layer)['_data']) == 0: # Check if user returned image without labels
        sys.exit("❌ Stopping program: User saved image without masking.")

    # Set the applied contrast to image for visualization purposes
    layer = viewer.layers['image']
    m, M = layer.contrast_limits
    rescaled = (layer.data - m) / (M - m)
    image = np.clip(rescaled, 0, 1)

    # Create mask from points layer
    mask,_, layer_type = points_to_labels(points_layer,image_layer) # Convert napari points layer to labels

    assert layer_type == 'Labels'

    points_layer_data = points_layer.data
    points_layer_label = points_layer.properties["label"]

    # Make sure indecies with negative coordinates and with coordinates not exceeding the image size are removed
    indices_to_remove = np.where((points_layer_data[:, 0] >= image_layer.data.shape[0]) | (points_layer_data[:, 1] >= image_layer.data.shape[1]) | (points_layer_data[:, 0] < 0) | (points_layer_data[:, 1] < 0))[0]

    points_layer_data = np.delete(points_layer_data, indices_to_remove, axis=0)
    points_layer_label = np.delete(points_layer_label, indices_to_remove)

    # Convert the binary mask to labelled mask using the centroid coordinates
    seeds = util.label_points(points_layer_data, image_layer.data.shape)
    distances = -ndi.distance_transform_edt(mask)
    mask = segmentation.watershed(distances, seeds, mask=mask)

    features = {'label': np.array(points_layer_label, dtype=int)}

    # update labelled mask label values to correspond the ID list
    updated_values = {0: 0}  
    updated_values.update({i+1: value for i, value in enumerate(features['label'])})
    mask = np.vectorize(updated_values.get, otypes=[mask.dtype])(mask)

    return image, mask


    
    
        











'''
            

def select_roi_extract_traces(filepath,filename, image_condition, image_replicate,image_id, first_frame, param_point_size, param_output_folder, param_annotated_images_df):

    label_text = {"string": "label",'size': 6,'color': '#eb1717'}
    roi_color = '#66abd9'
    all_labels = set() # reset all_labels variable


    image_mean, img_stack = read_and_project_image(filepath, first_frame)
    
    # Initialize napari viewer with the image
    viewer, image_layer = napari.imshow(image_mean, title=f'{filename}_{image_condition}_{image_replicate}')
    viewer.window._qt_viewer.layerButtons.hide() # Hide unnecessary buttons
    state = viewer.window.qt_viewer.view.camera.get_state()
    viewer.window.qt_viewer.view.camera.set_state(state)

    features = {'label': np.empty(0, dtype=int)}
    points_layer = viewer.add_points(features = features,name='points', size=param_point_size, face_color=roi_color, symbol='disc',edge_width=0,text=label_text) # Adding points layer to viewer

    @points_layer.events.data.connect
    def update_points_layer():
        current_labels = {label for label in points_layer.properties.get('label', [])}
        # Update the global set with the current labels
        all_labels.update(current_labels)
        # Assign the next label ID
        next_label = max(all_labels, default=0) + 1  # Use default=0 for an empty set
        points_layer.feature_defaults['label'] = next_label

    # fire function set_point_labels once at the start
    update_points_layer()

    # Start the event loop and show the viewer
    
    
    if len(vars(points_layer)['_data']) == 0: # Check if user returned image without labels
        sys.exit("❌ Stopping program: User saved image without masking.")
        #sys.exit(0)
        #print("❌ [ANNOTATE ERROR] User saved image without masking. Skipping downstream.", file=sys.stderr)
        #sys.exit(1)


    layer = viewer.layers['image_mean']
    m, M = layer.contrast_limits
    rescaled = (layer.data - m) / (M - m)
    image_mean = np.clip(rescaled, 0, 1)

    def create_mask(points_layer, image_layer, img_stack):
        img_mask,_, layer_type = points_to_labels(points_layer,image_layer) # Convert napari points layer to labels

        assert layer_type == 'Labels'

        points_layer_data = points_layer.data
        points_layer_label = points_layer.properties["label"]

        # Make sure indecies with negative coordinates and with coordinates not exceeding the image size are removed
        indices_to_remove = np.where((points_layer_data[:, 0] >= image_layer.data.shape[0]) | (points_layer_data[:, 1] >= image_layer.data.shape[1]) | (points_layer_data[:, 0] < 0) | (points_layer_data[:, 1] < 0))[0]

        points_layer_data = np.delete(points_layer_data, indices_to_remove, axis=0)
        points_layer_label = np.delete(points_layer_label, indices_to_remove)

        # Convert the binary mask to labelled mask using the centroid coordinates
        seeds = skimage.util.label_points(points_layer_data, image_layer.data.shape)
        distances = -ndi.distance_transform_edt(img_mask)
        img_mask = skimage.segmentation.watershed(distances, seeds, mask=img_mask)

        features = {'label': np.array(points_layer_label, dtype=int)}

        # update labelled mask label values to correspond the ID list
        updated_values = {0: 0}  
        updated_values.update({i+1: value for i, value in enumerate(features['label'])})
        img_mask = np.vectorize(updated_values.get, otypes=[img_mask.dtype])(img_mask)
        
        # Get mean intensity for each label 
        properties = regionprops_table_all_frames(img_stack[:,np.newaxis,:,:], np.full_like(img_stack[:,np.newaxis,:,:], img_mask, dtype=int), intensity=True,size=True )
        properties_df = pd.DataFrame(properties)
        # Extract raw traces from mean intensity
        intensity_df = properties_df[['label','mean_intensity','area']]
        intensity_df = intensity_df.groupby('label').agg({'mean_intensity': list,'area': list}).reset_index()
        raw_traces = pd.DataFrame(intensity_df['mean_intensity'].tolist()).T.to_numpy()
        #intensity_df['centroid-0'], intensity_df['centroid-1'] = points_layer_data[:, 0].tolist(),points_layer_data[:, 1].tolist()
        #features = {'label': np.array(points_layer.properties["label"], dtype=int)}
        intensity_df['label'] = features['label'] 

        #label_centroids = points_layer.data.tolist()
        return points_layer_data, features, raw_traces, img_mask, properties_df

    label_centroids, features, raw_traces, image_mask, properties_df = create_mask(points_layer, image_layer, img_stack)

    def generate_output(filename, image_condition, image_replicate,image_id, image_mean, image_mask, raw_traces, properties_df, label_centroids,common_labels, filepath, param_output_folder, param_annotated_images_df): 
        # Generate output folder, overwrite if exists
        output_folder = os.path.join(param_output_folder, f'{image_condition}/{image_replicate}/{filename}')

        if os.path.exists(output_folder):
            shutil.rmtree(output_folder)
        os.makedirs(output_folder)
        os.makedirs(os.path.join(output_folder, 'Plots/Annotation/1_Raw_traces'))
        # os.mkdir(f"{output_folder}/5_dff_traces/")

        # plot_cell_traces(dff_traces, output_folder, filename,common_labels)
        plot_individual_traces(raw_traces, os.path.join(output_folder, 'Plots/Annotation/1_Raw_traces'), filename, common_labels)


        # Plot raw and dff-normalized calcium signal and mean calcium signal of all ROI
        df_raw = pd.DataFrame(raw_traces)

        # Write traces and label properties to file
        df_raw.to_csv(f"{output_folder}/raw_traces.csv", header=False, index=False)
        properties_df.to_csv(f'{output_folder}/label_properties.csv')

        f = plt.figure(figsize=(15,10))
        for col in df_raw.columns:
            plt.plot(df_raw[col],linewidth=1)
        plt.savefig(f"{output_folder}/Plots/Annotation/2_Raw_traces.pdf",dpi=300,bbox_inches="tight")  
          
        # Saving image and mask to file
        skimage.io.imsave(f"{output_folder}/image_projected.tiff", image_mean)
        skimage.io.imsave(f"{output_folder}/image_mask.tiff", image_mask.astype('uint16')) 

        def overlay_roi(image_mean, image_mask, label_centroids, labels):
            # Show the selected ROIs on top of the image
            image_mean_visualize = cv.normalize(image_mean, None, alpha=0, beta=255, norm_type=cv.NORM_MINMAX, dtype=cv.CV_8U)
            image_mean_visualize = cv.cvtColor(image_mean_visualize, cv.COLOR_GRAY2RGB)  
            
            
            img_mask_visualize = np.float32(image_mask)
            img_mask_visualize = cv.cvtColor(img_mask_visualize, cv.COLOR_GRAY2RGB)
            image_mean_visualize[np.all(img_mask_visualize != (0, 0, 0), axis=-1)] = tuple(int(roi_color.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))

            fig, ax = plt.subplots()
            ax.imshow(image_mean_visualize)
            kwargs={'horizontalalignment':'center','fontsize':2}
            for i in range(len(label_centroids)): # Numbering the labels
                ax.text(label_centroids[i][1]+0.55, label_centroids[i][0]+3,  str(labels[i]),color='#eb1717',**kwargs)
            return image_mean_visualize
    


        overlay_img = overlay_roi(image_mean,image_mask,label_centroids, common_labels)
        plt.imshow(overlay_img)
        plt.axis('off')
        #os.path.join(output_folder, 'Plots/Annotation/aw_traces')
        plt.savefig(os.path.join(output_folder, f'Plots/Annotation/3_Overlay_{filename}.pdf' ), dpi=400)


        # Save output to csv file
        param_annotated_images_df.loc[len(param_annotated_images_df.index)] = [filename,image_condition, image_replicate,image_id, filepath, f"{output_folder}/image_{filename}.tiff", f"{output_folder}/mask_{filename}.tiff", label_centroids,common_labels,f'{output_folder}/Plots/overlay_{filename}.pdf']
        param_annotated_images_df.to_csv(f'{param_output_folder}/annotated_images.csv')

        return param_annotated_images_df

    # def generate_output(filename, image_condition, image_replicate, image_mean, image_mask, raw_traces, properties_df, label_centroids,common_labels, filepath, param_output_folder): 
    param_annotated_images_df = generate_output(filename, image_condition, image_replicate,image_id, image_mean, image_mask, raw_traces, properties_df, label_centroids,features['label'].tolist(), filepath, param_output_folder, param_annotated_images_df)

    print(f'\n{filename}.tif annotated.')

    return param_annotated_images_df

    #return image_mean, label_centroids, features, raw_traces, image_mask, properties_df



def launch_annotation(param_image_df, param_annotated_image_df, nuclear_stain, param_output_folder, segmentation_mode, segmentation_parameters):

    if nuclear_stain == 'None':
        first_frame = 0
    else: 
        first_frame = 1

    for i in range(0,len(param_image_df)):
        annotated = param_annotated_image_df['image_id'].values

        image_id = param_image_df.iloc[i].image_id
        filename = param_image_df.iloc[i].filename
        image_condition = param_image_df.iloc[i].condition
        image_replicate = param_image_df.iloc[i].biological_replicate

        if segmentation_mode == 'Cellpose':
            run_cellpose(param_image_df.iloc[i].filepath, first_frame, param_output_folder, segmentation_parameters[0], segmentation_parameters[1], segmentation_parameters[2], segmentation_parameters[3])

        else:
            # Check that image not annotated yet
            if image_id not in annotated:
                param_annotated_image_df = select_roi_extract_traces(param_image_df.iloc[i].filepath,filename, image_condition, image_replicate,image_id, first_frame, segmentation_parameters[0], param_output_folder, param_annotated_image_df)


            else: # if image is annotated, skip it and move to the next one
                print(f'\nOmitting {image_id}: Image is already annotated.\n')

'''