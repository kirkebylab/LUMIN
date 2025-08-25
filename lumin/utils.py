import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as mcolors
import os, shutil
import pandas as pd
from skimage.measure import regionprops_table
from skimage import io
from collections import Counter
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import traceback
import matplotlib.pyplot as plt



# Max projection
# read_and_project_image function reads the image stack from path,  computes max projection across frames, and returns both the max projected image and original stack
def read_and_project_image(filepath: str = None, first_frame: int = 0):
    print(f'Reading image stack from {filepath}...')
    # Read and project image stack to its max intensity
    img_stack = io.imread(filepath)
    image_projected = np.max(img_stack[first_frame:], axis = (0)) 

    return image_projected, img_stack # return projected image and stack



def parse_input_output(input_file: str = None, project_dir: str = None, selection_mode: str = None, co_stain: bool = False):
    try:
        image_df = pd.read_csv(input_file,  index_col=None, sep=None)
        image_df['image_id'] = image_df['filename'].astype(str) + '_' + image_df['plate_id'].astype(str) 


        # Parse output
        if not os.path.isdir(project_dir) or (selection_mode == 'Automated' or selection_mode is None):
            if os.path.exists(project_dir): shutil.rmtree(project_dir, ignore_errors=True)
            os.makedirs(project_dir)

        # Check if output csv file in output folder. If yes, read the images already annotated
        if os.path.exists(f'{project_dir}/annotated_images.csv'):
            annotated_image_df = pd.read_csv(f'{project_dir}/annotated_images.csv', index_col=0, sep=None) 
            
        else: # Otherwise create output DF
            if os.path.exists(project_dir): 
                shutil.rmtree(project_dir, ignore_errors=True)
            if co_stain == True:
                annotated_image_df = pd.DataFrame(columns=['image_id','filename','biological_replicate', 'stimulation','plate_id', 'marker_name','max_label','image_stack','mask'])
            else: 
                annotated_image_df = pd.DataFrame(columns=['image_id','filename','biological_replicate', 'stimulation','plate_id', 'max_label','image_stack','mask'])

    except Exception as e: 
        print("\nAn error occurred:", e)
        traceback.print_exc()

    return image_df, annotated_image_df

# Create DF contining each ovelrlapping nuclei
# Overlapping fraction
# ID of cell
# Size of cell

def nuclei_cell_intersection(mask_nuclear: np.ndarray, df_nuclear: pd.DataFrame, mask_cyto: np.ndarray, df_cyto: pd.DataFrame):
    # Step 1: Get binary masks and overlap
    mask_overlap_binary = (mask_nuclear > 0) & (mask_cyto > 0)

    # Step 2: Get overlapping label pairs
    nuclear_labels = mask_nuclear[mask_overlap_binary]
    cell_labels = mask_cyto[mask_overlap_binary]
    label_pairs = np.stack([nuclear_labels, cell_labels], axis=1)
    # Step 3: Count occurrences of each (nuc_id, cell_id) pair = overlap area
    pair_counts = Counter(map(tuple, label_pairs))  # {(nuc_id, cell_id): area}
    # Step 4: Build lookup dictionaries for nuclear and cyto areas
    nuc_area_dict = dict(zip(df_nuclear['label'], df_nuclear['area']))
    cell_area_dict = dict(zip(df_cyto['label'], df_cyto['area']))
    # Step 5: Assemble final results
    records = []

    for (nuc_id, cell_id), overlap_area in pair_counts.items():
        records.append({
            'nuclear_id': nuc_id,
            'nuclear_area': nuc_area_dict[nuc_id],
            'label': cell_id,
            'cell_area': cell_area_dict[cell_id],
            'overlap_area': overlap_area,
            'overlap_fraction_nuclear': overlap_area / nuc_area_dict[nuc_id],
            'overlap_fraction_cyto':  overlap_area / cell_area_dict[cell_id]
        })

    # Filtering the df
    overlap_df = pd.DataFrame(records)
    overlap_df['overlap_fraction_sum'] = overlap_df['overlap_fraction_nuclear'] + overlap_df['overlap_fraction_cyto']

    pivot_overlap = overlap_df.pivot(index='nuclear_id', columns='label', values='overlap_fraction_sum').fillna(0)
    # Create a cost matrix (negative overlap_area for maximizing the overlap)
    cost_matrix = -pivot_overlap.values  # Maximize overlap -> minimize the negative overlap

    # Solve the linear assignment problem (minimization of negative overlap)
    nuc_idx, cell_idx = linear_sum_assignment(cost_matrix)

    merged_ids_list = [f"{a}_{b}" for a, b in zip(pivot_overlap.index[nuc_idx], pivot_overlap.columns[cell_idx])]
    overlap_df['merged_ids'] = overlap_df['nuclear_id'].astype(str) + '_' + overlap_df['label'].astype(str)
    overlap_df = overlap_df[overlap_df.merged_ids.isin(merged_ids_list)]

    # Keep nuclei / cells that are in overlap_df
    keep_mask_nuclear = np.isin(mask_nuclear, overlap_df.nuclear_id.tolist())
    filtered_mask_nuclear = np.where(keep_mask_nuclear, mask_nuclear, 0) # Set all other labels to 0 (background)

    keep_mask_cyto = np.isin(mask_cyto, overlap_df.label.tolist())
    mask_cyto_temp = np.where(keep_mask_cyto, mask_cyto, 0)

    return overlap_df, mask_cyto_temp , filtered_mask_nuclear



def filter_labels(mask: np.ndarray, cell_properties_df: pd.DataFrame, mask_nuclear: np.ndarray = None, intensity_min: float = None, intensity_max: float = None, cell_area_min: float = None, cell_area_max: float = None, 
                  nuclear_overlap: float = None, nuclear_area_min: float = None, nuclear_area_max: float = None):

    mask_temp = mask.copy()

    valid = pd.Series(True, index=cell_properties_df.index)

    if cell_area_min is not None:
        valid &= cell_properties_df.cell_area >= cell_area_min
    if cell_area_max is not None:
        valid &= cell_properties_df.cell_area <= cell_area_max
    if nuclear_area_min is not None:
        valid &= cell_properties_df.nuclear_area >= nuclear_area_min
    if nuclear_area_max is not None:
        valid &= cell_properties_df.nuclear_area <= nuclear_area_max
    if intensity_min is not None:
        valid &= cell_properties_df.ca_intensity >= intensity_min
    if intensity_max is not None:
        valid &= cell_properties_df.ca_intensity <= intensity_max
    if nuclear_overlap is not None:
        valid &= cell_properties_df.overlap_fraction_nuclear > nuclear_overlap

    valid_cell_labels = cell_properties_df[valid].label
    mask_temp[~np.isin(mask_temp, valid_cell_labels)] = 0

    if 'nuclear_id' in cell_properties_df.columns and mask_nuclear is not None:
        mask_nuclear_temp = mask_nuclear.copy()

        valid_nuclear_labels = cell_properties_df[valid].nuclear_id

        mask_nuclear_temp[~np.isin(mask_nuclear_temp, valid_nuclear_labels)] = 0
        cell_properties_df = cell_properties_df[cell_properties_df.label.isin(valid_cell_labels)]
        return mask_temp, mask_nuclear_temp, cell_properties_df
    
    else:
        cell_properties_df = cell_properties_df[cell_properties_df.label.isin(valid_cell_labels)]

        return mask_temp, cell_properties_df


# Extract raw traces
def extract_raw_traces(image_stack: np.ndarray, mask: np.ndarray, cell_properties_df: pd.DataFrame):
    all_props = []
    for frame_index, frame in enumerate(image_stack):
        df_prop = pd.DataFrame(regionprops_table(mask,intensity_image=frame,  properties=('label', 'intensity_mean')))
        df_prop['frame'] = frame_index  
        all_props.append(df_prop)

    intensity_df = pd.concat(all_props)

    intensity_df = intensity_df.groupby('label').agg({'intensity_mean': list}).reset_index()
    intensity_df = intensity_df.rename(columns={'intensity_mean': 'raw'})

    cell_properties_df = cell_properties_df.merge(intensity_df[['label', 'raw']],on='label', how='left') # Merge traces to original label properties df

    return cell_properties_df


def get_cell_properties(mask: np.ndarray, image: np.ndarray):
    cell_properties_df = pd.DataFrame(regionprops_table(mask, intensity_image=image, properties=['label','centroid', 'intensity_mean', 'area']))
    cell_properties_df = cell_properties_df.rename(columns={'intensity_mean': 'ca_intensity'})
    return cell_properties_df
    

def compute_auc(cell_properties_df: pd.DataFrame, start_frame: int  = None, end_frame: int = None, column: str='AUC'):
    auc_list = []
    for index, cell in cell_properties_df.iterrows():
        auc_list.append(sum(cell['dff'][start_frame:end_frame]))

    cell_properties_df[column] = auc_list
    return cell_properties_df

# Check for the column names, very risky function
'''def percentage_responding_baseline(cell_properties_df: pd.DataFrame):

    columns = list(set(cell_properties_df.columns).difference(['label', 'centroid-0', 'centroid-1', 'ca_intensity',
    'overlap_fraction_nuclear', 'nuclear_area', 'cell_area', 'nuclear_id','area',
    'raw', 'mask_path', 'filepath', 'Unnamed: 0', 'dff', 'baseline', 'AUC','AUC_kcl', 'response']))

    response_perc_well_df = cell_properties_df[columns].drop_duplicates().reset_index(drop=True)

    if 'marker' in cell_properties_df.columns:
        grouped = cell_properties_df.groupby(['stimulation', 'image_id', 'marker'])
    else:
        grouped = cell_properties_df.groupby(['stimulation', 'image_id'])

    # Count how many are "above" and "below" per group
    counts = grouped['response'].value_counts().unstack(fill_value=0)

    # Make sure "above" and "below" columns exist
    for col in ['above', 'below']:
        if col not in counts.columns:
            counts[col] = 0

    # Calculate proportions
    counts['proportion_positive_cells'] = round(counts['above'] / counts.sum(axis=1) * 100, 2)
    counts['proportion_negative_cells'] = round(counts['below'] / counts.sum(axis=1) * 100, 2)

    if 'marker' in cell_properties_df:
        # Merge back to your response_perc_well_df
        response_perc_well_df = response_perc_well_df.merge(
            counts[['proportion_positive_cells', 'proportion_negative_cells']],
            on=['stimulation', 'image_id', 'marker'],
            how='left'
        )
    else:
        response_perc_well_df = response_perc_well_df.merge(
            counts[['proportion_positive_cells', 'proportion_negative_cells']],
            on=['stimulation', 'image_id'],
            how='left'
        )


    response_perc_well_df = response_perc_well_df.sort_values(['stimulation','proportion_positive_cells'], ascending = False)
    if 'marker' in cell_properties_df:
        response_perc_rep_df = response_perc_well_df.groupby(['biological_replicate', 'stimulation', 'marker'], observed=True)['proportion_positive_cells'].mean().reset_index()
    else: 
        response_perc_rep_df = response_perc_well_df.groupby(['biological_replicate', 'stimulation'], observed=True)['proportion_positive_cells'].mean().reset_index()


    return response_perc_well_df, response_perc_rep_df'''


''''response_perc_well_df = cell_properties_df[columns].drop_duplicates().reset_index(drop=True)

    if 'marker' in cell_properties_df.columns:
        grouped = cell_properties_df.groupby(['stimulation', 'image_id', 'marker'])
    else:
        grouped = cell_properties_df.groupby(['stimulation', 'image_id'])

    # Count how many are "above" and "below" per group
    counts = grouped['response'].value_counts().unstack(fill_value=0)

    # Make sure "above" and "below" columns exist
    for col in ['above', 'below']:
        if col not in counts.columns:
            counts[col] = 0
'''


def percentage_responding(cell_properties_df: pd.DataFrame, analysis_type = None):
    columns = list(set(cell_properties_df.columns).difference(['label', 'centroid-0', 'centroid-1', 'ca_intensity',
    'overlap_fraction_nuclear', 'nuclear_area', 'cell_area', 'nuclear_id', 'area','AUC_kcl','AUC','response',
    'raw', 'mask_path', 'filepath', 'Unnamed: 0', 'dff', 'baseline', 'peak_location','rise_time','amplitude', 'decay_time', 'low_quality_peaks', 'prominence', 'width', 'frequency']))

    response_perc_well_df = cell_properties_df[columns].drop_duplicates().reset_index(drop=True)

    if 'marker' in cell_properties_df.columns:
        group_cols = ['stimulation', 'image_id', 'marker']
    else:
        group_cols = ['stimulation', 'image_id']

    # Count active/inactive cells
    counts = cell_properties_df.groupby(group_cols)['response'].value_counts().unstack(fill_value=0)

    col_dict = {'spontaneous': ['active', 'inactive'], 'baseline': ['above', 'below']}
    # Ensure both columns exist
    for col in col_dict[analysis_type]:
        if col not in counts.columns:
            counts[col] = 0

    if analysis_type == 'spontaneous':
        # Calculate proportion
        counts['proportion_responding'] = round(counts['active'] / counts.sum(axis=1) * 100, 2)
        response_col = 'proportion_responding'

    elif analysis_type == 'baseline':
        # Calculate proportions
        counts['proportion_positive_cells'] = round(counts['above'] / counts.sum(axis=1) * 100, 2)
        counts['proportion_negative_cells'] = round(counts['below'] / counts.sum(axis=1) * 100, 2)
        response_col = 'proportion_positive_cells'

    counts = counts.reset_index()
    response_perc_well_df = response_perc_well_df.merge(counts, on=group_cols, how='left')


    # Compute replicate-level means
    if 'marker' in cell_properties_df.columns:
        response_perc_rep_df = response_perc_well_df.groupby(
            ['biological_replicate', 'stimulation', 'marker'], observed=True
        )[response_col].mean().reset_index()
    else:
        response_perc_rep_df = response_perc_well_df.groupby(
            ['biological_replicate', 'stimulation'], observed=True
        )[response_col].mean().reset_index()

    return response_perc_well_df, response_perc_rep_df




'''

def percentage_responding_spontaneous(cell_properties_df: pd.DataFrame):
    columns = list(set(cell_properties_df.columns).difference(['label', 'centroid-0', 'centroid-1', 'ca_intensity',
    'overlap_fraction_nuclear', 'nuclear_area', 'cell_area', 'nuclear_id', 'area',
    'raw', 'mask_path', 'filepath', 'Unnamed: 0', 'dff', 'baseline', 'peak_location','rise_time','amplitude', 'decay_time', 'low_quality_peaks', 'prominence', 'width', 'frequency']))


    response_perc_well_df = cell_properties_df[columns].drop_duplicates().reset_index(drop=True)

    prop_active_cells_list = []

    for img in response_perc_well_df['image_id']:
        inactive = sum(cell_properties_df.loc[cell_properties_df.image_id == img].frequency == 0)
        active = sum(cell_properties_df.loc[cell_properties_df.image_id == img].frequency > 0)
        
        prop_active_cells_list.append(round((active / (inactive + active)*100),2))

    response_perc_well_df['proportion_active_cells'] = prop_active_cells_list

    response_perc_rep_df = response_perc_well_df.groupby(['biological_replicate', 'stimulation'], observed=True)['proportion_active_cells'].mean().reset_index()

    
    return response_perc_well_df, response_perc_rep_df
'''
def scale_spike_properties(cell_properties_df: pd.DataFrame):

    features_df = cell_properties_df[['frequency','width','rise_time','decay_time','amplitude']]
    X_scaled = StandardScaler().fit_transform(features_df.to_numpy())
    features_df = pd.DataFrame(X_scaled, columns = [f'{col}_scaled' for col in features_df.columns])
    cell_properties_df = cell_properties_df.reset_index(drop=True)
    for col in features_df.columns:
        cell_properties_df[col] = features_df[col]
    return cell_properties_df


def k_means_clustering(cell_properties_df: pd.DataFrame, array:np.array, n_clusters: int = None):
    
    kmeans = KMeans(n_clusters=n_clusters, random_state=10)
    cell_properties_df['cluster'] = kmeans.fit_predict(array)
    cell_properties_df['cluster'] = cell_properties_df['cluster'].astype('category')
    cell_properties_df['cluster'] = cell_properties_df['cluster'].cat.reorder_categories(sorted(cell_properties_df['cluster'].cat.categories))

    colors = list(plt.cm.viridis(np.linspace(0, 1, len(sorted(cell_properties_df.cluster.unique())))))
    colors = [mcolors.to_hex(c[:3]) for c in colors]
    colors_dict = {}
    for cluster in sorted(cell_properties_df.cluster.unique()):
        colors_dict[cluster] = colors[cluster]

    cell_properties_df['cluster_color'] = (cell_properties_df['cluster'].astype(int).map(colors_dict))


    return cell_properties_df

def cluster_centroids(cell_properties_df:pd.DataFrame):
    # Computes cluster centroids
    cluster_dict = {}
    for cluster in sorted(cell_properties_df.cluster.unique()):
        cluster_dict[cluster] = np.mean(np.array(cell_properties_df[cell_properties_df.cluster == cluster].dff.to_list()), axis=0)

    # Make color dict
    colors_dict = dict(cell_properties_df[['cluster','cluster_color']].drop_duplicates().values.tolist())

    return  cluster_dict, colors_dict







