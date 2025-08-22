
import numpy as np
import statistics
import pandas as pd





def pre_stimulation(cell_properties_df: pd.DataFrame = None, stimulation_frame: int = None):

    dff_traces_list, f0_list = [], []

    for i, cell in cell_properties_df.iterrows():
        raw_trace = cell['raw']
        mean_baseline = np.mean(raw_trace[:stimulation_frame]) 
        dff_traces_list.append([(F-mean_baseline)/mean_baseline for F in raw_trace])
        f0_list.append([mean_baseline] * len(raw_trace))

    cell_properties_df['dff'] = dff_traces_list
    cell_properties_df['baseline'] = f0_list

    return cell_properties_df





def sliding_window(cell_properties_df: pd.DataFrame = None, sliding_window_size: int = None, percentile_threshold: int = None):

    dff_traces_list, f0_list = [], []

    for i, cell in cell_properties_df.iterrows():
        f0_trace = []
        raw_trace = cell['raw']
        sliding_window = raw_trace[:sliding_window_size]

        pctl_thresh = np.percentile(sliding_window,percentile_threshold)
        f0 = statistics.mean([i for i in sliding_window if i < pctl_thresh]) # Baseline fluorescence 
        f0_temp = [f0] * sliding_window_size 
        f0_trace.extend(f0_temp)
        dff = [(value - f0) / f0 for value in sliding_window]
        
        for location,value in enumerate(raw_trace[sliding_window_size:], start=sliding_window_size+1):
            sliding_window = raw_trace[location-sliding_window_size-1:location-1]
            pctl_thresh = np.percentile(sliding_window,10)
            f0 = statistics.mean([i for i in sliding_window if i < pctl_thresh])
            f0_trace.append(f0)
            dff.append((value - f0) / f0)

        dff_traces_list.append(dff)
        f0_list.append(f0_trace)

    cell_properties_df['baseline'] = f0_list
    cell_properties_df['dff'] = dff_traces_list

    return cell_properties_df

        

