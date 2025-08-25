import pandas as pd
import numpy as np
from scipy.signal import find_peaks



def spike(cell_properties_df: pd.DataFrame = None, prominence: float = None, amplitude_width_ratio: float = 0, imaging_interval: float = None, start_frame: int=None, end_frame: int=None):
    filtered_peak_list,  amplitude_list,  prominence_list, width_list, rise_time_list, decay_time_list, low_quality_peak_list = [],[],[], [], [], [],[]
    
    
    for index, cell in cell_properties_df.iterrows():
        trace = cell['dff']

        if end_frame is not None:
            if end_frame < 0: end_frame = None

        if start_frame is None:
            start_frame = 0
        if end_frame is None:
            end_frame = len(trace)  
  

        window_trace = trace[start_frame:end_frame]
        print('trace length spike calling:', len(window_trace))
        print(start_frame, end_frame)

        peak_location, peak_properties = find_peaks(trace, prominence = prominence, width=0, height=0)
        peak_properties.update({'peak_location':peak_location})
        peak_properties_df = pd.DataFrame(peak_properties)
        peak_properties_df = peak_properties_df[(peak_properties_df.peak_location >= start_frame) & (peak_properties_df.peak_location < end_frame)]
        peak_properties_df['amplitude'] = [trace[i] for i in peak_properties_df['peak_location']]
        peak_properties_df['amplitude_width_ratio'] = peak_properties_df['amplitude'] / peak_properties_df['widths']
        peak_properties_df['pass_qc']  = peak_properties_df['amplitude_width_ratio'] > amplitude_width_ratio
        peak_properties_df['rise_time'] = peak_properties_df['peak_location'] - peak_properties_df['left_ips']
        peak_properties_df['decay_time'] = peak_properties_df['right_ips'] - peak_properties_df['peak_location']
        low_quality_peaks = peak_properties_df[peak_properties_df['pass_qc'] == False].peak_location.tolist()
        
        low_quality_peak_list.append(low_quality_peaks)

        indices = [peak_properties_df['peak_location'].tolist().index(x) for x in peak_properties_df['peak_location'].tolist() if x in low_quality_peaks]
        
        filtered_peak_list.append(list(np.delete(peak_properties_df['peak_location'].to_numpy(), indices)))
        amplitude_list.append(list(np.delete(peak_properties_df['amplitude'].to_numpy(), indices)))
        prominence_list.append(list(np.delete(peak_properties_df['prominences'].to_numpy(), indices)))
        width_list.append(list(np.delete(peak_properties_df['widths'].to_numpy(), indices)))
        rise_time_list.append(list(np.delete(peak_properties_df['rise_time'].to_numpy(), indices)))
        decay_time_list.append(list(np.delete(peak_properties_df['decay_time'].to_numpy(), indices)))
    
    total_time_min = len(window_trace) * imaging_interval / 60

    cell_properties_df['peak_location'] = filtered_peak_list
    cell_properties_df['frequency'] = [len(x) / total_time_min for x in filtered_peak_list]
    cell_properties_df['amplitude'] = [np.mean(x) if len(x) != 0 else 0 for x in amplitude_list ]
    cell_properties_df['prominence'] = [np.mean(x) if len(x) != 0 else 0 for x in prominence_list ]
    cell_properties_df['width'] = [np.mean(x) * imaging_interval if len(x) != 0 else 0 for x in width_list ]
    cell_properties_df['rise_time'] = [np.mean(x) * imaging_interval if len(x) != 0 else 0 for x in rise_time_list ]
    cell_properties_df['decay_time'] = [np.mean(x) * imaging_interval if len(x) != 0 else 0 for x in decay_time_list ]
    cell_properties_df['low_quality_peaks'] = low_quality_peak_list

    cell_properties_df["response"] = np.where(cell_properties_df["frequency"] > 0, "active", "inactive")
    

    return cell_properties_df, start_frame, end_frame


def baseline_change(cell_properties_df: pd.DataFrame = None,  control_condition: str = None,  std_threshold: float = None ):

    std_dict, mean_dict, response_l = {},{}, []

    for exp_replicate in cell_properties_df.plate_id_biological_replicate.unique():
        mean_dict[exp_replicate] = cell_properties_df[(cell_properties_df.plate_id_biological_replicate == exp_replicate) & (cell_properties_df.stimulation == control_condition)].AUC.mean()
        std_dict[exp_replicate] = cell_properties_df[(cell_properties_df.plate_id_biological_replicate == exp_replicate) & (cell_properties_df.stimulation == control_condition)].AUC.std()

    for index, row in cell_properties_df.iterrows():
        if mean_dict[row.plate_id_biological_replicate] + std_threshold * std_dict[row.plate_id_biological_replicate] < row.AUC:
            response_l.append('above')

        elif mean_dict[row.plate_id_biological_replicate] - std_threshold * std_dict[row.plate_id_biological_replicate] > row.AUC:
            response_l.append('below')

        else: response_l.append('no response')

    cell_properties_df['response'] = response_l

    return cell_properties_df


'''
    cell_properties_df_list = []
    
    if len(cell_properties_df.stimulation.cat.categories.to_list()) == 2 and control_condition in cell_properties_df.stimulation.cat.categories.to_list():


        def classify_response(df, mean_auc_well, std_auc_well, exp_replicate):
            threshold_above = mean_auc_well + std_threshold * std_auc_well
            threshold_below = mean_auc_well - std_threshold * std_auc_well

            sample = df['plate_id_biological_replicate'] == exp_replicate

            df.loc[sample, 'response'] = np.where(
                df.loc[sample, 'AUC'] > threshold_above, "above",
                np.where(df.loc[sample, 'AUC'] < threshold_below, "below", "no response")
            )

        for exp_replicate in cell_properties_df['plate_id_biological_replicate'].unique():

            mean_auc_control = cell_properties_df[(cell_properties_df.plate_id_biological_replicate == exp_replicate) & (cell_properties_df.stimulation == control_condition)].AUC.mean()
            std_auc_control = cell_properties_df[(cell_properties_df.plate_id_biological_replicate == exp_replicate) & (cell_properties_df.stimulation == control_condition)].AUC.std()
            
            classify_response(cell_properties_df, mean_auc_control, std_auc_control, exp_replicate)

        #cell_properties_df_list.append(cell_properties_df)

        
        return cell_properties_df
    
    else: print('Input dataframe not valid')'''


'''
def baseline_change(cell_properties_df: pd.DataFrame = None,  control_condition: str = None, std_threshold: float = None ):

    
    response_df_control = cell_properties_df[cell_properties_df['condition'] == control_condition]
    response_df_treatment = cell_properties_df[cell_properties_df['condition'] != control_condition]
    
    def classify_response(df, mean_auc_well, std_auc_well, exp_replicate):
        threshold_above = mean_auc_well + std_threshold * std_auc_well
        threshold_below = mean_auc_well - std_threshold * std_auc_well

        sample = df['plate_id_biological_replicate'] == exp_replicate

        df.loc[sample, 'response'] = np.where(
            df.loc[sample, 'AUC'] > threshold_above, "above",
            np.where(df.loc[sample, 'AUC'] < threshold_below, "below", "no response")
        )

    for exp_replicate in response_df_control['plate_id_biological_replicate'].unique():
        mean_auc_well = response_df_control.loc[response_df_control.plate_id_biological_replicate == exp_replicate, 'AUC'].mean()
        std_auc_well = response_df_control.loc[response_df_control.plate_id_biological_replicate == exp_replicate, 'AUC'].std()

        classify_response(response_df_control, mean_auc_well, std_auc_well, exp_replicate)
        classify_response(response_df_treatment, mean_auc_well, std_auc_well, exp_replicate) 

    cell_properties_df = pd.concat([response_df_control,response_df_treatment])  
    return cell_properties_df

'''









def peak_calling(filename,  stimulation,biological_replicate, image_id, label_id, dff_traces, param_prominence_threshold, param_amplitude_width_ratio, output_folder):

    filtering = True
    peak_properties_dict = {}
    peak_properties_all_dict = {}
    print(dff_traces)
    print(type(dff_traces))
    
    # Call all peaks in the data and plot
    for cell_idx,ca_trace in enumerate(dff_traces.T):
        print('CELL index')
        print(cell_idx)
        #print(ca_trace)
        #print(type(ca_trace))
        #print(ca_trace.shape)
        #print(type(param_prominence_threshold))

        amplitude_list, low_quality_peak_list = [],[]

        peak_locations, peak_properties = find_peaks(ca_trace, prominence = param_prominence_threshold,width=0, height=0)
            
        peak_properties.update({'peak_locations':peak_locations})

        rise_time = [peak_locations[i] - peak_properties['left_ips'][i] for i in range(0,len(peak_locations))]
        decay_time = [peak_properties['right_ips'][i] - peak_locations[i] for i in range(0,len(peak_locations))]
        
        for i in peak_locations:
            amplitude_list.append(ca_trace[i])
    
        peak_properties.update({'peak_amplitudes':amplitude_list})
        peak_properties.update({'peak_rise_time':rise_time})
        peak_properties.update({'peak_decay_time':decay_time})
    
        for i in range(0, len(amplitude_list)):
            if float(peak_properties['peak_amplitudes'][i] / peak_properties['widths'][i]) < param_amplitude_width_ratio:
                low_quality_peak_list.append(peak_locations[i])

        peak_properties_all = peak_properties.copy()
        peak_properties_all.update({'low_quality_peaks':low_quality_peak_list})
        peak_properties_all_dict.update({cell_idx:peak_properties_all}) 
    
        # Filter low quality peaks
        if filtering == True:
            properties_to_remove_idx = np.intersect1d(peak_locations,np.array(low_quality_peak_list),return_indices=True)[1]
            if len(properties_to_remove_idx) > 0:

                for key in peak_properties.keys():
                    peak_properties[key] = np.delete(peak_properties[key],properties_to_remove_idx)
    
        peak_properties.update({'cell_id':label_id[cell_idx]})
        print(peak_properties)

        peak_properties_dict.update({cell_idx:peak_properties}) 
        print(peak_properties_dict)