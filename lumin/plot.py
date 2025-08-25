import matplotlib.pyplot as plt
from skimage.segmentation import mark_boundaries
import os
import numpy as np
from skimage import io
from cv2 import circle
from skimage.color import gray2rgb
from skimage.util import img_as_ubyte
import pandas as pd
from typing import Literal
from sklearn.preprocessing import minmax_scale
import seaborn as sns
import matplotlib.lines as mlines
from pca import pca


import warnings
warnings.filterwarnings("ignore", message="Clipping input data to the valid range for imshow with RGB data")



def create_palette(unique_stimulations):
    colors = sns.color_palette(n_colors=len(unique_stimulations))
    hex_colors = ["#{:02X}{:02X}{:02X}".format(int(color[0]*255), int(color[1]*255), int(color[2]*255)) for color in colors]
    return dict(zip(unique_stimulations, hex_colors))


def segmentation(image: np.ndarray, mask: np.ndarray,  title_image: str = 'Image',  title_mask: str = 'Mask Outline',  output_path: str = None,  file_name: str = 'segmentation_plot'):

    with plt.rc_context({"figure.dpi": (350), 'figure.figsize':(10, 10)}):
        fig, ax = plt.subplots()

        # Plot 1: Input image
        #axes[0].imshow(mark_boundaries(image,np.zeros_like(mask),  color = (1.0, 0.0, 1.0) , mode='thick'), cmap='gray')
        #axes[0].axis("off")
        #axes[0].set_title(title_image)

        # Plot 2: Overlay boundaries on input
        ax.imshow(mark_boundaries(image,mask,  color = (1.0, 0.0, 1.0) , mode='thick'), cmap='gray')
        ax.axis("off")
        #"axes[0].set_title(title_mask)

        plt.savefig(os.path.join(output_path, f'{file_name}.pdf'), bbox_inches = 'tight')
        plt.close('all')



'''def segmentation(image: np.ndarray, mask: np.ndarray,  title_image: str = 'Image',  title_mask: str = 'Mask Outline',  output_path: str = None,  file_name: str = 'segmentation_plot'):

    with plt.rc_context({"figure.dpi": (350)}):
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))

        # Plot 1: Input image
        axes[0].imshow(mark_boundaries(image,np.zeros_like(mask),  color = (1.0, 0.0, 1.0) , mode='thick'), cmap='gray')
        axes[0].axis("off")
        axes[0].set_title(title_image)

        # Plot 2: Overlay boundaries on input
        axes[1].imshow(mark_boundaries(image,mask,  color = (1.0, 0.0, 1.0) , mode='thick'), cmap='gray')
        axes[1].axis("off")
        axes[1].set_title(title_mask)

        plt.tight_layout()
        plt.savefig(os.path.join(output_path, f'{file_name}.pdf'), bbox_inches = 'tight')'''


def overlay_labels(image: np.ndarray, mask: np.ndarray,  cell_properties_df: pd.DataFrame, mask_nuclear: np.ndarray = None,   output_path: str = None,  file_name: str = 'overlay_plot'):

    image_rgb = gray2rgb(img_as_ubyte(image))

    with plt.rc_context({"figure.dpi": (350), 'figure.figsize':(10, 10)}):

        # Extract labels and their loaction
        centroids = cell_properties_df[['centroid-0','centroid-1']].values.tolist()
        labels = cell_properties_df['label'].values.tolist()

        fig, ax = plt.subplots()

        # Color the mask
        mask_bool = mask.astype(bool)
        roi_rgb = np.array([102, 171, 217], dtype=np.uint8)  
        image_rgb[mask_bool] = roi_rgb

        # Create figure
        if mask_nuclear is not None:
            img = mark_boundaries(mark_boundaries(image_rgb, mask,   color = (38, 86, 245) , mode='thick' ), mask_nuclear,  color = (1.0, 0.0, 1.0) , mode='thick')
            ax.imshow(img, cmap='gray')

        else:
            img = mark_boundaries(image_rgb, mask,  color = (38, 86, 245) , mode='thick' )
            ax.imshow(img, cmap='gray')
            
            kwargs={'horizontalalignment':'center','fontsize':6}
            for i in range(len(centroids)): # Numbering the labels
                ax.text(centroids[i][1], centroids[i][0]+6,  str(labels[i]),color='#eb1717',**kwargs)

        plt.axis('off')
        if output_path is not None:
            plt.savefig(os.path.join(output_path, f'{file_name}.pdf'), bbox_inches = 'tight')
            plt.close('all')
        else: 
            return img



def overlaid_traces(cell_properties_df: pd.DataFrame,trace: str = None, mean: bool= False):
    # dff traces with mean       
    filename = cell_properties_df['filename'].iloc[0]
    with plt.rc_context({"figure.dpi": 350}):
        fig, ax = plt.subplots(figsize=(15,10))

        for i, cell in cell_properties_df.iterrows():
             ax.plot(cell[trace],linewidth=1)
        
        if mean == True:
            ax.plot(np.array(cell_properties_df[trace].tolist()).mean(axis=1), color='black',linewidth = '3')

    return ax



def overlaid_traces_two_groups(cell_properties_df: pd.DataFrame = None, trace: str = None, output_path: str = None, mean: bool= False, control_condition: str = None, treatment_condition: str = None, start_frame: int = None, end_frame: int = None, stimulation_frame: int = None, ax = None, palette:dict = None, imaging_interval:float = None, kcl_frame:int=None):
    with plt.rc_context({"figure.dpi": 350}):
    
        if ax is None:
            fig, ax = plt.subplots()
        
        unique_stimulations = cell_properties_df["stimulation"].cat.categories

        if palette is None: 
            palette = create_palette(unique_stimulations)

        def darken_color(hex_color, factor=0.7):
            # Convert HEX to RGB
            hex_color = hex_color.lstrip("#")
            r, g, b = [int(hex_color[i:i+2], 16) for i in (0, 2, 4)]
        
            # Apply darkening factor
            r = int(r * factor)
            g = int(g * factor)
            b = int(b * factor)
        
            # Clamp values between 0 and 255
            r, g, b = [max(0, min(255, v)) for v in (r, g, b)]
        
            # Convert back to HEX
            return "#{:02X}{:02X}{:02X}".format(r, g, b)
    
        if ax is None:
            fig, ax = plt.subplots()

        

        ax.spines['right'].set_color(None)
        ax.spines['top'].set_color(None)

        #fig, ax = plt.subplots()
        ax.plot(np.array(cell_properties_df[cell_properties_df.stimulation==treatment_condition][trace].tolist()).T, color=palette[treatment_condition],lw=0.4 , alpha=0.5)
        ax.plot(np.array(cell_properties_df[cell_properties_df.stimulation==control_condition][trace].tolist()).T, color=palette[control_condition],lw=0.4 , alpha=0.5)
    
        ax.plot(np.mean(np.array(cell_properties_df[cell_properties_df.stimulation==treatment_condition][trace].tolist()).T, axis=1), color=darken_color(palette[treatment_condition]),lw=2)
        ax.plot(np.mean(np.array(cell_properties_df[cell_properties_df.stimulation==control_condition][trace].tolist()).T, axis=1), color=darken_color(palette[control_condition]),lw=2)

    
        # Dummy lines for legend
        ax.plot([], [], color=palette[control_condition], lw=2, label=control_condition)
        ax.plot([], [], color=palette[treatment_condition], lw=2, label=treatment_condition)
    
        ax.axvspan(start_frame, end_frame, color='#555555', alpha=0.48, label='Analysis window', lw=0)
        ax.axvspan(0, stimulation_frame, color='lightgray', alpha=0.5, label='Baseline recording', lw=0)

        if kcl_frame is not None:
            ax.axvline(x=kcl_frame,  color='red', linewidth=1.5, ls='--', label = 'KCl')
        
        ax.legend(loc="upper left", ncol=1, prop={'size': 14})
    
        ax.set_ylim(-1,25)
        ax.tick_params(labelsize=12)
    
        ax.set_ylabel(r"Ca$^{2+}$ signal ($\Delta F/F_0$)", fontsize='xx-large')

        # Set ticks
        num_ticks = 6
        rough_step = len(np.mean(np.array(cell_properties_df[cell_properties_df.stimulation==treatment_condition][trace].tolist()).T, axis=1)) / (num_ticks - 1)
        step = int(np.round(rough_step / 10) * 10) # Round to nearest 10 for clean tick spacing
        xticks = np.arange(0, len(np.mean(np.array(cell_properties_df[cell_properties_df.stimulation==treatment_condition][trace].tolist()).T, axis=1)) , step)
        tick_labels = (xticks * imaging_interval).astype(int)
        ax.set_xticks(xticks)
        ax.set_xticklabels(tick_labels, rotation=0)
        ax.set_xlabel('Time (s)', fontsize='xx-large')  # Adjust labels  
        
        return ax



def cellwise_traces(cell_properties_df: pd.DataFrame, trace: str = None, baseline: bool = False, spikes: bool = False,
                    spikes_mode: Literal["all", "filtered"] = "all", output_path: str = None):


    filename = cell_properties_df['filename'].iloc[0]
    plot_list = []

    cell_counter = 0
    total_cells = cell_properties_df.shape[0]

    while cell_counter < total_cells:
        # Determine how many cells to plot in this chunk
        cells_in_chunk = min(20, total_cells - cell_counter)

        # Subset the DataFrame for this chunk
        chunk_df = cell_properties_df.iloc[cell_counter:cell_counter + cells_in_chunk]

        # Create figure and axes
        fig, ax = plt.subplots(nrows=cells_in_chunk, ncols=1, figsize=(15, 20), sharex=True, sharey=True, constrained_layout=False)
        if cells_in_chunk == 1:
            ax = [ax]

        for i, (idx, cell) in enumerate(chunk_df.iterrows()):
            ax_plot = ax[i]
            ax_plot.plot(cell[trace], color= '#1c1fb0',  linewidth=2)

            if baseline and 'baseline' in cell:
                ax_plot.plot(cell['baseline'], color = '#fc4a2b', alpha=0.7,  linewidth=1)

            ax_plot.set_ylabel(f"{cell['label']}.", rotation=0, fontsize=12, labelpad=10)
            ax_plot.yaxis.tick_right()
            ax_plot.set_ylim(-0.2, cell_properties_df[trace].explode().max())

            if spikes and 'peak_location' in cell and isinstance(cell['peak_location'], (list, np.ndarray)):
                for peak in cell.peak_location:
                    ax_plot.axvline(x=peak, ymin=-0.2, ymax = cell_properties_df[trace].explode().max(), ls='--', color='#fca02b', linewidth=2)
                if spikes_mode == "all" and 'low_quality_peaks' in cell and isinstance(cell['low_quality_peaks'], (list, np.ndarray)):
                    for peak in cell.low_quality_peaks:
                        ax_plot.axvline(x=peak,ymin=-0.2, ymax = cell_properties_df[trace].explode().max(), ls='--', color="#e24438", linewidth=2)

        #fig.tight_layout()
        fig.subplots_adjust(bottom=0.05)
        color_dict = {}
        if trace == 'dff':
            color_dict['dff signal'] = '#1c1fb0'
        elif trace == 'raw':
            color_dict['raw signal'] = '#1c1fb0'

        if baseline:
            color_dict['baseline'] = '#e4746d'

        if spikes:
            if spikes_mode == "all":
                color_dict['spike'] = '#fca02b'
                color_dict['removed spike'] = "#e24438"
            elif spikes_mode == "filtered":
                color_dict['spike'] = "#fca02b"

        # Create legend handles manually
        legend_elements = [
            mlines.Line2D([0], [0], color=color, lw=2, label=label)
            for label, color in color_dict.items()
        ]

        # Add the custom legend
        fig.legend(handles=legend_elements, loc="lower center", ncol=len(legend_elements), prop={'size': 20})


        plot_list.append(fig)



        # Save if needed
        if output_path is not None:
            plt.savefig(os.path.join(output_path, f"{(cell_counter // 20) + 1}_{filename}.pdf"), dpi=350, bbox_inches="tight")
            plt.close('all')


        cell_counter += cells_in_chunk

    return plot_list


def overlay_events(cell_properties_df: pd.DataFrame):
    stack = io.imread(cell_properties_df['filepath'].values[0])
    mask = np.zeros_like(stack, dtype=np.uint16)
    for index, row in cell_properties_df.iterrows():
        if 'marker' in cell_properties_df.columns:
            locations, centroids, area = row['peak_location'], (int(row['centroid-1']), int(row['centroid-0'])), row['area']
        else:
            locations, centroids, area = row['peak_location'], (int(row['centroid-1']), int(row['centroid-0'])), row['cell_area']
            
        if len(locations) > 0:
            radius = int(np.sqrt(area / np.pi))
            n_frames = mask.shape[0]
            for frame_idx in locations:
                for offset in [-1, 0, 1, 2]:
                    idx = frame_idx + offset
                    if 0 <= idx < n_frames:   # only draw if within valid range
                        circle(mask[idx], centroids, int(radius*1.5), color=100, thickness=2)


    return stack, mask


'''def draw_vertical_brace(ax, yspan, xx, text):
    """Draws an annotated vertical brace on the axes."""
    ymin, ymax = yspan
    yspan = ymax - ymin
    #yspan = abs(ymax - ymin)
    ax_ymin, ax_ymax = ax.get_ylim()
    yax_span = ax_ymax - ax_ymin
    #yax_span = abs(ax_ymax - ax_ymin)

    xmin, xmax = ax.get_xlim() 
    xspan = xmax - xmin
    resolution = int(yspan / yax_span * 100) * 2 + 1  # guaranteed uneven
    beta = 300. / yax_span  # the higher this is, the smaller the radius

    y = np.linspace(ymin, ymax, resolution)
    y_half = y[:int(resolution / 2) + 1]
    x_half_brace = (1 / (1. + np.exp(-beta * (y_half - y_half[0])))
                    + 1 / (1. + np.exp(-beta * (y_half - y_half[-1]))))
    x = np.concatenate((x_half_brace, x_half_brace[-2::-1]))
    x = xx + (.05 * x - .01) * xspan   # adjust horizontal position

    ax.autoscale(False)
    ax.plot(x, y, color='black', lw=1, clip_on=False)  # <-- key fix here

    ax.text(xx + .09 * xspan, (ymax + ymin) / 2., text,
            ha='left', va='center', rotation='vertical', clip_on=False) 
    

def generate_beeswarm(distributions: list, tick_labels: list, max_plot_width: int = 1, alpha=0.7,
                            number_of_segments=12,
                            separation_between_plots=0.1,
                            separation_between_subplots=0.1,
                            vertical_limits=None,
                            grid=False,
                            remove_outlier_above_segment=None,
                            remove_outlier_below_segment=None,
                            y_label=None,
                            title=None, ax=None, palette = None):

    #unique_conditions = tick_labels
    
    #if palette is None: 
        #palette = create_palette(unique_conditions)
        
    #palette = {cond: palette[cond] for cond in unique_conditions}

    number_of_plots = len(distributions)

    ax.set_xlim(left=0, right=number_of_plots * (max_plot_width + separation_between_plots) + separation_between_plots)

    ticks = [separation_between_plots + max_plot_width / 2 + (max_plot_width + separation_between_plots) * i
             for i in range(0, number_of_plots)]
    
    max_counts = 0.0
    counts_filled_list = []
    segment_indices_list = []
    for i in range(len(distributions)):
        distribution = distributions[i]
        segments = np.linspace(np.min(distribution), np.max(distribution), number_of_segments + 1)[1:-1]
        
        segment_indices = number_of_segments - 1 - np.where(segments[:, None] >= distribution[None, :], 1, 0).sum(0)
        
        if remove_outlier_above_segment:
            a = remove_outlier_above_segment[i]
            distribution = distribution[segment_indices <= a]
            segment_indices = segment_indices[segment_indices <= a]

        if remove_outlier_below_segment:
            b = remove_outlier_below_segment[i]
            distribution = distribution[segment_indices >= b - 1]
            segment_indices = segment_indices[segment_indices >= b - 1]
        segment_indices_list.append(segment_indices)

        values, counts = np.unique(segment_indices, return_counts=True)
        if np.max(counts) > max_counts:
            max_counts = np.max(counts)
        counts_filled = []
        j = 0
        for k in range(number_of_segments):
            if k in values:
                counts_filled.append(counts[j])
                j += 1
            else:
                counts_filled.append(0)
        counts_filled_list.append(counts_filled)


    for i in range(len(distributions)):    
        variances = (max_plot_width / 2) * (counts_filled_list[i] / max_counts)
        jitter_unadjusted = np.random.uniform(-1, 1, len(distributions[i])) 
        jitter = np.take(variances, segment_indices_list[i]) * jitter_unadjusted

        ax.scatter(jitter + ticks[i], distributions[i], alpha=0.85, s=3, linewidth=0.2, edgecolor='black', c = list(palette.values())[i])


    ax.spines['right'].set_color(None)
    ax.spines['top'].set_color(None)
    ax.tick_params(axis='y', labelsize='x-small') 
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, size='large')
    #ax.set_yticks(fontsize='x-large')

    

    return ax, ticks

def beeswarm(cell_properties_df: pd.DataFrame, control_condition: str = None,  std_threshold: float = None, ax = None, palette: dict = None, brace: bool = False, ycolumn: str = None, control_condition_mean: bool = False):
    with plt.rc_context({"figure.dpi": 350}):
        cell_properties_df['stimulation'] = cell_properties_df['stimulation'].cat.remove_unused_categories()
        
        if ax is None:
            fig, ax = plt.subplots()

        if palette is None: 
            palette = create_palette(cell_properties_df["stimulation"].cat.categories)
            
        # Handle this to work with several conditions, check the colors
        distributions = [cell_properties_df[cell_properties_df.stimulation == cond][ycolumn].values for cond in cell_properties_df.stimulation.cat.categories]
        palette = {cond: palette[cond] for cond in cell_properties_df.stimulation.cat.categories}
        ax , ticks =  generate_beeswarm(distributions, tick_labels=cell_properties_df.stimulation.cat.categories, number_of_segments=200,
                                grid=False, ax=ax, palette = palette)
        
        if control_condition_mean == True:
            ax.axhline(np.mean(distributions[0]), color='#636363', linewidth = 1, label=f'Mean ({control_condition})')
        
        else:
            for i, dist in enumerate(distributions):
                center = ticks[i]
                ax.boxplot(
                    dist,
                    positions=[center],         # manual x-location
                    widths=0.5,
                    showfliers=False,
                    showcaps=True,
                    patch_artist=True,
                    boxprops={"facecolor": "none", "edgecolor": "black", "linewidth": 0.5},
                    whiskerprops={"color": "black", "linewidth": 0.5},
                    capprops={"color": "black", "linewidth": 0.5},
                    medianprops={"color": "black", "linewidth": 0.7}
                )

            
        ax.set_xticks(ticks)
        ax.set_xticklabels(cell_properties_df.stimulation.cat.categories)

        if std_threshold is not None:
            ax.axhline(np.std(distributions[0])*std_threshold, color='red', linewidth = 1, linestyle='--', label='std threshold')
            ax.legend(loc='upper left',bbox_to_anchor=(-0.03, 1), frameon=False, fontsize='small', handlelength=1.2, handletextpad=0.3)

        if brace == True:
            max_value = np.max([np.max(distributions[0]), np.max(distributions[1])])
            std_threshold = np.std(distributions[0]) * std_threshold
            if max_value > std_threshold:
                draw_vertical_brace(ax, (std_threshold, max_value),2.3, 'Responding')
        
        return ax  '''

def draw_vertical_brace(ax, yspan, xx, text):
    """Draws an annotated vertical brace on the axes."""
    ymin, ymax = yspan
    yspan = ymax - ymin
    ax_ymin, ax_ymax = ax.get_ylim()
    yax_span = ax_ymax - ax_ymin

    xmin, xmax = ax.get_xlim()
    xspan = xmax - xmin
    resolution = int(yspan / yax_span * 100) * 2 + 1
    beta = 300. / yax_span

    y = np.linspace(ymin, ymax, resolution)
    y_half = y[:int(resolution / 2) + 1]
    x_half_brace = (1 / (1. + np.exp(-beta * (y_half - y_half[0])))
                    + 1 / (1. + np.exp(-beta * (y_half - y_half[-1]))))
    x = np.concatenate((x_half_brace, x_half_brace[-2::-1]))
    x = xx + (.05 * x - .01) * xspan

    ax.autoscale(False)
    ax.plot(x, y, color='black', lw=1, clip_on=False)
    ax.text(xx + .09 * xspan, (ymax + ymin) / 2., text,
            ha='left', va='center', rotation='vertical', clip_on=False)
    


def generate_beeswarm(distributions, positions, ax, colors,
                      max_plot_width=1, number_of_segments=200):
    """Low-level beeswarm plotting given explicit positions + colors."""

    max_counts = 0.0
    counts_filled_list = []
    segment_indices_list = []

    for distribution in distributions:
        if len(distribution) == 0:
            segment_indices_list.append([])
            counts_filled_list.append([0]*number_of_segments)
            continue

        segments = np.linspace(np.min(distribution),
                               np.max(distribution),
                               number_of_segments + 1)[1:-1]
        segment_indices = number_of_segments - 1 - np.where(
            segments[:, None] >= distribution[None, :], 1, 0).sum(0)

        values, counts = np.unique(segment_indices, return_counts=True)
        max_counts = max(max_counts, np.max(counts))
        counts_filled = [counts[values.tolist().index(k)] if k in values else 0
                         for k in range(number_of_segments)]
        counts_filled_list.append(counts_filled)
        segment_indices_list.append(segment_indices)

    for i, dist in enumerate(distributions):
        if len(dist) == 0:
            continue
        variances = (max_plot_width / 2) * (np.array(counts_filled_list[i]) / np.max(counts_filled_list[i]))
        jitter_unadjusted = np.random.uniform(-1, 1, len(dist))
        jitter = np.take(variances, segment_indices_list[i]) * jitter_unadjusted

        ax.scatter(jitter + positions[i], dist, alpha=0.85, s=5,
                   linewidth=0.2, edgecolor='black', c=colors[i])

    return ax


def beeswarm(cell_properties_df: pd.DataFrame,
             x: str = "stimulation",
             y: str = None,
             hue: str = None,
             palette: dict = None,
             ax=None,
             max_plot_width=0.8,
             hue_separation=0.3,
             separation_between_plots=1,
             brace=False,
             std_threshold=None, 
             control_condition: str = None,
             control_condition_mean: bool = False):
    """Beeswarm plot with optional hue grouping."""

    if y is None:
        raise ValueError("You must provide the y-axis column name.")

    if ax is None:
        fig, ax = plt.subplots()

    # categories
    x_categories = (cell_properties_df[x].cat.categories
                    if hasattr(cell_properties_df[x], "cat")
                    else np.unique(cell_properties_df[x]))
    if hue:
        hue_categories = (cell_properties_df[hue].cat.categories
                          if hasattr(cell_properties_df[hue], "cat")
                          else np.unique(cell_properties_df[hue]))
    else:
        hue_categories = [None]

    # palette
    if palette is None:
        unique = hue_categories if hue else x_categories
        cmap = plt.get_cmap("tab10")
        palette = {u: cmap(i % 10) for i, u in enumerate(unique)}

    distributions = []
    colors = []
    positions = []

    # ticks for main x groups
    base_ticks = [separation_between_plots + i * separation_between_plots
                  for i in range(len(x_categories))]

    # offsets for hue groups
    n_hue = len(hue_categories)
    offset_range = np.linspace(-hue_separation,hue_separation, n_hue) if n_hue > 1 else [0]

    for i, xc in enumerate(x_categories):
        for j, hc in enumerate(hue_categories):
            if hue:
                mask = (cell_properties_df[x] == xc) & (cell_properties_df[hue] == hc)
                data = cell_properties_df.loc[mask, y].values
                color = palette[hc]
            else:
                mask = (cell_properties_df[x] == xc)
                data = cell_properties_df.loc[mask, y].values
                color = palette[xc]

            distributions.append(data)
            colors.append(color)
            positions.append(base_ticks[i] + offset_range[j])

    # draw beeswarm
    generate_beeswarm(distributions, positions, ax, colors,
                      max_plot_width=max_plot_width)
    if control_condition_mean == True:
            ax.axhline(np.mean(distributions[0]), color='#636363', linewidth = 1, label=f'Mean ({control_condition})')
        
    else:
        for i, dist in enumerate(distributions):
            center = positions[i]
            ax.boxplot(
                dist,
                positions=[center],         # manual x-location
                widths=max_plot_width-0.1,
                showfliers=False,
                showcaps=True,
                patch_artist=True,
                boxprops={"facecolor": "none", "edgecolor": "black", "linewidth": 0.5},
                whiskerprops={"color": "black", "linewidth": 0.5},
                capprops={"color": "black", "linewidth": 0.5},
                medianprops={"color": "black", "linewidth": 0.7}
            )

    # x labels only for main groups
    ax.set_xticks(base_ticks)
    ax.set_xticklabels(x_categories)

    # legend if hue is used
    if hue:
        handles = [plt.Line2D([0], [0], marker='o', color='w',
                              markerfacecolor=palette[hc], label=str(hc),
                              markersize=6)
                   for hc in hue_categories]
        ax.legend(handles=handles, title=hue, frameon=False)

    # optional std threshold line
    if std_threshold is not None:
        ref_vals = cell_properties_df.loc[cell_properties_df[x] == x_categories[0], y].values
        ax.axhline(np.std(ref_vals)*std_threshold, color='red',
                   linewidth=1, linestyle='--', label='std threshold')
        ax.legend(loc='upper left',bbox_to_anchor=(-0.03, 1), frameon=False, fontsize='small', handlelength=1.2, handletextpad=0.3)
        
        
    #if std_threshold is not None:
        #ax.axhline(np.std(distributions[0])*std_threshold, color='red', linewidth = 1, linestyle='--', label='std threshold')

    # optional brace
    if brace and len(x_categories) > 1:
        d0 = cell_properties_df.loc[cell_properties_df[x] == x_categories[0], y].values
        d1 = cell_properties_df.loc[cell_properties_df[x] == x_categories[1], y].values
        max_value = np.max([np.max(d0), np.max(d1)])
        thr = np.std(d0) * std_threshold if std_threshold else np.mean(d0)
        if max_value > thr:
            draw_vertical_brace(ax, (thr, max_value), 2.45, 'Responding')


    ax.spines['right'].set_color(None)
    ax.spines['top'].set_color(None)

    return ax
    

def heatmap(cell_properties_df: pd.DataFrame,  imaging_interval: float = None, palette: dict = None, cmap: str = 'plasma', minmax_bool: bool = True, vmin:int = None, vmax:int = None):
    cell_properties_df['stimulation'] = cell_properties_df['stimulation'].cat.remove_unused_categories()
    
    unique_stimulations = cell_properties_df["stimulation"].cat.categories

    if palette is None: 
        palette = create_palette(unique_stimulations)
    
    # Set colors
    color_list = cell_properties_df['stimulation'].astype(str).map(palette).values
    


    # Extract dff traces and min-max scale
    ca_traces = np.stack(cell_properties_df['dff'].values)
    if minmax_bool == True:
        ca_traces = minmax_scale(ca_traces, axis=1)
    
    if vmin is None:
        vmin = ca_traces.min()
    if vmax is None:
        vmax = ca_traces.max()

    kws = dict(cbar_kws=dict(ticks=[vmin, vmax], orientation='vertical'))

    #kws = dict(cbar_kws=dict(ticks=[ca_traces.min(), ca_traces.max()], orientation='vertical'))

    # Plotting
    with plt.rc_context({"figure.dpi": 350}):
        cm = sns.clustermap(ca_traces,
                            col_cluster=False,
                            row_cluster=False,
                            standard_scale=None,
                            figsize=(2, 3),
                            xticklabels=True,
                            yticklabels=False,              
                            cmap=cmap,vmin=vmin, vmax=vmax,
                            row_colors=color_list, **kws)

        cm.ax_row_dendrogram.set_visible(False)
        cm.ax_col_dendrogram.set_visible(False)
        cm.ax_heatmap.tick_params(labelsize=3)
        cm.gs.update(wspace=0.035)
        
        cm.ax_cbar.set_yticklabels(['Min',  'Max'])
        cm.ax_cbar.set_position([0.95, 0.159, 0.03, 0.15])
        cm.ax_cbar.tick_params(labelsize=6,axis='y', width=0, length=0,pad=1)

        # Set ticks
        num_ticks = 6
        rough_step = ca_traces.shape[1] / (num_ticks - 1)
        step = int(np.round(rough_step / 10) * 10) # Round to nearest 10 for clean tick spacing
        xticks = np.arange(0, ca_traces.shape[1] , step)
        tick_labels = (xticks * imaging_interval).astype(int)
        cm.ax_heatmap.set_xticks(xticks)
        cm.ax_heatmap.set_xticklabels(tick_labels, rotation=0, fontsize=6)
        cm.ax_heatmap.tick_params(pad=1, width=0.5, length=2)

        handles2 = [mlines.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=5, label=label)
                              for label, color in palette.items()]


        # Add the second legend
        legend2 = cm.ax_heatmap.legend(
            handles=handles2, loc='upper right', bbox_to_anchor=(1.45, 0.8), frameon=False,
            handletextpad=0, prop={'size': 6}
        ).get_title().set_position((-10,0))

        

        # Set xlabel
        ax = cm.ax_heatmap
        ax.set_xlabel('Time (s)', fontsize=8, labelpad=2)

        return ax

    
def all_conditions_barplot(dataframe: pd.DataFrame, palette: dict = None, ycolumn: str = None, xcolumn: str = None, hue:str = None, alpha: float = 1):

    if hue is None: color_cat = xcolumn
    else: color_cat = hue
    dataframe[color_cat] = dataframe[color_cat].cat.remove_unused_categories()

    unique_stimulations = dataframe[color_cat].cat.categories
    
    if palette is None: palette = create_palette(unique_stimulations)

    with plt.rc_context({"figure.dpi": 350}):
        fig, ax = plt.subplots()  # Create a dedicated figure
        sns.barplot(
            x=xcolumn, 
            y=ycolumn, 
            hue=hue,
            data=dataframe, 
            capsize=.2,
            errwidth=1,
            errcolor='black',
            edgecolor='black',
            palette=palette,
            ax=ax,
            alpha=alpha
        )

        sns.stripplot(
            data=dataframe, 
            x=xcolumn, 
            y=ycolumn, 
            color='black', 
            dodge=True, 
            hue=hue,
            ax=ax,
            alpha=0.7,
            legend=False, 
            marker='o'
        )

        ax.spines['right'].set_color(None)
        ax.spines['top'].set_color(None)

        handles = [
            mlines.Line2D([], [], color=palette[cond], marker="o", linestyle="None", markersize=8, label=cond)
            for cond in unique_stimulations
        ]
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1, 0.5),
                frameon=False, handletextpad=-0.1, fontsize='x-large')

    return  ax

def two_conditions_barplot(response_perc_mean_df: pd.DataFrame, palette: dict = None, x: str = None, y: str = None, hue: str = None):
    response_perc_mean_df['stimulation'] = response_perc_mean_df['stimulation'].cat.remove_unused_categories()

    if palette is None: 
        palette = create_palette(response_perc_mean_df["stimulation"].cat.categories)

    with plt.rc_context({"figure.dpi": 350}):
    
        # Create the figure and axis
        fig, ax = plt.subplots()
        
        # Plot lineplot first (on top)
        if hue is None:
            sns.lineplot(
                data=response_perc_mean_df, 
                style='biological_replicate', 
                x=x, # stim
                y=y, # "proportion_positive_cells"
                estimator=None,  
                markers=['o'], 
                markeredgewidth=0, 
                dashes={line: (2, 2) for line in response_perc_mean_df['biological_replicate'].unique()}, 
                markersize=3.5, 
                lw=0.5, 
                ax=ax, 
                color='black',             # Set the color to black
                alpha=0.7
            )
        
        # Create bar plot without error bars (plotted below)
        sns.barplot(
            data=response_perc_mean_df, 
            x= x, #"stimulation", 
            y=y, #"proportion_positive_cells", 
            hue=hue,
            palette=palette,
            ci=None, ax=ax  # Disables error bars
        )

        if hue is not None:
            sns.stripplot(
                data=response_perc_mean_df, 
                x=x, 
                y=y, 
                hue=hue, 
                dodge=True,           
                jitter=True,           
                marker="o", 
                size=3.5, 
                linewidth=0, 
                color='black',
                ax=ax,
                alpha=0.8, legend=False
            )
        
        # Customize plot
        plt.xlabel("")
        plt.xticks(fontsize='large')
        
        plt.yticks(fontsize='x-small')
    
        plt.legend(loc='center left',bbox_to_anchor=(1,0.5),frameon=False)
        
        if hue is None:
            ax.get_legend().remove()
    
        ax.spines['right'].set_color(None)
        ax.spines['top'].set_color(None)
    
        return ax
    
'''def cell_count_barplot(df: pd.DataFrame , palette: dict, x:str=None, y:str=None, hue:str=None ):
    with plt.rc_context({"figure.dpi": 350}):
        fig, ax = plt.subplots()
        
        sns.barplot(df, y = y, x = x, hue=hue, palette=palette,
                        capsize=.3, errwidth=1, ax = ax, edgecolor='black')
        sns.stripplot(
            x=x, 
            y=y, hue=hue, color='black',
            data=df, dodge=True, alpha=0.7, ax=ax, legend=None
        )
        plt.xticks(rotation=90)
        ax.spines['right'].set_color(None)
        ax.spines['top'].set_color(None)
        
        plt.ylabel('Cells per image')
        plt.xlabel('')
        
        handles = [
            mlines.Line2D([], [], color=palette[cond], marker="o", linestyle="None", markersize=8, label=cond)
            for cond in df.stimulation.cat.categories
        ]
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1, 0.5),
                frameon=False, handletextpad=-0.1, fontsize='x-large')
        
        plt.title('')

        return ax'''
    
def biplot(cell_properties_df: pd.DataFrame, arrow_scale_factor: float = 1, palette: dict = None):

    # Scale
    features_df = cell_properties_df[['frequency_scaled','width_scaled','rise_time_scaled','decay_time_scaled','amplitude_scaled']]
    features_df.columns = ['frequency','width','rise_time','decay_time','amplitude']
    
    if palette is None: 
            palette = create_palette(cell_properties_df["stimulation"].cat.categories)

    # PCA
    pca_model = pca(n_components=2, verbose=0)
    pca_results = pca_model.fit_transform(X=features_df, verbose=0)
    cell_properties_df[['PC1','PC2']] = pca_results['PC']

    with plt.rc_context({"figure.dpi": 350}):

        # Create the main figure and axis
        fig = plt.figure(figsize=(12, 12))
        gs = fig.add_gridspec(5, 5, wspace=0.01, hspace=0.01)
        
        # Main biplot axis
        ax_main = fig.add_subplot(gs[1:, :-1])
        
        # Top marginal axis for x-axis density
        ax_marg_x = fig.add_subplot(gs[0, :-1], sharex=ax_main)
        
        # Right marginal axis for y-axis density
        ax_marg_y = fig.add_subplot(gs[1:, -1], sharey=ax_main)
        
        # Completely hide ticks, labels, and gridlines on marginal plots
        for ax in [ax_marg_x, ax_marg_y]:
            ax.tick_params(axis="both", which="both", length=0, labelbottom=False, labelleft=False)
            ax.grid(False)  # Turn off gridlines
            ax.set_facecolor("none")  # Ensure no background color
            for spine in ax.spines.values():  # Remove all spines
                spine.set_visible(False)
        
        # Plot the biplot
        pca_model.biplot(
            n_feat=6,title='',
            c=cell_properties_df['colors'],
            s=120,
            color_arrow='black',
            fontsize=0,
            arrowdict={'fontsize': 0, 'color_text': 'white', 'weight': 'bold', 'alpha': 1, 'scale_factor' : arrow_scale_factor},
            ax=ax_main, legend=False, verbose=0
        )
        
        # Plot the marginal densities
        for treatment, group_data in cell_properties_df.groupby("stimulation"):
            # X-axis density
            sns.kdeplot(
                data=group_data,
                x="PC1",
                ax=ax_marg_x,
                color=palette[treatment],
                linewidth=2.5, fill=True
            )
            # Y-axis density
            sns.kdeplot(
                data=group_data,
                y="PC2",
                ax=ax_marg_y,
                color=palette[treatment],
                linewidth=2.5, fill=True
            )
        
        # Add labels and titles
        ax_main.set_xlabel(f"PC1 ({round(pca_results['variance_ratio'][0]*100,1)}%)", fontsize=38)
        ax_main.set_ylabel(f"PC2 ({round(pca_results['variance_ratio'][1]*100,1)}%)", fontsize=38)

        handles = [plt.Line2D([0], [0], color='white',markerfacecolor=color,marker='o', markersize=24, label=treatment) for treatment, color in palette.items()]
        ax_main.legend(
            handles=handles, 
            #title="Treatment", 
            fontsize=30, 
            #title_fontsize=14, 
            loc="upper left", 
            bbox_to_anchor=(0.99, 1),  # Position the legend outside the plot
            frameon=False, handletextpad=-0.2
        )
        
        text_pos={}
        coeff = pca_results['loadings'].iloc[[0, 1], :]
        for col in coeff.columns:
            max_axis = np.max(np.abs(pca_results['PC'].iloc[:, [0, 1]]).min(axis=1))
            max_arrow = np.abs(coeff).max().max() 
            scale = (np.max([1, np.round(max_axis / max_arrow, 3)])) * 2 + 1
            text_pos[col] = [(coeff[col].values[0] * scale), (coeff[col].values[1] * scale)]

        for col, pos in text_pos.items():
            ax_main.text(pos[0], pos[1], col, fontsize=32,weight='bold', color='black', ha='center', va='center',zorder=10 )

        
        return  ax_main, cell_properties_df


def pca_property(cell_properties_df: pd.DataFrame, color_by: str = 'cluster', cmap:str='plasma'):
    with plt.rc_context({"figure.dpi": 350}):
        fig, ax = plt.subplots(figsize=(2, 2))

        if color_by in ['frequency', 'width','rise_time','decay_time','amplitude']:
            values = np.log(cell_properties_df[color_by].values)
        else:
            values = cell_properties_df[color_by].values
        
        scatter = ax.scatter(
            cell_properties_df['PC1'],
            cell_properties_df['PC2'],
            c=values,
            cmap=cmap,
            edgecolors='black',
            linewidths=0.2,
            s=5
        )

        ax.set_xlabel('PC1', fontsize=12)
        ax.set_ylabel('PC2', fontsize=12)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        if color_by not in ['frequency', 'width','rise_time','decay_time','amplitude']:
            unique_vals = sorted(cell_properties_df[color_by].unique())
            handles = [
                plt.Line2D([0], [0], marker='o', color='w',
                           label=str(val), markerfacecolor=scatter.cmap(scatter.norm(val)),
                           markersize=5)
                for val in unique_vals
            ]
            ax.legend(handles=handles,  fontsize=6, markerscale=1.5,
                      loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)

        plt.title(color_by)

        return ax
    

def cluster_heatmap(cell_properties_df: pd.DataFrame, vmax: int = None, vmin: int = None, cbar:bool = False):
    mean_df = cell_properties_df[['frequency_scaled','width_scaled','rise_time_scaled','decay_time_scaled','amplitude_scaled','cluster']].groupby(['cluster']).mean()
    
    with plt.rc_context({"figure.dpi": 350}):
        fig, ax = plt.subplots()
        ax = sns.heatmap(mean_df,  cmap='plasma',vmin = vmin, vmax = vmax, cbar = cbar)
        ax.set_ylabel('Cluster', size=13)
        ax.set_xticklabels([col.rsplit('_', 1)[0] for col in mean_df.columns])
        plt.yticks(rotation=0, size=11)
        plt.xticks( size=13)

    return ax

def cluster_centroids( cluster_dict:dict,palette:dict,  imaging_interval:float = None, alpha: float = 1):


    with plt.rc_context({"figure.dpi": 350}):

        fig, ax = plt.subplots()  # create the figure once

        for cluster, centroid in cluster_dict.items():
            ax.plot(centroid, label=cluster, lw=3, alpha=alpha, color=palette[cluster])
        ax.legend(title='Cluster', frameon=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.title('Cluster centroids', fontsize=14)
        ax.set_ylabel(r"Ca$^{2+}$ signal ($\Delta F/F_0$)", fontsize=14)
        ax.set_xlabel(r"Time (s)", fontsize=14)

        # Set ticks
        num_ticks = 5
        rough_step = len(centroid) / (num_ticks - 1)
        step = int(np.round(rough_step / 10) * 10) # Round to nearest 10 for clean tick spacing
        xticks = np.arange(0, len(centroid) , step)
        tick_labels = (xticks * imaging_interval).astype(int)
        ax.set_xticks(xticks)
        ax.set_xticklabels(tick_labels)

    return ax


'''def cluster_barplot(cluster_percentages_df: pd.DataFrame, palette: dict = None):
    with plt.rc_context({"figure.dpi": 350}):
        fig, ax = plt.subplots()
        ax = sns.barplot(data=cluster_percentages_df, x='percentage', y='cluster', hue='stimulation', legend=None,palette=palette, ci=None, edgecolor='black', linewidth=0.5)
        sns.stripplot(data=cluster_percentages_df, x='percentage', y='cluster', hue='stimulation',  legend=None,ax=ax, color='black', dodge=True, size=2)
        
        #ax.set_ylabel('Cluster', size=13)
        #ax.set_xlabel('% of cells', size=13)
        
        plt.yticks(rotation=0)
        #plt.xticks( size=13)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    return ax'''