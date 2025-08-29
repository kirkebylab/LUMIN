# LUMIN (Live-cell User Module for Imaging and analysis of Neuronal activity)

Lumin is a ....

This repository contains ....

## Documentation
To be updated ... 


## Installation
Lumin has been tested on a MacBook Pro with M2 chip and Windows 11 Pro x64 with NVIDIA GPU. Separate installation guides are available for macOS and Windows operating systems. It is recommended to install LUMIN into a conda environment. LUMIN requires a computer equipped with a GPU for efficient data processing. 

If Miniconda (or Anaconda) is not installed on the system, one can do so by following the instructions on [Miniconda installation page](https://www.anaconda.com/docs/getting-started/miniconda/install). 


#### macOS
1. Clone LUMIN GitHub repository:
```bash
git clone https://github.com/kirkebylab/LUMIN.git
```

3. Navigate to the cloned folder:
```bash
cd LUMIN
```

5. Create conda environment:
```bash
conda create -n lumin_env python=3.10 -c conda-forge -y
```

6. Activate conda environment:
```bash
conda activate lumin_env
```

7. Install LUMIN in the conda environment:
```bash
pip install -e .
```

8. Install Jupyter kernel (optional, enables downstream analysis or figure polishing of LUMIN-generated tabular output):
```bash
python -m ipykernel install --user --name=lumin_env
```

#### Windows
1. Clone LUMIN GitHub repository:
```bash
git clone https://github.com/kirkebylab/LUMIN.git
```

3. Navigate to the cloned folder:
```bash
cd LUMIN
```

5. Create conda environment:
```bash
conda create -n lumin_env python=3.10 cudatoolkit=11.2 cudnn=8.1.0 -c conda-forge
```

6. Activate conda environment:
```bash
conda activate lumin_env
```

7. Install LUMIN in the conda environment:
```bash
pip install -e .
```

8. Install jupyter kernel (optional, enables downstream analysis or figure polishing of LUMIN-generated tabular output):
```bash
python -m ipykernel install --user --name=ca_env
```

9. Enable developer mode in windows settings. This is done to avoid [error](https://github.com/stardist/stardist/issues/287) caused by StarDist python package. 



## Data loading
We provide a small dataset to test the pipeline. It's a subset of data used in the manuscript.

Testing data can be loaded from this dropbox [link](https://www.dropbox.com/scl/fo/z1gg916e09zk5gmb6yrqn/ALrZjemEdyoNE28Tss1efxs?rlkey=1ljkizlzhdwxpeoqggp40zqpm&dl=0).

Download it to `{local_path}/LUMIN` folder and unzip

To run LUMIN with proper data, user can use the provided .csv files as templates to generate input input data file for LUMIN. The pipelines requires the following columns to be present in the input `plate_id`, `filename`, `filepath`, `biological_replicate` and `stimulation`. User can provide any additional sample associated metadata to the input data file (Note `marker` column is reserved for the co-expression and shouldn't be used in input file). 

## Example analysis workflow


LUMIN is built around Napari and uses its data [layers](https://napari.org/dev/howtos/layers/image.html) widgets to display data during parameter fine-tuning process. The calcium imaging data analysis happens through Segmentation and signal extraction and Single-cell data analysis pipelines, which can be configured using custom made Napari widgets. To get started with LUMIN here we provide an example pipeline configuration for two different analysis setups: (1) spontanoeusly active neurons (Transient activity analysis) and (2) quiescent stimulus evoked neurons (Baseline shift analysis).

### Transient activity analysis
####
1. Launch napari:
```bash
napari
```

1. Select Segmentation and signal extraction -pipeline (Plugins > LUMIN > Segmentation and signal extraction) and apply the following settings through the GUI input fields:
```
Input file: {Select ca_spontanoeus_input_data.csv from test_data folder}
Project directory: {Select project folder, for instance, generate Output/Spontaneous_project to LUMIN folder}
ROI segmentation mode: Automated
Nuclear stain: None
Stain to segment: Cytoplasmic (Cellpose)
Model: cyto2
Diameter: 30
Cell probability threshold: 0.0
Flow threshold: 0.4
```

3. Press `Test settings on random image` -button to sample a random image from the input. Once the image appears in the canvas, apply the following post-processing settings and evaluate the outcome of the configured settings:

```
Cell area: 350 - 5600
Calcium intensity: 350 - 8000
```

The scale of post-processing settings is determined based on the maximum value of sampled images * 2 (by default). To increase the scale, you might need to sample another image.

4. The user can play around with the segmentation and post-processing settings and continue sampling random images until satisfied with the results. (This process will not save any output.)
5. Press `Run` -button to execute the Segmentation and signal extraction pipeline. This will process all files indicated in the `Input file` and create a `Segmentation` folder in the specified `Project directory`

The pipeline process can be followed from the terminal. Once the pipeline finishes, the user can move to the quantification step.

7. To further process and quantify the extracted signal, open Single-cell data analysis -pipeline (Plugins > LUMIN > Single-cell data analysis) and apply the following settings through GUI input fields:

```
Project directory: {Select same project directory ({local_path}/Spontanoeus_project) than in Segmentation and signal extraction -pipeline}
Analysis mode: Spontaneous activity
Control condition: Control
Normalization: Sliding window
Sliding window size: 75
Percentile threshold: 25
Prominence threshold: 0.2
Amplitude width ratio: 0.003
Imaging interval: 0.2
KCl stimulation frame: -1
Number of clusters: 6
```

7. Press `Test settings on random image` -button to sample a random image from input. The user can explore the baseline estimation and spike detection using line plots, or play the calcium video to visualize detected spikes overlaid with the video. The user can adjust the peak detection settings and continue sampling random recordings until content with the results.

8. Press `Run` -button to execute the Single-cell data analysis pipeline. This will process all files indicated in the `Input file` and create a `Quantification` output folder specified in `Project directory`.



### Baseline shift analysis
1. Launch napari (if not open):
```bash
napari
```


1. Select Segmentation and signal extraction -pipeline (Plugins > LUMIN > Segmentation and signal extraction) and apply the following settings through the GUI input fields:
```
Input file: {Select ca_evoked_input_data.csv from test_data folder}
Project directory: {Select project folder, for instance, generate Output/Evoked_project to LUMIN folder}
ROI segmentation mode: Automated
Nuclear stain: First frame
Stain to segment: Nuclear (StarDist) and cytoplasmic (Cellpose)
Probability/Score Threshold: 0.6
Overlap threshold: 0.3
Model: cyto2
Diameter: 35
Cell probability threshold: 0.0
Flow threshold: 3.0
```

3. Press `Test settings on random image` -button to sample a random image from the input. Once the image appears in the canvas, apply the following post-processing settings and evaluate the outcome of the configured settings:

```
Filtering - Nuclear overlap: 0.7
Filtering - Nuclear area: 30 - 3442
Filtering - Cell area: 350 - 6554
Filtering - Calcium intensity: 350 - 10470
```

The scale of post-processing settings is determined based on the maximum value of sampled images * 2 (by default). To increase the scale, you might need to sample another image. With the testing data, you might not be able to reach the values used in the manuscript.

4. The user can play around with the settings and continue sampling random images until satisfied with the results. (This process will not save any output.)
5. Press `Run` -button to execute the Segmentation and signal extraction pipeline. This will process all files indicated in the `Input file` and create a `Segmentation` folder in the specified `Project directory`
6. To further process and quantify the extracted signal, open Single-cell data analysis -pipeline (Plugins > LUMIN > Single-cell data analysis) and apply the following settings through GUI input fields:

```
Project directory: {Select same project directory ({local_path}/Evoked_project) than in Segmentation and signal extraction -pipeline}
Analysis mode: Compound-evoked activity
Activity type: Baseline change
Control condition: Control
Normalization: Pre-stimulus window
Stimulation frame: 20
Analysis window start: 50
Analysis window end: 70
Standard deviation threshold: 3
Imaging interval: 0.5
KCl stimulation frame: 100
Number of clusters: 3
```

7. Press `Test settings on random image` -button to sample a random image from input. The user can explore the normalization and analysis settings using line plots and activity classification using swarm plots. The user can adjust the settings and continue sampling random recordings until content with the analysis setup.

8. Press `Run` -button to execute the Single-cell data analysis pipeline. This will process all files indicated in the `Input file` and create a `Quantification` output folder specified by `Project directory`.




