# LUMIN (Live-cell User Module for Imaging and analysis of Neuronal activity)

Lumin is a ....

This repository contains ....

## Documentation
To be updated ... 


## Installation
Lumin has been tested on MacBook Pro with M2 chip and Windows 11 Pro x64 with NVIDIA GPU. Separate installation guides should be followed for macOS and Windows operating systems. It is recommended to install LUMIN into a conda environment hosting all software dependencies of LUMIN. LUMIN requires computer equipped with GPU for efficient data processing. 

If Miniconda (or Anaconda) is not installed to the system one can do so by following the instructions on [Miniconda installation page](https://www.anaconda.com/docs/getting-started/miniconda/install). 


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

7. Install LUMIN to the conda environment:
```bash
pip install -e .
```
8. Install jupyter kernel (optional, enables downstream analysis or figure polishing of LUMIN generated tabular output):
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
``````

6. Activate conda environment:
```bash
conda activate lumin_env
```

7. Install LUMIN to the conda environment:
```bash
pip install -e .
```
8. Install jupyter kernel (optional, enables downstream analysis or figure polishing of LUMIN generated tabular output):
```bash
python -m ipykernel install --user --name=ca_env
```

https://napari.org/dev/howtos/layers/image.html

## Data loading
We provide small dataset to test the pipeline. It's a subset of data used in the manuscript.

Testing data can be loaded from https://www.dropbox.com/scl/fo/z1gg916e09zk5gmb6yrqn/ALrZjemEdyoNE28Tss1efxs?rlkey=1ljkizlzhdwxpeoqggp40zqpm&dl=0

Download it to `{local_path}/LUMIN` folder and unzip

## Example analysis workflow


Settings indicated below are the ones used during in the manuscript, and provided here as examples.

### Transient activity analysis
####
1. Launch napari:
```bash
napari
```

1. Select Segmentation and signal extraction -pipeline (Plugins > LUMIN > Segmentation and signal extraction) and apply following settings thorugh GUI input fields:
```
Input file: {Select ca_spontanoeus_input_data.csv from test_data folder}
Project directory: {Select project folder for instance generate Output/Spontaneous_project to LUMIN folder}
ROI segmentation mode: Automated
Nuclear stain: None
Stain to segment: Cytoplasmic (Cellpose)
Model: cyto2
Diameter: 30
Cell probability threshold: 0.0
Flow threshold: 0.4
```

3. Press `Test settings on random image` -button to sample random image from input. Once the image appears in the canvas apply following post-processing settings and evaluate impact of configured settings:

```
Cell area: 350 - 5600
Calcium intensity: 350 - 8000
```

The scale of post-processing settings are determined based on maximum value of sampled images * 2 (by default). To increase the scale you might need to sample another image.

4. User can play around with the segmentation and post-processing setting and continue sampling random image until satisfied with the results. (This process will not save any output)
5. Press `Run` -button to execute the Segmentation and signal extraction -pipeline. This will process all files indicated in the `Input file` and create a `Segmentation` folder to specificed `Project directory`

The pipeline process can be followed from terminal. Once the pipeline finnishes user can move to the quantification step.

7. To furhter process and quantify the extracted signal open Single-cell data analysis -pipeline (Plugins > LUMIN > Single-cell data analysis) and apply following settings thorugh GUI input fields:

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
NUmber of clusters: 6
```

7. Press `Test settings on random image` -button to sample random image from input. User can explore the baseline estimation and spike detection using line plots, or play the calcium video to visualize detected spikes overlaid with the video. User can adjust the peak detection settings and continue sampling random recordings until content with the results.

8. Press `Run` -button to execute the Single-cell data analysis -pipeline. This will process all files indicated in the `Input file` and create a `Quantification` output folder specificed by `Project directory`.



### Baseline shift analysis
1. Launch napari (if not open):
```bash
napari
```


1. Select Segmentation and signal extraction -pipeline (Plugins > LUMIN > Segmentation and signal extraction) and apply following settings thorugh GUI input fields:
```
Input file: {Select ca_evoked_input_data.csv from test_data folder}
Project directory: {Select project folder for instance generate Output/Evoked_project to LUMIN folder}
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

3. Press `Test settings on random image` -button to sample random image from input. Once the image appears in the canvas apply following post-processing settings and evaluate impact of configured settings:

```
Filtering - Nuclear overlap: 0.7
Filtering - Nuclear area: 30 - 3442
Filtering - Cell area: 350 - 6554
Filtering - Calcium intensity: 350 - 10470
```

The scale of post-processing settings are determined based on maximum value of sampled images * 2 (by default). To increase the scale you might need to sample another image. With the testing data you might not be able to reach the values used in the manuscript.

4. User can play around with the setting and continue sampling random image until satisfied with the results. (This process will not save any output)
5. Press `Run` -button to execute the Segmentation and signal extraction -pipeline. This will process all files indicated in the `Input file` and create a `Segmentation` folder to specificed `Project directory`
6. To furhter process and quantify the extracted signal open Single-cell data analysis -pipeline (Plugins > LUMIN > Single-cell data analysis) and apply following settings thorugh GUI input fields:

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

7. Press `Test settings on random image` -button to sample random image from input. User can explore the normalization and analysis settings using line plots and activity classification using swarm plots. User can adjust the settings and continue sampling random recordings until content with the analysis setup.

8. Press `Run` -button to execute the Single-cell data analysis -pipeline. This will process all files indicated in the `Input file` and create a `Quantification` output folder specificed by `Project directory`.





























.csv from test_data folder}
Project directory: {Select project folder for instance generate Spontanoeus_project to desired location}
ROI selection mode: Automated
Nuclear stain: None
Stain to segment: Cytoplasmic (Cellpose)
Model: cyto2
Diameter: 30
Cell probability threshold: 0.0
Flow threshold: 0.4
```

3. Press `Test settings on random image` -button to sample random image from input. Once the image appears in the canvas apply following post-processing settings and evaluate impact of configured settings:

```
Cell area: 350 - 5600
Calcium intensity: 350 - 8000
```

The scale of post-processing settings are determined based on maximum value of sampled images * 2 (by default). To increase the scale you might need to sample another image.

4. User can play around with the setting and continue sampling random image until satisfied with the results. (This process will not save any output)
5. Press `Run` -button to execute the Segmentation and signal extraction -pipeline. This will process all files indicated in the `Input file` and create a `Segmentation` folder to specificed `Project directory`
6. To furhter process and quantify the extracted signal open Single-cell data analysis -pipeline (Plugins > LUMIN > Single-cell data analysis) and apply following settings thorugh GUI input fields:

```Project directory: /Volumes/T9/Ca_data/exp_pharmacology/Analysis_new_env
Analysis mode: Spontaneous activity
Control condition: Control
Normalization mode: Sliding window
Sliding window size: 75
Percentile threshold: 25
Prominence threshold: 0.2
Amplitude width ratio: 0.003
Imaging interval: 0.2
KCl stimulation frame: -1
NUmber of clusters: 6
```

7. Press `Test settings on random image` -button to sample random image from input. User can explore the baseline estimation and spike detection using line plots, or play the calcium video to visualize detected spikes overlaid with the video. User can adjust the peak detection settings and continue sampling random recordings until content with the results.

8. Press `Run` -button to execute the Single-cell data analysis -pipeline. This will process all files indicated in the `Input file` and create a `Quantification` output folder specificed by `Project directory`.



Install conda (if not available):
https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions


Clone the repository:
`git clone https://github.com/kirkebylab/LUMIN.git`

Navigate 
cd LUMIN



Installation:
conda create -n lumin_env python -c conda-forge

Activate:
conda activate ca_env

Install plugin:
pip install -e .

Jupyter kernel (optional):
python -m ipykernel install --user --name=ca_env

Launching napari:
napari


new project dir for each segmentation



conda create -n lumin_env python=3.10 -c conda-forge

conda create -n lumin_env python=3.10 cudatoolkit=11.2 cudnn=8.1.0 -c conda-forge


tested with x64 windows 11 pro














