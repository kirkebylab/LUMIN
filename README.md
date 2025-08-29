# LUMIN (Live-cell User Module for Imaging and analysis of Neuronal activity)

Lumin is a software that integrates GUI-based pipeline configuration with automated calcium imaging data analysis, containing steps for deep-learning ROI segmentation, signal extraction, ΔF/F normalization, and response quantification.

This repository contains a guide on how to get started with LUMIN, as well as notebooks used to recreate figures for the manuscript with correct dimensions.


## Installation
Lumin has been tested on a MacBook Pro with M2 chip and Windows 11 Pro x64 with NVIDIA GPU. Separate installation guides are available for macOS and Windows operating systems. LUMIN requires a computer equipped with a GPU for efficient data processing. It is recommended to install LUMIN into a [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/getting-started.html) environment. This requires Miniconda or Anaconda to be installed.

If not installed, one can do so by following the instructions on [Miniconda installation page](https://www.anaconda.com/docs/getting-started/miniconda/install). 


#### macOS

1. Open Terminal.

2. Clone LUMIN GitHub repository:
```bash
git clone https://github.com/kirkebylab/LUMIN.git
```

3. Navigate to the cloned folder:
```bash
cd LUMIN
```

4. Create conda environment:
```bash
conda create -n lumin_env python=3.10 -c conda-forge -y
```

5. Activate conda environment:
```bash
conda activate lumin_env
```

6. Install LUMIN in the conda environment:
```bash
pip install -e .
```

7. Install Jupyter kernel (optional, enables custom downstream analysis or figure polishing of LUMIN-generated tabular output):
```bash
python -m ipykernel install --user --name=lumin_env
```

#### Windows
1. Open Anaconda Prompt.

2. Clone LUMIN GitHub repository:
```bash
git clone https://github.com/kirkebylab/LUMIN.git
```

3. Navigate to the cloned folder:
```bash
cd LUMIN
```

4. Create conda environment:
```bash
conda create -n lumin_env python=3.10 cudatoolkit=11.2 cudnn=8.1.0 -c conda-forge -y
```

5. Activate conda environment:
```bash
conda activate lumin_env
```

6. Install LUMIN in the conda environment:
```bash
pip install -e .
```

7. Install jupyter kernel (optional, enables custom downstream analysis or figure polishing of LUMIN-generated tabular output):
```bash
python -m ipykernel install --user --name=lumin_env
```

8. Enable developer mode in Windows settings. This is done to avoid an [error](https://github.com/stardist/stardist/issues/287) caused by StarDist python package.


#### Linux
Note, installation for Linux is not tested. It's anticipated that the GPU is detected if user follows these steps:

1. Follow installation guide for macOS. 

2. In active lumin_env run:
```
pip install tensorflow[and-cuda]
```
If this command returns non-empty list GPU is detected:
```
python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"
```



## Data loading
A small dataset is provided to test the pipeline. It's a subset of data used in the manuscript.

Testing data can be loaded from this dropbox [link](https://www.dropbox.com/scl/fo/z1gg916e09zk5gmb6yrqn/ALrZjemEdyoNE28Tss1efxs?rlkey=1ljkizlzhdwxpeoqggp40zqpm&dl=0).

Download it to `{local_path}/LUMIN` folder and unzip

To run LUMIN with proper data, use the provided `.csv` files as template to generate an input data file for LUMIN. The pipelines require the following columns to be present in the input file `plate_id`, `filename`, `filepath`, `biological_replicate`, and `stimulation`. User can provide any additional sample-associated metadata to the input data file. 

## Example analysis workflow


LUMIN is built around Napari and uses its data [layer](https://napari.org/dev/howtos/layers/image.html) widgets to display data during the parameter fine-tuning process. 

The calcium imaging data analysis happens through (1) Segmentation and signal extraction and (2) Single-cell data analysis pipelines, which can be configured using two custom-made Napari widgets. After the configuration fine-tuning process both pipelines process the data autonomously. To get started with LUMIN, we provide an example pipeline configuration for two different analysis setups: (1) spontaneously active neurons (Transient activity analysis) and (2) quiescent stimulus-evoked neurons (Baseline shift analysis).

The process of the pipeline execution can be followed from terminal window. The Napari window is unresponsible while the analysis is running.

When using the testing data, LUMIN should be executed from `{local_path}/LUMIN` because the test input file do not contain absolute data paths. Lumin is launched from a Terminal (macOS) or Anaconda Prompt (Windows) window.

### Transient activity analysis
####
1. Activate conda environment (if not active):
```bash
conda activate lumin_env
```

2. Launch napari window (if not open):
```bash
napari
```

3. Select Segmentation and signal extraction -pipeline (Plugins > LUMIN > Segmentation and signal extraction) and apply the following settings using the GUI input fields:
```
Input file: {Select ca_spontanoeus_input_data.csv from test_data folder}
Project directory: {Select project folder, for instance, generate Output/Spontaneous_project inside LUMIN folder}
ROI segmentation mode: Automated
Nuclear stain: None
Stain to segment: Cytoplasmic (Cellpose)
Model: cyto2
Diameter: 30
Cell probability threshold: 0.0
Flow threshold: 0.4
```

4. Press `Test settings on random image` -button to sample a random image from the input. Once the image appears in the canvas, apply the following post-processing settings and evaluate the outcome of the parameter configuration:

```
Cell area: 350 - 5600
Fluorescence intensity: 350 - 8000
```

The scale of post-processing settings is determined based on the maximum value of sampled images * 2 (by default). To increase the scale, you might need to sample another image. With the testing data, you might not be able to reach the values indicated above.

5. The user can play around with the segmentation and post-processing settings and continue sampling random images until satisfied with the results (this process will not save any output).

6. Press `Run` -button to execute the Segmentation and signal extraction pipeline. This will process all files indicated in the `Input file` and create a `Segmentation` folder in the specified `Project directory`.

Once the pipeline finishes, the user can move to the quantification step.

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

8. Press `Test settings on random image` -button to sample a random image from input. The user can explore the baseline estimation and spike detection using line plots, or play the calcium video to visualize detected spikes overlaid with the video. The user can adjust the peak detection settings and continue sampling random recordings until content with the results (this process will not save any output).

9. Press `Run` -button to execute the Single-cell data analysis pipeline. This will create a `quantification` folder in the specified `Project directory`



### Baseline shift analysis
1. Activate conda environment (if not active):
```bash
conda activate lumin_env
```

2. Launch napari window (if not open):
```bash
napari
```


3. Select Segmentation and signal extraction pipeline (Plugins > LUMIN > Segmentation and signal extraction) and apply the following settings through the GUI input fields:
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

4. Press `Test settings on random image` -button to sample a random image from the input. Once the image appears in the canvas, apply the following post-processing settings and evaluate the outcome of the configured settings:

```
Nuclear overlap: 0.7
Nuclear area: 30 - 3442
Cell area: 350 - 6554
Fluorescence intensity: 350 - 10470
```

The scale of post-processing settings is determined based on the maximum value of sampled images * 2 (by default). To increase the scale, you might need to sample another image. With the testing data, you might not be able to reach the values indicated above.

5. The user can play around with the settings and continue sampling random images until satisfied with the results (this process will not save any output).
6. Press `Run` -button to execute the Segmentation and signal extraction pipeline. This will process all files indicated in the `Input file` and create a `Segmentation` folder in the specified `Project directory`.
7. To further process and quantify the extracted signal, open Single-cell data analysis -pipeline (Plugins > LUMIN > Single-cell data analysis) and apply the following settings through GUI input fields:

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

8. Press `Test settings on random image` -button to sample a random image from input. The user can explore the normalization and activity classification settings using line and swarm plots. The user can adjust the settings and continue sampling random recordings until content with the analysis setup (this process will not save any output).

9. Press `Run` -button to execute the Single-cell data analysis pipeline. This will process all files indicated in the `Input file` and create a `Quantification` output folder specified by `Project directory`.







