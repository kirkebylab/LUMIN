# LUMIN (Live-cell User Module for Imaging and analysis of Neuronal activity)

Lumin is a ....


## Documentation
To be updated ... 


## Installation
Lumin has been tested on MacBook Pro with M2 chip and Windows 11 Pro x64 with NVIDIA GPU. Separate installation guides should be followed for macOS and Windows operating systems. It is recommended to install LUMIN into a conda environment hosting all software dependencies of LUMIN.

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
conda create -n lumin_env python -c conda-forge -y
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

## Example analysis workflow
1. Launch napari:
```bash
napari
```

2. Select desiread pipeline from Napari plugins:
   - Plugins > LUMIN > Segmentation and signal extraction
   - Plugins > LUMIN > Single-cell data analysis
  
Settings indicated below are the ones used during in the manuscript, and provided here as examples.

### Transient activity analysis
####
1. Select Segmentation and signal extraction -pipeline and apply following settings thorugh GUI input fields:
```
Input file: /Volumes/T9/Ca_data/exp_pharmacology/input_data.csv
Project directory: /Volumes/T9/Ca_data/exp_pharmacology/Analysis_new_env
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
Clusters:6
```

7. Press `Test settings on random image` -button to sample random image from input. User can explore the baseline estimation and spike detection using line plots, or play the calcium video to visualize detected spikes overlaid with the video. User can adjust the peak detection settings and continue sampling random recordings until content with the results.

8. Press `Run` -button to execute the Single-cell data analysis -pipeline. This will process all files indicated in the `Input file` and create a `Quantification` output folder specificed by `Project directory`.



### Baseline shift analysis



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









