Install conda (if not available):
https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions


Clone the repository:
git clone https://github.com/kirkebylab/LUMIN.git

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


