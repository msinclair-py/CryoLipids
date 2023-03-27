# CryoLipids
A workflow for modeling incomplete lipids from cryo-em experiments.

## Installation Guide
To install CryoLipids simply clone this repo or download the code manually. In
order to run CryoLipids I have provided a yaml file for generating a conda environment.
It is strongly recommended to use an environment unless you are comfortable with managing
multiple packages with conflicting dependencies yourself. For the uninitiated you can
do so with the following command:
`conda env create --file environment.yaml`

To enter the environment after installing you can use the following:
`conda activate cryo`
If you get package missing errors when executing the code you are likely still in the
base environment or a different env altogether.

If you want to build this environment yourself you can follow the below workflow. Note that
because scikit-tda depends on pip to be installed it MUST be installed last due to conda's
inability to successfully manage any additional dependencies installed via pip (e.g. sklearn,
numpy, etc). sklearn specifically will break the code if you install mdanalysis after
scikit-tda. Additionally, in the future python 3.11 will presumably be in the main
anaconda channel but due to mdanalysis being in conda-forge you will likely still need
to include conda-forge as a channel for anaconda.

```conda create -c conda-forge python=3.11 -n cryo
conda activate cryo
conda install mdanalysis
pip install scikit-tda```

## Running CryoLipids
Simply running the main script `cryolipids.py` will perform all necessary operations. To
clarify the required positional and keyword arguments you can enter into your terminal
the following command:
`python cryolipids.py -h`
which will provide a succinct description of all user inputs.
