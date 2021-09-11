# PytomGUI

PyTom is a toolbox developed for interpreting cryo electron tomography data. All steps from reconstruction, localization, alignment and classification are covered with standard and improved methods.

## Getting Started



### Prerequisites

PyTomGUI is designed to run on linux systems, but can also be installed on MacOSX. It requires the following software package to be installed:

```
# General packages 
- python (>= 3.7 )
- openmpi 
- fftw3
- gcc (version 5-7) 
- libxml2
- swig (>= 3.0.12)

# Python Packages
- numpy
- scipy
- boost
- lxml 
- mrcfile
- tqdm
- scikit-image
- matplotlib
- cupy (Optional for GPU acceleration)

# Optional Software Packages used by GUI
- PyQt5
- pyqtgraph
- motioncor2 ( >=1.2.1)
- imod (=4.10.25)

```

### Installing

Before you can install PyTom you need to have an account on github, and to the account you need to link a token for online identification. For more information click on the following link.

Furthermore, the software packages git needs to be install. Git can be installed by sudo apt install git or yum install git. After git has been installed, run the following lines:

```
git clone --recursive https://github.com/FridoF/PyTomPrivate.git
cd PyTomPrivate
bash installMiniconda.sh
```

Please remember the location where you decide to install conda (CONDA_INSTALL_DIR). 

```
conda env create -f pytom_env.yml
conda activate pytom_env
python3.8 setup.py install --prefix [CONDA_INSTALL_DIR]/envs/pytom_env
```

## Versioning

For the versions available, see the [tags on this repository]. 

## Authors

* **Gijs van der Schot** - *PyTomGUI* 
* **Thomas Hrabe**       - *PyTom* 
* **Yuxiang Chen**       - *PyTom*
* **Friedrich Forster**  - *PyTom* 

See also the list of [contributors] who participated in this project.

## License

Copyright (c) 2021

Utrecht University

http://pytom.sites.uu.nl

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The complete license can be obtained from 
http://www.gnu.org/licenses/gpl-2.0.html.
