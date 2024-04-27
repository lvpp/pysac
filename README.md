# Pairwise interacting surface activity coefficient, reference implementation in Python
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6962631.svg)](https://doi.org/10.5281/zenodo.6962631)
[![DOI:10.1016/j.fluid.2022.113611](http://img.shields.io/badge/DOI-10.1016/j.fluid.2022.113611-B31B1B.svg)](https://doi.org/10.1016/j.fluid.2022.113611)
[![DOI:10.1016/j.fluid.2024.114113](http://img.shields.io/badge/DOI-10.1016/j.fluid.2024.114113-B31B1B.svg)](https://doi.org/10.1016/j.fluid.2024.114113)

This repository contains a reference implementation of pairwise interacting surface
activity coefficient (SAC) equations.
There are two accompanying papers showing the derivation of all equations used:
 - [Ref. 1 - Beyond activity coefficients with pairwise interacting surface (COSMO-type) models](https://doi.org/10.1016/j.fluid.2022.113611).
 - [Ref. 2 - Unraveling Order and Entropy with Modern Quasi-Chemical Models](https://doi.org/10.1016/j.fluid.2024.114113).

All the relevant SAC equations are implemented in the [pysac.py module](./pysac.py)
favoring readability over efficiency.
The base code is not dependent on how the surface interaction energies are actually calculated and
then it can be adapted to work as a COSMO-RS/SAC, F-SAC, COSMO-SAC-Phi, etc. implementation.

## Interactive notebooks

In addition to the SAC equations implementation, there are many [jupyter](https://jupyter.org/)
interactive notebooks.

Notebooks related to Ref. 1:
 - [List all microstates for a given number of segments](./notebook/microstates.ipynb)
 - [Nonrandom factors for a five compund mixture](./notebook/five_compounds.ipynb)
 - [Excess properties for a size-symmetric nonideal mixture](./notebook/ue_symmetric_direct.ipynb)
 - [Excess properties for a size-symmetric nonideal mixture, adjusted UNIQUAC parameters](./notebook/ue_symmetric.ipynb)
- [Excess properties for a size-asymmetric nonideal mixture](./notebook/ue_asymmetric.ipynb)
- [Study on the interaction probabylities of a size-symmetric mixture](./notebook/prob_symmetric.ipynb)
- [Study on the interaction probabylities of a size-asymmetric mixture](./notebook/prob_asymetric.ipynb)
- [Study on the interaction probabylities of a nearly ideal mixture](./notebook/prob_nearly_ideal.ipynb)
- [Excess properties with temperature dependent interaction energies](./notebook/fsac_ft.ipynb)
- [Numerical check for excess energy with temperature dependent interaction energies](./notebook/fsac_u_check.ipynb)

Notebooks related to Ref. 2:
 - [Model polar/inert system](./notebook/egner_polar_inert.ipynb)
 - [Model hydrogen bonding/inert system](./notebook/egner_hb_inert.ipynb)
 - [Excess properties for the chloroform/methanol mixture](./notebook/chloroform_methanol.ipynb)
 - [Excess properties for the n-butanol/n-hexane mixture](./notebook/butanol_hexane.ipynb)
 - [Excess properties for the tetrahydrofuran/water mixture](./notebook/tetrahydrofuran_water.ipynb)
 
## Running on your machine

### Requirements

You will need a [jupyter](https://jupyter.org/) installation on your machine.

### Clone or download the repository

Clone this repository using git or [download the code](https://github.com/lvpp/pysac/archive/refs/heads/main.zip) and extract the zip contents.

### Running

Enter the root directory of the project and run:
```console
jupter-notebook
```

This should start the notebook app and open your web browser so that you can run and edit the notebooks.

Another option is simply to open the project directory with [Visual Studio Code](https://code.visualstudio.com/) and install the suggested Python extensions.

## Other related projects

The reader might be also interested in the following related projects:
 - [LVPP sigma profile database](https://github.com/lvpp/sigma)
 - [JCOSMO package](https://www.ufrgs.br/lvpp/2023/09/22/jcosmo-download)
 - [Benchmark COSMO-SAC implementation](https://github.com/usnistgov/COSMOSAC)
 - [openCOSMO-RS_cpp](https://github.com/TUHH-TVT/openCOSMO-RS_cpp)
 - [openCOSMO-RS_py](https://github.com/TUHH-TVT/openCOSMO-RS_py)
