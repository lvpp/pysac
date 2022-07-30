# Pairwise interacting surface activity coefficient, reference implementation in Python

This repository contains a reference implementation of pairwise interacting surface
activity coefficient (SAC) equations.
There is an accompanying paper showing the derivation of all equations used (link to be provided after publication).

All the relevant SAC equations are implemented in the [pysac.py module](./pysac.py)
favoring readability over efficiency.
The base code is not dependent on how the surface interaction energies are actually calculated and
then it can be adapted to work as a COSMO-RS/SAC, F-SAC, COSMO-SAC-Phi, etc. implementation.

## Interactive notebooks

In addition to the SAC equations implementation, there are many [jupyter](https://jupyter.org/)
interactive notebooks:
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
 - [Benchmark COSMO-SAC implementation](https://github.com/usnistgov/COSMOSAC)
 - [openCOSMO-RS_cpp](https://github.com/TUHH-TVT/openCOSMO-RS_cpp)
 - [openCOSMO-RS_py](https://github.com/TUHH-TVT/openCOSMO-RS_py)