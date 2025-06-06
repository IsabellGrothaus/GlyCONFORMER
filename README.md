# GlyCONFORMER

**GlyCONFORMERS** is a Python package that assigns labels to glycan conformers, based on their torsion angle values, in order to differentiate various 3D-structures of one single glycan. It enables the automated assignment of the GlyCONFORMER string that was introduced in:

 **Grothaus et al. 2022, Exploration, Representation, and Rationalization of the Conformational Phase Space of N-Glycans, J. Chem. Inf. Model. 2022, 62, 20, 4992–5008** https://pubs.acs.org/doi/full/10.1021/acs.jcim.2c01049 for N-glycans. 

![](https://github.com/IsabellGrothaus/GlyCONFORMER/blob/v1.0.0-alpha/TUTORIAL/Conformer_string.png?raw=true)

Check the paper or tutorial for a detailed explanation of the GlyCONFORMER string generation. The workflow is not exclusivly designed for N-glycans, but can also be applyied to any other glycan type there is. 

## Installation

GlyCONFORMER is available via pip and docker:

## pip

```
pip install GlyCONFORMER
```
Stable performance was only tested and verified up to python version 3.10.
Pandas<=1.3.5 is required!

Recommendation - build via conda environment with following commands:

```
conda create -n name python=3.10 --no-default-packages
conda activate name
pip install GlyCONFORMER
``` 

## Docker Image

This image includes:
- GlyCONFORMER pre-installed via `pip install .`
- Jupyter Notebook launched on container start
- Required dependencies like `plumed`, `numpy`, `scikit-learn`, etc.

```bash
docker pull grothaus/glyconformer:latest

docker run -p 8888:8888 grothaus/glyconformer:latest
```

## Tutorial

The tutorial includes different N-glycan types and different complexity levels of how to obtain a GlyCONFORMER label string for custom glycan types and their recorded torsion angle values. The minimum example is given by the high-mannose type N-glycan M5, where all necessary information are read from the LIBRARY_GLYCANS folder by specifying the **glycantype = "Man5**":

```
M5 = Glyconformer(glycantype = "Man5")
```
Raw torsion angle values get converted to letters corresponding to their values and associate them to a certain minima of the free energy profile along that torsion angle. 

The raw input data can be accessed via:

```
M5.colvar
```

![](https://github.com/IsabellGrothaus/GlyCONFORMER/blob/v1.0.0-alpha/TUTORIAL/Input.png?raw=true)

and its translated variant via:

```
M5.binary
```

![](https://github.com/IsabellGrothaus/GlyCONFORMER/blob/v1.0.0-alpha/TUTORIAL/Output.png?raw=true)

Additionally, the occurance of each conformer string is counted and outputted via:

```
M5.count
```

![](https://github.com/IsabellGrothaus/GlyCONFORMER/blob/v1.0.0-alpha/TUTORIAL/Count.png?raw=true)

All this is performed in the background and final results can be outputted using various plotting functions, e.g.:

``` 
M5.distribution
```

Conformer labels are given on the x-axis and deviations from the most populated conformer indicated by explicit letters, where dots are used when no change in that torsion angle could be detected. 

![](https://github.com/IsabellGrothaus/GlyCONFORMER/blob/v1.0.0-alpha/TUTORIAL/M5_example/Conformer_distribution.png?raw=true)

For more elaborate examples and a detailed explanation of how free energy profiles of each torsion angle are classified, check out the **Tutorial_GlyCONFORMER.ipynb** notebook. 

## Documentation

See documentation https://glyconformer.readthedocs.io/en/latest/index.html
