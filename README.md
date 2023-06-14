# GlyCONFORMER

**GlyCONFORMERS** is a Python package that assigns labels to glycan conformers, based on their torsion angle values, in order to differentiate various 3D-structures of one single glycan. It enables the automated assignment of the GlyCONFORMER string that was introduced in:

 **Grothaus et al. 2022, Exploration, Representation, and Rationalization of the Conformational Phase Space of N-Glycans, J. Chem. Inf. Model. 2022, 62, 20, 4992–5008** https://pubs.acs.org/doi/full/10.1021/acs.jcim.2c01049 for N-glycans. 

![](https://github.com/IsabellGrothaus/GlyCONFORMER/blob/v1.0.0-alpha/TUTORIAL/Conformer_string.png?raw=true)

Check the paper or tutorial for a detailed explanation of the GlyCONFORMER string generation. The workflow is not exclusivly designed for N-glycans, but can also be applyied to any other glycan type there is. 

## Installation

To use GlyCONFORMER, first install it using pip:
```
pip install GlyCONFORMER
```
Stable performance was only tested and verified with python version 3.8.
Pandas<=1.3.5 is required!
 
## Tutorial

The tutorial juypter notebook should be run from within the GlyCONFORMER package folder or you have to change the path directing to the TUTORIAL folder.

The tutorial includes different N-glycan types and different complexity levels of how to obtain a GlyCONFORMER label string for custom glycan types and their recorded torsion angle values. The minimum example is given by the high-mannose type N-glycan M5, where only the file **M5_angles.dat** with torsion angle values is used as input:

![](https://github.com/IsabellGrothaus/GlyCONFORMER/blob/v1.0.0-alpha/TUTORIAL/Input.png?raw=true)

It is used by the glyconformer package, whereas remaining necessary information are read from the LIBRARY_GLYCANS folder by specifying the **glycantype = "M5**": 

```
conformer = glyconformer(inputfile = "TUTORIAL/M5_example/M5_angles.dat", glycantype = "M5")
```

When executing the run command:

```
binary, population = conformer.run()
```

a **binary** dataframe is produced, where the torsion angles have been converted to letters corresponding to their values and associate them to a certain minima of the free energy profile along that torsion angle. 

![](https://github.com/IsabellGrothaus/GlyCONFORMER/blob/v1.0.0-alpha/TUTORIAL/Output.png?raw=true)

Additionally, the occurance of each conformer string is counted and outputted to the **population** dataframe:

![](https://github.com/IsabellGrothaus/GlyCONFORMER/blob/v1.0.0-alpha/TUTORIAL/Count.png?raw=true)

The obtained information can also be used to plot a histogram, displaying the conformer distribution by:

``` 
conformer.plot()
```

Conformer labels are given on the x-axis and deviations from the most populated conformer indicated by explicit letters, where dots are used when no change in that torsion angle could be detected. 

![](https://github.com/IsabellGrothaus/GlyCONFORMER/blob/v1.0.0-alpha/TUTORIAL/M5_example/Conformer_distribution.png?raw=true)

For more elaborate examples and a detailed explanation of how free energy profiles of each torsion angle are classified, check out the **Tutorial_GlyCONFORMER.ipynb** notebook. 

## Documentation

See documentation https://glyconformer.readthedocs.io/en/latest/index.html
