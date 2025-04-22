# Para-Equilibrium Thermodynamic Modeling

## Description
This open-source program calculates para-equilibrium states for thermodynamic data sets, advancing the thermodynamic modeling of multicomponent phases in hydrogen-para-equilibrium. It was developed to support the research presented in the paper "Advancing the thermodynamic modeling of multicomponent phases in hydrogen-para-equilibrium."

## Installation
To install and set up this program, please refer to the detailed instructions provided in the supplementary materials of the publication. You can access the paper and supplementary materials via the following DOI: [10.1016/j.actamat.2024.120529](https://doi.org/10.1016/j.actamat.2024.120529).

## Usage
The program consists of five main Python files:

1. **Calculate_PCI.py**: Calculates the plateau pressures and pressure values corresponding to mole fractions in single-phase regions.
2. **CONSTANTS.py**: Contains constants used for calculating the common tangent.
3. **Equations.py**: Minimizes the Gibbs energy curve of the solid phase by finding the global minimum for discrete compositions of xH.
4. **database.py**: Contains the thermodynamic variables to describe different phases.
5. **Execute.py**: Provides example usage of the files to generate PCI diagrams.

## Running the Code
To run the program, the file `Execute.py` gives example usage. The script will import the necessary files and prompt the user to input the following variables:

- **solids**: A list of n lists representing the composition of each sublattice.
- **multiplicity**: A list of integers defining the multiplicity of the sublattices.
- **site_fractions**: A list of n lists representing the site fractions of the metal elements in each sublattice.

## Outputs
Running Execute.py will generate the PCI diagrams as described in Sections S1-S3 of the supplementary material.

## Citation
If you use this work, please cite it as:

@article{HANNAPPEL2025120529,
title = {Advancing the thermodynamic modeling of multicomponent phases in hydrogen-para-equilibrium},
journal = {Acta Materialia},
volume = {284},
pages = {120529},
year = {2025},
issn = {1359-6454},
doi = {https://doi.org/10.1016/j.actamat.2024.120529},
url = {https://www.sciencedirect.com/science/article/pii/S1359645424008784},
author = {Peter Hannappel and Felix Heubner and Mateusz Balcerzak and Thomas Weißgärber}
}
