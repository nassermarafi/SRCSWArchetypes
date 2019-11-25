# SRCSWArchetypes
 
This repository generates OpenSees Scripts for 69 special reinforced concrete shear wall archetypes ranging from 
4- to 40-stories. The archetypes were designed using either ASCE 7-10 or ASCE 7-16 for Seattle and 
are described in further detail in the following journal publications:

Marafi, N. A., Makdisi, A. J., Eberhard, M. O., and Berman, J. W. (2020) Impacts of M9 Cascadia 
Subduction Zone Earthquake and Seattle Basin on Performance of RC Core-Wall Buildings. Journal of Structural Engineering,
 [in-press](https://nassermarafi.github.io/papers/Marafietal2020-RCWallPerformanceM9.pdf)

Marafi, N. A., Makdisi, A. J., Berman, J. W., and Eberhard, M. O. (2020) Impacts of Design Strategies 
to Account for Effects of Sedimentary Basins on Reinforced Concrete Walls. Earthquake Spectra, 
[in-press](https://nassermarafi.github.io/papers/Marafietal2020-DesignStrategies.pdf)

The modelling methodology used in these wall archetypes are described in more detail in:

Marafi, N. A., Ahmed, K. A., Lowes, L. N., and Lehman, D. E. (2019) Sensitivity of Collapse Analysis 
for Reinforced Concrete Walls to Modelling Parameters. Journal of Structural Engineering, 
doi: [10.1061/(ASCE)ST.1943-541X.0002311](https://doi.org/10.1061/(ASCE)ST.1943-541X.0002311).

## To Get Started

To start using these archetypes, make sure to do the following first:

1. Download the latest OpenSees executable file and make sure it is in your path environment.

* Make sure uniaxial material models [Steel02Fatigue](https://github.com/OpenSees/OpenSees/blob/master/SRC/material/uniaxial/Steel02Fatigue.cpp) 
and [Concrete02IS](https://github.com/OpenSees/OpenSees/blob/master/SRC/material/uniaxial/Concrete02IS.cpp) 
are included in the OpenSees executable.
     
2. Install [OpenSeesAPI](https://github.com/nassermarafi/OpenSeesAPI).

3. Install [Jupyter Lab](https://github.com/jupyterlab/jupyterlab) (or Jupyter Notebook).

4. Open WallAnalysisResults.ipynb using Jupyter to get started.