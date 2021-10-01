A Quadratic Programming Approach to Microwave Imaging
=====================================================

[André Batista](mailto:andre-costa@ufmg.br) and Lucas Batista
Operations Research and Complex Systems Laboratory - ORCS Lab  
Universidade Federal de Minas Gerais  
Belo Horizonte, Brazil

Ricardo Adriano
Electrical Engineering Deparment
Universidade Federal de Minas Gerais  
Belo Horizonte, Brazil

***

This repository contains the implementation of a Quadratic Programming (QP) approach for solving a 2D formulation of a Microwave Imaging (MWI) problem. The developed approach is embedded into the Born Iterative Method as the linear regularizer. Our goal was to implement such approach and to evaluated it with case studies based in geometric figures and breast phantoms, exploring different scenarios.

This folder includes the following directories:

* *main functions*: here is the implementation of the main functions of the project. So it contains not only the algorithm itself, but a bunch of *.m* files which implement common computations such as Green's function calculation, Finite-Difference Time-Domain Method, among other. The names of the files are self-inductive.
* *basic_experiments*: here is the scripts which implements the experiments runned with traditional geometries (triangles, squares, etc) exploring different scenarios such as high contrast ones, noisy ones, multiple materials, among others. You may find not only the scripts but the results stored in the *.mat* files and plotted in figures.
* *breast_phantoms*: here is the scripts that runs the experiments with breast phantoms. The phantoms are divided according to class and the slices removed from the 3D models. The models are the ones provided by the [University of Wisconsin Cross-Disciplinary Electromagnetics Laboratory](https://uwcem.ece.wisc.edu/phantomRepository.html). The results were not upload here due GitHub restrictions (they are +100MB)

Besides of having **MATLAB**, you will need [**Gurobi Optimizer**](https://www.gurobi.com/documentation/9.1/refman/matlab_api_overview.html) as well to run the inverse solver.

This work was published in IEEE Transactions on Antennas and Propagation and you will find the paper clicking [here](https://ieeexplore.ieee.org/abstract/document/9362200). If you that these codes will contribute to your project, please, **cite us**:

A. C. Batista, L. S. Batista and R. Adriano, "A Quadratic Programming Approach for Microwave Imaging," in *IEEE Transactions on Antennas and Propagation*, vol. 69, no. 8, pp. 4923-4934, Aug. 2021, doi: 10.1109/TAP.2021.3060092.

If you have any question, let me know by e-mail!

Have fun!
André Batista