## Numerical modeling of resistive switching in RRAM device
__Comprehensive numerical modeling of filamentary RRAM devices including voltage ramp-rate and cycle-to-cycle variations__

The set of files are the product of my PhD research. These matlab files need COMSOL and livelink to MATLAB to operate. 
Executing them will reproduce RRAM device operation including SET and RESET process. 

The three main files, named Complete_I_V, generates average device characteristics, ramp-date dependency, and cycle-to-cycle variation. 
Please read the Catalog for a detailed description including how to execute, found under Section 3. Implementation.  

Note: One needs LiveLink to MATLAB to run these files. LiveLink to Matlab is a COMSOL add-on product that enables MATLAB to access COMSOL Solver.

LiveLink to MATLAB comes with the COMSOL Multiphysics with MATALB executable. Clicking this will automatically open a MATLAB linked to COMSOL.

[Catalog](./CATALOG_git.pdf)

To Do: Replace repeated minimization codes (Coarse minimization + Brents Minimization) with a function

## Related Work
### Theory of Resistive Switching
1. [V. G. Karpov, D. Niraula, I. V. Karpov, and R. Kotlyar, “Thermodynamics of Phase Transitions and Bipolar Filamentary Switching in Resistive Random-Access Memory,” Phys. Rev. Appl., 8, 024028 (2017).](https://doi.org/10.1103/PhysRevApplied.8.024028)

2. [V. Karpov, D. Niraula, and I. Karpov, “Thermodynamic analysis of conductive filament,” Appl. Phys. Letts. 109, 093501 (2016).](https://doi.org/10.1063/1.4962136)

3. [V. G. Karpov and D. Niraula, “Resistive Switching in nanostructures,” Scientific Reports, 8, 12212 (2018).](https://doi.org/10.1038/s41598-018-30700-6)
   
### Numerical Modeling
4. [Dipesh Niraula and Victor Karpov, “Comprehensive numerical modeling of filamentary RRAM devices including voltage ramp-rate and cycle-to-cycle variations,” J. Appl. Phys., 124, 174502 (2018).](https://doi.org/10.1063/1.5042789)

5. [D. Niraula and V. Karpov, “Comprehensive Numerical Modeling of Filamentary RRAM device,” Proceeding of the 2018 COMSOL conference, Boston, MA, Oct 3-5, 2018](https://www.comsol.com/paper/download/566882/niraula_paper.pdf)

6. [D. Niraula and V. Karpov, “Numerical Modeling of Resistive Switching in a RRAM device,” Proceeding of the 2017 COMSOL conference, Boston, MA, Oct 2017.](https://www.comsol.com/paper/download/437382/niraula_paper.pdf) [Best paper award]

### Quantum effects due to extremely small dimension
7. [D. Niraula, C. R. Grice, and V. G. Karpov, “Dimensional quantization effects in the thermodynamics of conductive filaments,” Nanotechnology, 29, 265202 (2018).](https://doi.org/10.1088/1361-6528/aabdcb)

### OFF State Conduction
8. [V. G. Karpov and D. Niraula, “OFF state conduction in filamentary RRAM,” IEEE El. Dev. Letts., 40, 550 (2019).](https://doi.org/10.1109/LED.2019.2899258)

### Heat Transfer Analysis
9. [Dipesh Niraula and Victor G. Karpov, “Heat transfer in filamentary RRAM devices,” IEEE Trans. El. Dev., 64, 4160 (2017).](https://doi.org/10.1109/TED.2017.2741782)

### Noise Analysis
10.	[V.G. Karpov and D. Niraula, “Log-Normal Statistics in Filamentary RRAM Devices and Related Systems,” IEEE El. Dev. Letts. 38, 1240 (2017).](https://doi.org/10.1109/LED.2017.2734961)

### PhD Desertation
11.	[Dipesh Niraula, "Physics and applications of conductive filaments in electronic structures: from metal whiskers to solid state memory." Doctoral dissertation, University of Toledo, 2019.](http://rave.ohiolink.edu/etdc/view?acc_num=toledo1561471348406944)
