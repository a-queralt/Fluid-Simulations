# Basic Navier Stokes simulator
## Description
C++ that simulates basic computational fluid dynamics (CFD) problems.
Boundary conditions are hard-coded, so for wach new problem the cod emust be modified and recompiled. Results will be saved in .txt files on the same folder as the .exe
## Post process
Post process of the results can be done with many different tools: MatLab or Python libraries or GnuPlot for isntance. 
A post processing tutorial with [paraView](https://www.paraview.org/) tutorial will be available in the wiki shortly.
## Future improvements
1. Addition of UDS, SUDS and QUICK convective schemes.
2. Use of Hypervolic concentrated meshes
3. Coupling of the heat equation
4. BOundary condition introduction via .txt file
