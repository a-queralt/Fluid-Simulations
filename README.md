# Basic Navier Stokes simulator
## Description
C++ that simulates basic computational fluid dynamics (CFD) problems. See Wiki for examples.
Boundary conditions are hard-coded, so for wach new problem the cod emust be modified and recompiled. Results will be saved in .txt files on the same folder as the .exe
## Post process
Post process of the results can be done with many different tools: MatLab or Python libraries or GnuPlot for isntance. 
In the Wiki a post processing with [paraView](https://www.paraview.org/) tutorial is available.
## Future improvements
1. Addition of UDS, SUDS and QUICK convective schemes.
2. Use of Hypervolic concentrated meshes
3. Coupling of the heat equation
4. BOundary condition introduction via .txt file
