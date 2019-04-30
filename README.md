# Description
Data Processing and Plotting for Electrochemical Mesasurements
This tools processes electrochemical data exported using NOVA software. The tool perform the following analyses:
- Cyclic and Linear Sweep Voltomograms
- KL Plot inlcuding reaction number determination
- Tafal Plot
- Electrochemical Surface Area (ECSA) Determination

# Requirements
* Python 2.7 or higher
* matplotlib
* Scipy
* Pandas
* Numpy


# To Use
Use the input file in the Example folder as a template for an input file then simply
```bash
$ elechemcat.py input.inp
```
Note: it is recommened to define the location of the elechemcat.py file in the global environmental variables file so the function can be called from any location


# Examples
Cyclic Voltomogram
![CV](Example/plots/CV.png)

KL Plot
![KLR](Example/plots/KL-SSvsLSV.png)

Electrochemical Surface Area (ECSA) Determination
![ECSA](Example/plots/Sample1-ECSA.png)
