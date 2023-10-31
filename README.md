# kappa_SMpeptide_deisgn
Necessary code to reproduce work found in a soon to be submitted manuscript.

![alt text](https://github.com/kdeibler/kappa_SMpeptide_deisgn/blob/main/figure1.png)

# General requirements
## Rosetta and Pyrosetta Software
Rosetta software is free for academic use. You can get the license from the following link:     
https://www.rosettacommons.org/software/license-and-download

You may download and build pyrosetta after obtaining academic license by following protocols in this page:
http://www.pyrosetta.org/dow
 
To install and build, follow the protocols in the link below: https://www.rosettacommons.org/demos/latest/tutorials/install_build/install_build
Please note that with 1-4 cores, it will take *hours* to fully build Rosetta on a desktop.

The code provided in this paper was run on a 2018 Rosetta version. If you want to run on versions after that, you need to update resfile and .xml usage to use `PackerPalette`.

## Other packages and softwares
1. Python (version 3.6 works for most of the code. Code in method 1 is older and was tested on python 2.7)
2. Numpy, Pandas and Seaborn
3. Open Babel (http://openbabel.org/wiki/Main_Page)
4. Jupyter Notebook (only required to run the demo for Docking analysis)

# Rosetta Design Protocols
This folder contains subdirectories with code and detailed isntruction on how to run each step design methods. Below is a brief description of what each folder contains

### src
Contains main scripts for backbone generation and design

### inputs
Includes all the files that are used in more than one folder:
1. resfiles
2. params file for CVV and CYY
3. composition files

### submission files
slurm based 

## Output
Because the processes of design and packing are stochastic, the output files can vary from run to run. Each step as an output of either a PDB file or a silent file.

# Molecular Dynamic methods
This folder contains subdirectories with code for performing the molecular dynamic simulations.
