# MergeSplit
##Codes and Data for *"Fission-fusion dynamics and group-size dependent composition in heterogeneous populations"*, 2018

##Codes:
1) C++11 source-code (merge_split_13.cpp), a Monte-Carlo algorithm to simulate the Coagulation-Fragmentation model for heterogeneous populations. Written by Gokul G. Nair
To compile the code:
> g++ -std=c++11 merge_split_13.cpp -o MergeSplit13

**It is recommended to create a folder in the same directory as the executable by the name of '*Output*'**
This creates the executable, MergeSplit13 which can be executed with:
> ./MergeSplit13

The Program asks the user to input: *'s'* (number of sites), *'N'* (total population), *'N1'* (type-I population), *'Ps0'* (denoted as *p_0'* in the manuscript, the base split rate), *'Pm'* (denoted as *'q'* in the manuscript, merge rate), *'d'* (denoted as *'\delta'* in the manuscript, excess split rate), *'a'* (an aditional parameter that is not mentioned in the manuscript. **Always give 1**), and the number of events (typically 10^7-10^9).

For a more detailed explanation about the source-code, dependencies, and output, refer to *MergeSplit_Documentation.pdf*.

2) Python script (FPI.py), to find the steady-state solution to the iterative equations (Fixed point iteration) derived in the manuscript.
To run it:
> python FPI.py

**The file input100_f has to be in the same directory**
Inputs for the program can be given inside the source-code. It runs for a sequence of delta values and creates a file called *FPI_N=xx_N1=xx_P0=xx_d=xx_Pm=xx.csv* for each value of delta. It also creates a folder for each of these files containing the plots of the distribution for different slices.

3) Data files
The data files used to create the figures in the manuscript are given in this .zip file. Filenames containing *MC-* are from the Monte-Carlo simulation, while those containing *FPI-* are from the fixed point iteration.
