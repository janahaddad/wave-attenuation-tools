These are data from the RBRsolos deployed at PKSAq in october 2019 for preliminary data. 

make_array_PKSAq.m is the script that reads text files outputted by Ruskin and makes a mat file with 
pressure, time, and height above bed. The output of this is LS_PKSAq_prelim.mat 

make_spectra_PKSAq.m reads the mat file & met information located in /met and does the spectral analysis.
It also plots and saves the output spectrum. 

wavenumber.m is a matlab function that solves the wave dispersion relation. make_spectra uses this. 

*note: pretty sure I set one of these sensors to record in UTC time by accident, so the output may 
look wonky as a result (Jana)
