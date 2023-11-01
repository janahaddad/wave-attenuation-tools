# Wave Analysis Workshop at UNC Chapel Hill Institute of Marine Sciences (IMS) 
These files were compiled for a 2021 workshop/tutorial on wave analysis in natural marshes and living shorelines, for collaborators at UNC-IMS. 

## make_array & make_spectra
These scripts use data from the RBRsolos deployed at PKSAq in october 2019 for preliminary data. 

make_array_PKSAq.m is the script that reads text files outputted by the RBR solftware and makes a mat file with 
pressure, time, and height above bed. The output of this is LS_PKSAq_prelim.mat 

make_spectra_PKSAq.m reads the mat file & met information located in /met and does the spectral analysis.
It also plots and saves the output spectrum. 

wavenumber.m is a matlab function that solves the wave dispersion relation. make_spectra uses this. 
