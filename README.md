# SN_data_analysis
Tools for managing, analyzing, and plotting supernova (SN) data:

1. Supernova Dataset Manager (dir: sn_data_manager):

	- Parses .json files from The Open Supernova Catalog (https://sne.space))
	  for data including light curves to process and write to plot files 
	- I compile with: g++-10 -std=c++17 -lstdc++fs sn_data.cxx
		- you will need a compiler that has experimental/filesystem

2. Analyze Spectra (dir: analyze_spectra):
	
	- Reads spectral data from a csv file, computes the signal-to-noise (SNR)
	  ratio, plots spectra (according to user specifications) as well as SNR.
	- I run as: python analyze_spectra.py data1.csv 
	 
