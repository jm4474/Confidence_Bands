This README.txt file was generated on 08/08/2018 by 

José Luis Montiel Olea and Mikkel Plagborg-Moller


----------------------
i) GENERAL INFORMATION
----------------------

The folders

1SimInferenceClass
2Gertler_Karadi_application
3Head_Mayer_Ries_application
4Additional_Figures
Reg_Sens
VAR_IRF

contain .csv files, Matlab scripts/functions/classes, and STATA do files to generate the figures reported in the paper "Simultaneous Confidence Bands: Theory, Implementation, and an Application to SVARs" by José Luis Montiel Olea and Mikkel Plagborg-Moller.  


--------------------------
ii) HARDWARE/SOFTWARE 
(specifications and requirements)
-------------------------- 

All the files have been tested on both:
 
* A MacBook Pro @2.4 GHz Intel Core i7 (8 GB 1600 MHz DDR3) running Matlab 2016b and Stata 13.

* A Lenovo Thinkpad @2.3 GHz Intel Core i5 (8 GB RAM) running Matlab R2017a and Stata 15.


--------------------------
iii) RECOMMENDED CITATION
-------------------------- 

When using this code please cite:

"Simultaneous Confidence Bands: Theory, Implementation, and an Application to SVARs", Montiel Olea, J.L. and Plagborg-Moller, Journal of Applied Econometrics, 2018.


---------------------
iv) DATA & MAIN FILE OVERVIEW
---------------------

* 1SimInferenceClass
    
This folder contains the "SimInference.m" Matlab class file, which collects different Matlab functions that are used to create the sup-t band, and other popular simultaneous bands (such as Bonferroni, Sidak, and Projection). This Matlab class also contains a simple algorithm to implement the "calibrated" Bootstrap/Bayes sup-t band.    

NOTE: Both applications call the SimInference.m class.      

* 2Gertler_Karadi_application

This folder contains the .csv files and Matlab scripts to replicate the figures related to the Gertler-Karadi Structural VAR application. The two main files for replication (both in the /Script folder) are:

run_gk_iv.m
run_gk_chol.m

The first file replicates Figure 2 and the second file Figure 3 in the paper (simply run the files on the Matlab command window or section by section). 

NOTE: To generate Figure 6, simply change line 43 and 50 in run_gk_iv.m. To generate Figure 7, simply change line 44 and 51 in run_gk_chol.m 

* 3Head_Mayer_Ries_application 

This folder contains the Stata file and Matlab scripts to replicate the figures related to our sensitivity analysis for the Head-Mayer-Ries application. The main file for replication (in the /Script folder) is:

run_hmr.m

This file generates Figure 8 in the paper.

NOTE: To run the Matlab file, you must perform the following three steps first:
a) Download the following zip file: http://econ.sciences-po.fr/sites/default/files/file/tmayer/data/col_regfile09.zip
b) Unzip the Stata data file "col_regfile09.dta" and place it in the subfolder 3Head_Mayer_Ries_application/Data
c) Run the Stata do-file create.do in the subfolder 3Head_Mayer_Ries_application/Data
These steps will create a large .csv file used by the above-mentioned Matlab script "run_hmr.m". The latter file is currently set to draw only 100 bootstrap and Bayes draws, which takes a couple of hours on a personal laptop. To increase the number of draws to 2,000 as in the paper, simply change lines 13 and 14 in "run_hmr.m".

* 4Additional_Figures

This folder contains Matlab scripts to replicate Figures 1, 4, and 5.


---------------------
iv) Additional Folders
---------------------

The folders Reg_Sens and VAR_IRF contain application-specific functions for the regression sensitivity analysis and for the VAR application.

