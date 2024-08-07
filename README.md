# EPSILOMETER
aleboyer@ucsd.edu | ncouto@ucsd.edu

This is the epsilometer matlab library.
It can read the .ascii files generated from the SOM acquisition system or directly the file store on the SD card.

- In matlab add this library to your path. 
- change directory to the location of your data.
- type something like : "data=epsi_class";

The library will create an epis_class which 
- read the data files,  
- convert counts to volts if needed
- contains methods to plot and analyze data. 

The derivation of epsilon and chi is one of the methods. 

Logic flow of the derivations:
- split data into "profiles"
- split profiles in segments
- compute spectra over these segments
- compute epsi and chi over these segments. 


 -------

Epsi processing manual can be found at https://www.overleaf.com/read/qbpchcnywdqp
Additional how-tos can be found at https://nicolecouto.notion.site/MOD-Fish-Wiki-0c1a74e213a4447eae5c8f826a1c34cd
Contact Nicole Couto for questions or suggestions - ncouto@ucsd.edu


 -------



