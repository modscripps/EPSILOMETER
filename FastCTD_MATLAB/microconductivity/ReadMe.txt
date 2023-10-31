FCTD microconductivity processing

Scripts: 
B_FCTD_addchi_grid_template.m
B_FastCTD_GridData.m
FCTD_ucon_cal.m
add_chi_microMHA_v2.m
FCTD_DefaultChiParam.m
+ function dependencies are in the other_necessary_scripts folder

To run the processing:
The main script is B_FCTD_addchi_grid_template. All the processing steps happen in this script.
Edit the file locations and the steps that you're interested in doing in
B_FCTD_addchi_grid_template (First section). Also specify the time period you
want files from (Second section). Check that the parameters in
FCTD_DefaultChiParam.m are correct for this deployment.
If additional variables are required (e.g. fluorometer, bottom depth), go
to the data gridding section (fifth section) and include those in VarsToGrid and
add in any additional steps needed for those (BLT steps are currently commented out).
Then the script can be run as a whole. 
