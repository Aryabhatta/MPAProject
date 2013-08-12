MPAProject
==========

GENERAL INSTRUCTIONS:

1. To avoid long running times, I had taken just a single point in file 
 'foutr1p1082f180'. I'm attaching the file that I used for running my project.

2. Right now, the project runs for single mode 'HALPHA'. I had to that for
 debugging. Also, I didn't have other .FITS files (grid) in the correct format.

3. The inputs to the project are:
       i. File of grid points  e.g. 'foutr1p1082f180'
       ii. Fits file containing observed flux e.g. 'spPlate-1962-53321.fits'
       iii. Definition file e.g. 'lf_grid6520.def'
       iii. The corresponding grid file in .Fits format e.g. 'lf_grid6520.fits'

4. The output of project are below files:
       i. 'thrSpectrum.log' containing thereotical spectrum
       ii. 'ObsSpectrum.log' containing observed spectrum
      iii. logfile as per fits file name e.g. 'spPlate-1960-53287_500.log'


INSTRUCTIONS FOR RUNNING THE PROJECT:

1. Change the file paths in file getchi.cpp & lfp.cpp. I have marked the places where the paths need to by changed by below block

/************************************************
 * CHANGE PATH BELOW 
 ***********************************************/

Please change the path for the line just below this comment block.

2. Build the project by doing 'make'. Creates executable 'project'

3. Run the project by './project'

4. To clean the project, do 'make cleanAll'
