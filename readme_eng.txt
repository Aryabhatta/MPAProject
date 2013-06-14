 1990-2000 J. Reetz, Universitaets-Sternwarte LMU, Muenchen 
 2012      M. Bergemann, Max-Planck Institute for Astrophysics
--------------------------------------------------------------------

System to construct a grid of theoretical spectra.

Components:


 lf_grid.def    = Basic definitions of parameters and parameter space
                  lf_grid<n>.def is located under /SIU/mod/grid/
                  and depends on the wavelength region, e.g.
                  lf_grid6520.def was computed for the wavelength
                  region 6420 - 6620 Angstrem.

 gen_grid.pro   = creates a grid of theoretical spectra
                  usign the input parameters from lf_grid<n>.def
                  (see definitions in the program). Important:
		  lf_grid<n>.log contains the filenames of each 
		  pre-computed spectrum. Spectra which failed, are
                  marked with 'F'.
		  The data can be read with 
                  read_bin,<filename>,x,y,id=id,obj=obj,com=com

 lfp.pro       =  performs interpolation in the grids of spectra using the
 	          the from lf_grid<n>.def and lf_grid<n>.grid.
		  'extrapol' parameter indicates extrapolation.
 Example:
 lfp, range=<Spectral region, e.g.."Halpha">,
      grid =<Grid, e.g. "high-res">,
      w,f,teff=teff,logg=logg,logz=logz,xi=xi,eps_dev=eps_dev,
      gauss=gauss,gamma=gamma,vsini=vsini,expo=expo,rt=rt,extrapol=extrapol

 getchi.pro    =  basis input and output for comparison of the
                  observed and theoretical spectra. 
		  calls: lfp.pro, ecorr.pro

 ecorr.pro     =  compute G-test statistics for the theoretical
                  spectrum (tx - wavelength array, ty - flux array) 
                  and the observed spectrum (ox,oy).
