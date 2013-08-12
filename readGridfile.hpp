/*************************************************************************
*
* Copyright:		Max Planck Institute for Astrophysics (MPA)
* 
* File:				readGridFile.hpp
* 
* Routine Info:		Declaration of methods for reading from grid file
* 					( which is in FITS format )
*
* Author:
*
* Modification 
* Log:	    		
*		        	
*
**************************************************************************/

#include <iostream>
#include <string>
#include <fitsio.h>
using namespace std;

#ifndef READGRID_H
#define READGRID_H

void readGrid( string GridFileSpecs , float * fWavelen, float * fFlux);

// Read flux paramaters belonging to flux #iIdx
bool readGridParams( string strGridFileSpecs, float * pars, int iIdx);

// returns true if the row belonging to iIdx exists in the grid file (FITS HDU 1)
bool GridDescrExists( string strGridFileSpecs, int iIdx);

// Read wavelength array, returns status , #wavelengths in fWave array
bool readGridWave( string strGridFileSpecs, float ** fWave,  int * iWaveCnt);

// Read flux belonging to certain iIdx
bool readGridFlux( string strGridFileSpecs, float * flux, int iIdx);

// Read Dimensions
bool readGridDim( string strGridFileSpecs, int iHdu,  int * naxis1, int * naxis2 );

#endif