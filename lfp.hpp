/*************************************************************************
*
* Copyright:		Max Planck Institute for Astrophysics (MPA)
* 
* File:				lfp.hpp
* 
* Routine Info:		Declaration of methods for lfp.cpp
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
#include "read_lf_grid.hpp"

using namespace std;

#ifndef LFP_H
#define LFP_H

bool lfp( float ** fSx, float **fSy, int * iThrElem, float fTeff, float fLogg, float fLogz, float fXi, float fEps_dev[2], float fGauss, 
		float fGamma, int iExtrapol, string strGrid, string strRange, int iNomessage,int iCnvl, int iNomc);

#endif