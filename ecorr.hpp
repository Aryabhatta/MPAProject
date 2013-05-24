#include <iostream>
#include <math.h>
#include "idlFuncn.hpp"
#include "modules.hpp"
using namespace std;

typedef unsigned char BYTE;

#ifndef ECORR_H
#define ECORR_H

double * ecorr( float * fThrWave, float * fThrData, int iThrElem, float * fObsWave, float * fObsdata, int iObsElem,
		        int iXsigma, float * fRf = 0, int iSzrf = 0 , BYTE * mask = 0, int iSzmask = 0, int iNomessage = 0, 
		        double dPn = 0, double dPp =0, float * fWr = 0, float fSigma = 0 );

#endif
		       
// RULE: You cannot have a non-defualt argument after your defualt arguments begin !!!
		        
		        
		        
		        
		        
		        
		        
		        
		        
		        
		        
		        
		        
		        
		        
		        
		        
		        
		        