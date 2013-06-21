#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "idlFuncn.hpp"
#include "modules.hpp"

using namespace std;

#ifndef GETRV_H
#define GETRV_H

void FFT_Prep( double * dNX_SP, float * fYArr, double * dNY_SP , int iSzNXSP);

float get_rv( float * ObsWave, float * ObsFlux, int iObsElem, float * ThrWave, float * ThrFlux, int iThrElem, float * fXr,  
		      int iNomessage = 0,int iAbsolute = 0, float eps = 0.0, int iLog = -1);

bool Log_Lin_Corr( float * fObsW, float * fObsF, int iObsSz, float * fThrW, float * fThrF, int iThrSz, float * fConvX, 
		           float * fConvY, int * iConvSz, bool bComplete = true , int iNomc = 0, float eps =8*pow(10,-7));

#endif