
#include <iostream>
#include <math.h>
#include "idlFuncn.hpp"
#include "modules.hpp"
using namespace std;

#ifndef CONTRSCL_H
#define CONTRSCL_H

float * cont_rscl( float * fWavelen, float * fData, float * fSx, float * fSy, int iElements, float * fRy,
		           float * fSig = 0, int iElemSig = 0, float * fRf = 0, int iElemRf = 0,
		           int iMessage = 0, int iStrong = 0, int iConst = 0, float fSigma = 0.03f, int iUpcut = -1,
		           float * fW_highprio = 0, int iElemhighprio = 0, int iIter = -1, int iLin = -1, int iVdop = -1 );

#endif