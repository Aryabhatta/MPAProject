
#include <iostream>
#include <math.h>
#include "idlFuncn.hpp"
#include "modules.hpp"
using namespace std;

#ifndef CONTRSCL_H
#define CONTRSCL_H

void cont_rscl( float * fWavelen, float * fData, int iObsElem, float * fSx, float * fSy, int iThrElem, float * fRy,
		           float * fSig = NULLPTR, int iElemSig = NotDefined, float * fRf = NULLPTR, int iElemRf = NotDefined,
		           int iMessage = 0, int iStrong = 0, int iConst = 0, float fSigma = 0.03f, int iUpcut = NotDefined,
		           float * fW_highprio = NULLPTR, int iElemhighprio = NotDefined, int iIter = NotDefined, 
		           int iLin = NotDefined, int iVdop = NotDefined );

#endif