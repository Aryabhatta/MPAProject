#include <iostream>
#include <string>
#include "read_lf_grid.hpp"

using namespace std;

#ifndef LFP_H
#define LFP_H

bool lfp( float * fSx, float *fSy, float fTeff, float fLogg, float fLogz, float fXi, float fEps_dev[2], float fGauss, float fGamma, int iExtrapol, string strGrid, string strRange, int iNomessage,int iCnvl, int iNomc);

#endif