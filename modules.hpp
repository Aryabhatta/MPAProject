#include <iostream>
#include <string>
#include "idlFuncn.hpp"

using namespace std;

#ifndef MODULES_H
#define MODULES_H

// regular grids
//double * integ( double * dNX_SP, double * dNY_SP ); // to pick the routine from Numerical Recipes site
void integ( double * X, double * Y, int iArraySz, double * dRes );

void gaussFold( float * fWavelen, float * fFlux, int iSize, float fSigma, int iRadius = 2 ); 

//done with writing below functions
void convol( float * fObs, float * fData, float * fObs1, float * fData1, int iElements, int iVdop);

void poly( float * fX, int iSizefx, float Coeff1, float Coeff2 , float * fOut);
void poly( float * fX, int iSizefx, float Coeff1, float Coeff2, float Coeff3, float * fOut );

void interpol( float * Y, float * X, int iSize, double * X1 , float * Y1out, int iOutSz );
void interpol( float * Y, float * X, int iSize, float * X1 , float * Y1out, int iOutSz );

void maxMedianFilter( float * fData, int numElements, int iRange );

void strTrim( string &strS , int iMode );

void strUpper( string &strS );

string strNexttoken( string &strS, char cDelim );

#endif