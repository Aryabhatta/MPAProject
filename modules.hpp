/*************************************************************************
*
* Copyright:		Max Planck Institute for Astrophysics (MPA)
* 
* File:				modules.hpp
* 
* Routine Info:		Declaration of necessary numerical methodss
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
#include <fstream>
#include "idlFuncn.hpp"
using namespace std;

#ifndef MODULES_H
#define MODULES_H

#define NotDefined -32766 // -32767 + 1
#define NULLPTR 0

string ReadInput( string name );

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

bool keyword_set(int iNumber);

#endif