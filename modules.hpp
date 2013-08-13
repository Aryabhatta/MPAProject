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

// Routine to log array for plotting it afterwards with gnuplot
// NOTE: USE bPlot AS SWITCH TO ALLOW LOGGING OR NOT
template<class T>
void createLog( string FileName, bool bPlot, T * Arr1, int iArrSz1, T * Arr2=0, int iArrSz2=0 )
{
	if( bPlot == false )
		return;
	
	ofstream logFile;
	
	// Logs the observed spectra at standard log directory
	string strLogDir = ReadInput( "DIR:LOGDIR" );
	
	string strFileSpecs = strLogDir + FileName;
	
	logFile.open( strFileSpecs.data(), ios::out );
	
	if( iArrSz2 != 0) // need to log both arrays
	{	
		for( int i=0; i< iArrSz1; i++)
		{
			logFile << Arr1[i] << "\t" << Arr2[i] << endl;
		}
	}
	else // logging only one array
	{
		for( int i=0; i< iArrSz1; i++)
		{
			logFile << Arr1[i] << endl;
		}
	}
	
	logFile.close();
	logFile.clear();
}

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