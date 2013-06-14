/*************************************************************************
*
* Copyright:		Max Planck Institute for Astrophysics (MPA)
* 
* File:				modules.cpp
*
* Routine Info:		Implementation of supplementary functions for 
* 					operations on strings
* 					- also implements functions from Numerical
* 					recipes C++
* 					       			
* Author:
*
* Modification 
* Log:	    		Added routine main
*		        	Shrikant V	Dec 11, 2012	16:20
*
**************************************************************************/

#include "modules.hpp"
#include <algorithm>
#include <vector>
#include <math.h>

using namespace std;

// NOTE: IT WOULD BE A GOOD IDEA TO SEPARATE THE IDL IMPLEMENTATIONS IN DIFFERENT FILES ALTOGETHER

// Implemtation of IDL poly with 2 coefficients
void poly( float * fX, int iSizefx, float Coeff1, float Coeff2 , float * fOut)
{
	// POLYNOMIAL EVALUATION
	// fX - variable array (can be scalar, vector or array)
	// Coeff* - Coefficients	
	// Computes c0 + x * c1 + x^2 * c2 +...
	// degree of poly = #Coeff* - 1
	float x = 0;
	for( int i=0; i< iSizefx; i++)
	{
		x =  fX[i]; 
		fOut[i] = Coeff1 + x * Coeff2;
	}	
}

// Implemtation of IDL poly with 3 coefficients
void poly( float * fX, int iSizefx, float Coeff1, float Coeff2, float Coeff3, float * fOut )
{
	// POLYNOMIAL EVALUATION
	// fX - variable array (can be scalar, vector or array)
	// Coeff* - Coefficients	
	// Computes c0 + x * c1 + x^2 * c2 +...
	// degree of poly = #Coeff* - 1
	float x = 0;
	for( int i=0; i< iSizefx; i++)
	{
		x =  fX[i]; 
		fOut[i] = Coeff1 + x * ( Coeff2 + x * Coeff3);
	}
}
// Implementation of IDL convol
void convol( float * fObs, float * fData, float * fObs1, float * fData1, int iElements, int iVdop)
{
	// The CONVOL function convolves an array with a kernel, and returns the result. Convolution is a general process 
	// that can be used for various types of smoothing, signal processing, shifting, differentiation, edge detection, etc.
	// The CENTER keyword controls the alignment of the kernel with the array and the ordering of the kernel elements. 
	// If CENTER is explicitly set to 0, convolution is performed in the strict mathematical sense, otherwise the kernel 
	// is centered over each data point.
	
	
	// ###################  NOTE   ################################
	// Control never enters the place where convol is written in cont_rscl.cpp.
	// so skipping the implementation of convol
	// Here, this function is just a dummy to avoid compilation errors
}

// Implementation of IDL integ for regular grids ( two input parameters )
// Integrated a function provided as an arry of points
void integ( double * X, double * Y, int iArraySz, double * dRes )
{
    // THEORY FROM IDL 
    /*
    PROCEDURE - Simpson Integration, where the mid-interval points are obtained
    from cubic interpolation using Neville's Algorithm
    */
	
	// Y -> vector containing y coordinates of data
	// X -> vector containing x coordinates
	// iArraySize -> size of array (X or Y )
	// Result -> resultant vector containing integration of Y coordinates from X(0) to X(i)
	
	// Setting all values in dRes = 0
	for(int i=0; i < iArraySz; i++)
	{
		dRes[i] = 0.0d;
	}
	
	// kernel to do integration ( using trapeziodal rule )
    int ilo = 0; // low of absicca
    int ihi = 0; // highest sample
    
    int n = 0;
    double z[iArraySz];
    double sum = 0.0d;
    /******************************************************
     * As per this logic, dRes will contain only iArraySz-2 
     * valid values, dRes[iArraySz-1] would be 0
     * ****************************************************/
	    
    // not taking into consideration imin & imax logic
    // just perforing a simple integration'
    // code help taken from www.astro.washington + integ + idl    
    for(int i=0; i < iArraySz-1; i++)
    {
    	ihi = i;
    	n = ihi - ilo;
    	sum = 0.0d;
    	
    	for( int j=0; j< i; j++)// to compute integral upto 'i'
	    {
    		z[j] = 0; //resetting value
	    	z[j] = Y[j+1] + Y[j]; // computing value
	    }
    	
    	for( int j=0; j< i; j++) // to compute integral upto 'i'
	    {
	    	sum += ( (z[j]/2) * (X[j+1]-X[j]) );
	    }
    	
    	dRes[i] = sum; 
    }

/*    
    // kernel to do integration ( using trapeziodal rule )
    int ilo = 0; // low of absicca
    int ihi = iArraysz-1; // highest sample
    
    int n = ihi - ilo;
    double z[iArraySz];
    double sum = 0.0d;
    
    for( int i=0; i< iArraysz-1; i++)
    {
    	z[i] = Y[i+1] + y[i];
    }
    
    for( int i=0; i< iArraySz-1; i++)
    {
    	sum += ( (z[i]/2) * (X[i+1]-X[i]) );
    }
 */
    
}


// Implementation of IDL interpol
// Interpol function performs linear, quadratic or spline interpolation on 
// vectors with a regular or irregular grid
void interpol( float * Y, float * X, int iSize, double * X1 , float * Y1out, int iOutSz )
{
	// As in the parent IDL code, everywhere linear interpolation is used
	// Implementation of linear interpolation here
	
	// Code from Numerical recipes is not available
	// Writing own version of interpol routing
	
	// Y - Input vector
	// X - Abcissa values for Y (assuming X monotonically increasing)
	// X1 - Absicca values for Y1out
	// Y1out- result
	int i,j;
	float x = 0;
	float x1, x2,y1,y2;
	for( i=0; i< iOutSz; i++ ) //outer loop for every value in X1
	{
		x = (float) X1[i];
		
		x1 = 0; y1 = 0;
		x2 = 0; y2 = 0;
		
		// find values x1, x2 that are closest to x, s.t. x1 < x < x2
		for( j = 1; j < iSize; j++ )
		{
			if( X[j] > x )
			{
				x1 = X[j-1];
				x2 = X[j];
				y1 = Y[j-1];
				y2 = Y[j];
				break;
			}
		}
		
		if( j == iSize ) //loop didn't broke 
		{
			Y1out[i] = -1; // invalid value, shouldn't happen
		}
		else
		{
			// Linear interpolation
			Y1out[i] = y1 + (y2-y1) * (x-x1) / (x2-x1);
		}
	}
}

// Implementation of IDL interpol
// Interpol function performs linear, quadratic or spline interpolation on 
// vectors with a regular or irregular grid
void interpol( float * Y, float * X, int iSize, float * X1 , float * Y1out, int iOutSz )
{
	// As in the parent IDL code, everywhere linear interpolation is used
	// Implementation of linear interpolation here
	
	// Code from Numerical recipes is not available
	// Writing own version of interpol routing
	
	// Y - Input vector
	// X - Abcissa values for Y (assuming X monotonically increasing)
	// X1 - Absicca values for Y1out
	// Y1out- result
	int i,j;
	float x = 0;
	float x1, x2,y1,y2;
	for( i=0; i< iOutSz; i++ ) //outer loop for every value in X1
	{
		x = X1[i];
		
		x1 = 0; y1 = 0;
		x2 = 0; y2 = 0;
		// find values x1, x2 that are closest to x, s.t. x1 < x < x2
		for( j = 1; j < iSize; j++ )
		{
			if( X[j] > x )
			{
				x1 = X[j-1];
				x2 = X[j];
				y1 = Y[j-1];
				y2 = Y[j];
				break;
			}
		}
		
		if( j == iSize ) //loop didn't broke 
		{
			Y1out[i] = -1; // invalid value, shouldn't happen
		}
		else
		{
			// Linear interpolation
			Y1out[i] = y1 + (y2-y1) * (x-x1) / (x2-x1);
		}
	}
}

// Implementation of IDL gaussfold
// Smoothing of plot by convolving with a Gaussian profile
void gaussFold( float * fWavelen, float * fFlux, int iSize, float fSigma, int iRadius )
{
    // SMOOTHING OF FLUX WITH A GAUSSIAN PROFILE
	// IMPLEMENTATION
	
	// ****************** NOTE ********************
	// Assuming a radius of 2 by default, need to check
	// Content taken from
	// http://homepages.inf.ed.ac.uk/rbf/HIPR2/gsmooth.htm
	
	int kSize = iRadius*2 + 1;
	
	// Compute Gaussian Profile oder Kernel
	float kernel[kSize];
	float x=0, expo= 0;
	for(int i =0; i < kSize; i++)
	{
		x = i - iRadius;
		expo = -1 * (x*x) / (fSigma*fSigma); 
		kernel[i] = (1/ (sqrt(2*M_PI) * fSigma)) * exp(expo); 
	}
	
	// temporary array to hold smoothed flux
	float smflux[iSize];
	smflux[0] = fFlux[0];
	smflux[1] = fFlux[1];
	
	// convoluting the spectrum with kernel that we have calculated
	for( int i=2; i < iSize; i++)
	{
		for(int j=0; j< kSize; j++)
		{
			smflux[i] = fFlux[(i-iRadius)+j] * kernel[j];
		}
	}
	
	//copying back the values to original flux array
	for( int i=0; i< iSize; i++)
	{
		fFlux[i] = smflux[i];
	}
}

// Implementation of median filtering from getchi.cpp
void maxMedianFilter( float * fData, int numElements, int iRange )
{
    // Implementation of fData = fData/ max( median(fData,3))
    float fTemp[ numElements ];

    float fRange[ iRange ];
    float fMaxMedian = - 1000000.00;

    // calculating median
    for( int i = 0; i < numElements; i++ )
    {
        if( i >= iRange - 1 )
        {
            for( int j =0; j < iRange; j++ )
            {
                // copy to fRange
                fRange[j] = fData[i-(iRange/2)+j];
            }
        
            // sort
            std::vector<float> fVector( fRange, fRange + iRange);
    
            std::sort( fVector.begin(), fVector.end() );

            // find median (assuming odd iRange) & put in array
            fTemp[i] = fVector[ iRange/2 ];

            // calculate max median
            if( fTemp[i] > fMaxMedian )
            {
                fMaxMedian = fTemp[i];
            }             
        }
        else
        {
            fTemp[i] = fData[i];

            // calculate max median
            if( fTemp[i] > fMaxMedian )
            {
                fMaxMedian = fTemp[i];
            }
        }
    } 

    // transforming array
    for( int i = 0; i < numElements; i++ )
    {
        fData[i] /= fMaxMedian;
    }
}

// function for trimming the standard C++ string
// On the lines of IDL, leading blanks are removed if mode is 1
// and both leading & trailing are removed if mode is 2
void strTrim( string &strS , int iMode )
{
    if( iMode == 1 || iMode == 2 ) // remove leading spaces
    {
        strS.erase(0, strS.find_first_not_of(" \n\r\t"));
    }
    if( iMode == 2 ) // remove trailing spaces
    {
        strS.erase( strS.find_last_not_of(" \n\r\t") + 1 );
    }
    
    return;
}

// function to change to the Upper case of string
void strUpper( string &strS )
{
    for( int iCntr = 0; iCntr < strS.length() ; iCntr++ )
    {
        if( strS[iCntr] >= 97 && strS[iCntr] <= 122 )
        {
            strS[iCntr] = 65 + (strS[iCntr]-97);
        }
    }
}


// function for tokenising the string based on Delim character
// need generalised function as the string can be tokenised many times
string strNexttoken( string &strS, char cDelim )
{
    int iIndex = 0;
    string strToken("");

    iIndex = strS.find_first_of(cDelim);

    if( iIndex != 0 )   // token exists
    {
        strToken = strS.substr(0, iIndex);

        strS.erase(0,iIndex+1 );  // trimming the parent string
        strTrim( strS, 2 );

        strTrim( strToken, 2 );
        strUpper( strToken );
    }

    if( iIndex == 0  /*|| strS->empty()*/ )
    {
        // search for characters other than cDelim which can be delimiters
        iIndex = strS.find_first_of("\0\n\r");  

        if( iIndex != 0 ) // last token exists
        {
            strToken = strS.substr(0, iIndex-1);
            strTrim( strToken, 2 );
            strUpper( strToken );
        }

        // At this point no more tokens exists
        strS.erase(0);  // erase remaining unwanted characters (if any)
    }

    return strToken;
}


