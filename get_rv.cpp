/*************************************************************************
*
* Copyright:		Max Planck Institute for Astrophysics (MPA)
* 
* File:				get_rv.cpp
*
* Routine Info:		Calculates the radial velocity
*       			
* Author:
*
* Modification 
* Log:	    		Added routine main
*		        	Shrikant V	Dec 11, 2012	16:20
*
**************************************************************************/

#include <iostream>
#include <fitsio.h>
#include "get_rv.hpp"

using namespace std;

// Routine to calculate the radial velocity
// copy header part from get_rv.pro
float get_rv( float * ObsWave, float * ObsFlux, int iObsElem, float * ThrWave, float * ThrFlux, int iThrElem, float * fXr,  int iNomessage,int iAbsolute, float fEps, int iLog)	
{
    float fXrv = 0.0;
    int i  = 0;
    
    if( fEps = 0.0 ) {     fEps = 8 * pow(10,-7);    } // Min Floating accuracy
    // iAbsolute is default 0
    // iNomessage is default 0
    
    int IndexArr[ iObsElem ];
    int iCnt = 0;
    fXr[0] = 6457.48;// dummy
    fXr[1] = 6500;// dummy

    iCnt=IdlWhere( ObsWave, ">=", fXr[0], "<=", fXr[1], iObsElem, IndexArr);
    
    // filtered observed spectrum
    float fObsW[ iCnt ];
    float fObsF[ iCnt ];
    int iObsSz = iCnt;
    
    for( i=0; i< iCnt; i++ )
    {
        fObsW[i] = ObsWave[ IndexArr[i] ];
        fObsF[i] = ObsFlux[ IndexArr[i] ];
    }
    
    int IndexArr1[ iThrElem ];
    int iCnt1 = 0;
    
    iCnt1=IdlWhere( ThrWave, ">=", fXr[0], "<=", fXr[1], iThrElem, IndexArr1);
    
//    cout << "Xr0= " << fXr[0] << " Xr1=" << fXr[1] << " Size: " << iThrElem << endl;
    
    // filtered thereotical spectrum
    float fThrW[ iCnt1 ];
    float fThrF[ iCnt1 ];
    int iThrSz = iCnt1;

    for( i=0; i< iCnt1; i++ )
    {
//    	cout << "Index " << i << ":" << IndexArr1[i] << " ";
        fThrW[i] = ThrWave[ IndexArr1[i] ];
        fThrF[i] = ThrFlux[ IndexArr1[i] ];
    }
    
    cout << "size 1 = " << iCnt1 << endl;
    
    /* IF NO OF PARAMETERS < 2, LOAD SOLAR FLUX SPECRUM AS REFERENCE
       SOURCE FILE read_kpno MISSING
    */
    if( iThrElem==0 )
    {
	    string strSolarFlux("/home/shrikant/Desktop/MPA/Files/solarflux.fits");
	    
	    cout << endl << "Loading solar flux spectrum from fits file as reference" << endl;
	    
	    fitsfile * fptr;
	    int iStatus = 0; // Initialise status
	    
	    // Opening the FITS file
	    fits_open_file( &fptr, strSolarFlux.data(), READONLY, &iStatus );
	    
	    int naxis1;
	    char * comment = new char [100];
	    
	    // Since we know the format of this file, just reading naxis1 which is
	    // #of rows in the table
	    fits_read_key( fptr, TINT, "NAXIS1" , &naxis1, comment, &iStatus);
	    cout << endl << "Naxis1 = #cols = " << naxis1 << "  Comment = " << \
	    comment<< endl;
	    
	    // Creating a large float array to read elements from flux file
	    // NOTE: however, creating such a large array is discouraged & 
	    // the possibility to access small chunks of it should be seen
	    	      
	    int hdutype;
	    fits_get_hdu_type(fptr, &hdutype, &iStatus);
	    switch(hdutype)
	    {
	    case IMAGE_HDU: cout << "Image HDU " << endl; break;
	    case ASCII_TBL: cout <<  "Ascii Table" << endl; break;
	    case BINARY_TBL: cout << "Binary Table" << endl; break;
	    }
	    
//	    int hdunum =0;
//	    fits_get_num_hdus(fptr, &hdunum, &iStatus );
//	    cout << "No of HDU's:" << hdunum << endl;
	    
//	    long nrows = 0;
//	    fits_get_num_rows( fptr, &nrows, &iStatus );
//	    cout << "No of rows:" << nrows << endl;
	    
//	    fits_read_col( fptr, TFLOAT, 2, 1, 1, nelements, &nullval, fSolarflux, &anynull, &iStatus);
	    
	    //reading a small chunk of 10 elements
	    int iNoData = 10;
	    	    
	    float fSolarflux[iNoData]; //just contains flux, wavelength to be stored in separate array
	    
	    float nullval = 0;
	    int anynull = 0;
	    int nelements = iNoData;
	    int firstelem = 1; // reading from start, can be changed
	    
	    fits_read_img(fptr, TFLOAT, firstelem, nelements, &nullval, fSolarflux, &anynull, &iStatus);
	    	    
	    //printing fsolarflux
	    for(i=0;i<iNoData;i++)
	    {
	    	cout<< fSolarflux[i] << " ";
	    }

	    fits_close_file( fptr, &iStatus );

	    if( iStatus ) /* print any error messages */
	    {
	        fits_report_error( stderr, iStatus );
	    }

	    delete [] comment;
    }
    
    float * fConvX = 0, * fConvY = 0;
    bool bError = false;
    int iConvSz = 0;
    
    bError = Log_Lin_Corr(fObsW, fObsF, iObsSz, fThrW, fThrF, iThrSz, fConvX, fConvY, iConvSz);
    
    if( bError)
    {
        cout << "Error in LOG_LIN_CORR" << endl;
        cout << "Aborting from get_rv.." << endl;
        return -1;
    }
  
    double dC = 0;// find out what is C -CHANGE- 

    for( i = 0; i < iConvSz; i++ )
    {
        fConvX[i] = ( exp(fConvX[i]) - 1.00 ) * dC * pow(10,-5);
    }
   
    float fCoef =  0;
    float fConvYA[ iConvSz ];
    
    float fMaxCY = IdlMax( fConvY, iConvSz);
    if( abs( IdlMin( fConvY, iConvSz ) ) > fMaxCY )
        fMaxCY = abs( IdlMin( fConvY, iConvSz ) );

    for( i =0; i < iConvSz; i++ )
    {
        fConvYA[i] = fConvY[i] / fMaxCY;
    } 

    float fTemp = IdlMax( fConvY, iConvSz );

    fXrv = fConvX[ (int)dC ]; // corresponding to rv=cx(!C) in get_rv.cpp
                        // in Idl , ! is used to reference system variabbles
    
    int iSt = 1;

    if( ! iAbsolute )
    {
        for( i = 0; i< iConvSz; i++ )
        {
            fConvY[i] = fConvYA[i];
        }
    }

    if( iLog != -1 ) // iLog defined
    {
        fEps = pow(10,-4);
        int iIndex[ iConvSz ];
        int iCount = IdlWhere( fConvY, ">", fEps, iConvSz, iIndex );

        if( iCount > 0 )
        {
            for( i = 0; i < iCount; i++)
            {
                fConvY[ iIndex[i]] = log10( fConvY[ iIndex[i]] ) - log10(fEps);
            }
        }

        iCount = IdlWhere( fConvY, "<", -fEps , iConvSz, iIndex );
    
        if( iCount > 0 )
        {
            for( i = 0; i < iCount; i++ )
            {
                fConvY[ iIndex[i]] = (-1 * log10( -1* fConvY[iIndex[i]])) - \
                                    log10(fEps);
            }
        }

        iCount = IdlWhere( fConvY, ">=", (-fEps), "<=", fEps, iConvSz, iIndex);
        
        if( iCount > 0 )
        {
            for( i=0; i<iCount; i++ )
            {
                fConvY[ iIndex[i] ] = 0;
            }
        }
    }

    float fXe = 0, fYe = 0;
    // Call to EXTREMUM_B line 215, get_rv.pro  -CHANGE-

    if( iSt == 1 )
    {
        fXrv = fXe;
        fCoef = fYe;        
    }
    

    // delete memory allocated in Log_Lin_Corr
    delete [] fConvX;
    delete [] fConvY;

    return fXrv;
}


bool Log_Lin_Corr( float * fObsWave, float * fObsF, int iObsSz, float * fThrWave, float * fThrF, int iThrSz, float * fConvX, float * fConvY, int &iConvSz, bool bComplete, int iNomc, float eps)
{
    // defaults for bComplete, iNomc & eps set in header file
    int i =0;
    bool bError = false;
    
    // containing log of values from Observed & thereotical wavelengths
    float fObsW[iObsSz], fThrW[iThrSz];

    for( i = 0; i < iObsSz; i++ )
        fObsW[i] = log( fObsWave[i] );
    for( i = 0; i < iThrSz; i++ )
        fThrW[i] = log( fThrWave[i] );

    // Spectrum must have same widths
    float fOW_start = fObsW[0];
    float fOW_last = fObsW[iObsSz-1];
    float fTW_start = fThrW[0];
    float fTW_last = fThrW[iThrSz-1];
    float fOWdiff = fOW_last-fOW_start;
    float fTWdiff = fTW_last - fTW_start;

    if( fOWdiff != fTWdiff )
    {
        if(fOWdiff > fTWdiff)  
        {
            fOW_start = fObsW[0] + ( fOWdiff - fTWdiff)/2;
            fOW_last = fOW_start + fTWdiff;
        }
        else
        {
            fTW_start = fThrW[0] + ( fTWdiff - fOWdiff)/2;
            fTW_last = fTW_start + fOWdiff;
        }
    }

    // Finding min resampling distance
    float fRdistOW = 10000000, fRdistTW = 1000000; // max dummy values
    for( i =1; i <= iObsSz; i++ )
        if( fRdistOW < (fObsW[i] - fObsW[i-1]))
            fRdistOW = fObsW[i] - fObsW[i-1];

    for( i =1; i <= iThrSz; i++ )
        if( fRdistTW < (fThrW[i] - fThrW[i-1]))
            fRdistTW = fThrW[i] - fThrW[i-1];
    
    float fRdist = std::max( fRdistOW, fRdistTW ) * 0.5f ;
    
    // if resampling distance < than eps, eps taken as new resampling distance
    fRdist = std::max(fRdist, eps);

    long int lNrPix = abs((fOW_last - fOW_start)/fRdist) + 1L;
    long int lNrPixChk = abs((fTW_last - fTW_start)/fRdist) + 1L;
    
    if( lNrPix != lNrPixChk )
    {
        cout << "PIX.NR(SPECTRUM) - PIX.NR(RV-STD): " <<(lNrPix - lNrPixChk);
        cout << endl;
        if( lNrPix > lNrPixChk )
            lNrPix = lNrPixChk;
    } 

    // line 92 logic get_rv.pro, need to understand again
    if( lNrPix % 2 != 1 ) // if not odd number then, make one
    {
        lNrPix += 1;
        fRdist = (fOW_last - fOW_start)/ (float) (lNrPix -1 );
    }

    double dNX_SP[ lNrPix ];
    double dNX_RV[ lNrPix ];

    for( i = 0; i < lNrPix ; i++)
    {
        dNX_SP[i] = i * fRdist + fOW_start;
        dNX_RV[i] = i * fRdist + fTW_start; 
    }

    // Resampling
    // No need to compare the size of Arrays fObsW & fThrW as it is of the 
    // same size as fObsWave & fThrWave
    if( iNomc == 0 )
    {
        cout << "Points for FFT: " << lNrPix << endl;
        cout << "Interpolation:" << endl;
    } 
    //float fYArr[iObsSz], fBArr[iThrSz]; // initialising to dummy size as of now
    float fYAarr[lNrPix], fYBarr[lNrPix];
    
    interpol( fObsF, fObsW, iObsSz, dNX_SP, fYAarr, lNrPix );
    interpol( fThrF, fThrW, iThrSz, dNX_SP, fYBarr, lNrPix );

    double dNY_SP[ lNrPix ];
    double dNY_RV[ lNrPix ];
    
    FFT_Prep( dNX_SP, fYAarr, dNY_SP, lNrPix);
    FFT_Prep( dNX_RV, fYBarr, dNY_RV, lNrPix);

    // CORRELATION
    long lShift = lNrPix/2;
    fConvX = new float [ lNrPix ];
    iConvSz = lNrPix;

    double dMaxNxrv = IdlMax( dNX_RV, lNrPix );
    double dTemp = (dMaxNxrv - dNX_RV[0])/2.0d + dNX_RV[0];

    // ConvX
    for( i = 0; i < lNrPix; i++ )
    {
        fConvX[i] = dNX_SP[i] - dTemp;
    }

    double dNorm_SP = 0.0, dNorm_RV = 0.0;
    if( bComplete )
    {
        if( iNomc == 0 )
            cout << "Normalization of correlation function" << endl << endl;
    
        //double dMean_SP[ lNrPix], dMean_RV[ lNrPix ];
        double * dMean_SP = new double [ lNrPix ];
        double * dMean_RV = new double [ lNrPix ];
        
        integ( dNX_SP, dNY_SP, lNrPix, dMean_SP);
        integ( dNX_RV, dNY_RV, lNrPix, dMean_RV);
        
        double dDiv1 = IdlMax( dNX_SP, lNrPix) - dNX_SP[0];
        double dDiv2 = IdlMax( dNX_RV, lNrPix) - dNX_RV[0];
        
        for( i = 0; i < lNrPix; i ++ )
        {
            dMean_SP[i] = dMean_SP[i] / dDiv1;
            dNY_SP[i] = dNY_SP[i] - dMean_SP[i];

            dMean_RV[i] = dMean_RV[i] / dDiv2;
            dNY_RV[i] = dNY_RV[i] - dMean_RV[i];
        }
        delete [] dMean_SP;
        delete [] dMean_RV;
        
        double dTotal_NYSP = 0.0d, dTotal_NYRV = 0.0d;

        for( i = 0; i < lNrPix; i++ )
        {   
            dTotal_NYSP += (dNY_SP[i] * dNY_SP[i]);
            dTotal_NYRV += (dNY_RV[i] * dNY_RV[i]);
        }

        dNorm_SP = sqrt( fRdist * \
                     (dTotal_NYSP - 0.5 * (dNY_SP[0] * dNY_SP[0] + \
                      dNY_SP[lNrPix-1] * dNY_SP[lNrPix-1] )));

        dNorm_RV = sqrt( fRdist * \
                     (dTotal_NYRV - 0.5 * (dNY_RV[0] * dNY_RV[0] + \
                      dNY_RV[lNrPix-1] * dNY_RV[lNrPix-1] )));
    
        if( iNomc == 0 )
        {
            //cout << "Mean sp = " << 
            cout << "Norm sp = " << dNorm_SP ;
            cout << "  Norm rv = " << dNorm_RV << endl;
        }
    }
    else
    {
        if( iNomc == 0 )
            cout << "Correlation function is not normalized..." << endl;

        dNorm_SP = 1.0;
        dNorm_RV = 1.0;
    }

    if( iNomc == 0 )
        cout << "FFT (RV-standard and spectrum)" << endl;

    // FFT functions to calculate the fourier & inverse fourier transforms
    // need to zero down on which library to use
    // functions calculate the value of fConvY - important
    // TODO -CHANGE-    

    // allocate memory for fConvY
    // fConvY = new float [ lNrPix ];
    // FFT
    
    return bError;
}

void FFT_Prep( double * dNX_SP, float * fYArr, double * dNY_SP , int iSzNYSP)
{
    // Determine the DC- part of spectrum
    // Spectrum is assumed to be equidistant
    // DC = inted( SX, SY ) / (MAX(SX) - MIN(SX))
    // OUT_SY = SY - DC

    for( int i = 0; i < iSzNYSP; i++ )
    {
        dNY_SP[i] = fYArr[i];
    } 

   /*
    Apodisation by multiplication of both Functions with Hanning Function

    Shrikant: Procedure to multiply by Hanning function now commented in 
    original code get_rv.pro
    */
}
