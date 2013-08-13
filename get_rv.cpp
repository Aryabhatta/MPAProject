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
#include <fstream>
#include <string>
#include <fitsio.h>
#include <fftw3.h>
#include "get_rv.hpp"

using namespace std;

// Routine to calculate the radial velocity
float get_rv( float * ObsWave, float * ObsFlux, int iObsElem, float * ThrWave, float * ThrFlux, int iThrElem, float * fXr,  
		      int iNomessage,int iAbsolute, float fEps, int iLog)	
{
    float fXrv = 0.0;
    int i  = 0;
    
    string strVerbose = ReadInput( "ARG:VERBOSE" );
    int iVerbose = atoi( strVerbose.c_str() );
    
    if( iVerbose > 0 )
    cout << endl << "===================== in GETRV===================================" << endl;
    
    if( fEps == 0.0 ) 
    {     
    	fEps = 8 * pow(10,-7);
    	if( iVerbose > 0 )
    	cout << "fEps = " << fEps << endl;
    } // Min Floating accuracy
    
    // iAbsolute is default 0
    // iNomessage is default 0
    
//    fXr[0] = 6457.48;// dummy 
//    fXr[1] = 6500;// dummy
    
    float * fObsW;
    float * fObsF;
    float * fThrW;
    float * fThrF;
    int iObsSz;
    int iThrSz;
    
    if( fXr[0] != 0.00 && fXr[1]!= 0.00 )
    {
    	int IndexArr[ iObsElem ];
	    int iCnt = 0;
	    
    	iCnt=IdlWhere( ObsWave, ">=", fXr[0], "<=", fXr[1], iObsElem, IndexArr);
    
    	// filtered observed spectrum
    	fObsW = new float[ iCnt ];
    	fObsF = new float[ iCnt ];
    	iObsSz = iCnt;
    
    	for( i=0; i< iCnt; i++ )
    	{
    		fObsW[i] = ObsWave[ IndexArr[i] ];
    		fObsF[i] = ObsFlux[ IndexArr[i] ];
    	}
    	
        if( iThrElem != 0 ) // if thereotical flux is given
        {
	    	int IndexArr1[ iThrElem ];
	    	int iCnt1 = 0;
	    
	    	iCnt1=IdlWhere( ThrWave, ">=", fXr[0], "<=", fXr[1], iThrElem, IndexArr1);
	    
	    	// filtered thereotical spectrum
	    	fThrW = new float[ iCnt1 ];
	    	fThrF = new float[ iCnt1 ];
	    	iThrSz = iCnt1;
	
	    	for( i=0; i< iCnt1; i++ )
	    	{	
	    		fThrW[i] = ThrWave[ IndexArr1[i] ];
	    		fThrF[i] = ThrFlux[ IndexArr1[i] ];
	    	}
        }
    }
    else // if fXr not declared
    {
    	fObsW = new float[ iObsElem ];
    	fObsF = new float[ iObsElem ];
    	iObsSz = iObsElem;
    	
    	for( i=0; i< iObsElem; i++)
    	{
    		fObsW[i] = ObsWave[i];
    		fObsF[i] = ObsFlux[i];
    	}
    	
    	if( iThrElem != 0 ) // if thereotical flux is given
        {
    		fThrW = new float[ iThrElem ];
	    	fThrF = new float[ iThrElem ];
	    	iThrSz = iThrElem;
	    	for( i=0; i<iThrElem; i++)
	    	{
	    		fThrW[i] = ThrWave[i];
	    		fThrF[i] = ThrFlux[i];
	    	}
        }
    }
    
	/* 
     * IF NO OF PARAMETERS < 2, LOAD SOLAR FLUX SPECRUM AS REFERENCE
     */       

    if( iThrElem==0 )
    {
    
    	string strInputDir = ReadInput( "DIR:INPUTDIR" );
    	
	    //string strSolarFlux("/home/shrikant/Desktop/MPA/Files/solarflux.fits");
    	string strSolarFlux = strInputDir + "solarflux.fits";
	    
	    cout << endl << "Loading solar flux spectrum from fits file as reference" << endl;
	    
	    fitsfile * fptr;
	    int iStatus = 0; // Initialise status
	    
	    // Opening the FITS file
	    fits_open_file( &fptr, strSolarFlux.data(), READONLY, &iStatus );
	    
	    int naxis1;
	    char * comment = new char [100];
	    
	    // Since we know the format of this file, just reading naxis1 which is
	    // #of cols in the table
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
    
	    //reading a small chunk of 10 elements TODO can change
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
	    
	    // Copy the flux to fThrW array
		iThrSz = iNoData;
		fThrW = new float[ iThrSz ];
		fThrF = new float[ iThrSz ];
    	for( i=0; i<iThrSz; i++)
    	{
    		//fThrW[i] = ThrWave[i];
    		fThrF[i] = fSolarflux[i];
    	}
    	
    	// How to calculate wavelength values here ? TODO

	    fits_close_file( fptr, &iStatus );

	    if( iStatus ) /* print any error messages */
	    {
	        fits_report_error( stderr, iStatus );
	    }

	    delete [] comment;
    }    
        
    float * fConvX, * fConvY;
    bool bError = false;
    int iConvSz = 0;
      
    // Passing fConvX & fConvY as an double pointer since memory is allocated in Log_Lin_Corr
    // Need to pass as double pointer, no other option
    bError = Log_Lin_Corr( fObsW, fObsF, iObsSz, fThrW, fThrF, iThrSz, &fConvX, &fConvY, &iConvSz);
        
    if( bError)
    {
        cout << "Error in LOG_LIN_CORR" << endl;
        cout << "Aborting from get_rv.." << endl;
        return -1;
    }
    
    double dC =299792.458;	// light velocity, km/s
    
    if( iVerbose > 0 )
    {
    	cout << "iCOnvSz:" << iConvSz;
    	cout << "ALL IS WELL !!" << endl;
    }
    
    for( i = 0; i < iConvSz; i++ )
    {
        fConvX[i] = ( exp(fConvX[i]) - 1.00 ) * dC * pow(10,-5);
    }
   
    float fCoef =  0;
    float fConvYA[ iConvSz ];
    
    float fMaxCY = IdlMax( fConvY, iConvSz);
    
    if( iVerbose > 0 )
    cout << "fMaxCY" << fMaxCY << endl;

    if( abs( IdlMin( fConvY, iConvSz ) ) > fMaxCY )
        fMaxCY = abs( IdlMin( fConvY, iConvSz ) );
    
    if( iVerbose > 0 )
    cout << "fMaxCY" << fMaxCY << endl;

    for( i=0; i < iConvSz; i++ )
    {
        fConvYA[i] = fConvY[i] / fMaxCY;
    } 

    float fTemp = IdlMax( fConvY, iConvSz );

    //fXrv = fConvX[ (int)dC ]; // corresponding to rv=cx(!C) in get_rv.cpp
                        // in Idl , ! is used to reference system variabbles
    
    // Above Idl syntax materialises to fConvX[0]
    fXrv = fConvX[ 0 ]; // corresponding to rv=cx(!C) in get_rv.cpp
    
    if( iVerbose > 0 )
    cout << "fXrv" << fXrv << endl;
    
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
    iSt=2;
    if( iSt == 1 )
    {
        fXrv = fXe;
        fCoef = fYe;        
    }
    
    if( iVerbose > 0 )
    cout << endl << "fXrv:" << fXrv << endl;
    
    // delete memory allocated in Log_Lin_Corr
    delete [] fConvX;
    delete [] fConvY;
    
    delete [] fObsW;
	delete [] fObsF;
	delete [] fThrW;
	delete [] fThrF;

    return fXrv;
}

bool Log_Lin_Corr( float * fObsWave, float * fObsF, int iObsSz, float * fThrWave, float * fThrF, int iThrSz, 
		           float ** fConvX, float ** fConvY, int * iConvSz, bool bComplete, int iNomc, float eps)
{
	// LOG-LIN-CONVOLUTION: X -> LOG(X) -> RESAMPLING OF (LOG(X),Y) -> SCALE FOR RV-DET
	
	string strVerbose = ReadInput( "ARG:VERBOSE" );
	int iVerbose = atoi( strVerbose.c_str() );
	
	if( iVerbose > 0 )
	{
		cout << "In log lin corr" << endl;
		cout << "ObsSz: " << iObsSz << " ThrSz: " << iThrSz << endl;
	}

    // defaults for bComplete, iNomc & eps set in header file
    int i =0;
    bool bError = false;
    
    // containing log of values from Observed & thereotical wavelengths
    float fObsW[iObsSz], fThrW[iThrSz];
    
    // Taking log
    for( i = 0; i < iObsSz; i++ )
        fObsW[i] = log( fObsWave[i] );
    for( i = 0; i < iThrSz; i++ )
        fThrW[i] = log( fThrWave[i] );
   
    // Spectrum must have same widths
    // Algorithm for that
    float fOW_start = fObsW[0];
    float fOW_last = fObsW[iObsSz-1];
    float fTW_start = fThrW[0];
    float fTW_last = fThrW[iThrSz-1];
    float fOWdiff = fOW_last - fOW_start;
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
    float fRdistOW = 10000000, fRdistTW = 10000000; // initial max values
    
    // calculate min resampling distance for observed spectra
    for( i=1; i < iObsSz; i++ )
    {
        if( fRdistOW > (fObsW[i] - fObsW[i-1]))
        {
            fRdistOW = fObsW[i] - fObsW[i-1];
        }
    }

    for( i=1; i < iThrSz; i++ )
    {
        if( fRdistTW > (fThrW[i] - fThrW[i-1]))
        {
            fRdistTW = fThrW[i] - fThrW[i-1];
        }
    }
    
    float fRdist = std::max( fRdistOW, fRdistTW ) * 0.5f ;
    
    // if resampling distance < than eps, eps taken as new resampling distance
    fRdist = std::max(fRdist, eps);
    
    if( iVerbose > 0 )
    cout << "fRdistOW:" << fRdistOW << " fRdistTW:" << fRdistTW << " max:" << fRdist << endl;

    long int lNrPix = abs((fOW_last - fOW_start)/fRdist) + 1L;
    long int lNrPixChk = abs((fTW_last - fTW_start)/fRdist) + 1L;
    
    if( iVerbose > 0 )
    cout << "lNrPix:" <<lNrPix << " lNrPixChk:" <<lNrPixChk << endl;
    
    if( lNrPix != lNrPixChk )
    {
        cout << "PIX.NR(SPECTRUM) - PIX.NR(RV-STD): " <<(lNrPix - lNrPixChk);
        cout << endl;
        if( lNrPix > lNrPixChk )
            lNrPix = lNrPixChk;
    }
    
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
    // same size as fObsWave & fThrWave // Assuming, recheck TODO
    if( iNomc == 0 )
    {
        cout << endl << "Points for FFT: " << lNrPix << endl;
        
        if( iVerbose > 0 )
        cout << "Interpolation:" << endl;
    } 
    
    ofstream logFile;
    
    // log for plotting - 2nd argument determines to log or not
    createLog( "interpolb4.log", false, fThrW, iObsSz, fThrF, iObsSz );
    
    float fYAarr[lNrPix], fYBarr[lNrPix];
      
    // calling interpolation
    // Input - flux, log of wavelengths(input absicca), dNX_**( output absicca)
    // Result to be stored in fY*arr
    interpol( fObsF, fObsW, iObsSz, dNX_SP, fYAarr, lNrPix );
    interpol( fThrF, fThrW, iThrSz, dNX_RV, fYBarr, lNrPix );
    
    string stra4 = ReadInput("DIR:LOGDIR") + "interpola4.log";
    if( false )
    {
	    logFile.open( stra4.data(), ios::out );
	    for( i=0; i< lNrPix; i++)
	    {
	    	logFile << dNX_RV[i] << "\t" << fYBarr[i] << endl;
	    }
	    logFile.close();
	    logFile.clear();
    }

    double dNY_SP[ lNrPix ];
    double dNY_RV[ lNrPix ];
    
    // In current implementation, just copies the fYarr to dNY_SP
    FFT_Prep( dNX_SP, fYAarr, dNY_SP, lNrPix);
    FFT_Prep( dNX_RV, fYBarr, dNY_RV, lNrPix);

    // CORRELATION
    long lShift = lNrPix/2;
    *fConvX = new float [ lNrPix ];
    *iConvSz = lNrPix;
    
    if( iVerbose > 0 )
    cout << "lShift" << lShift << endl;

    double dMaxNxrv = IdlMax( dNX_RV, lNrPix );
    double dTemp = (dMaxNxrv - dNX_RV[0])/2.0d + dNX_RV[0];
    
    if( iVerbose > 0 )
    cout << "dMaxNxrv:" << dMaxNxrv << " dTemp:" << dTemp << " dNX_RV[0]: " << dNX_RV[0] << endl;

    // ConvX
    for( i = 0; i < lNrPix; i++ )
    {
        *(*fConvX + i) = dNX_SP[i] - dTemp;
    }

    double dNorm_SP = 0.0, dNorm_RV = 0.0;
    if( bComplete )
    {
        if( iNomc == 0 )
            cout << "Normalization of correlation function" << endl << endl;
    
        double * dMean_SP = new double [ lNrPix ];
        double * dMean_RV = new double [ lNrPix ];
        
        integ( dNX_SP, dNY_SP, lNrPix, dMean_SP);
        integ( dNX_RV, dNY_RV, lNrPix, dMean_RV);
        
        if( iVerbose > 0 )
        cout << "First elem: " << dMean_SP[0] << " Last element:" << dMean_SP[lNrPix-1] << endl; 
        
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
    
        if( iVerbose > 0 )
        {	
            cout << " Norm sp = " << dNorm_SP ;
            cout << " Norm rv = " << dNorm_RV << endl;
        }
    }
    else
    {
        if( iNomc == 0 )
            cout << "Correlation function is not normalized..." << endl;

        dNorm_SP = 1.0;
        dNorm_RV = 1.0;
    }

    if( iVerbose > 0 )
        cout << "FFT (RV-standard and spectrum)" << endl;

    // FFT functions to calculate the fourier & inverse fourier transforms
    // need to zero down on which library to use
    // functions calculate the value of fConvY - important
    
    if( iVerbose > 0 )
    cout << "lnrpix" << lNrPix << endl;
    
    // prepare arrays for FFT
    for( i=0; i < lNrPix; i++ )
    {
    	dNY_RV[i] = dNY_RV[i]/ dNorm_RV;
    	dNY_SP[i] = dNY_SP[i]/ dNorm_SP;
    }
    
    if( false )
    {
//    	string strNyrv("/home/shrikant/Desktop/MPA/Log/RV.log");
//    	string strNysp("/home/shrikant/Desktop/MPA/Log/SP.log");
    	
    	string strNyrv = ReadInput("DIR:LOGDIR") + "RV.log";
    	string strNysp = ReadInput("DIR:LOGDIR") + "SP.log";
    	
    	logFile.open(strNysp.data(), ios::out );
    	for( i=0; i<lNrPix; i++)
    	{
    		logFile << dNX_SP[i] << "\t" << dNY_SP[i] << endl;
    	}
    	logFile.close();
    	logFile.clear();
    	
    	logFile.open(strNyrv.data(), ios::out );
    	for( i=0; i<lNrPix; i++)
    	{
    		logFile << dNX_RV[i] << "\t" << dNY_RV[i] << endl;
    	}
    	logFile.close();
    	logFile.clear();
    }
    
    // Taking the FFTW transforms
    fftw_complex *in, *out;
    fftw_plan p;
    
    in  = (fftw_complex*)malloc(sizeof(fftw_complex) * lNrPix);
    out = (fftw_complex*)malloc(sizeof(fftw_complex) * lNrPix);
    
    // Copy NY_SP to in
    for( i=0; i< lNrPix; i++)
    {
    	// just copying the real part
    	in[i][0] = dNY_RV[i];
    }
    
    // take a Inverse Fourier Transform of NY_RV
    p = fftw_plan_dft_1d( lNrPix, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    
    fftw_destroy_plan(p);
    
//    // print
//    cout << endl << "Printing dNY_RV" << endl;
//    for( i=0; i< lNrPix; i++)
//    {
//    	cout << dNY_RV[i] << " " ;
//    }
//    
//    cout << endl << "Print inverse transform of dNY_RV" << endl;
//    for( i=0; i< lNrPix; i++)
//    {
//    	cout << out[i][0] << "," << out[i][1] << " " ;
//    }
    
         
    // take a Forward Fourier Transform of NY_SP
    fftw_complex *in1, *out1, *temp, *out2;
    fftw_plan p1, p2;
    
    in1  = (fftw_complex*)malloc(sizeof(fftw_complex) * lNrPix);
    out1 = (fftw_complex*)malloc(sizeof(fftw_complex) * lNrPix);
    temp = (fftw_complex*)malloc(sizeof(fftw_complex) * lNrPix);
    out2 = (fftw_complex*)malloc(sizeof(fftw_complex) * lNrPix);
    
    // Copy NY_SP to in
    for( i=0; i< lNrPix; i++)
    {
    	// just copying the real part
    	in1[i][0] = dNY_SP[i];
    }
    
    p1 = fftw_plan_dft_1d( lNrPix, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p1);    
    fftw_destroy_plan(p1);
    
    // Multiply & combine in single temp array
    for( i=0; i< lNrPix; i++)
    {
    	// Manually copying the real & complex part, Not sure however if 
    	// this is the right way to do so
    	temp[i][0] = out[i][0] * out1[i][0];
    	temp[i][1] = out[i][1] * out1[i][1];
    }
    
    // take inverse fourier transform of temp array
    p2 = fftw_plan_dft_1d( lNrPix, temp, out2, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p2);    
    fftw_destroy_plan(p2);
    
    // Now out2 contains the final transformed array which needs to be shifted & stored in ConvY
    *fConvY = new float[ lNrPix ];
    float * fConvYTemp = new float[ lNrPix];
    
    for( i=0; i< lNrPix; i++ )
    {
    	// Just copying the real part of result
    	// Need to cross check if that's the RIGHTA way
    	fConvYTemp[i] = out2[i][0];
    }
    
    // Unallocate fftw arrays
    fftw_free(in);
    fftw_free(out);
    fftw_free(in1);
    fftw_free(out1);
    fftw_free(temp);
    
    // Shift by lshift & store in ConvY //CIRCULAR SHIFT
    for( i=0; i<lNrPix; i++)
    {
    	*(*fConvY + i) = fConvYTemp[(i+lShift)%lNrPix];
    }
    
    // deleting the memory allocated for temporary array for shifts
    delete [] fConvYTemp;
    
    // Update ConvY depending on ConvX
    for( i=0; i< lNrPix; i++)
    {
    	*(*fConvY + i) = (( *(*fConvX + (lNrPix-1))-*(*fConvX + 0) )/lNrPix) * ( *(*fConvY + i) );
    }
    
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
