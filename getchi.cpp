/*************************************************************************
*
* Copyright:		Max Planck Institute for Astrophysics (MPA)
* 
* File:				getchi.cpp
*
* Routine Info:		MC generator using grids of pre-computed synthetic
*       			spectra , evalutation of chi square
*
* Author:
*
* Miscellaneous:	We follow Hungarian notation for variable naming
*
* Modification 
* Log:	    		Added routine main
*		        	Shrikant V	Dec 11, 2012	16:20
*
**************************************************************************/
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <fitsio.h>
#include <unistd.h>
#include "lfp.hpp"
#include "get_rv.hpp"
#include "cont_rscl.hpp"
#include "ecorr.hpp"
#include "modules.hpp"

using namespace std;

/*------------------------------------------------------------------------
ROUTINE
	getchi.cpp
INPUTS
	fin: input file with grid points: T,g,m (.def file)
	fid: fiber number						( in {1,640})
	
OUTPUT

-------------------------------------------------------------------------*/

int main(int iArgCnt, char * pArg[])
{

//int 	iFid  = 500;
string Fid = ReadInput( "ARG:FID"); 
int iFid = atoi(Fid.c_str());

int	    iConv = 1;
int 	iNomc = 0;

float	fRes  = 2000.0;
float	fPlot_it = 0.0;
float 	fNumplot = 1.0;

float	fLimCont = 1.01;	// error in spectrum interpolation

int 	i = 0; // iteration counter

enum eMode {Halpha, Mgb, CaT, Hbeta}; // any update here, pls also 
                                      // update iNoModes below
int iNoModes = 4;

// Directory paths
//string strSpdir( "/afs/mpa/project/stellar/mbergema/sspp/output1/" );
//string strSpobs( "/afs/mps/project/stellar/mbergema/sspp/"         );
//string strLocal( "/afs/mpa/data/mbergema/SIU/spec/output1"         );

//string strSpecs( "spPlate-1960-53287.fits"                         );
string strSpecs = ReadInput( "FILE:FITSFILE");

//string strFin( "Input_File_Name" );
string strFin = ReadInput("FILE:GRIDPOINTS"); // file of grid points Rlogz, Rteff, Rlogg

/*
 * Here in all References:
 * Rlogz - determines the metallicity of star
 * Rteff - effective temperature of star
 * Rlogg - specific gravity of star
 */									

string strLog( strSpecs );

string strTemp( "_" );
strTemp.append( Fid ); 
strTemp.append( ".log" );

// Creatung a name for a log file
strLog.replace( strLog.find(".fits"), 5, strTemp);

string strGrid( "noconv" );	// explicit convolution

float fRgauss = 	5;
float fRgamma = 	0;

int iNomessage= 1;		// 1 = skip all warnings
int iUse_rf   =	1;		// 1 = reliability function should be used
int iUse_cont_rscl = 1;		// if continnum of observed spectrum is to
			            	// be optimised

double dCCC = 299792.458;	// light velocity, km/s

ofstream logFile; 	    // for log file
ifstream inputFile; 	// for input file

/************************************************
 * Changing path for log directory:
 * Go to ProgramInputs.txt. Change the path in
 * section #DIR against LOGDIR
 ***********************************************/
//string strFileSpecs( "/home/shrikant/Desktop/MPA/Log/" );
string strLogDir = ReadInput( "DIR:LOGDIR");
string strFileSpecs = strLogDir;

// Open file for logging
if ( !strLog.empty() )
{
    strFileSpecs.append( strLog );

    logFile.open( strFileSpecs.data() , ios::out | ios::app );

    if( ! logFile.is_open() )
    {
        cout << "The path to logfile is wrong !! Aborting !!! ";
        return 0;
    }

    // Status Msg
    cout << endl << "Logfile = " << strFileSpecs << endl; 
    
    // logFile << "This is to check whether \n the code is \n working";
}

// get time to calculate the program runtime.
timeval time1;
gettimeofday( &time1, NULL );

int iValid = 0;

// Input file specification
string strInputDir = ReadInput("DIR:INPUTDIR");

/************************************************
 * Changing path for input directory:
 * Go to ProgramInputs.txt. Change the path in
 * section #DIR against INPUTDIR 
 ***********************************************/
//string strInputFilePath( "/home/shrikant/Desktop/MPA/Files/foutr1p1082f180" ); // sample file name
string strInputFilePath = strInputDir + strFin; // contains list of stars with Rlogz, Rteff, Rlogg

// Status Msg
cout << "Reading grid points from " << strInputFilePath << endl;

// Open input file
inputFile.open( strInputFilePath.data(), ifstream::in );

if( ! inputFile.is_open() )
{
   cout << "Wrong path to Input file !!! Aborting !!!";
   return 0;
}

//
// Set microturbulence
//
int iT0 = 5500 ;
float fG0 = 4.0 ;
float fRxi = 2.0 ;
iValid = 1;

int iCnt = 0;
int iRf = 0; // reliablity function
   
string strRange("");
string strWcen("");

int arrWr [2];
int arrWrc[2];

/*
 * The FITS file contains the observations for a particular plate, MJD such that we can read 
 * the flux from the file corresponding to 640 fibers & 3857 wavelengths for each of them
 */

/************************************************
 * Changing path for input directory:
 * Go to ProgramInputs.txt. Change the path in
 * section #DIR against INPUTDIR 
 ***********************************************/
//string strFitsFilePath( "/home/shrikant/Desktop/MPA/Files/spPlate-1962-53321.fits" ); 
string strFitsFilePath = strInputDir + strSpecs;

// Status Msg
cout << "Reading from FITS file " << strFitsFilePath << endl;

// ----------------------------------------------------------------
// Calling the below routines using CFITSIO library
// Output Wavelengths array and data array correspondin to fid
// ----------------------------------------------------------------

fitsfile *fptr;
int iStatus = 0;    // Must initialise 'status'

// Opening the FITS file
fits_open_file( &fptr, strFitsFilePath.data(), READONLY, &iStatus );

int naxis1,naxis2;
char * comment = new char [100];

// Read no of coloumns in the table
fits_read_key( fptr, TINT, "NAXIS1" , &naxis1, comment, &iStatus);
cout << endl << "Naxis1 = #cols = " << naxis1 << "  Comment = " << \
comment<< endl;

// Read no of rows in the table
fits_read_key( fptr, TINT, "NAXIS2" , &naxis2, comment, &iStatus);
cout << "Naxis2 = #rows = " << naxis2 << "  Comment = " << \
comment<< endl;

float fCrval1, fCoeff1;

// Read CRVAL1
fits_read_key( fptr, TFLOAT, "CRVAL1" , &fCrval1, comment, &iStatus);
cout << "Crval1 = " << fCrval1 << "  Comment = " << \
comment<< endl;

// Read COEFF1
fits_read_key( fptr, TFLOAT, "COEFF1" , &fCoeff1, comment, &iStatus);
cout << "Coeff1 = " << fCoeff1 << "  Comment = " << \
comment<< endl << endl;

// * creating Wavelength array *
float * fWavelengths = new float [ naxis1 ];

int iCntr1 = 0;

// cases on exceptional values of dCrval1 and dCoeff1 not handled
// Filling the wavelengths array with wavelengths (corresponding to CRVAL1, COEFF1)
while( iCntr1 < naxis1 )
{
    fWavelengths[ iCntr1 ] = powf( 10.0, fCrval1 );
    fCrval1 += fCoeff1;
    iCntr1++;
}

// Status
cout << "Wavelength array created from values of "; 
cout << "CRVAL1, COEFF1 & NAXIS1 !" << endl;

// Reading a row of image/table/data  
long firstelem = 0, nelements = 0;
float nullval = 0;
float * data = new float[ naxis1 ];
int anynull;

// Reading data from file corresponding to particular fiber for which we read Rlogz,RTeff, Rlogg 
// Status Msg
cout << "Reading Data for FID = " << iFid << endl; 

firstelem = ( ( iFid - 1 ) * naxis1 ) + 1;
nelements = naxis1;

// routine to read the row corresponding to iFid (values in data)
fits_read_img( fptr, TFLOAT, firstelem, nelements, &nullval, data, &anynull, &iStatus);

// Status
cout << "Row corresponding to FID read !" << endl;

fits_close_file( fptr, &iStatus );

if( iStatus ) /* print any error messages */
{
    fits_report_error( stderr, iStatus );
}

delete [] comment;

//------------------------------------------------------------------
// Can convert above logic from CFITSIO to CCFITS as it needs to
// be C++ compatible supporting object orinted functionality
// but no need as it fulfills the functionality
//------------------------------------------------------------------

// Loop for number of elements in mode (checking for every possible mode listed)
// Procedure - choosing a range & then processing it 
//for( int iCntr = 0; iCntr < iNoModes; iCntr++ ) /* TODO:Uncomment this line & comment next line to run the program for all modes */
for( int iCntr = 0; iCntr < 1; iCntr++ )
{
    iCnt = 0;	// re-initialized for each wav segment
    iRf = 0;	// reliability function

    switch( iCntr )
    {
        case Mgb:   strRange = "Mgb";
                    strWcen = "W5200";
                    arrWr[0] = 5000;
                    arrWr[1] = 5350;
                    arrWrc[0] = arrWr[0];
                    arrWrc[1] = arrWr[1];
                    break;
                    
        case CaT:   strRange = "CaT";
                    strWcen = "W8600";
                    arrWr[0] = 8400;
                    arrWr[1] = 8800;
                    arrWrc[0] = arrWr[0];
                    arrWrc[1] = arrWr[1];
                    break;
                    
        case Halpha:strRange = "Halpha";
                    strWcen = "W6520";
                    arrWr[0] = 6400;
                    arrWr[1] = 6640;
                    arrWrc[0] = arrWr[0];
                    arrWrc[1] = arrWr[1];
                    break;
                    
        case Hbeta: strRange = "Hbeta";
                    strWcen = "W4850";
                    arrWr[0] = 4600;
                    arrWr[1] = 5000;
                    arrWrc[0] = arrWr[0];
                    arrWrc[1] = arrWr[1];
                    break;
                    
        default:    break;
    }
       
    // Status
    cout << endl << "Range = " << strRange << "\t";
    cout << "Wcen = " << strWcen << "\t";
    cout << "wr = [ " << arrWr[0] << " , " << arrWr[1] << " ]" << endl;
    
    // Status
    cout << "Sampling the only wavelengths that lie between ";
    cout << arrWr[0] << " and " << arrWr[1] << endl;

    // Sampling the wavelengths that lie between arrWr[0] and arrWr[1]
    int iElements = 0;
    for( iCntr1 = 0; iCntr1 < naxis1; iCntr1++ )
    {
        if( fWavelengths[ iCntr1 ] > arrWr[0] &&
            fWavelengths[ iCntr1 ] < arrWr[1] )
        {
            iElements++;    // getting the count for # valid elements
        }

        // since elements in fWavelengths are sorted, we can break once-
        // -it is greater than arrWr[1]
        if( fWavelengths[ iCntr1 ] > arrWr[1] )
        {
            break;
        }
    }

    float * fWavelen = new float [ iElements ];
    float * fData = new float [ iElements ];

    iElements = 0;
    
    // getting relevant data by chopping unnecessary data
    // copying only wavelength & flux that lies in the given range
    for( iCntr1 = 0; iCntr1 < naxis1; iCntr1++ )
    {
        if( fWavelengths[ iCntr1 ] > arrWr[0] &&
            fWavelengths[ iCntr1 ] < arrWr[1] )
        {
            fWavelen[ iElements ] = fWavelengths[ iCntr1 ];
            fData[ iElements ] = data[ iCntr1 ];
        
            iElements++;
        }

        // since elements in fWavelengths are sorted, we can break once
        // it is greater than arrWr[1]
        if( fWavelengths[ iCntr1 ] > arrWr[1] )
        {
            break;
        }
    }

    // Status
    cout << "No of observed data points: " << iElements << endl;

    // Implementing the median filtering part (Smoothing the flux)
    maxMedianFilter( fData, iElements, 3 );

    // Extreme initial values for DataMin and DataMax
    float fDataMin = 100000000, fDataMax = -1000000;

    // Integrating logic of finding fDataMin & fDataMax
    for( iCntr1 = 0; iCntr1 < iElements; iCntr1++ )
    {
        if( fData[ iCntr1 ] > fDataMax )
        {
            fDataMax = fData[ iCntr1 ];  // finding max // max flux
        }

        if( fData[ iCntr1 ] < fDataMin )
        {
            fDataMin = fData[ iCntr1 ];  // finding min // min flux
        }
    }

    // Status
    cout << "DataMax = " << fDataMax << "\t" ;
    cout << "DataMin = " << fDataMin << endl;

    // Byte array for masking   
    unsigned char * arrMask = new unsigned char[ iElements ];
    
    for( iCntr1=0; iCntr1 < iElements; iCntr1++)
    {
    	arrMask[iCntr1] = '1';
    }

    // For every range, mask the unwanted wavelengths with NULL
    // corresponding to IF mode(m) eq 'Halpha' then begin
    if( iCntr == Halpha )
    {
        for( int iC = 0; iC < iElements; iC++ )
        {
            if( fWavelen[ iC ] < 6480.0     ||
                fWavelen[ iC ] > 6600.0     ||
              ( fWavelen[ iC ] > 6510.0 && fWavelen[ iC ] < 6543.0 ) )
            {
                arrMask[ iC ] = '\0'; // NULL
            } 
        }
    }
    else if( iCntr == CaT )
    {
        for( int iC = 0; iC < iElements; iC++ )
        {
            if( fWavelen[ iC ] < 8480.0     ||
                fWavelen[ iC ] > 8680.0     ||
              ( fWavelen[ iC ] > 8510.0 && fWavelen[ iC ] < 8535.0 )||                        
              ( fWavelen[ iC ] > 8550.0 && fWavelen[ iC ] < 8640.0 ) )
            {
                arrMask[ iC ] = '\0'; //  NULL
            } 
        }
    }
    else if( iCntr == Mgb )
    {
        for( int iC = 0; iC < iElements; iC++ )
        {
            if( fWavelen[ iC ] < 5115.0     ||
                fWavelen[ iC ] > 5285.0     ||
              ( fWavelen[ iC ] > 5150.0 && fWavelen[ iC ] < 5160.0 )||                        
              ( fWavelen[ iC ] > 5195.0 && fWavelen[ iC ] < 5215.0 )||                        
              ( fWavelen[ iC ] > 5235.0 && fWavelen[ iC ] < 5258.0) )
            {
                arrMask[ iC ] = '\0'; // NULL
            } 
        }
    }
    
    // TODO Hbeta missing here in the above logic

    float fGauss = fRgauss;
    float fGamma = fRgamma;

    float fRlogz = 0.0 ;
    float fRteff = 0.0 ;
    float fRlogg = 0.0 ;

    // Make inputfile point to beginning of file for every mode 
    inputFile.clear();
    inputFile.seekg(0,inputFile.beg);

    float fXrv = 0.0;
    float fXr[2] = { 0.0, 0.0 };

    // Variables for reading theoretical spectrum
    float * fSx, * fSy;
    int iThrElem;

    // read till end of file is reached     
    while( ! inputFile.eof() ) // loop that reads no of points
    {	
        string strRow;
        getline( inputFile, strRow );   
         
        // parsing the string to give float values for Rlogz, Rteff, Rlogg
        if( ! strRow.empty())
        {
            // tokenising strRow
            char cDelim = ' ';

            strTemp.clear();
            strTemp = strNexttoken( strRow, cDelim );
            fRlogz = strtof( strTemp.c_str(), NULL );

            strTemp.clear();
            strTemp = strNexttoken( strRow, cDelim );
            fRteff = strtof( strTemp.c_str(), NULL );

            strTemp.clear();
            strTemp = strNexttoken( strRow, cDelim );
            fRlogg = strtof( strTemp.c_str(), NULL );

            // Status Msg
            cout << "Rlogz = " << fRlogz << "\t";
            cout << "Rteff = " << fRteff << "\t";
            cout << "Rlogg = " << fRlogg << endl;
        }
        else // EOF reached
        {
            break;
        }
    
        float fTeff = 0.0, fLogg = 0.0, fLogz = 0.0;

        fTeff = fRteff;
        fLogg = fRlogg;
        fLogz = fRlogz;      

        float fXi = 0;
        
        // Calculate fXi depending on values of fTeff,fLogg
        if( fTeff > 5250.0 && fLogg >= 3.5 ) // MS and RGB
                                             // Teff >= 5250, logg >= 3.5
        {
            fXi = 1.15 + ( 2.0 * pow10( -4.0 ) ) * ( fTeff - iT0 ) + \
                  ( 3.95 * pow10( -7 ) ) * powf( (fTeff - iT0), 2 ) - \
                  0.13 * ( fLogg - fG0 ) + 0.13 * powf( (fLogg - fG0), 2 );
        }
        else if( fTeff < 5250.0 && fLogg >= 3.5 ) // MS, Teff <= 5250
        {
            fXi = 1.15 + ( 2.0 * pow10( -4.0 ) ) * ( fTeff - iT0 ) + \
                  ( 3.95 * pow10( -7 ) ) * powf( (fTeff - iT0), 2 ) - \
                  0.13 * ( fLogg - fG0 ) + 0.13 * powf( (fLogg - fG0), 2 );
        }
        else if( fTeff <= 5250.0 )     // RGB / AGB
        {
            fXi = 0.94 + (2.2 * pow10(-5)) * ( fTeff - iT0)  \
                  - ( 0.5 * pow10( -7 )) * powf( ( fTeff - iT0 ), 2 ) \
                  - 0.1 * ( fLogg - fG0 ) + 0.04 * powf( ( fLogg - fG0 ), 2)\
                  - 0.37 * fLogz - 0.07 * powf( fLogz, 2 );
        }
        
        float fMg_eps = 0.0;
        if( fLogz < -0.6 )
        {
            fMg_eps = 0.4f;
        }

        bool bError;
		float fEps_dev[2];
		
		fEps_dev[0] = 12.00;
		fEps_dev[1] = fMg_eps;

		int iCnvl = 0, iExtrapol = 0;

        strUpper( strRange );
        
          // Thereotical Spectra
         /*
          * Uncomment below region of code if you want to skip lfp
          * by making thereotical spectra = observed spectra
          */        

/*
        fSx = new float[ iElements ];
        fSy = new float[ iElements ];
        iThrElem = iElements;
        
        // Need to change when LFP starts working perfectly
        for( iCntr1 = 0; iCntr1 < iElements; iCntr1 ++ )
        {
            fSx[ iCntr1 ] = fWavelen[ iCntr1 ];
            fSy[ iCntr1 ] = fData[ iCntr1];
        }
        bError = false;
*/
        
        // PLOT Observed Spectrum
        int iPlot = 0; // change iPlot to 1 if you want to plot
        if( iPlot == 1)
        {
			ofstream logFile;
			// Logs the observed spectra at standard log directory
			//string strObsSpectra("/home/shrikant/Desktop/MPA/Log/ObsSpectrum.log");
			string strObsSpectra = strLogDir + "ObsSpectrum.log";
			
			logFile.open( strObsSpectra.data(), ios::out );
			for( int i=0; i< iElements; i++)
			{
				logFile << fWavelen[i] << "\t" << fData[i] << endl;
			}
			logFile.close();
			logFile.clear();        	
        }        

        /**************************************************************************
         * !!! Important !!!
         * If you are commenting/uncommenting LFP, remember to reflect the changes
         * in the intialisation of Sx , Sy above !!!
         * ***********************************************************************/
        // Call to lfp (with thereotical spectra & values of Teff, Logg, Logz & fXi         		
		bError = lfp( &fSx, &fSy, &iThrElem, fTeff, fLogg, fLogz, fXi, fEps_dev, fGauss, fGamma, iExtrapol, strGrid, strRange, iNomessage, iCnvl, iNomc );
		
		iPlot = 1;
		if( iPlot == 1 )
		{	
			ofstream logFile;
			
			// Logging thereotical spectra in log directory
			
			//string strThrSpectra("/home/shrikant/Desktop/MPA/Log/thrSpectrum.log");
			string strThrSpectra = strLogDir + "thrSpectrum.log";
			
			logFile.open( strThrSpectra.data(), ios::out );
			for( int i=0; i< iThrElem; i++)
			{
				logFile << fSx[i] << "\t" << fSy[i] << endl;
			}
			logFile.close();
			logFile.clear();
		}
		
		if( bError || iThrElem == 0 ) 
        {
            // goto NEXT iteration
            iCnt++; // iCnt initialised at starting of loop over no of modes
            continue;
        }

        // max value in thereotical spectrum
        float fSyMax = 0.0;
        fSyMax = IdlMax<float>( fSy, iThrElem);
        
        int iFlag;
        if( fSyMax > fLimCont )
        {
            iFlag = 0;
        }
        else
        {
            iFlag = 1;
        }
        
        iPlot = 0;
        if( iPlot == 1 )
        {
        	logFile.close();
        	logFile.clear();
        	string strFileSpecsb4 = strFileSpecs + "before.log";
        	logFile.open( strFileSpecsb4.data() , ios::out );
        	logFile << "Before convolution" << endl;
        	for(iCntr1=0;iCntr1<iElements;iCntr1++)
        	{
        		logFile << fSx[iCntr1] << "\t" << fSy[iCntr1] << endl;
        	}
        	logFile.close();
        	logFile.clear();
        }
        
        if( iConv )
        {
            gaussFold( fSx, fSy, iElements, ( IdlMean(arrWr,2)/fRes) );
        }
        
         iPlot = 0;
         if( iPlot == 1 )
         {
        	string strFileSpecsa4 = strFileSpecs + "after.log";
            logFile.open( strFileSpecsa4.data()  , ios::out );
            logFile << "After convolution" << endl;
        	for(iCntr1=0;iCntr1<iElements;iCntr1++)
        	{
        		logFile << fSx[iCntr1] << "\t" << fSy[iCntr1] << endl;
        	}
        	logFile.close();
        	logFile.clear();
            logFile.open( strFileSpecs.data() , ios::out  );
        }
         
        float  fRy[iElements];
        
        //
        //  Initial estimate of radial velocity fXrv
        //
        if( iCnt == 0 || fXrv == 0)
        {
            fXrv = get_rv( fWavelen, fData, iElements, fSx, fSy, iThrElem, fXr,  iNomessage);
        
            cout << "Radial Velocity in getchi: " << fXrv << endl;
	        if( abs( fXrv ) < 1000.0 )
	        {
	            double fDiv = double ( fXrv ) / dCCC  + 1.0;
	            for( iCntr1 = 0; iCntr1 < iElements; iCntr1++)
	            {
	                fWavelen[ iCntr1 ] =  fWavelen[ iCntr1 ] / fDiv;
	            }
	        }
	        else
	        {
	            fXrv = 0.0;
	        }

		    //
		    // Rescale continnum
		    //		    
		    if( iUse_cont_rscl > 0 )
		    {
		        cont_rscl( fWavelen, fData, iElements, fSx, fSy, iThrElem, fRy);        
		    }
		    else
		    {
		    	for( iCntr1 =0; iCntr1<iElements; i++ )
		    	{
		    		fRy[iCntr1] = 1.0f; 
		    	}
		    }
	    
		    //
		    // Final estimate of radial velocity
		    //
		    float fXrv2 = 0.0f;
	
	        int iIndex1[ iElements ], iIndex2[ iThrElem];
	        
	        int iCount1 = IdlWhere( fWavelen, ">=", (float) arrWrc[0], "<=", (float)arrWrc[1], iElements, iIndex1 );
	        int iCount2 = IdlWhere( fSx, ">=", (float) arrWrc[0], "<=", (float) arrWrc[1], iThrElem, iIndex2 );
	
		    // Creating temporary array to contain values of fWavelen, fData...
		    // according to indices iii1 = iCount1 & iii2 = iCount2
		    float tmp_fWavelen[iCount1];
		    float tmp_fDataRy[iCount1];
		    
		    float tmp_fSx[iCount2];
		    float tmp_fSy[iCount2];
		    
		    for( i=0; i<iCount1; i++)
		    {
		    	tmp_fWavelen[i] = fWavelen[iIndex1[i]];
		    	tmp_fDataRy[i] = fData[iIndex1[i]] * fRy[iIndex1[i]];
		    }
		    for( i=0; i< iCount2; i++)
		    {
		    	tmp_fSx[i] = fSx[iIndex2[i]];
		    	tmp_fSy[i] = fSy[iIndex2[i]];
		    }
		    
		    // Final estimate of Radial Velocity
		    fXrv2 = get_rv( tmp_fWavelen, tmp_fDataRy, iCount1, tmp_fSx, tmp_fSy, iCount2, fXr, iNomessage);
	    
		    if( abs( fXrv2 ) < 1000.0 )
	        {
	            double fDiv = double ( fXrv2 ) / dCCC  + 1.0;
	            for( iCntr1 = 0; iCntr1 < iElements; iCntr1++)
	            {
	                fWavelen[ iCntr1 ] =  fWavelen[ iCntr1 ] / fDiv;
	            }
	        }
	        else
	        {
	            fXrv2 = 0.0;
	        }
	    
		    // Adding the velocities
		    fXrv = fXrv + fXrv2;
		    
		    cout << endl << "Final estimate of radial velocity in km/s: " << fXrv << endl;
	    
		    // converting float to string
		    std::ostringstream ss;
		    ss << fXrv;
		    string strfXrv( ss.str());
		    
		    string Text1 = "RV: ";
		    Text1.append(strfXrv);
		    Text1.append(" km/s");
        }
	    
	    //
	    // xsig -- number of flux points to be considered: +-1 sigma in the error
	    //
	    double * dCoef = 0;
	    
	    // temporary array to hold y * ry = fData * fRy
	    float fDataRy[iElements];
	    for( i=0; i< iElements; i++)
	    {
	    	fDataRy[i] = fData[i] * fRy[i];
	    }
	    
	    float * fRf = (float *)malloc(sizeof(float) * iElements );
	    for( i=0; i< iElements; i++)
	    {
	    	fRf[i] = 1.00f;
	    }
	    
	    int iXsigma = 1;
	    
	    dCoef = ecorr( fSx, fSy, iThrElem, fWavelen, fDataRy, iElements, iXsigma, fRf, iElements, arrMask, iElements);
	    
	    delete [] fRf;
	    
	    // converting double to string
	    std::ostringstream ss;
	    ss << dCoef;
	    string strdCoef( ss.str());
	    
	    string Text2 = "cf: ";
	    Text2.append( strdCoef );
	     
	    if( iCnt % 50 == 0 )
	    {
	    	cout << endl << endl ;
	    	cout << "IterNo: " << iCnt << " Teff: " << fTeff << " Logg: " << fLogg << " Xi: " << fXi;
	    	cout << " Gauss: " << fGauss << " Mg_eps: " << fMg_eps << " Wcen: " << strWcen ;
	    	cout << " Coef: "  << dCoef[0] << " Xrv: " << fXrv << endl;
	    }
	    
	    // Time has come to print to the log file.. Yay !
	    logFile<< "Cnt: " << iCnt << " Teff: " << fTeff << " Logg: " << fLogg << " Xi: " << fXi;
	    logFile<< " Gauss: " << fGauss << " Mg_eps: " << fMg_eps << " Wcen: " << strWcen ;
	    logFile<< " Coef: "  << dCoef[0] << " Xrv: " << fXrv << endl;
	    
	    iCnt++;
	    
	    // Moved inside as memory allocated in every iteration
	    delete [] fSx;
        delete [] fSy;

    } // end of iterations over grid points | counter - while loop
    
    delete [] fWavelen; //**can replace dynamic array allocation with static**
    delete [] fData;
    delete [] arrMask;

} // end of loop for no of modes | counter - iCntr
    
// deleting arrays 'fWavelengths' and 'fData' 
delete [] fWavelengths;
delete [] data;

// get time to calculate the program runtime.
timeval time2;
gettimeofday( &time2, NULL );

// calculate time spend
double dExecTime = 0.0;
dExecTime =  ( time2.tv_sec * 1000000 + time2.tv_usec ) - \
             ( time1.tv_sec * 1000000 + time1.tv_usec ) ;

// Status
cout << endl << "Execution time in micro sec:" << dExecTime << endl << endl;

logFile.close();
inputFile.close();

return 0;
}












