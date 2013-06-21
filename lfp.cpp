/*************************************************************************
*
* Copyright:		Max Planck Institute for Astrophysics (MPA)
* 
* File:				lfp.cpp
*
* Routine Info:		 					
*       			
* Author:
*
* Modification 
* Log:	    		Added routine main
*		        	Shrikant V	Dec 11, 2012	16:20
*
**************************************************************************/

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include "lfp.hpp"
#include "readGridfile.hpp"

using namespace std;

bool lfp( float * fSx, float * fSy, float fTeff, float fLogg, float fLogz, float fXi, float fEps_dev[2], float fGauss, float fGamma, 
		int iExtrapol, string strGrid, string strRange, int iNomessage,int iCnvl, int iNomc )
{

iNomessage = 0;
iCnvl = 0;

long int lH_cnt = 0;

// printing information when the program executes
cout << endl << "Interpolation within precomputed grid of precomputed ";
cout << "synthetic spectra." << endl;/*
cout << "Please refer to /SIU/code/readme.txt" << endl;
cout << "lfp, w, f, range = , teff = ,logg=, logz=, xi =, eps_dev =, expo=, ";
cout << "gauss=, gamma =, rt =, vsini=, grid =, extrapol=, nomessage=,";
cout << "nomc=" << endl;
cout << "Available ranges: 'HALPHA', 'HBETA', 'GBAND', 'MGB'" << endl;
cout << "(params = teff, logg, logz, xi, eps_dev=[12,x] ) " << endl;
cout << "Available grids : 'lores', 'hires', 'noconv' "<< endl;
cout << "(gauss = 0 - 10, 90 - 150     0   ) " << endl;
cout << "(gamma = 0     , 0  - 500     0   ) " << endl;*/

// elements
string elements[] = { "", "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", \
                      "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", \
                      "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", \
                      "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR"\
                       , "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU",   \
                      "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I",\
                       "XE", "CS", "BA", "LA", "CE", "PR", "ND", "PM",    \
                      "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB"\
                      , "LU", };

// directory paths
string strGrid_dir( "/afs/mpa/data/mbergema/SIU/mod/grid/" );
string strGrid_dir_cnv( "/afs/mpa/data/mbergema/SIU/mod/grid/cnv/" );

// Initialising uninitialiased variables
if( strRange.empty())   {   strRange.append( "HALPHA" );    }
if( strGrid.empty() )   {   strGrid.append(  "hires"  );    }
if( fTeff == 0.0    )   {   fTeff = 5777.0;                 }
if( fLogg == 0.0    )   {   fLogg = 4.44;                   }
if( fLogz == 0.0    )   {   fLogz = 0.0;                    }
if( fXi   == 0.0    )   {   fXi = 1.0;                      }

float fVsini = 0.0;
float fRt    = 0.0;
float fExpo  = 0.0;

// redundant code to make sure the variables are initialised properly
if( fVsini == 0.0    )   {   fVsini = 0.0;                   }
if( fGauss == 0.0    )   {   fGauss = 0.0;                   }
if( fGamma == 0.0    )   {   fGamma = 0.0;                   }
if( fRt    == 0.0    )   {   fRt = 0.0;                      }
if( fExpo  == 0.0    )   {   fExpo = 0.0;                    }

// Initialising the array, if not already initialised
if( fEps_dev[0] == 0.0 ) { fEps_dev[0] = 0.0;                }
if( fEps_dev[1] == 0.0 ) { fEps_dev[1] = 0.0;                }

// below fields are EDITABLE
string strSran[] = { "HALPHA", "MGB", "HBETA", "CAT", "HAL", "HAR"};

// please never change fWcen_std !!!
float fWcen_std[5] = { 6520.0, 5200.0, 4850.0, 4300.0, 8600.0 }; // Central standard wavelenths
float fWran_std[5] = { 120.0, 100.0, 150.0, 100.0, 200.0 };     // range values around central wavelengths

float fWcen[5] = { 0.0 };   // Initialising all elements to 0.0
float fWran[5] = { 0.0 };   // Initialising all elements to 0.0

int iCntr = 0;  // Loop counter/iterator

for( iCntr = 0; iCntr < 5; iCntr++ )
{
    fWcen[ iCntr ] = fWcen_std[ iCntr ];
    fWran[ iCntr ] = fWran_std[ iCntr ];
}

int iIdx = 0;   

// may replace expression strSran[x] with "HALPHA" // more robust -TODO-
if( strRange == strSran[0] )        // H_alpha
{
    iIdx = 0;                       // corresponds to fWcen[0]
}
else if( strRange == strSran[1] )   // Mg_b
{
    iIdx = 1;                       // corresponds to fWcen[1]
}
else if( strRange == strSran[2] )   // H_beta
{
    iIdx = 2;                       // corresponds to fWcen[2]
}
else if( strRange == strSran[3] )   // H_gamma
{
    iIdx = 4;                       // corresponds to fWcen[4];
}
else if( strRange == strSran[4] )   // CaT
{
    iIdx = 0;                       // corresponds to fWcen[0]
                                    // define grid to be read
}
else if( strRange == strSran[5] )   // H_alpha
{
    iIdx = 0;                       // corresponds to fWcen[0]
}
else
{
    iIdx = -1;
}

if( iIdx < 0 )  // error, wavelength range not found
{
    cout << endl << "Error: wavelength range not found! Current ranges:";
    cout << endl << "[ ";
    for( iCntr = 0; iCntr< 6; iCntr++ )
    {
        cout << "\t" << strSran[ iCntr ];  
    }
    cout << " ]" << endl;
    
    return true;
}

string strGriddef( "" );

if( strGrid == "noconv" )
{
    switch(iIdx)
    {
        case 0: strGriddef = strGrid_dir + "lf_grid6520.def";
                break;
        case 1: strGriddef = strGrid_dir + "lf_grid5200.def";
                break;
        case 2: strGriddef = strGrid_dir + "lf_grid4850.def";
                break;
        case 4: strGriddef = strGrid_dir + "lf_grid8600.def";
                break;
        case 5: strGriddef = strGrid_dir + "lf_grid6520.def";
                break;
        default:strGriddef = strGrid_dir + "lf_grid.def";
                break;
    }
}
else
{
    strGriddef = strGrid_dir_cnv + "lf_grid_cnv.def";

    if( iCnvl != 0 )
    {
        switch( iIdx )
        {
            case 1: strGriddef = strGrid_dir + "lf_grid5200.def";
                    break;
            default:strGriddef = strGrid_dir + "lf_grid.def";
                    break;
        }
    }
}

// Vector which will be filled in read_lf_grid routine 
std::vector<container> gArray;

// Call to Read_LF_Grid
Read_LF_Grid( strGriddef, gArray , iNomc);

// Status
cout << "No of elements in vector:" << gArray.size() << endl;
cout << "Snapshot of vector array:" << endl;

// Iterating over all the elements in vector
cout << "gArray.p" << " gArray.def" << " gArray.unit" << " gArray.ion" << " gArray.min" << " gArray.delta" << " gArray.n" << endl;

for( int i = 0; i< gArray.size() ; i++ )
{
    cout<<gArray[i].p<<"\t"<<gArray[i].def << "\t" << gArray[i].unit << "\t" << gArray[i].ion <<"\t\t" << gArray[i].min<< "\t"<< gArray[i].delta << "\t"<<gArray[i].n<< endl;
    
}

/* COMMENTS FROM ORIGINAL IDL FILE lfp.pro
; wcen may be superseded
; ----------------------------------------------
; idx = WHERE(wcen EQ wcen_std(idx))
; idx = idx(0)
;----------------------------------------------
*/ 

strGrid_dir = strGrid_dir_cnv;

string strGridFile( "" );
string strCnv_log( "lf_grid" );

// converting float to string
std::ostringstream ss;
ss << fWcen[iIdx];
string strWave( ss.str());
strWave = strWave.substr( 0, strWave.find_first_of('.'));
strTrim(strWave, 2 );

// specification for grid file, fits for us
strGridFile = strCnv_log + strWave + ".fits" ;

cout << endl << "Interpolation within grid file:" << strGridFile << endl \
     << endl;

float fWmin = 0.0, fWmax = 0.0 ;
int iNpar;
long int nModel;

/*
 * Calculate min wavelength & max wavelength depending on central wavelength
 * and range around it
 */
fWmin = fWcen[iIdx] - fWran[iIdx]; 
fWmax = fWcen[iIdx] + fWran[iIdx];

// No of parameters
iNpar = gArray.size();
nModel = 1; // dnt knw why it is defined Long in .pro

iIdx = 0;
int iNoElements = 0;

// Parameter definition
string strPdef[ gArray.size() ];// for copying gArray.def here (except elements)

for( int i = 0; i < gArray.size(); i++ )
{
    if( gArray[i].def == "ABUND" )
    {
    	// if 'ABUND', get the element no
        strPdef[i].append( elements[gArray[i].ion] );
    }
    else
    {
        strPdef[i] = gArray[i].def;
    }
}

long int nFac1[ iNpar ];
long int nFac2[ iNpar ];

// Initialising to 1
for( int i = 0; i< iNpar; i++ )
{
    nFac1[ i ] = 1;
    nFac2[ i ] = 1;
}

iExtrapol = 0;
double dP[ iNpar ];

// counting no of elements in fEps_dev
int iSizefEps = 2; // Hardcoding here, for use further

for( int i = 0; i< iNpar-1 ; i++)
{
    nModel = nModel * gArray[ i ].n;

    for ( int j = i; j < iNpar-1; j++ )
    {
        nFac1[ i ] = nFac1[ i ] * gArray[ j ].n;
    }
    for( int j = i+1; j < iNpar-1; j++ )
    {
        nFac2[ i ] = nFac2[ i ] * gArray[ j ].n;
    }

    // Comparing strPdef
    if      ( strPdef[i] == "TEFF"  )  { dP[i] = fTeff; } 
    else if ( strPdef[i] == "LOGG"  )  { dP[i] = fLogg; } 
    else if ( strPdef[i] == "LOGZ"  )  { dP[i] = fLogz; } 
    else if ( strPdef[i] == "XI"    )  { dP[i] = fXi; }
    else if ( strPdef[i] == "VSINI" )  { dP[i] = fVsini; } 
    else if ( strPdef[i] == "GAUSS" )  { dP[i] = fGauss; } 
    else if ( strPdef[i] == "EXPO"  )  { dP[i] = fExpo; } 
    else if ( strPdef[i] == "RT"    )  { dP[i] = fRt; } 
    else if ( strPdef[i] == "GAMMA" )  { dP[i] = fGamma; } 
    else
    {
        // RECHECK THIS LOGIC, MAY BE INCONSISTENT
        if( fEps_dev[0] != 0 )
        {
            // lindgen intricacies
            long int lindgenArr [iSizefEps/2] ; // declaring size of array // evaluates to 1, since iSizeofEps=2
            int k = 0;
            for( k =0; k < (iSizefEps/2); k++ )
            {
                lindgenArr[k] = k;              // lindgen property
            }
            string elems[k];
            for( int l=0; l<k; l++ )
            {
                elems[l] = elements[ int ( fEps_dev[ 2 * lindgenArr[l]])]; // accesses elements[0]//'Mg'
            } // logic will never give OutofBound array exception
              // explicit type casting for array subscript

            int iCntr = 0;
            int iIdx = 0; // counter as iIdx[0] in original code lfp.pro   
            for( int m = 0; m < k ; m ++ )
            {
                if( strPdef[i] == elems[m] )
                {   
                    iCntr ++;
                    if( iCntr == 1 )
                    {
                        iIdx = m; // only set one time
                    }
                } 
            }

            if( iCntr == 1 )
            {
                dP[i] = fEps_dev[2*iIdx+1];
                cout << endl << "Calculating dP[i] depeding on fEpd_dev !!!" << endl;
            }
        }
    } // end of else

    // lindgen intricacies
    int iSizelArray = gArray[i].n;
    
    float fV[ iSizelArray ];
    int lindgenArr1[ iSizelArray ];

    for( int k =0; k < iSizelArray ;k ++ )
    {
        lindgenArr1[k] = k;
        
        fV[k] = gArray[i].min * lindgenArr1[k] * gArray[i].delta;
    }
    
    // mp
    int iCnt1 = 0; // no of elements satisfying condition
    for( int k =0; k < iSizelArray; k++ )
    {
        if( fV[k] < dP[i] ) 
        {      
            iCnt1 ++;        
        }
    }
    
    int iMp [ iCnt1 ];
    int iCntr = 0;

    for( int k =0; k< iSizelArray; k++ )
    {
        if( fV[k] < dP[i] ) // Redundant logic , can change TODO
        {
            iMp[iCntr] = k;
            iCntr ++;
        }
    }

    int iSizeMp = iCntr;
    
    if( iSizeMp > 0 )
    {    
        int iMpMax =0;
        iMpMax =  IdlMax <int> ( iMp, iSizeMp );
    
        if( iMpMax >= gArray[i].n -1 )
        {
            gArray[i].i = gArray[i].n - 2;
            cout << "extrapolation (up) of " << gArray[i].def << endl;
            iExtrapol = 1;
        }   
        else
        {
            gArray[i].i = iMpMax;
        }
    }
    else
    {
        gArray[i].i = 0;
        cout << "extrapolation (low) of " << gArray[i].def << endl;
        iExtrapol = 1;  
    } 
}

cout << endl << "No of Models = " << nModel << endl;

/*
 * Working Perfectly till here
 */


 /*****************************************************
 // Routine to read the grid file
 *****************************************************/
// Temporarily pointing to file that we have, change TODO
 string strGridFileSpecs("/home/shrikant/Desktop/MPA/Files/lf_grid4300.fits");

// Read the grid file with specification in arguments
readGrid( strGridFileSpecs , fSx, fSy );
//readGrid( strGridFile);

/*
 * BLACKBOX THAT READS THE FITS FILES AND READS THE 
 * CORRESPONDING WAVELENGTH/FLUX
 */

// The remaining logic of lfp.pro starts from here

long nr_new = 0L;
int iHcnt = 0;

new_start:
iHcnt = iHcnt + 1;

if( iHcnt > 6 )
{
	//FREE LUN
	return true;
}

// Array to hold properties (Teff,logg,logz,Metallicity..) for all models
double xpar[nModel][iNpar]; // can it be directly read from fits file ? TODO

for( int i=0; i< nModel; i++)
{
	for( int j=0; j< iNpar; j++)
	{
		// nFac1 & nFac2 are arrays of size iNpar
		iIdx = (( i % nFac1[j])/nFac2[j]);
		
		// setting the ith row of xpar
		xpar[i][j] = gArray[i].min + (iIdx * gArray[i].delta); 
	}
}

iIdx = 0;
for( int i=0; i< iNpar; i++)
{
	 iIdx += (gArray[i].i * nFac2[i]);
	 cout << "I:" << i << "  gArray[i].i:" << gArray[i].i << " nFac1[i]:"<< nFac1[i] 
	 << " nFac2[i]:"<< nFac2[i] << " iIdx: " << iIdx << endl; 
}

// Copying a row from gArray to p0
double p0[iNpar];
for( int i=0; i<iNpar; i++)
{
	p0[i] = xpar[iIdx][i]; //iIdx th coloumn copied
	cout << "P0[" << i << "]: " << p0[i] << endl;
}

double dDp[iNpar];
for( int i=0; i<iNpar; i++)
{
	// TODO // size of dP is iNpar?, changed to iNpar
	dDp[i] = dP[i] - p0[i];
}

double x[2][iNpar];

for(int i=0; i<iNpar; i++)
{
	// setting first row
	x[1][i] = dDp[i]/gArray[i].delta;
	
	// setting 0th row
	x[0][i] = 1.0 - x[1][i];
}

if( nr_new > 0 )
{
	if( !iNomessage)
	{
		cout << "printing p0" << endl;
		for( int i=0; i<iNpar; i++)
		{
			cout << p0[i] << " "; 
		}
	}
}

int iNpar2 = pow(2,iNpar);

int Lii[iNpar];
for(int i=0; i< iNpar; i++)
{
	Lii[i] = (iNpar-1) - i; // REVERSE(lindgen())
}

for(int i=0; i< iNpar2; i++)
{
	long addx[iNpar];
	
	iIdx = 0;
	for( int j=0; j<iNpar; j++)
	{
		addx[j] = (i % (int)pow(2,(Lii[j]+1)))/ pow(2,Lii[j]);
		iIdx += ((gArray[j].i + addx[j]) * nFac2[j]);
	}
	
	if( iIdx >= nModel )
	{
		cout << endl << "iIdx GE nModel.. Aborting" << endl;
		return true;
	}	
	
	long gi[iNpar];
	string strSrc("");
	
	// replacing exist definition by followin	
	if( iIdx == 0 ) // exist is to check in the byte array exist // TODO change
	{
		if( !iNomessage )
		{
			cout << endl << "Gridpoint missing !" << endl;
			for( int i=0; i < iNpar; i++)
			{
				cout << xpar[iIdx][i] << " "; //print iIdx row
			}				
		}
		
		if( nr_new == 0L )
		{
			for( int j=0; j<iNpar;j++)
			{
				gi[j] = gArray[j].i;
			}
			strSrc = "LOGG"; 
		}
		
next1:
		if( nr_new == 1 )
		{
			strSrc = "TEFF";
		}
		int ipa;
		for(int j=0; j< iNpar; j++ )
		{
			if( gArray[j].def == strSrc )
			{
				ipa = j; //get index
				break;
			}
		}
		
		for( int u = 1; u<=2; u++)
		{
			gArray[ipa].i = gi[ipa] + u; 
			for(int j=0; j<iNpar; j++)
			{
				iIdx += (gArray[ipa].i * nFac2[j]);
			}
			
			if( gi[ipa] < gArray[ipa].n -2  )
			{
				goto next2;
			}
			
			if( iIdx >= nModel )
			{
				cout << endl << "iIdx GE nModel.. Aborting" << endl;
				return true;
			}
			
			if( iIdx != 0)// exist is to check in the byte array exist // TODO change
			{
				if( !iNomessage )
				{
					cout << endl << "extrapolating (low) " ;
					cout << strSrc << endl;
				}
				
				goto nexta;
			}
		}
		
next2:
		
		for( int u=1; u<=2; u++)
		{
			gArray[ipa].i = gi[ipa] - u;
			iIdx = 0;
			for(int j=0; j<iNpar; j++)
			{
				iIdx += (gArray[ipa].i * nFac2[j]);
			}
			if( gi[ipa] < 0 )
			{
				goto next3;
			}
			if( iIdx >= nModel )
			{
				cout << endl << "iIdx GE nModel.. Aborting" << endl;
				return true; // return error true
			}
			
			if( iIdx != 0)// exist is to check in the byte array exist // TODO change
			{
				if( !iNomessage )
				{
					cout << endl << "extrapolating (up) " ;
					cout << strSrc << endl;
				}
				
				goto nexta;
			}
		}

next3:
		if( nr_new == 0 )
		{
			nr_new = nr_new + 1;
			goto next1;
		}
		// FREE_LUN, ut
		// f = fltarr(dim(0))
		if( !iNomessage )
		{
			cout << endl << "No Spectrum !!!" << endl;
			bool bError = true;
			return bError;
		}

nexta:
		nr_new = nr_new + 1;		
		goto new_start;
		
	} // endif
	
	// READ FROM FITS FILE
	//y = flux( descr(iIdx));
	
	double dCoef[iNpar];
	
	for( int m=0; m< iNpar; m++)
	{
		dCoef[m] = 1.0D;
		dCoef[m] = dCoef[m] * x[addx[m]][m];
	}
	
	if( i == 0)
	{
		// f=y*coef
	}
	else
	{
		//f = f+y*coef
	}
} // endfor

// free_lun

if( iCnvl != 0 && (fExpo > 0 || fGauss > 0 || fRt > 0 || fGamma > 0 || fVsini > 0))
{
//	
// Skipping the convolution part.. since its not needed
//
	float res = std::max( fGauss, fExpo);
	res = std::max(res, fRt);
	res = std::max(res, fVsini);
	res = res * 0.1;
}

if( !iNomessage )
{
	cout << endl << strPdef << "=";
	for( int j=0; j<iNpar; j++)
	{
		cout << dP[j] << " ";
	}
	cout << "  ";
	for( int j=0; j<iNpar; j++)
	{
		cout << gArray[j].unit << " ";
	}
	cout << ";";	
}

bool bError = false;
return false; // return error // oder returning bError 
}
















