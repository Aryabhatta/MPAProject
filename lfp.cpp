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
#include "lfp.hpp"
#include "readGridfile.hpp"

using namespace std;

bool lfp( float * fSx, float *fSy, float fTeff, float fLogg, float fLogz, float fXi, float fEps_dev[2], float fGauss, float fGamma, int iExtrapol, string strGrid, string strRange, int iNomessage,int iCnvl, int iNomc )
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
float fWcen_std[5] = { 6520.0, 5200.0, 4850.0, 4300.0, 8600.0 };
float fWran_std[5] = { 120.0, 100.0, 150.0, 100.0, 200.0 }; 

float fWcen[5] = { 0.0 };   // Initialising all elements to 0.0
float fWran[5] = { 0.0 };   // Initialising all elements to 0.0

int iCntr = 0;  // Loop counter/iterator

for( iCntr = 0; iCntr < 5; iCntr++ )
{
    fWcen[ iCntr ] = fWcen_std[ iCntr ];
    fWran[ iCntr ] = fWran_std[ iCntr ];
}

int iIdx = 0;   

// may replace expression strSran[x] with "HALPHA" // more robust -CHANGE-
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
    cout<<"\t"<<gArray[i].p<<"\t"<<gArray[i].def << "\t" << gArray[i].unit << "\t" << gArray[i].ion <<"\t\t" << gArray[i].min<< "\t"<< gArray[i].delta << "\t"<<gArray[i].n<< endl;
    
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

// specification for grid file
strGridFile = strCnv_log + strWave + ".grid" ;

cout << endl << "Interpolation within grid file:" << strGridFile << endl \
     << endl;

float fWmin = 0.0, fWmax = 0.0 ;
int iNpar;
long int nModel;

fWmin = fWcen[iIdx] - fWran[iIdx];
fWmax = fWcen[iIdx] + fWran[iIdx];

iNpar = gArray.size();
nModel = 1; // dnt knw why it is defined Long in .pro

iIdx = 0;
int iNoElements = 0;
string strPdef[ gArray.size() ];// for copying gArray.def here (except elements)

for( int i = 0; i < gArray.size(); i++ )
{
    if( gArray[i].def == "ABUND" )
    {
        strPdef[i].append( elements[gArray[i].ion] );
    }
    else
    {
        strPdef[i] = gArray[i].def;
    }
}

long int nFac1[ iNpar ];
long int nFac2[ iNpar ];

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
            long int lindgenArr [iSizefEps/2] ; // declaring size of array
            int k = 0;
            for( k =0; k < (iSizefEps/2); k++ )
            {
                lindgenArr[k] = k;              // lindgen property
            }
            string elems[k];
            for( int l=0; l<k; l++ )
            {
                elems[l] = elements[ int ( fEps_dev[ 2 * lindgenArr[l]])];
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
                        iIdx = m;
                    }
                } 
            }

            if( iCntr == 1 )
            {
                dP[i] = fEps_dev[2*iIdx+1];
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
    // ------ where intricacies ------
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
        if( fV[k] < dP[i] )
        {
            iMp[iCntr] = k;
            iCntr ++;
        }
    }

    int iSizeMp = iCntr;
    // ------------ end of where ---------------------------

/*    cout << endl << "dP[i]=" << dP[i] << endl << "fV[k] " << endl;
    for( int k = 0; k < iSizelArray ; k++ )
    {
        cout << fV[k] << " ";
    }

    cout << endl << "mp" << endl;
    for( int k = 0; k < iSizeMp ; k++ )
    {
        cout << iMp[k] << "  ";
    }

    cout << endl << "Count Mp = " << iSizeMp << endl; */

    if( iSizeMp > 0 )
    {    
        int iMpMax =0;
        iMpMax =  IdlMax <int> ( iMp, iSizeMp );
        //  cout << "MpMax = " << iMpMax << endl << endl;
    
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

    /*****************************************************
    // Routine to read the grid file
    *****************************************************/
    //LOGIC HERE MOVED BELOW FOR SIMPLCITY TO DEBUG( JUST 1 CALL)
 
}  

string strGridFileSpecs("/home/shrikant/Desktop/MPA/Files/lf_grid4300.fits");

// Read the grid file with specification in arguments
readGrid( strGridFileSpecs );

return false; // return error 
}
