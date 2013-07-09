#include <iostream>
#include <stdlib.h>
#include "ecorr.hpp"

using namespace std;

//------------------------------------------------------------------------------------
// Fit Merit Function
//------------------------------------------------------------------------------------
//int iNomessage = 0, double dPn = 0, double dPp =0, float fWr=0,float fSigma = 0 );
double * ecorr( float * fThrWave, float * fThrData, int iThrElem, float * fObsWave, float * fObsData, int iObsElem,
		        int iXsigma, float * fRf, int iSzrf, BYTE * mask, int iSzmask, int iNoMessage, double dPn, 
		        double dPp, float * fWr,float fSigma)
{
	float fSmall = 0.00f;
	int i = 0; // iterators for 'for' loop
	
	// mask for evaluation of parameter specific merits, mask scales with ox !!!
	
	// KEYWORD SETTING
//	if( iSzrf == 0) { fRf = (float *)malloc(sizeof(float) * (iObsElem + 1) );}
//	if( iSzmask == 0 ) { mask = (BYTE * )malloc(sizeof(BYTE) * (iObsElem + 1));}
	
	// Some defualt settings
	if( dPp == 0 ) { dPp = 5.0d; }
	if( dPn == 0 ) { dPn = 5.0d; }
	// iNomessage set to 0 by default
	
	bool delFwr = false;
	if( fWr == 0)
	{
		delFwr = true;
		fWr = (float *) malloc( sizeof(float) * 2 );
		fWr[0] = std::max( fThrWave[0], fObsWave[0] );
		fWr[1] = std::min( fThrWave[ iThrElem-1], fObsWave[iObsElem-1]); // was going out of bound, why the hell it did not give any error?		
	}
	
	if( fSigma == 0) { fSigma = 0.02; } // 2% tolerance
	
	//
	// Scale consistence  ---------------------------------------------------------------------------
	//
	
	// How come the size of mp defined 1 in original code ?
	/*
	 * Since the following statement evaluates to 1
	 * "N_ELEMENTS(mask(0,*)"
	 * we are declaring variable 'dMp' as double array of size 1 i.e a double
	 *
	double * dMp = (double *)malloc( sizeof(double) * (iObsElem)); // parameter specific merits
	 															  //  size of mask from original code
	 */
	double * dMp = (double *)malloc( sizeof(double) * 1 );
	 															   
	
	int iIndex1[iObsElem];
	int iCount1 = IdlWhere( fObsWave, ">=", fWr[0], "<=", fWr[1], iObsElem, iIndex1 );

	if( iCount1 <= 0 )
	{
		if( iNoMessage != 0 )
		{
			cout << endl << " t and o have no overlapping regions !!" << endl;
		}
		return 0; // aborting from the function
	}
	
	float fObsData1[iCount1];
	float fRf1[iCount1];
	float fTxx[iCount1];	
	
	for( i=0; i< iCount1; i++)
	{
		fObsData1[i] = fObsData[iIndex1[i]];
		fRf1[i] = fRf[iIndex1[i]];
		fTxx[i] = fObsWave[iIndex1[i]];
	}
	
	
	float fTyy[iCount1];
	
	// Interpolation
	interpol( fThrData, fThrWave, iThrElem, fTxx, fTyy, iCount1 );
	
	// Plot before interpol
	int iPlot = 0; // results okay
	if( iPlot == 1)
	{
		ofstream logFile;
		string strBefore("/home/shrikant/Desktop/MPA/Log/interpolb4Ecorr.log");
		
		logFile.open( strBefore.data(), ios::out );
		for( i=0; i < iCount1; i++ )
		{
			logFile << fThrWave[i] << "\t" << fThrData[i] << endl;			
		}
		logFile.close();
		logFile.clear();
		
		string strAfter("/home/shrikant/Desktop/MPA/Log/interpola4Ecorr.log");
		logFile.open( strAfter.data(), ios::out | ios::app );
		for( i=0; i < iCount1; i++ )
		{
			logFile << fTxx[i] << "\t" << fTyy[i] << endl;			
		}
		logFile.close();
		logFile.clear();
	}
	
	
	float fSigmaY = fSigma * 2; // CHANGE:Assuming that sig is not an array & why the hell is it indexed in ecorr.pro?
	BYTE mmm[iCount1];
	
	cout << "Count1: " << iCount1 << endl;

	for( i=0; i<iCount1; i++)
	{
		mmm[i] = mask[iIndex1[i]];
	}
	
	float fEps = 0.001f * IdlTotal(fObsData, iObsElem) / iCount1;
	
	float fDyy[iObsElem];
	
	//following calculation is done assuming that fObsData or fTyy is never 0!
	for( i=0; i< iObsElem; i++)
	{
		if( fObsData[i] > fEps )
		{
			fDyy[i] = log( fObsData[i]/ fTyy[i] );			
		}
		else
		{
			fDyy[i] = log( fEps/fTyy[i]);
		}
	}
	
	long lMinPoints = long( iCount1 * 0.1 ); 	
	//------------------------------------------------------------------------------------------------
	
	float fDln_sig[iObsElem];
	
	for( i=0; i < iObsElem; i++ )
	{
		fDln_sig[i] = fSigma / std::max( fObsData[i], fEps);
	}
	/*
	 * Since , the following statement always evaluates to 1, we will not enclose following
	 * statements in for
	 * "FOR i=0,N_ELEMENTS(mmm(0,*))-1 DO BEGIN" ( line 42, ecorr.pro )
	 */ 
		
	// get no of nonzero elements from array mmm
	int iCountNZ = 0;
	int iIndexNZ[iCount1];
	
	iCountNZ = IdlWhere(mmm, ">", (unsigned char)'0', iCount1, iIndexNZ);
	
	if( iCountNZ <= lMinPoints ) // Is it time to exit ?
	{
		cout << endl << "Mask 0" << endl;
		dMp[0] = -1;
		if( iNoMessage != 0 )
		{
			cout << "ECORR: No of available points (" << iCountNZ << ") too small" << endl;
			cout << "Minimum no of points" << lMinPoints << endl;
		}
	}
	
	// Copying rf to phi & then operating 
	float fPhi[iObsElem];
	for( i=0; i< iObsElem; i++)
	{
		fPhi[i] = fRf[i];
		//cout << fPhi[i] << " ";
	}
	
	int iIndex2[iObsElem];
	int iCount2 = IdlWhere( fPhi, "<", 0.1f, iObsElem, iIndex2 );
	
	if( iCount2 > 0 ) // control never actually  goes here !
	{
		for( i=0; i<iCount2; i++)
		{
			fPhi[iIndex2[i]] = fSmall;		
		}
	}
	
	// fPd & fNd are arrays since fDln_sig is array
	float fPd[iObsElem], fNd[iObsElem];
	for( i=0; i< iObsElem; i++)
	{
		fPd[i] = iXsigma * fDln_sig[i];
		fNd[i] = iXsigma * fDln_sig[i];
	}

	int iPind[iCountNZ]; // idxp in orig code
	int iPcnt=0; // pcnt
	
	for( i=0; i<iCountNZ; i++)
	{
		if( fDyy[iIndexNZ[i]] > fPd[iIndexNZ[i]] )
		{
			iPind[iPcnt++] = i;
		}		
	}
	
	if( iPcnt > 0)
	{
		for( i=0; i< iPcnt; i++)
		{
			iPind[i] = iIndexNZ[iPind[i]];
		}
	}
	
	int iNind[iCountNZ]; //idxn
	int iNcnt=0; // ncnt
	
	for( i=0; i<iCountNZ; i++)
	{
		if( fDyy[iIndexNZ[i]] < (-fNd[ iIndexNZ[i] ]))
		{
			iNind[iNcnt++] = i;
		}
	}

	if( iNcnt > 0 )
	{
		for( i=0;i< iNcnt; i++)
		{
			iNind[i] = iIndexNZ[ iNind[i]];	
		}
	}
	
	float fSum = 0.0d;
	float fds = 1.0d/ IdlTotal( fPhi, iObsElem);
	cout << "fds: " << fds << endl;
	
	if( iPcnt > 0 )
	{
		float fTemp[iPcnt];
		for( i=0; i<iPcnt; i++)
		{
			float temp = ( fDyy[iPind[i]] - fPd[i]) / fPd[i];
			temp = std::min( pow( temp, dPp ), 1.00000d );  // dPp should be high
			
			fTemp[i] = fPhi[ iPind[i]] * temp;			
		}
		fSum = fSum + IdlTotal( fTemp, iPcnt );
	}
	
	if( iNcnt > 0 )
	{
		float fTemp[iNcnt];
		for( i=0; i<iNcnt; i++)
		{
			float temp = (-fNd[i] - fDyy[iNind[i]]) / fNd[i];
			temp = std::min( pow( temp, dPn ), 1.00000d ); // dPn should be small
			
			fTemp[i] = fPhi[ iNind[i]] * temp;			
		}
		fSum = fSum + IdlTotal( fTemp, iNcnt );
	}
	cout << "fSum: " << fSum << endl;

	dMp[0] = fSum * fds; 
		
	if( iNoMessage != 0 )
	{
		cout << endl << "pcnt = " << iPcnt << ", ncnt = " << iNcnt << ", pd[0] = " << fPd[0] << endl;
	}
	
	if( iNoMessage != 0 )
	{
		cout << endl << "ECORR: " << endl;
		for( i=0; i<1 ;i++)
		{
			cout << dMp[i] << endl;
		}
	}
	
	if( delFwr )
	{
		delete [] fWr;
	}

	return dMp;
}
















