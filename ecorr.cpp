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
	if( iSzrf == 0) { fRf = (float *)malloc(sizeof(float) * (iObsElem + 1) );}
	if( iSzmask == 0 ) { mask = (BYTE * )malloc(sizeof(BYTE) * (iObsElem + 1));}
	if( dPp == 0 ) { dPp = 5.0d; }
	if( dPn == 0 ) { dPn = 5.0d; }
	// iNomessage set to 0 by default
	
	if( fWr == 0)
	{
		fWr = (float *) malloc( sizeof(float) * 2 );
		fWr[1] = std::max( fThrWave[0], fObsWave[0] );
		fWr[2] = std::min( fThrWave[ iThrElem-1], fObsWave[iObsElem-1]);		
	}
	
	if( fSigma == 0) { fSigma = 0.02; } // 2% tolerance
	
	//
	// Scale consistence  ---------------------------------------------------------------------------
	//
	double * dMp = (double *)malloc( sizeof(double) * (iObsElem+1)); //parameter specific merits
																	  // size of mask from original code
	
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
	
	float fSigmaY = fSigma * 2; // CHANGE:Assuming that sig is not an array & why the hell is it indexed in ecorr.pro?
	BYTE mmm[iCount1];
	
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
	
	float fDln_sig = fSigma / std::max( fObsData[0], fEps); // CHANGE: assuming oyy refers to first element
														 // what the hell, oyy refers to. It is an array and not an
														// element, Oh God !
	// get no of nonzero elements from array mmm
	int iMIndex[iCount1];
	int iMcnt = 0;
	for( i=0; i< iCount1; i++)
	{
		if( mmm[i] > 0)
		{
			iMIndex[iMcnt++] = i;
		}
			
	}
	
	// Copying rf to phi & then operating 
	float fPhi[(iObsElem + 1)];
	for( i=0; i< iObsElem+1; i++)
	{
		fPhi[i] = fRf[i];
	}
	
	int iIndex2[iObsElem+1];
	int iCount2 = IdlWhere( fPhi, "<", 0.1f, (iObsElem+1), iIndex2 );
	
	if( iCount2 > 0 )
	{
		for( i=0; i<iCount2; i++)
		{
			fPhi[iIndex2[i]] = fSmall;
		}
	}
	
	float fPd = iXsigma * fDln_sig;
	float fNd = iXsigma * fDln_sig;
	
	int iPind[iMcnt]; // idxp in orig code
	int iPcnt; // pcnt
	
	for( i=0; i<iMcnt; i++)
	{
		if( fDyy[iMIndex[i]] > fPd )
		{
			iPind[iPcnt++] = i;
		}		
	}
	
	if( iPcnt > 0)
	{
		for( i=0; i< iPcnt; i++)
		{
			iPind[i] = iMIndex[iPind[i]];
		}
	}
	
	int iNind[iCount1]; //idxn
	int iNcnt; // ncnt
	
	for( i=0; i<iCount1; i++)
	{
		if( fDyy[iMIndex[i]] < (-fNd))
		{
			iNind[iNcnt++] = i;
		}
	}
	
	if( iNcnt > 0 )
	{
		for( i=0;i< iNcnt; i++)
		{
			iNind[i] = iMIndex[ iNind[i]]; 
		}
	}
	
	float fSum = 0.0d;
	float fds = 1.0d/ IdlTotal( fPhi, (iObsElem + 1));
	
	if( iPcnt > 0 )
	{
		float fTemp[iPcnt];
		for( i=0; i<iPcnt; i++)
		{
			float temp = (fDyy[iPind[i]] - fPd) / fPd;
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
			float temp = (-fNd - fDyy[iNind[i]]) / fNd;
			temp = std::min( pow( temp, dPn ), 1.00000d ); // dPn should be small
			
			fTemp[i] = fPhi[ iNind[i]] * temp;			
		}
		fSum = fSum + IdlTotal( fTemp, iNcnt );
	}
	
	for( i=0; i < iCount1; i++ )
	{
		if( iMcnt <= lMinPoints )
		{
			cout << endl << "Mask 0";
			dMp[i] = -1.0d;
			if( iNoMessage != 0 )
			{
				cout << endl << "ECORR: number of available points " << iMcnt << \
				" too small. Min = " << lMinPoints << endl;
			}
			break;
		}
		
		dMp[i] = fSum * fds; 
		
		if( iNoMessage != 0 )
		{
			cout << endl << "pcnt = " << iPcnt << ", ncnt = " << iNcnt << ", pd = " << fPd << endl;
		}
	}// end for
	
	if( iNoMessage != 0 )
	{
		cout << endl << "ECORR: " << endl;
		for( i=0; i<(iObsElem+1) ;i++)
		{
			cout << dMp[i] << endl;
		}
	}
	
	delete[] fRf; // should delete ? oder?
	return dMp;
}
















