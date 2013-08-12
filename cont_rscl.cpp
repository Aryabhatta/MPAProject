/*************************************************************************
*
* Copyright:		Max Planck Institute for Astrophysics (MPA)
* 
* File:				cont_rscl.cpp
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

#include "cont_rscl.hpp"
#include <iostream>

using namespace std;

void cont_rscl( float *  fObsWave, float * fObsData, int iObsElem, float * fThrWave, float * fThrData, int iThrElem, float * fRy,
		       float * fSig, int iElemSig, float * fRf, int iElemRf, int iMessage, int iStrong, int iConst, 
		       float fSigma, int iUpcut, float * fW_highprio, int iElemhighprio, int iIter, int iLin, int iVdop )
{
	float fWaveMax = std::max( fObsWave[0], fThrWave[0] );
	
	float fOper = std::min(IdlMax( fThrWave, iThrElem ), IdlMax( fObsWave, iObsElem ));
	
	int iIndex1[iObsElem], iIndex2[iThrElem];
	
	// loop iterators
	int i = 0;
	
	int iCount1 = IdlWhere( fObsWave, ">=", fWaveMax, "<=", fOper, iObsElem, iIndex1);
	int iCount2 = IdlWhere( fThrWave, ">=", fWaveMax, "<=", fOper, iThrElem, iIndex2);
	
	// initialising the arrays as per count from above IdlWhere
	float fObsW[iCount1];
	float fObsD[iCount1]; 
	float fThrW[iCount2];
	float fThrD[iCount2];
	
	// Fill values in the arrays
	for( i=0; i < iCount1; i++ )
	{
		fObsW[i] = fObsWave[ iIndex1[i] ];
		fObsD[i] = fObsData[ iIndex1[i] ];
	}
	
	for( i=0; i< iCount2; i++)
	{
		fThrW[i] = fThrWave[ iIndex2[i] ];
		fThrD[i] = fThrData[ iIndex2[i] ];
	}
	
	float fRf1 = NotDefined;
	
	if( fSig != 0 && iElemSig == iCount1 )
	{
		fSigma = fSig[ iCount1 ];
	}
	
	if( fRf != 0 && iElemRf == iCount1 )
	{
		fRf1 = fRf[ iCount1 ];
	}
	
	float fSmall = 0.0;
	
	// NOTE: all keyword set messages at line 24 in cont_rscl.pro are set as default arguments
	
	double fS[ iCount1];
	float fPhi[ iCount1 ];
	
	if( keyword_set( iUpcut ) ) // Upcut set to some value
	{
		int iIndex3[ iCount1 ];
		
		int iCount3 = IdlWhere( fObsD, ">", (float)iUpcut, iCount1, iIndex3 );
		
		for( i = 0; i < iCount1; i++)
		{
			if( fObsD[i] < iUpcut )
			{
				fS[i] =(double) fObsD[i];
			}
			else
			{
				fS[i] =(double) iUpcut;
			}
		}
		
		float fSMax = IdlMax( fS, iCount1 );
		
		for( i=0; i<iCount1;i++)
		{
			fPhi[i] = fS[i]/ fSMax;
		}
		
		if( iCount3 > 0)
		{
			for( i =0; i<iCount3; i++)
			{
				fPhi[iIndex3[i]] = fSmall;
			}
		}
	}
	else
	{
		for( i=0; i<iCount1; i++)
		{
			fS[i] = double( fObsD[i]);
		}
		
		float fSMax = IdlMax( fS, iCount1 );
		for( i=0; i<iCount1;i++)
		{
			fPhi[i] = fS[i]/ fSMax;
		}
	}	
	
	if( fW_highprio != NULLPTR ) // the array fW_highprio has some values
	{	
		int iInd = (iElemhighprio/3) - 1;
		
		for( int ki=0; ki <= iInd; ki++)
		{
			int iIndex4[iCount1];
			int iCount4 = IdlWhere( fObsW, ">=", fW_highprio[3*iInd], "<=", fW_highprio[(3*iInd) + 1], iCount1, iIndex4 );
			
			int iIndex5[iCount1];
			int iCount5 =0;
			
			for( i =0; i< iCount1; i++)
			{
				if( fPhi[ iIndex4[ i ] ] > 0.5)
				{
					iIndex5[iCount5] = iIndex4[i];
					iCount5 ++;
				}
			}
			
			if( iCount5 > 0)
			{
				for( i = 0; i<iCount5; i++)
				{
					fPhi[iIndex5[i]] = fW_highprio[3*iInd + 2]; // will not overwrite memory 
				}
			}
		}
	}
	
	if( fRf != NULLPTR ) // keyword set fRf
	{
		int iIndex6[iCount1];		
		int iCount6 = IdlWhere( fPhi, ">=", 0.1f, iCount1, iIndex6 );
		
		if( iCount6 > 0 )
		{
			for( i=0; i<iCount6; i++)
			{
				if( fRf[ iIndex6[i]]  > fSmall )
				{
					fPhi[ iIndex6[i] ] = fRf[ iIndex6[i] ];
				}
				else
				{
					fPhi[ iIndex6[i] ] = fSmall;
				}                     
			}
		}
	}
	
	int iIT = NotDefined;
	
	// keyword set iIter
	if( keyword_set( iIter )) 	{	iIT = 1;	}
	else						{	iIT = 0;	}
	
	// keyword set iLin
	if( keyword_set( iLin )) 	{	iLin = 1;	}
	else						{	iLin = 0;	}	// Linear LSQ
	
	if( keyword_set( iVdop )) // keyword set Vdop	 	
	{	
		// NOTE: Control never enters here, so skipping the implementation of convol
		// in Modules.cpp, there's just a dummy to avoid compilation errors
		
		// Make copy of fObsW, fObsD, fThrW, fThrD before convolution // for methods ahead 
		float fObsW_old[iCount1];		// put unnecessary things inside the block
		float fObsD_old[iCount1];
		float fThrW_old[iCount2];
		float fThrD_old[iCount2];
		
		for( i=0; i<iCount1; i++)
		{
			fObsW_old[i] = fObsW[i];
			fObsD_old[i] = fObsD[i];
		}
		for( i=0; i<iCount2; i++)
		{
			fThrW_old[i] = fThrW[i];
			fThrD_old[i] = fThrD[i];
		}
		
		float fObsWlog[iCount1], fObsW_oldlog[iCount1]; 
		for( i=0; i< iCount1; i++)
		{
			fObsWlog[i] = log( fObsW[i] );
			fObsW_oldlog[i] = log( fObsW_old[i] );
		}

		// call to convolution
		convol( fObsW, fObsD, fObsW, fObsD, iCount1, iVdop);
		convol( fThrW, fThrD, fThrW, fThrD, iCount2, iVdop);
		
		if( iMessage != 0 )
		{
			cout << endl << "Gaussian !" << endl;
		}
				
		interpol( fObsD, fObsWlog, iCount1, fObsW_oldlog, fObsD, iCount1 ); 
		
		for( i=0; i < iCount1; i++)
		{
			fObsW[i] = fObsW_old[i];
		}
	}
	
	float fThrWlog[iCount2];
	float fT[iCount2];
	
	for( i=0; i < iCount2; i++)
	{
		fThrWlog[i] = log( fThrW[i] );
	}
	
	float fObsWlog[iCount1]; 
	for( i=0; i< iCount1; i++)
	{
		fObsWlog[i] = log( fObsW[i] );		
	}
	
	interpol( fThrD, fThrWlog, iCount2, fObsWlog, fT, iCount2 ); // ERROR: ambiguity here in iCount1 & iCount2 TODO
																 // Currently keeping iCount2 as it matches with original for interpolation
																 // Confirm if the no of points in theoretical & observed spectra 
																 // will always be same or not
	
	int iPlot = 0;
	
	if( iPlot == 1 )
	{
		ofstream logFile;
		string strThr = ReadInput("DIR:LOGDIR") + "thr.log";
		//string strObs("/home/shrikant/Desktop/MPA/Log/obs.log");
		string strObs = ReadInput("DIR:LOGDIR") + "obs.log";
		
		logFile.open( strThr.data(), ios::out );
		for( i=0; i < iCount2; i++)
		{
			logFile << fThrWlog[i] << "\t" << fThrD[i] << endl;
		}
		logFile.close();
		logFile.clear();
		
		logFile.open( strObs.data(), ios::out );
		for( i=0; i < iCount2; i++)
		{
			logFile << fObsWlog[i] << "\t" << fT[i] << endl;
		}
		logFile.close();
		logFile.clear();
	}
	
		
	float fXm = ( IdlMax( fObsWlog, iCount1 ) + IdlMin( fObsWlog, iCount1 ) ) * 0.5f;
	
	for( i=0; i < iCount1; i++)	
	{ 				
		fObsWlog[i] -= fXm;
	}
	
	
	//-----------------------------------------------------------------------------
	//iterate:
	
	// fRy is calculated from polynomial evaluation of array fX(size iElements) 
	// fX is variable array for polynomial evaluation
	// So, allocating space to fRy equal to elements in fX
	// fRy = new float[iElements];

iterate:
	double fS2[ iCount1];
	
	for( i=0;i< iCount1; i++)
	{
		fS2[i] = fS[i] * fPhi[i];
	}
	
	if( iStrong==1 &&  !iIT  )
	{
		if( iMessage != 0 )
		{
			cout << endl << "Strong coupling to highpoints" << endl;
		}
 
		int iIndex7[iCount1];
		int iCount7 = IdlWhere( fPhi, "<=", 0.001f, iCount1, iIndex7 );
 
		if( iCount7 > 0 )
		{
			for( i=0; i < iCount7; i++ )
			{
				fPhi[ iIndex7[i] ] = fSmall;
			}
		}
		
		iCount7=0;
		iCount7 = IdlWhere( fPhi, "<=", 0.8f, iCount1, iIndex7 );
	 
		if( iCount7 > 0 )
		{
			for( i=0;i < iCount7; i++)
			{
				fS2[ iIndex7[i] ] = fS2[ iIndex7[i] ] * pow( fPhi[ iIndex7[i] ], 30 ) ;
			}
		}	 
	}		
	
	if( iLin != 0 )
	{
		float fS3[ iCount1 ];
		float fX2[ iCount1 ];
		int iSizeS2s =(iCount1<iCount2?iCount1:iCount2); // if throreotical = observed, then iCount1=iCount2
		int iSizeS2g =(iCount1>iCount2?iCount1:iCount2); // this logic can change is found that #points in obs spectra
														 // is always equal to #pts in thr spectra
		float fS2T[ iSizeS2g ];
		
		for( i=0; i < iCount1; i++)
		{
			fS3[i] = fS2[i] * fS[i];
			fX2[i] = fObsWlog[i] * fObsWlog[i]; 
		}
		
		// there's discrepancies in the size of fST and fT
		// so, first multiply till both arrays have valid indexes
		for( i=0; i< iSizeS2s; i++)
		{
			fS2T[i] = fS2[i] * fT[i]; 
		}		
		
		// Now Copy remaining elements
		if( iCount1 > iSizeS2s )
		{
			for( i = iSizeS2s; i<iCount1; i++)
			{
				fS2T[i] = fS2[i];
			}
		}
		else
		{
			for( i = iSizeS2s; i<iCount2; i++)
			{
				fS2T[i] = fT[i];
			}
		}
		
		float fA0 = IdlTotal( fS2T, iSizeS2g );
		
		// to compute the sum of s2t(iSizeS2g) * x(iCount1), need to do 
		// hell indexing again & again
		iSizeS2s = (iSizeS2g < iCount1)? iSizeS2g : iCount1;
		iSizeS2g = (iSizeS2g > iCount1)? iSizeS2g : iCount1;
		float temp1[iSizeS2g];
		for( i=0; i< iSizeS2s; i++)
		{
			temp1[i] = fS2T[i] * fObsWlog[i];
		}
		// Copy remaining
		if( iSizeS2s == iCount1 )
		{
			for( i= iSizeS2s; i< iSizeS2g; i++ )
			{
				temp1[i] = fS2T[i];
			}
		}
		else
		{
			for( i= iSizeS2s; i< iSizeS2g; i++ )
			{
				temp1[i] = fObsWlog[i];
			}
		}
		
		float fA1 = IdlTotal( temp1, iSizeS2g);
		
		float fB0 = IdlTotal( fS3, iCount1 );
		
		float temp[iCount1];
		for(i=0; i< iCount1; i++)
		{
			temp[i] = fS3[i] * fObsWlog[i]; 
		}
		
		float fB1 = IdlTotal( temp, iCount1 );
		
		for(i=0; i< iCount1; i++)
		{
			temp[i] = fS3[i] * fX2[i]; 
		}
		float fB2 = IdlTotal( temp, iCount1 );
		
		float fDetA = (fB0*fB2)-(fB1*fB1);
		float fDdetA = 1.0d / fDetA;
		float fA = fDdetA * ((fA0*fB2) - (fA1*fB1));
		float fB = fDdetA * ((fB0*fA1) - (fB1*fA0));
		
		if( iIT == 0 )
		{
			float fX[iObsElem];
			for( i=0; i< iObsElem; i++)
			{
				fX[i] = log(fObsWave[i]) - fXm;
			}
			// input to the poly function is a vector
			// hence the output of poly will be a vector
			poly( fX, iObsElem, fA , fB, fRy );
		}
	} // end of lin prakar
	else
	{
		float fS3[ iCount1 ];
		float fX2[ iCount1 ];
		float fX3[ iCount1 ];
		float fX4[ iCount1 ];
		
		for( i=0; i < iCount1; i++)
		{
			fS3[i] = fS2[i] * fS[i];
			fX2[i] = fObsWlog[i] * fObsWlog[i];
			fX3[i] = fX2[i] * fObsWlog[i];
			fX4[i] = fX3[i] * fObsWlog[i];
		}
		
		int iSizeS2s =(iCount1<iCount2?iCount1:iCount2);
		int iSizeS2g =(iCount1>iCount2?iCount1:iCount2);
		float fS2T[ iSizeS2g ];
		
		// there's discrepancies in the size of fST and fT
		// so, first multiply till both arrays have valid indexes
		for( i=0; i< iSizeS2s; i++)
		{
			fS2T[i] = fS2[i] * fT[i]; 
		}		
		
		// Now Copy remaining elements
		if( iCount1 > iSizeS2s )
		{
			for( i = iSizeS2s; i<iCount1; i++)
			{
				fS2T[i] = fS2[i];
			}
		}
		else
		{
			for( i = iSizeS2s; i<iCount2; i++)
			{
				fS2T[i] = fT[i];
			}
		}
		
		float fA0 = IdlTotal( fS2T, iSizeS2g );

		// to compute the sum of s2t(iSizeS2g) * x(iCount1), need to do 
		// hell indexing again & again
		iSizeS2s = (iSizeS2g < iCount1)? iSizeS2g : iCount1;
		iSizeS2g = (iSizeS2g > iCount1)? iSizeS2g : iCount1;
		float temp1[iSizeS2g];
		for( i=0; i< iSizeS2s; i++)
		{
			temp1[i] = fS2T[i] * fObsWlog[i];
		}
		// Copy remaining
		if( iSizeS2s == iCount1 )
		{
			for( i= iSizeS2s; i< iSizeS2g; i++ )
			{
				temp1[i] = fS2T[i];
			}
		}
		else
		{
			for( i= iSizeS2s; i< iSizeS2g; i++ )
			{
				temp1[i] = fObsWlog[i];
			}
		}
		
		float fA1 = IdlTotal( temp1, iSizeS2g);
		
		// Computation of fA2
		for( i=0; i< iSizeS2s; i++)
		{
			temp1[i] = fS2T[i] * fX2[i];
		}
		
		// Copy remaining
		if( iSizeS2s == iCount1 )
		{
			for( i= iSizeS2s; i< iSizeS2g; i++ )
			{
				temp1[i] = fS2T[i];
			}
		}
		else
		{
			for( i= iSizeS2s; i< iSizeS2g; i++ )
			{
				temp1[i] = fX2[i];
			}
		}
		
		float fA2 = IdlTotal( temp1, iSizeS2g);
		
		// Computation of b0, b1, b2, b3 & b4
		float fB0 = IdlTotal( fS3, iCount1 );
		
		float temp[iCount1];
		for(i=0; i< iCount1; i++)
		{
			temp[i] = fS3[i] * fObsWlog[i]; 
		}
		
		float fB1 = IdlTotal( temp, iCount1 );
		
		for(i=0; i< iCount1; i++)
		{
			temp[i] = fS3[i] * fX2[i]; 
		}
		float fB2 = IdlTotal( temp, iCount1 );

		for(i=0; i< iCount1; i++)
		{
			temp[i] = fS3[i] * fX3[i]; 
		}
		float fB3 = IdlTotal( temp, iCount1 );

		for(i=0; i< iCount1; i++)
		{
			temp[i] = fS3[i] * fX4[i]; 
		}
		float fB4 = IdlTotal( temp, iCount1 );
		
		float fDetA = fB0*((fB2*fB4)-(fB3*fB3)) + fB1*((fB2*fB3)-(fB1*fB4)) + fB2*((fB1*fB3)-(fB2*fB2));
		float fDdetA = 1.0d / fDetA;
		
		float fA = fDdetA * (	fA0*((fB2*fB4)-(fB3*fB3)) +
								fB1*((fA2*fB3)-(fA1*fB4)) +
								fB2*((fA1*fB3)-(fA2*fB2))	);
		float fB = fDdetA * (	fB0*((fA1*fB4)-(fA2*fB3)) +
								fA0*((fB3*fB2)-(fB1*fB4)) +
								fB2*((fA2*fB1)-(fA1*fB2))	);
		float fC = fDdetA * (	fB0*((fA2*fB2)-(fA1*fB3)) +
								fB1*((fA1*fB2)-(fA2*fB1)) +
								fA0*((fB1*fB3)-(fB2*fB2))	);
		if( iIT == 0 )
		{
			float fX[iObsElem];
			for( i=0; i< iObsElem; i++)
			{
				fX[i] = log(fObsWave[i]) - fXm;
			}
			
			poly( fX, iObsElem, fA , fB, fC, fRy );
		}
	}
	
	if( iIT == 1 )
	{
		iIT = 0;
		float temp[iCount1];
		
		for( i=0; i < iCount1 ; i++)
		{
			temp[i] = fRy[i] * (fS[i] + 2.0 * fSigma ); // check fRy 
		}
		
		int iIndex8[ iCount1 ];
		
		// CHANGE 'ERROR': Currently writing fT[0] instead of some particular value
		// as in the idl version 't' is passed to "Where", but the value to be compared with
		// cannot be an array !!! O God  !
		int iCount8 = IdlWhere( temp, "<", fT[0], iCount1, iIndex8 );

		for( i=0; i < iCount1 ; i++)
		{
			temp[i] = fRy[i] * (fS[i] - 2.0 * fSigma ); // CHANGE: check 2. is 2.0 or not 
		}

		int iIndex9[ iCount1 ];
		int iCount9 = IdlWhere( temp, ">", fT[0], iCount1, iIndex9 );
		
		if( iCount8 > 0 )
		{
			if( iMessage != 0 )
			{
				cout << endl << "Ignoring " << iCount8 << "points:" << endl;
			}
			for( i=0; i<iCount8; i++)
			{
				cout << "Points:" << fObsW[iIndex8[i]];
				fS[iIndex8[i]] = fT[iIndex8[i]]/ fRy[i];	
			}
			
			float fPhi[iCount1];
			float fSmax = IdlMax(fS,iCount1);
			
			for( i=0; i< iCount1; i++)
			{
				fPhi[i] = fS[i] / fSmax;
			}
			for(i=0;i<iCount8;i++)
			{
				fPhi[iIndex8[i]] = 0.1; // low priority
			}
			
			if( fW_highprio != 0 )
			{
				int iIndex9[iCount1];
				int iIndex10[iCount1];
				int iIndex11[iCount1];
				
				int iCount9 = 0, iCount10 = 0, iCount11 =0;
				// Custom IDL where for satisfying 3 conditions on different arrays
				for( int ki=0; ki < (iElemhighprio/3)-1; ki++)
				{
					float fValue = fW_highprio[3*ki];
					float fValue1 = fW_highprio[3*ki+1];
					
					iCount9 = IdlWhere(fPhi, ">=", 0.5f, iCount1, iIndex9 );
					iCount10 = IdlWhere(fObsW, ">=", fValue, "<=", fValue1, iCount1, iIndex10);
					
					// we have to take intersection of iIndex9 array und iIndex10 array
					for( int x=0; x < iCount1; x++ )
					{
						float tempa = iIndex9[x];
						//search the index in other array
						for( int y=0; y<iCount1; y++)
						{
							if( iIndex10[y] == tempa )
							{
								iIndex11[iCount11++] = y;
								break;
							}
							else if(iIndex11[y]>tempa)
							{
								break;
							}
						}
					}// @ end of this loop, iIndex11 will have final indices
					
					if( iCount11 > 0)
					{
						for( i=0; i<iCount11; i++)
						{
							fPhi[iIndex11[i]] = fW_highprio[3*ki+2];
						}
					}
				}// endfor ki
			}// endif keywordset (highprio)
			
			if( fRf != 0 )
			{
				int iIndex12[iCount1];
				int iCount12 = IdlWhere( fPhi, ">=", 0.1f, iCount1, iIndex12 );
				
				if( iCount12 > 0 )
				{
					for( i=0; i< iCount12; i++)
					{
						if( fRf[iIndex12[i]] > fSmall )
							fPhi[iIndex12[i]] = fRf[iIndex12[i]];
						else
							fPhi[iIndex12[i]] = fSmall;
					}
				}
			}// endif rf
		}// endif count8 >0
		
		if( iCount9 > 0  && iStrong== 1)			
		{
			if( iMessage != 0 )
			{
				cout << endl << "Bounding " << iCount9 << "points:" << endl;
			}
			// ??????
			for( i=0; i<iCount9; i++)
			{
				cout << "Points:" << fObsW[iIndex9[i]];				
			}
			if(fRf != 0 )
			{
				int iCount13 = 0;
				int iIndex13[iCount9];
				
				// cannot use IdlWhere as indexes for phi are difference htere
				for( i=0; i< iCount9; i++)
				{
					if(fPhi[iIndex9[i]] >= 0.1)
					{
						iIndex13[iCount13++] = i;
					}
				}
				if( iCount13 > 0 )
				{
					for( i=0; i< iCount13; i++)
					{
						float temp = 20.0 * (fRf[iIndex9[iIndex13[i]]]);
						
						if( temp > 0.9 && temp > fSmall ) // COnfirm the logic once by running IDL script
						{
							fPhi[iIndex9[iIndex13[i]]] = temp;						
						}
						else
						{
							fPhi[iIndex9[i]] = 20.0;
						}
					}
				}
			}// endif fRf!=0			
		}// endif (count9 > 0 & iStrong ==1)
		
		if( iCount8  < (iCount1/3)) // '!<' = '>='
		{
			goto iterate; // iteration only if enough points
		}
	}
	
	iPlot = 0;
		
	if( iPlot == 1 )
	{
		ofstream logFile;
		
		//string strRy("/home/shrikant/Desktop/MPA/Log/ry.log");
		string strRy = ReadInput("DIR:LOGDIR") + "ry.log";
				
		logFile.open( strRy.data(), ios::out );
		for( i=0; i < iCount1; i++)
		{
			logFile << fRy[i] << endl;
		}
		logFile.close();
		logFile.clear();
	}
	
/**/
	return; 
}
















