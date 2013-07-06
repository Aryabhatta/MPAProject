/*************************************************************************
*
* Copyright:		Max Planck Institute for Astrophysics (MPA)
* 
* File:				read_lf_grid.cpp
*
* Routine Info:		reads the file *.def and tabulates it in a vector 
* 					array to be used afterwards
* 					       			
* Author:
*
* Modification 
* Log:	    		Added routine main
*		        	Shrikant V	Dec 11, 2012	16:20
*
**************************************************************************/

#include <iostream>
#include <stdlib.h>
#include <fstream>

#include "read_lf_grid.hpp"

using namespace std;


void Read_LF_Grid( string strGridDef, std::vector<container> &gArray, int iNomc)
{

// elements
string elements[] = { "", "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", \
                      "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", \
                      "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", \
                      "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR"\
                       , "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU",   \
                      "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I",\
                       "XE", "CS", "BA", "LA", "CE", "PR", "ND", "PM",    \
                      "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB"\
                      , "LU", }; // 72 ELEMENTS, 


// Manipulating strGridDef file for testing TODO-CHANGE-
//strGridDef = "/home/shrikant/Desktop/MPA/Files/lf_grid4850.def";
//strGridDef = "/home/shrikant/Desktop/MPA/Files/lf_grid6520.def";

// Status
cout << "Reading from grid definition file " << strGridDef << endl;

ifstream inputFile;     // for reading grid file

inputFile.open( strGridDef.data(),  ios::in );

if( ! inputFile.is_open() )
{
    cout << "Grid definition file cannot be opened !!! ERROR !!!" << endl;
    return;
}

string strRec;
string strDef_par( "" );  // the strings in the array can be appended to make  
                          // one large string separated by delimiters
string strPar(""); // same reason as strDef_par

int iMode = 0;

int iCntr = 0; // counter for strDef_par
int jCntr = 0; // counter for strPar

while( ! inputFile.eof() )//EOF
{
    // reading record from grid file
    getline( inputFile, strRec );
    
    // to get away with carriage return
    if( strRec.length()==1 && strRec.substr(0,1)=="\r")
    	continue;

    // tokenising the line read to find where it fits !
    if( !strRec.empty( ) && strRec.substr( 0, 1 ) != "#")
    {
        if( strRec.substr( 0, 3 ) == ">>>" )
        {
            strRec = strRec.substr(4, strRec.length() );

            // trimming on both sides
            strTrim( strRec, 2);

            if( strRec == "parameter definition" )
            {
                iMode = 1;
            }   
            else if( strRec == "parameter description" )
            {
                iMode = 2;
            }
            else if( strRec == "pathnames" )
            {
                iMode = 3;
            }
            else if( strRec == "spectral ranges" )
            {
                iMode = 4;
            }
            else
            {
                cout << endl << "Unknown Meta-Description: ";
                cout << strRec << endl ;
                iMode = 0;
            }
        }
        else
        {
            switch( iMode )
            {
                case 1: // "parameter definition"
                        // append to strDef_par                        
                        strDef_par.append( strRec );
                        strDef_par.append( "," ); // append delimiter
                        iCntr++;
                        break;

                case 2: // "parameter description"
                        // add to str_par array                    
                        strPar.append( strRec );
                        strPar.append( "," );   // delimiter
                        jCntr++;
                        break;
                case 3:// EXECUTE  (how to execute from C++ )
                        break;
                case 4:// EXECUTE
                        break;
                default: break;
            }
        }
    }

}// WHILE EOF , RUN THIS LOOP 

inputFile.close();  

// make the par array consisting of various fields
// below is the format in which the variables are stored in gArray
// as a class
// { p:'', def:'', unit:'', i:0L, min:0.D0, delta:0.D0, n:0L, ion:0L}

string strToken( "" );
string strTemp(  "" );

char cDelim = ',';
int iIndex = 0;

for( int i = 0; i < iCntr; i++ )
{
    container gObject;

    // get the token
    strToken = strNexttoken( strDef_par, cDelim );

    if( ! strToken.empty() )
    {
        strTemp = strNexttoken( strToken, '=' );

        gObject.p = strTemp;

        iIndex = 0;
        iIndex = strToken.find_first_of( '=' );

        if( iIndex == string::npos )
        {
            strTrim( strToken, 2 );
            strUpper( strToken );

            gObject.def = strToken;
        }
        else
        {
            strTemp = strToken.substr(0, iIndex);
            strToken.erase( 0, iIndex+1 );
        
            strTrim( strTemp, 2); // trimming on both ends
            strUpper( strTemp );    // change to Upper case

            gObject.def = strTemp;

            strTrim( strToken, 2 );
            strUpper( strToken );
            gObject.unit = strToken;
        }
    }
        
    if( gObject.def.substr( 0, 5 ) == "ABUND")
    {
        iIndex = 0;
        strTemp =  gObject.def;
        iIndex = strTemp.find_first_of( ':' );

        if( iIndex != 0)
        {
            gObject.def = "ABUND";
            strTemp.erase( 0, iIndex+1 );
            
            strTrim( strTemp, 2 ); // trimming on both ends
            strUpper( strTemp );    // change to Upper case

            int iElementIndex = 0;
            // 72 hard coded as # elements in elements array
            for( int iElemIndex = 0; iElemIndex < 72 ; iElemIndex ++ )
            {
                if( elements[iElemIndex] == strTemp )
                {
                    gObject.ion = iElemIndex;
                }
            }
        }
    }

    // Push the object thus created on the vector
    gArray.push_back( gObject );
}

cDelim = ',';

// Analysing strPar array
for( int i = 0; i< gArray.size(); i++ )
{
    strTemp = strNexttoken( strPar, cDelim );

    strToken = strNexttoken( strTemp, '=' );
    
    int iNoOccur = 0;
    iIndex = 0;

    for(int j = 0; j< gArray.size() ; j++ )
    {
        if( gArray[j].p == strToken )
        {
            if( iNoOccur == 0 )
            {
                iIndex = j;
            }
            iNoOccur ++;
        }
    }

    if( iNoOccur == 1 )
    {
        strToken.erase(0);

        strToken = strNexttoken( strTemp, ' ' );            
        gArray[iIndex].min = atof( strToken.data() );

        strToken.erase(0);

        strToken = strNexttoken( strTemp, ' ' );            
        gArray[iIndex].delta = atof( strToken.data() );

        gArray[iIndex].n = atof( strTemp.data() );      
    }
}
    // NOTE
    // from code it seems that there can be at the most one match
    // confirm this and to optimise, we can break out of the loop
    // once the match is found
return;
}

