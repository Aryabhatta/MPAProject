#include <iostream>

#include "readGridfile.hpp"
using namespace std;

#define print 0

bool readGridParams( string strGridFileSpecs, float * pars, int iIdx)
{
	fitsfile * fptr;
	int iStatus = 0;
	
	// Opening the FITS file
    fits_open_file( &fptr, strGridFileSpecs.data(), READONLY, &iStatus );
    if( iStatus ) 
    {  
    	fits_report_error( stderr, iStatus );
    	cout << endl << "Error with reading the file !" << "Aborting.." << endl;
    	return false;
    }
    
    int hdunum =0;
    int hdutype;
    int naxis1,naxis2;
    char * comment = new char [100];
    fits_get_num_hdus(fptr, &hdunum, &iStatus );
       
    // confirm if the no of HDU's are 3, if not, print erro
    if( hdunum != 3 )
    {
    	cout << "No of HDU's in FITS file different than 3.." << endl;
    	cout << "Aborting from reading FITS file" << endl;
    	fits_close_file( fptr, &iStatus );
    	return false; // routine unsuccessful    	
    }
	
   // Move to 2nd HDU - reading flux parameters for specific flux number
   hdunum = 2;    
   fits_movabs_hdu( fptr, hdunum, &hdutype, &iStatus );
   fits_get_hdu_num( fptr, &hdunum );
   
   if( print )
   switch(hdutype)
   {   
	   case IMAGE_HDU: cout << "Image HDU " << endl; break;
	   case ASCII_TBL: cout <<  "Ascii Table" << endl; break;
	   case BINARY_TBL: cout << "Binary Table" << endl; break;   
   }
   
   // Read no of cols in the table = different wavelengths
   fits_read_key( fptr, TINT, "NAXIS1" , &naxis1, comment, &iStatus);
   if( print )
   cout << endl << "Naxis1 = #cols = " << naxis1 << "  Comment = " << comment<< endl;
   
   // Read no of rows in the table = different fluxes
   fits_read_key( fptr, TINT, "NAXIS2" , &naxis2, comment, &iStatus);
   if( print )
   cout << "Naxis2 = #rows = " << naxis2 << "  Comment = " << comment<< endl;
   
   // sanitise input array
   for( int i=0; i< naxis2; i++)
   {
	   pars[i] = 0.00f;
   }

   int colnum = iIdx;
   int firstrow=1;
   // Starting from
   long firstelem;
   int nelements = 1;
   float nullval = 0;
   int anynull;
   
      
   for( int i=0; i<naxis2; i++ )
   {
	   firstelem = iIdx + i * naxis1;
	   fits_read_img( fptr, TFLOAT, firstelem, nelements, &nullval, (pars+i), &anynull, &iStatus);
	   if( iStatus ) 
	   {  
		   fits_report_error( stderr, iStatus );
		   cout << endl << "Error with reading the coloumn !" << "Aborting.." << endl;
		   return false;
	   }
	   if( print )
	   cout << " " << pars[i];
   }
   
   // control here => data read successfully
   
   fits_close_file( fptr, &iStatus );
   return true;
	
}

bool readGridDim( string strGridFileSpecs, int iHdu,  int * naxis1,  int * naxis2 )
{
	fitsfile * fptr;
	int iStatus = 0;
	
	// Opening the FITS file
    fits_open_file( &fptr, strGridFileSpecs.data(), READONLY, &iStatus );
    if( iStatus ) 
    {  
    	fits_report_error( stderr, iStatus );
    	cout << endl << "Error with reading the file !" << "Aborting.." << endl;
    	return false;
    }
    
    int hdunum =0;
    int hdutype;
    int naxis11,naxis22;
    char * comment = new char [100];

    fits_get_num_hdus(fptr, &hdunum, &iStatus );
    
    // confirm if the no of HDU's are 3, if not, print erro
    if( hdunum != 3 )
    {
    	cout << "No of HDU's in FITS file different than 3.." << endl;
    	cout << "Aborting from reading FITS file" << endl;
    	fits_close_file( fptr, &iStatus );
    	return false; // routine unsuccessful    	
    }
	
   // Move to 3nd HDU
   hdunum = iHdu;    
   fits_movabs_hdu( fptr, hdunum, &hdutype, &iStatus );
   fits_get_hdu_num( fptr, &hdunum );
   if( print )
   cout << endl << "Current HDU:" << hdunum << endl;
   
   if( print )
   switch(hdutype)
   {   
	   case IMAGE_HDU: cout << "Image HDU " << endl; break;
	   case ASCII_TBL: cout <<  "Ascii Table" << endl; break;
	   case BINARY_TBL: cout << "Binary Table" << endl; break;   
   }
   
   // Read no of cols in the table = different wavelengths
   fits_read_key( fptr, TINT, "NAXIS1" , &naxis11, comment, &iStatus);
   if( print )
   cout << endl << "Naxis1 = #cols = " << naxis11 << "  Comment = " << comment<< endl;
   
   // Read no of rows in the table = different fluxes
   fits_read_key( fptr, TINT, "NAXIS2" , &naxis22, comment, &iStatus);
   if( print )
   cout << "Naxis2 = #rows = " << naxis22 << "  Comment = " << comment<< endl;
   
   *naxis1 = naxis11;
   *naxis2 = naxis22;

   fits_close_file( fptr, &iStatus );
   
   if( iStatus ) /* print any error messages */
	{
	    fits_report_error( stderr, iStatus );
	}
   
   // if control here, then no error
   return true;   
}

bool readGridFlux( string strGridFileSpecs, float * flux, int iIdx)
{
	fitsfile * fptr;
	int iStatus = 0;
	
	// Opening the FITS file
    fits_open_file( &fptr, strGridFileSpecs.data(), READONLY, &iStatus );
    if( iStatus ) 
    {  
    	fits_report_error( stderr, iStatus );
    	cout << endl << "Error with reading the file !" << "Aborting.." << endl;
    	return false;
    }
    
    int hdunum =0;
    int hdutype;
    int naxis1,naxis2;
    char * comment = new char [100];
    fits_get_num_hdus(fptr, &hdunum, &iStatus );
       
    // confirm if the no of HDU's are 3, if not, print erro
    if( hdunum != 3 )
    {
    	cout << "No of HDU's in FITS file different than 3.." << endl;
    	cout << "Aborting from reading FITS file" << endl;
    	fits_close_file( fptr, &iStatus );
    	return false; // routine unsuccessful    	
    }
	
   // Move to 3nd HDU - reading flux for specific wavelength
   hdunum = 3;    
   fits_movabs_hdu( fptr, hdunum, &hdutype, &iStatus );
   fits_get_hdu_num( fptr, &hdunum );
   
   if( print )
   switch(hdutype)
   {   
	   case IMAGE_HDU: cout << "Image HDU " << endl; break;
	   case ASCII_TBL: cout <<  "Ascii Table" << endl; break;
	   case BINARY_TBL: cout << "Binary Table" << endl; break;   
   }
   
   // Read no of cols in the table = different wavelengths
   fits_read_key( fptr, TINT, "NAXIS1" , &naxis1, comment, &iStatus);
   if( print )
   cout << endl << "Naxis1 = #cols = " << naxis1 << "  Comment = " << comment<< endl;
   
   // Read no of rows in the table = different fluxes
   fits_read_key( fptr, TINT, "NAXIS2" , &naxis2, comment, &iStatus);
   if( print )
   cout << "Naxis2 = #rows = " << naxis2 << "  Comment = " << comment<< endl;
   
   // Sanitise the flux array before reading into it
   for( int i=0; i< naxis1; i++ )
   {
	   flux[i] = 0.00f;
   }
   
   // therefore naxis1 = diff wave.., naxis2 = diff fluxes
      
   // reading flux for diff wavelengths by reading row belonging to iIdx

   // Starting from
   long firstelem = 1 + naxis1 * iIdx ;
   
   // # Wavelengths to read
   int nelements = naxis1;
   float nullval = 0;
   int anynull;

   // routine to read the row corresponding to iFid (values in data)
   fits_read_img( fptr, TFLOAT, firstelem, nelements, &nullval, flux, &anynull, &iStatus);
   if( iStatus ) /* print any error messages */
   {
	   cout << "Error in reading the wavelengths !!! Aborting ..." << endl;
	   fits_report_error( stderr, iStatus );
	   return false;
   }
   
   fits_close_file( fptr, &iStatus );
   
   if( iStatus ) /* print any error messages */
	{
	    fits_report_error( stderr, iStatus );
	}
   
   // if control here, then no error
   return true;
}
bool readGridWave( string strGridFileSpecs, float ** fWave,  int * iWaveCnt)
{
	fitsfile * fptr;
	int iStatus = 0;
	
	// Opening the FITS file
    fits_open_file( &fptr, strGridFileSpecs.data(), READONLY, &iStatus );
    if( iStatus ) 
    {  
    	fits_report_error( stderr, iStatus );
    	cout << endl << "Error with reading the file !" << "Aborting.." << endl;
    	return false;
    }
    
    int hdunum =0;
    int hdutype;
    int naxis1,naxis2;
    char * comment = new char [100];
    fits_get_num_hdus(fptr, &hdunum, &iStatus );
    
    if( print )
    cout << "No of HDU's:" << hdunum << endl;
    
    // confirm if the no of HDU's are 3, if not, print erro
    if( hdunum != 3 )
    {
    	cout << "No of HDU's in FITS file different than 3.." << endl;
    	cout << "Aborting from reading FITS file" << endl;
    	fits_close_file( fptr, &iStatus );
    	return false; // routine unsuccessful    	
    }
    // Confirm whether it is the first HDU
    fits_get_hdu_num( fptr, &hdunum );
    if( print )
    cout << endl << "Current HDU:" << hdunum << endl;
    
    fits_get_hdu_type(fptr, &hdutype, &iStatus);
    
    if( print )
    switch(hdutype)
    {
    case IMAGE_HDU: cout << "Image HDU " << endl; break;
    case ASCII_TBL: cout <<  "Ascii Table" << endl; break;
    case BINARY_TBL: cout << "Binary Table" << endl; break;
    }
    
    // THE FILE IS WAVELENGTH FILE WITH 'Naxis1' ENTRIES
    
    // Read no of rows in the table
    fits_read_key( fptr, TINT, "NAXIS1" , &naxis1, comment, &iStatus);
    if( print )
    cout << endl << "Naxis1 = #rows = " << naxis1 << "  Comment = " << comment<< endl;
    
    // allocating size of original array
    *(iWaveCnt) = naxis1;
    
    // allocating memory in the original array
    *fWave = new float[ naxis1 ];

    // Starting from
    int firstelem = 1;
    
    // # Wavelengths to read
    int nelements = naxis1;
    float nullval = 0;
    int anynull;

    // routine to read the row corresponding to iFid (values in data)
    fits_read_img( fptr, TFLOAT, firstelem, nelements, &nullval, *fWave, &anynull, &iStatus);
    if( iStatus ) /* print any error messages */
	{
    	cout << "Error in reading the wavelengths !!! Aborting ..." << endl;
	    fits_report_error( stderr, iStatus );
	    return false;
	}
    
    fits_close_file( fptr, &iStatus );
    
    if( iStatus ) /* print any error messages */
	{
	    fits_report_error( stderr, iStatus );
	}
    
    // if control here, then no error
    return true;	
}

bool GridDescrExists( string strGridFileSpecs, int iIdx)
{
	fitsfile * fptr;
	int iStatus = 0;
	
	// Opening the FITS file
    fits_open_file( &fptr, strGridFileSpecs.data(), READONLY, &iStatus );
    if( iStatus ) 
    {  
    	fits_report_error( stderr, iStatus );
    	cout << endl << "Error with reading the file !" << "Aborting.." << endl;
    	return false;
    }
    
    int hdunum =0;
    int hdutype;
    int naxis1,naxis2;
    char * comment = new char [100];
    fits_get_num_hdus(fptr, &hdunum, &iStatus );
    
    if( print )
    cout << "No of HDU's:" << hdunum << endl;
    
    // confirm if the no of HDU's are 3, if not, print erro
    if( hdunum != 3 )
    {
    	cout << "No of HDU's in FITS file different than 3.." << endl;
    	cout << "Aborting from reading FITS file" << endl;
    	fits_close_file( fptr, &iStatus );
    	return false; // routine unsuccessful    	
    }
    
    // Move to 2nd HDU
   hdunum = 2;    
   fits_movabs_hdu( fptr, hdunum, &hdutype, &iStatus );
   fits_get_hdu_num( fptr, &hdunum );
   
   if( print )
   cout << endl << "Current HDU:" << hdunum << endl;
   
   if( print )
   switch(hdutype)
   {   
	   case IMAGE_HDU: cout << "Image HDU " << endl; break;
	   case ASCII_TBL: cout <<  "Ascii Table" << endl; break;
	   case BINARY_TBL: cout << "Binary Table" << endl; break;   
   }
   
   // Read no of cols in the table
   fits_read_key( fptr, TINT, "NAXIS1" , &naxis1, comment, &iStatus);
   
   if( print )
   cout << endl << "Naxis1 = #cols = " << naxis1 << "  Comment = " << comment<< endl;
   
   // Read no of rows in the table
   fits_read_key( fptr, TINT, "NAXIS2" , &naxis2, comment, &iStatus);
   
   if( print )
   cout << "Naxis2 = #rows = " << naxis2 << "  Comment = " << comment<< endl;
   
   // Reading data for iIdx
   if( print )
   cout << "Reading data belonging to coloumn: " << iIdx << endl;
   
   float fData[ naxis2 ]; // reading all rows
   int colnum = iIdx;
   int firstrow=1;
   // Starting from
   long firstelem;
   int nelements = 1;
   float nullval = 0;
   int anynull;
   
      
   for( int i=0; i<naxis2; i++ )
   {
	   firstelem = iIdx + i * naxis1;
	   fits_read_img( fptr, TFLOAT, firstelem, nelements, &nullval, (fData+i), &anynull, &iStatus);
	   if( iStatus ) 
	   {  
		   fits_report_error( stderr, iStatus );
		   cout << endl << "Error with reading the coloumn !" << "Aborting.." << endl;
		   return false;
	   }
	   if( print )
	   cout << " " << fData[i];
   }
   
   // control here => data read successfully
   
   fits_close_file( fptr, &iStatus );
   return true;
}

// Read the grid file & store the wavelength & flux in the arrays
void readGrid( string GridFileSpecs, float * fWavelen, float * fFlux )
{
	fitsfile *fptr;
	int iStatus = 0;    // Must initialise 'status'
	
	// Opening the FITS file
    fits_open_file( &fptr, GridFileSpecs.data(), READONLY, &iStatus );
    if( iStatus ) 
    {  
    	fits_report_error( stderr, iStatus );
    	cout << endl << "Error with reading the file !" << "Aborting.." << endl;
    	return;
    }
    
    int hdunum =0;
    fits_get_num_hdus(fptr, &hdunum, &iStatus );
    cout << "No of HDU's:" << hdunum << endl;
    
    int i=0;
    int hdutype;
    int naxis1,naxis2;
    char * comment = new char [100];
    // Reading a row of image/table/data  
    long firstelem = 0, nelements = 0;
    float nullval = 0;
    int anynull;
    
    // Confirm whether it is the first HDU
    fits_get_hdu_num( fptr, &hdunum );
    cout << endl << "Current HDU:" << hdunum << endl;
    
    fits_get_hdu_type(fptr, &hdutype, &iStatus);
    switch(hdutype)
    {
    case IMAGE_HDU: cout << "Image HDU " << endl; break;
    case ASCII_TBL: cout <<  "Ascii Table" << endl; break;
    case BINARY_TBL: cout << "Binary Table" << endl; break;
    }
    
    // THE FILE IS WAVELENGTH FILE WITH 102 ENTRIES
    
    // Read no of rows in the table
    fits_read_key( fptr, TINT, "NAXIS1" , &naxis1, comment, &iStatus);
    cout << endl << "Naxis1 = #rows = " << naxis1 << "  Comment = " << comment<< endl;
    
    int iNoWave = 10;
    
    float * fWave = new float[ iNoWave ];

    // Starting from
    firstelem = 1;
    
    // # Wavelengths to read
    nelements = iNoWave;

    // routine to read the row corresponding to iFid (values in data)
    fits_read_img( fptr, TFLOAT, firstelem, nelements, &nullval, fWave, &anynull, &iStatus);
    
    //print out the wavelengths read
    for( i=0; i<iNoWave; i++)
    {
    	cout << fWave[i] << "  ";
    }
    
    delete[] fWave;
    
    /*************************************************************
     * For HDU's 2nd & 3rd, current logic is to read them row wise
     * As reading by row is good for memory access with no jumps
     * But the files currently are in wrong order, rows now should be 
     * made coloumns & col as rows, waiting for that...
     * Leaving the function, keeping this logic intact
     **************************************************************/
   
        
    // Move to 2nd HDU
    hdunum = 2;    
    fits_movabs_hdu( fptr, hdunum, &hdutype, &iStatus );
    fits_get_hdu_num( fptr, &hdunum );
    cout << endl << "Current HDU:" << hdunum << endl;
    switch(hdutype)
    {
    case IMAGE_HDU: cout << "Image HDU " << endl; break;
    case ASCII_TBL: cout <<  "Ascii Table" << endl; break;
    case BINARY_TBL: cout << "Binary Table" << endl; break;
    }
    
    // THIS HDU IS STORED AS 14144 COL * 5 ROWS, REVERSE WAY WOULD HAVE BEEN GOOD
    
    // Read no of cols in the table
    fits_read_key( fptr, TINT, "NAXIS1" , &naxis1, comment, &iStatus);
    cout << endl << "Naxis1 = #cols = " << naxis1 << "  Comment = " << comment<< endl;
    
    // Read no of rows in the table
    fits_read_key( fptr, TINT, "NAXIS2" , &naxis2, comment, &iStatus);
    cout << "Naxis2 = #rows = " << naxis2 << "  Comment = " << comment<< endl;
    
    // reading a col of data corresponding to certain flux
    
    float * fData = new float[ naxis2 ];

    int colnum = 3, firstrow=1;
    // Starting from
    firstelem = 1;
    int dntknw=1;
    
    // # Wavelengths to read
    nelements = naxis2;

    // ERROR: says there are 2 col in this table.. need to find out col & rowsize through fits routine
    // request to change the format from col,row to row, col
    
    // routine to read the row corresponding to iFid (values in data)
    fits_read_img( fptr, TFLOAT, firstelem, nelements, &nullval, fData, &anynull, &iStatus);
    //fits_read_col( fptr, TFLOAT, colnum, firstrow, dntknw, nelements, &nullval, fData, &anynull, &iStatus);
   
    //print out the wavelengths read
    for( i=0; i<naxis2; i++)
    {
    	cout << fData[i] << "  ";
    }
    
    delete[] fData;


    // Move to 3rd HDU
    hdunum = 3;    
    fits_movabs_hdu( fptr, hdunum, &hdutype, &iStatus );
    fits_get_hdu_num( fptr, &hdunum );
    cout << endl << "Current HDU:" << hdunum << endl;
    switch(hdutype)
    {
    case IMAGE_HDU: cout << "Image HDU " << endl; break;
    case ASCII_TBL: cout <<  "Ascii Table" << endl; break;
    case BINARY_TBL: cout << "Binary Table" << endl; break;
    }

    // Read no of cols in the table
    fits_read_key( fptr, TINT, "NAXIS1" , &naxis1, comment, &iStatus);
    cout << endl << "Naxis1 = #cols = " << naxis1 << "  Comment = " << comment<< endl;
    
    // Read no of rows in the table
    fits_read_key( fptr, TINT, "NAXIS2" , &naxis2, comment, &iStatus);
    cout << "Naxis2 = #rows = " << naxis2 << "  Comment = " << comment<< endl;
    
    // Read elements
    float * fData1 = new float[ 10 ];
    
    // Starting from wavelength
    firstelem = 1;
    
    // # fluxes to read
    nelements = 10;
    
    // routine to read the row corresponding to iFid (values in data)
    fits_read_img( fptr, TFLOAT, firstelem, nelements, &nullval, fData1, &anynull, &iStatus);
   
    //print out the wavelengths read
    for( i=0; i<10; i++)
    {
    	cout << fData1[i] << "  ";
    }
    
    delete[] fData1;
    
    fits_close_file( fptr, &iStatus );
    
    if( iStatus ) /* print any error messages */
	{
	    fits_report_error( stderr, iStatus );
	}
    cout << endl << endl;
}