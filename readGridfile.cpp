#include <iostream>

#include "readGridfile.hpp"

using namespace std;

void readGrid( string GridFileSpecs )
{
	fitsfile *fptr;
	int iStatus = 0;    // Must initialise 'status'
	
	// Opening the FITS file
    fits_open_file( &fptr, GridFileSpecs.data(), READONLY, &iStatus );
    if( iStatus ) {  fits_report_error( stderr, iStatus ); }
    
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
    firstelem = 14145;
    int dntknw=1;
    
    // # Wavelengths to read
    nelements = naxis2;

    // ERROR: says there are 2 col in this table.. need to find out col & rowsize through fits routine
    // request to change the format from col,row to row, col
    
    // routine to read the row corresponding to iFid (values in data)
    //fits_read_img( fptr, TFLOAT, firstelem, nelements, &nullval, fData, &anynull, &iStatus);
    fits_read_col( fptr, TFLOAT, colnum, firstrow, dntknw, nelements, &nullval, fData, &anynull, &iStatus);
   
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
    
    fits_close_file( fptr, &iStatus );
    
    if( iStatus ) /* print any error messages */
	{
	    fits_report_error( stderr, iStatus );
	}
    cout << endl << endl;
}