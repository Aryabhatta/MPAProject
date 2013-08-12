/*************************************************************************
*
* Copyright:		Max Planck Institute for Astrophysics (MPA)
* 
* File:				read_lf_grid.hpp
* 
* Routine Info:		Declaration of methods & data structures for
* 					reading from .def file
*
* Author:
*
* Modification 
* Log:	    		
*		        	
*
**************************************************************************/

#include <iostream>
#include <string>
#include <vector>
#include "modules.hpp"

using namespace std;

#ifndef READLFGRID_H
#define READLFGRID_H

class container
{
public:
    string p;
    string def;
    string unit;
    
    long i;
    long n;
    long ion;
    double min;
    double delta;
//    float min;
//    float delta;

    container( void ) 
    {
        p = "";
        def = "";
        unit = "";
        i = 0;
        n = 0;
        ion = 0;
//        min = 0.0f;
//        delta = 0.0f;
        min = 0.0d;
        delta = 0.0d;
    }
};

void Read_LF_Grid( string strGriddef, std::vector<container> &gArray, int iNomc);

#endif