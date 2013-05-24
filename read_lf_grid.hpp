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

    container( void ) 
    {
        p = "";
        def = "";
        unit = "";
        i = 0;
        n = 0;
        ion = 0;
        min = 0.0D;
        delta = 0.0D;
    }
};

void Read_LF_Grid( string strGriddef, std::vector<container> &gArray, int iNomc);

#endif