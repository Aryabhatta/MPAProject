#include <iostream>
#include <string>

using namespace std;

#ifndef IDLFUNC_H
#define IDLFUNC_H

// Implementation of IDL total
template<class T>
T IdlTotal( T Arr[], int iSizeArr )
{
    T Sum = 0;
    for( int i = 0; i < iSizeArr; i++ )
    {
        Sum = Sum + Arr[i];
    }
    return Sum;
}

// very important
// implementation of IdlWhere- only one condition
// returns cnt - no of valid indices in IndArr
template<class T>
int IdlWhere( T Arr[], string strOper, T value, int iSizeArr, int * IndArr )
{
	int i = 0, cntr = 0;
	
    if( strOper.empty())
    {
        cout << "ERROR - Invalid comparison operators in IdlWhere.." << endl;
        cout << "Aborting IdlWhere.." << endl;
        return 0;
    }
	// Size of index array equal to size of Arr[] & space allocated in main calling program
	// cleaning the index array
	for ( i = 0; i < iSizeArr; i++ )
	{
		IndArr[i] = -1; // -1 implies undefined (not valid index)
	}
	
    if( strOper == "<")
    {
    	for( i = 0; i < iSizeArr; i++  )
    	{
    		if( Arr[i] < value)	{	IndArr[ cntr ] = i; 	cntr ++;   		}
    	}
    }
    else if( strOper == "<=")
    {
    	for( i = 0; i < iSizeArr; i++  )
    	{
    		if( Arr[i] <= value){	IndArr[ cntr ] = i; 	cntr ++;   		}
    	}
    }
    else if( strOper == ">")
    {
      	for( i = 0; i < iSizeArr; i++  )
       	{
       		if( Arr[i] > value){	IndArr[ cntr ] = i; 	cntr ++;   		}
       	}
    }
    else if( strOper == ">=")
    {
      	for( i = 0; i < iSizeArr; i++  )
       	{
       		if( Arr[i] >= value){	IndArr[ cntr ] = i; 	cntr ++;   		}
       	}
    }
    else if( strOper == "=")
    {
      	for( i = 0; i < iSizeArr; i++  )
       	{
       		if( Arr[i] == value){	IndArr[ cntr ] = i; 	cntr ++;   		}
       	}
    }
    else if( strOper == "!=")
    {
      	for( i = 0; i < iSizeArr; i++  )
       	{
       		if( Arr[i] != value){	IndArr[ cntr ] = i; 	cntr ++;   		}
       	}
    }
    
    return cntr;
};

// Function overloading in C++
// implementation of IDLwhere - two conditions (indexes of elements in range)
// returns cnt - no of valid indices in IndArr
template<class T>
int IdlWhere( T Arr[], string strOper1, T value, string strOper2, T value1, int iSizeArr, int * IndArr)
{
	int i = 0, cntr = 0, fctr  = 0;

    if( strOper1.empty() || strOper2.empty() )
    {
        cout << "ERROR - Invalid comparison operators in IdlWhere.." << endl;
        cout << "Aborting IdlWhere.." << endl;
        return 0;
    }
	
	// Size of index array equal to size of Arr[] & space allocated in main calling program

    // temporary array for holding indices of elements satisafying 1st conditn
    int tempInd[ iSizeArr ];

	// cleaning the index array
	for ( i = 0; i < iSizeArr; i++ )
	{
		IndArr[i] = -1; // -1 implies undefined (not valid index)
        tempInd[i] = -1;
	}
	
    // first condition checking
    if( strOper1 == "<")
    {
    	for( i = 0; i < iSizeArr; i++  )
    	{
    		if( Arr[i] < value)	{	tempInd[ cntr ] = i; 	cntr ++;   		}
    	}
    }
    else if( strOper1 == "<=")
    {
    	for( i = 0; i < iSizeArr; i++  )
    	{
    		if( Arr[i] <= value){	tempInd[ cntr ] = i; 	cntr ++;   		}
    	}
    }
    else if( strOper1 == ">")
    {
      	for( i = 0; i < iSizeArr; i++  )
       	{
       		if( Arr[i] > value){	tempInd[ cntr ] = i; 	cntr ++;   		}
       	}
    }
    else if( strOper1 == ">=")
    {
      	for( i = 0; i < iSizeArr; i++  )
       	{
       		if( Arr[i] >= value){	tempInd[ cntr ] = i; 	cntr ++;   		}
       	}
    }
    else if( strOper1 == "=")
    {
      	for( i = 0; i < iSizeArr; i++  )
       	{
       		if( Arr[i] == value){	tempInd[ cntr ] = i; 	cntr ++;   		}
       	}
    }
    else if( strOper1 == "!=")
    {
      	for( i = 0; i < iSizeArr; i++  )
       	{
       		if( Arr[i] != value){	tempInd[ cntr ] = i; 	cntr ++;   		}
       	}
    }

    // checking for 2nd condition
    if( strOper2 == "<" )
    {
        for( i = 0; i < cntr; i++ )
        {
            if(Arr[tempInd[i]] < value1){ IndArr[fctr] = tempInd[i]; fctr++; }
        }
    }
    else if( strOper2 == "<=")
    {
        for( i = 0; i < cntr; i++ )
        {
            if(Arr[tempInd[i]] <= value1){ IndArr[fctr] = tempInd[i]; fctr++; }
        }
    }
    else if( strOper2 == ">")
    {
        for( i = 0; i < cntr; i++ )
        {
            if(Arr[tempInd[i]] > value1){ IndArr[fctr] = tempInd[i]; fctr++; }
        }
    }
    else if( strOper2 == ">=")
    {
        for( i = 0; i < cntr; i++ )
        {
            if(Arr[tempInd[i]] >= value1){ IndArr[fctr] = tempInd[i]; fctr++; }
        }
    }
    else if( strOper2 == "=")
    {
        for( i = 0; i < cntr; i++ )
        {
            if(Arr[tempInd[i]] == value1){ IndArr[fctr] = tempInd[i]; fctr++; }
        }
    }
    else if( strOper2 == "!=")
    {
        for( i = 0; i < cntr; i++ )
        {
            if(Arr[tempInd[i]] != value1){ IndArr[fctr] = tempInd[i]; fctr++; }
        }
    }

    return fctr;    
};

template < class T >
T IdlMax( T Arr[], int iSizeArr )
{
    T max = -10000;
    for( int i = 0; i < iSizeArr ; i ++ )
    {
        if( Arr[i] > max )
        {
            max = Arr[i];
        }
    }

    return max;
};

template < class T >
T IdlMin( T Arr[], int iSizeArr )
{
    T min = 1000000;
    for( int i = 0; i < iSizeArr ; i ++ )
    {
        if( Arr[i] < min )
        {
            min = Arr[i];
        }
    }

    return min;
};

template < class T >
T IdlMean( T Arr[], int iSizeArr )
{
    T mean = 0;
    for( int i = 0; i < iSizeArr; i++ )
    {
        mean += Arr[i];
    }

    mean = mean / iSizeArr;
    return mean;
};


#endif
