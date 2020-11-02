#include <iostream>
#include "library.h"  // for Calculations
#include "DataType.h" // for Types
#include "cmath"


int main() {

    /*
     * my library is templated it means you can use (any) data type with defined operations!
     */

    //std::cout<<Differentiation<double,double>::Forward1D([](double x){return x;},0.0001,2);
    std::cout<<Algebric::Matrix<int>(3,3,0);
    std::cout<<Algebric::MultiDimPoint<int>((const int[]){1,2,3},3);
return 0;
}

