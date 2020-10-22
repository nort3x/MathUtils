#include <iostream>
#include "library.h"
#include "cmath"
#include "iomanip"

double h(double x){
    return (pow(x,3) - 25);
}
int main() {

    // How to use Differentiation lib
    std::cout<<std::setprecision(16)<<"sin(0.5): "<<sin(0.5)<<
    " \ncos(0.5): "<<cos(0.5)<<
    "\n\nWe Differentiate Sin function:\n3point:"
    <<Differentiation::ThreePointMidPoint1D(sin,0.001,0.5)<<"\n5point:"<<Differentiation::FivePointMidPoint1D(sin,0.001,0.5)<<std::endl;
    std::cout<<Differentiation::NthDiff(Differentiation::FivePointMidPoint1D,sin,0.01,0,2);

    std::cout<<"\n\nWe now differentiate sin function 3 times! \nsin'''(0):";
    std::cout<<Differentiation::NthDiff(Differentiation::Forward1D,sin,0.0001,0,3); //Nth Diff!
    return 0;
}

