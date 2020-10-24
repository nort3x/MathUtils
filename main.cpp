#include <iostream>
#include "library.h"
#include "cmath"
#include "iomanip"

double h(double x){
    return (pow(x,3) - 25);
}
int main() {

    // How to use Differentiation lib
//    std::cout<<std::setprecision(16)<<"sin(0.5): "<<sin(0.5)<<
//    " \ncos(0.5): "<<cos(0.5)<<
//    "\n\nWe Differentiate Sin function:\n3point:"
//    <<Differentiation::ThreePointMidPoint1D(sin,0.001,0.5)<<"\n5point:"<<Differentiation::FivePointMidPoint1D(sin,0.001,0.5)<<std::endl;
//    std::cout<<Differentiation::NthDiff(Differentiation::FivePointMidPoint1D,sin,0.01,0,2);
//
//    std::cout<<"\n\nWe now differentiate sin function 3 times! \nsin'''(0):";
//    std::cout<<Differentiation::NthDiff(Differentiation::Forward1D,sin,0.0001,0,3); //Nth Diff!




    //How to use Interpolation methods


//    double x[] = {-1,0,1};
//    double y[] = {1,0,1};
//    RealFunction1D func = Interpolation::LagrangePolynomialFromDataSetWithNPoint(3,x,y);
//    std::cout<<func(0.5)<<std::endl;


    // sorry for this unpreferable design but u have to allocate Array before using FromFunctionInterpolation
//    Point x[3];
//    Point y[3];
//    auto func1 = Interpolation::LagrangePolynomialFromFunction(sin,3,0.5,+1,x,y);
//    std::cout<<func1(0.1)<<std::endl;


//integration:

std::cout<<Integration::SimpsonIndefiniteClosed(sin,0,1)(M_PI);

    return 0;
}

