#include "library.h"
#include <iostream>



double Differentiation::Forward1D(RealFunction1D function, StepSize s, Point p) {
    return (function(p + s) - function(p))/s;
}
double Differentiation::Backward1D(RealFunction1D function, StepSize s, Point p) {
    return Differentiation::Forward1D(function,-s,p);
}
double Differentiation::ThreePointMidPoint1D(RealFunction1D function, StepSize s, Point p) {
    return (function(p+s)-function(p-s))/(2*s);
}
double Differentiation::FivePointMidPoint1D(RealFunction1D function, StepSize s, Point p) {

    return (double(4)*ThreePointMidPoint1D(function,s,p) - ThreePointMidPoint1D(function,2*s,p) )/double(3);
}
double Differentiation::NthDiff(double (*differentiator)(RealFunction1D, StepSize, Point), RealFunction1D rf,
                                StepSize s, Point p, int n) {

    if(n>0){
        auto func = [differentiator,s,rf](double x){
            return differentiator(rf,s,x);
        };
        return NthDiff(differentiator,func,s,p,n-1);
    } else{
        return rf(p);
    }
}