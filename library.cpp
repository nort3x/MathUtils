#include "library.h"
#include <iostream>
#include <utility>



double Differentiation::Forward1D(const RealFunction1D& function, StepSize s, Point p) {
    return (function(p + s) - function(p))/s;
}
double Differentiation::Backward1D(const RealFunction1D& function, StepSize s, Point p) {
    return Differentiation::Forward1D(function,-s,p);
}
double Differentiation::ThreePointMidPoint1D(const RealFunction1D& function, StepSize s, Point p) {
    return (function(p+s)-function(p-s))/(2*s);
}
double Differentiation::FivePointMidPoint1D(const RealFunction1D& function, StepSize s, Point p) {

    return (double(4)*ThreePointMidPoint1D(function,s,p) - ThreePointMidPoint1D(function,2*s,p) )/double(3);
}
double Differentiation::NthDiff(double (*differentiator)(RealFunction1D, StepSize, Point), const RealFunction1D& rf,
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


RealFunction1D Interpolation::LagrangeCoefficient(int n,const  Point *a,int k_point) {

    return [n,a,k_point](double x){
        Point answer = 1;
        Point k = a[k_point];
        Point w = a[1];
        for(int i = 0 ; i<n;i++){
            if(i!=k_point){
                answer=answer*((x - a[i])/(a[k_point]-a[i]));
            }
        }
        return answer;
    };
}
RealFunction1D Interpolation::LagrangePolynomialFromDataSetWithNPoint(int n, const Point *a,const  Point *y) {
    auto func = [n,a,y](double x){
        double answer = 0;
        for(int i = 0 ; i<n;i++){
            answer += (y[i])*(LagrangeCoefficient(n,a,i)(x));
        }
        return answer;
    };
    return func;
}
RealFunction1D Interpolation::LagrangePolynomialFromFunction(const RealFunction1D& function1D, int samples, Point start,
                                                             Point end,Point x[],Point y[]) {
    if(start>end){  // just to make sure a is after b :)
        Point temp = start;
        start = end;
        end = temp;
    }

    double steps = (end - start)/(double)(samples);
    for(int i=0;i<samples;i++){
        x[i] = start+(steps*i);
        y[i] = function1D(x[i]);
    }
    return LagrangePolynomialFromDataSetWithNPoint(samples,x,y);
}

RealFunction1D Integration::SimpleIndefinite(const RealFunction1D& function1D, Point lowerbound, double stepSize) {
    return [function1D,lowerbound,stepSize](double x){
        int chunks = std::abs(x-lowerbound)/stepSize;
        double answer =0 ;
        for(int i =0;i<chunks;i++){
            answer += function1D(lowerbound + (i+0.5)*stepSize)*(stepSize);
        }
        return answer;
    };
}
double Integration::SimpleDefinite(const RealFunction1D& function1D, Point LowerBound, Point UpperBound, double stepsize) {
    return Integration::SimpleIndefinite(function1D,LowerBound,stepsize)(UpperBound);
}

RealFunction1D Integration::SimpsonIndefiniteClosed(const RealFunction1D& function1D, Point LowerBound, double stepSize) {
    return [function1D,LowerBound,stepSize](double x){
        int n = (x-LowerBound)/stepSize;
        double answer = 0;
        for(int i = 0 ;i<n;i++){
            answer += (stepSize/(double)6)*(function1D(LowerBound + i*stepSize)
                    + 4*function1D(LowerBound+(i+0.5)*stepSize)
                    + function1D(LowerBound+(i+1)*stepSize));
        }
        return answer;
    };
}

double Integration::SimpsonDefiniteClosed(const RealFunction1D& function1D, Point LowerBound, Point UpperBound,
                                          double stepSize) {
    return SimpleIndefinite(function1D,LowerBound,stepSize)(UpperBound);
}