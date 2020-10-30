#include "library.h"
#include <iostream>
#include <utility>
#include <cmath>
#include <vector>


void Misc::StartBeforeEnd(Point &start, Point &end) {
    if(start>end){  // just to make sure a is after b :)
        Point temp = start;
        start = end;
        end = temp;
    }
}

Point EquationSolver::BisectionZeroFinder(RealFunction1D &function1D, Point a, Point b, StepSize accuracy) {

    if(function1D(a)*function1D(b)<0){
        Point answer = (a+b)*0.5;
        while(std::abs((function1D(answer)))>accuracy){
            if(function1D(answer)*function1D(a)<0){
                b = answer;
            } else{
                a = answer;
            }
            answer = (a+b)*0.5;
        }
        return answer;
    } else{
        throw ;
    }

}
Point EquationSolver::Newton_Raphson(RealFunction1D &function1D, Point init_p, StepSize accuracy,int MaxTry, double(*diifer)(const RealFunction1D& function, StepSize s, Point p)) {
    RealFunction1D f_prime = [&function1D,&accuracy,&diifer](double x){
        return diifer(function1D,accuracy,x);
    };
    int n = 0;
    while (std::abs(function1D(init_p))>accuracy){
        if(n>=MaxTry){
            break;
        }
        init_p -= function1D(init_p)/f_prime(init_p);
        n++;
    }
    return init_p;

}


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

Point Interpolation::ChebyshevNode(int n, Point a, Point b,int samples) {
    Misc::StartBeforeEnd(a,b);
    return ( 0.5*(a+b) + 0.5*(b-a)*cos(((2*n)-1)*M_PI/(Point)2*samples));
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
RealFunction1D Interpolation::LagrangePolynomialFromFunctionSimpleStep(const RealFunction1D &function1D, int samples,
                                                                       Point start, Point end, Point *x, Point *y) {
    Misc::StartBeforeEnd(start,end);
    double steps = (end - start)/(double)(samples);
    for(int i=0;i<samples;i++){
        x[i] = start+(steps*i);
        y[i] = function1D(x[i]);
    }
    return LagrangePolynomialFromDataSetWithNPoint(samples,x,y);
}
RealFunction1D  Interpolation::LagrangePolynomialFromFunctionChebyshevNodes( RealFunction1D &function1D,
                                                                           int samples, Point start, Point end,
                                                                           Point *x, Point *y) {
    Misc::StartBeforeEnd(start,end);
    for(int i=0;i<samples;i++){
        x[i] = ChebyshevNode(i,start,end,samples);
        y[i] = function1D(x[i]);
    }
    return LagrangePolynomialFromDataSetWithNPoint(samples,x,y);
}
MDPoint Interpolation::LinearRegression2D_parameters(std::vector<MDPoint> arr,Point accuracy,double(*diifer)(const RealFunction1D& function, StepSize s, Point p),int MaxTry) {
    if(arr.at(0).getDim() == 2){ //its 2d reg
        int NumberOfPoints = arr.size();
        Point X =0;
        Point Y =0;
        Point MostHigh = 1 ;
        Point MostLow = 1;
        for(int i=0;i<NumberOfPoints;i++){
            X += arr.at(i).getValue(0); // x values
            Y += arr.at(i).getValue(1); // y values
            if(arr.at(i).getValue(1)/arr.at(i).getValue(0) > MostHigh){
                MostHigh = arr.at(i).getValue(1)/arr.at(i).getValue(0);
            }else if(arr.at(i).getValue(1)/arr.at(i).getValue(0)<MostLow){
                MostLow = arr.at(i).getValue(1)/arr.at(i).getValue(0);
            }
        }

        RealFunction1D finalFunction_a = [&arr,&NumberOfPoints,&X,&Y](Point a){
            Point FirstSigma =0;
            Point SecondSigma = 0;
            Point temp;
            for(int i=0;i<NumberOfPoints;i++){
                temp = (arr.at(i).getValue(1) - (a*arr.at(i).getValue(0) + (1/double(NumberOfPoints))*(Y-a*X)));
                FirstSigma += pow(temp,2);
                SecondSigma+= arr.at(i).getValue(0)*temp;
            }
            return((a/(1+pow(a,2)))*FirstSigma + SecondSigma);
        };
        MDPoint answer = MDPoint(2);
        answer.setValue( EquationSolver::Newton_Raphson(finalFunction_a,MostHigh,accuracy,MaxTry,diifer),0);
        answer.setValue((1/(Point)(NumberOfPoints))*(Y-answer.getValue(0)*X),1);
        return answer;
    }
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
double Integration::SimpleDefinite(const RealFunction1D& function1D, Point LowerBound, Point UpperBound, double stepsize){
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