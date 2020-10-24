#ifndef MATHUTILS_LIBRARY_H
#define MATHUTILS_LIBRARY_H
#include <iostream>
#include <functional>

typedef double  StepSize;
typedef StepSize Point;
typedef std::function<double(double)> RealFunction1D;

class Differentiation{
public:
    static double Forward1D(const RealFunction1D& rf,StepSize s,Point p);
    static double Backward1D(const RealFunction1D& rf,StepSize s,Point p);
    static double ThreePointMidPoint1D(const RealFunction1D& rf,StepSize s,Point p);
    static double FivePointMidPoint1D(const RealFunction1D& rf,StepSize s,Point p);
    static double NthDiff(double (*differentiator)(RealFunction1D ,StepSize s,Point p),const RealFunction1D& rf,StepSize s,Point p,int n);
}
;

class Interpolation{
private:
    static RealFunction1D LagrangeCoefficient(int n, const Point x[], int K_point);
public:
    static RealFunction1D LagrangePolynomialFromDataSetWithNPoint(int n,const Point x[] ,const Point y[]);
    static RealFunction1D LagrangePolynomialFromFunction(const RealFunction1D& function1D,int samples,Point start,Point end,Point x[],Point y[]);
};

class Integration{
public:
    static RealFunction1D SimpleIndefinite(const RealFunction1D& function1D,Point LowerBound,double stepSize);
    static double SimpleDefinite(const RealFunction1D& function1D,Point LowerBound,Point UpperBound,double stepsize);

    static RealFunction1D SimpsonIndefiniteClosed(const RealFunction1D& function1D,Point LowerBound,double stepSize);
    static double SimpsonDefiniteClosed(const RealFunction1D& function1D,Point LowerBound,Point UpperBound,double stepSize);
};


#endif //MATHUTILS_LIBRARY_H
