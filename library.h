#ifndef MATHUTILS_LIBRARY_H
#define MATHUTILS_LIBRARY_H
#include <iostream>
#include <functional>

typedef double  StepSize;
typedef StepSize Point;
typedef std::function<double(double)> RealFunction1D;

class Differentiation{
public:
    static double Forward1D(RealFunction1D rf,StepSize s,Point p);
    static double Backward1D(RealFunction1D rf,StepSize s,Point p);
    static double ThreePointMidPoint1D(RealFunction1D rf,StepSize s,Point p);
    static double FivePointMidPoint1D(RealFunction1D rf,StepSize s,Point p);
    static double NthDiff(double (*differentiator)(RealFunction1D ,StepSize s,Point p),RealFunction1D rf,StepSize s,Point p,int n);
}
;


#endif //MATHUTILS_LIBRARY_H
