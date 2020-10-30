#ifndef MATHUTILS_LIBRARY_H
#define MATHUTILS_LIBRARY_H
#include <iostream>
#include <functional>
#include <vector>

typedef double  StepSize;
typedef StepSize Point;
typedef std::function<double(double)> RealFunction1D;


class MDPoint{
    std::vector<Point> mdp;
public:

    MDPoint(const Point p[],int dim){
        mdp.reserve(dim);
        for(int i =0 ;i<dim;i++){
            setValue(p[i],i);
        }
    }
    MDPoint(int dim){
        mdp.reserve(dim);
    }
    MDPoint(){}

    MDPoint& setValue(Point p,int index){
        mdp.insert(mdp.begin()+index,p);
        return *this;
    }
    Point getValue(int index){
        return mdp.at(index);
    }
    int getDim(){
        return mdp.capacity();
    }

    friend std::ostream& operator<<(std::ostream& o ,MDPoint& mdPoint){
        o<<"(";
        for(int i=0 ; i<mdPoint.getDim();i++){
            if(i!=mdPoint.getDim()-1){
                o<<mdPoint.getValue(i) <<" ,";
            } else{
                o<<mdPoint.getValue(i)<<")";
            }
        }
        return o;
    }
};

class Differentiation{
public:
    static double Forward1D(const RealFunction1D& rf,StepSize s,Point p);
    static double Backward1D(const RealFunction1D& rf,StepSize s,Point p);
    static double ThreePointMidPoint1D(const RealFunction1D& rf,StepSize s,Point p);
    static double FivePointMidPoint1D(const RealFunction1D& rf,StepSize s,Point p);
    static double NthDiff(double (*differentiator)(RealFunction1D ,StepSize s,Point p),const RealFunction1D& rf,StepSize s,Point p,int n);
};

class Interpolation{
private:
    static RealFunction1D LagrangeCoefficient(int n, const Point x[], int K_point);
public:
    static Point ChebyshevNode(int n,Point a ,Point b,int samples);
    static RealFunction1D LagrangePolynomialFromDataSetWithNPoint(int n,const Point x[] ,const Point y[]);
    static RealFunction1D LagrangePolynomialFromFunctionSimpleStep(const RealFunction1D& function1D,int samples,Point start,Point end,Point x[],Point y[]);
    static RealFunction1D LagrangePolynomialFromFunctionChebyshevNodes(RealFunction1D& realFunction1D,int samples,Point start,Point end,Point x[],Point y[]);
    static MDPoint  LinearRegression2D_parameters(std::vector<MDPoint> arr,Point accuracy,double(*diifer)(const RealFunction1D& function, StepSize s, Point p),int MaxTry);
};

class Integration{
public:
    static RealFunction1D SimpleIndefinite(const RealFunction1D& function1D,Point LowerBound,double stepSize);
    static double SimpleDefinite(const RealFunction1D& function1D,Point LowerBound,Point UpperBound,double stepsize);

    static RealFunction1D SimpsonIndefiniteClosed(const RealFunction1D& function1D,Point LowerBound,double stepSize);
    static double SimpsonDefiniteClosed(const RealFunction1D& function1D,Point LowerBound,Point UpperBound,double stepSize);
};

class EquationSolver{
public:
    static Point BisectionZeroFinder(RealFunction1D& function1D, Point a, Point b, StepSize accuracy);
    static Point Newton_Raphson(RealFunction1D &function1D, Point init_p, StepSize accuracy,int maxtry,double(*diifer)(const RealFunction1D& function, StepSize s, Point p));

};

class Misc{
public:
    static void StartBeforeEnd(Point &start,Point &end);
};

#endif //MATHUTILS_LIBRARY_H
