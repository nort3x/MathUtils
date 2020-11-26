#include "../../MathUtils.h"
#include "../../gnuplot.h"


inline double LogisticEquation(double xn, double r);
inline double LogisticIterator(int n,double x0,double r,double (*LogisticEquation)(double xn,double n));
int main(){

    Calculus::SingleVar::Sampler<double,double>::SampledFunction s;;
    for(double r =0 ;r<=4;r+=0.0001){
        for (int i = 100; i < 170; i++) {
            s.y.push_back(LogisticIterator(i,0.01,r,LogisticEquation));
            s.x.push_back(r);
        }
    }
    plt::TwoDim::DataPlot2D(s.x,s.y,"w p pt 7 ps 0.1 ","plot ");
    return 0;
}


double LogisticEquation(double xn, double r){
    return r*xn*(1-xn);
};

inline double LogisticIterator(int n,double x0,double r,double (*LogisticEquation)(double xn,double r)){
    double xn = x0;
    double xn_p1;
    for (int i = 1; i <= n; ++i) {
        xn_p1 = LogisticEquation(xn,r);
        xn = xn_p1;
    }
    return xn_p1;
}