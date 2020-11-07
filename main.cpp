#include <iostream>
#include "library.h"  // for Calculations
#include "DataType.h" // for Types
#include "cmath"
#include "gnuplot.h"
#include "Utils.h"

int main() {
    using namespace Calculus::SingleVar;
    auto l2 = [](int x)->void{
        Interpolation<double,double>::LagrangePolynomialFromFunctionSimpleStep(sin,x,-4,4)(1);
    };
    auto l3 = Utils::BenchMark<int>::methodDurationFunction(l2);
    plt::TwoDim<double,int>::FunctionPlot("",l3,0,100,100,"with l");



    return 0;
}
