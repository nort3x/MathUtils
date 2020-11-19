#include "../../MathUtils.h"
#include "../../gnuplot.h"


int main() {

    using namespace Calculus::SingleVar;

    std::function<double(double,double,double)> Field = [](const double yp,const double y,double t){  //define field
        return -y;
    };
    auto l = ODE<double,double>::EulerMethodSecondOrder(Field,0,0,1,0.001);

    plt::TwoDim::Gnuplot2D<double,double> g;
    g.addPlotCommand("cos(x)");
    g.addFunctionPlot(l,0,10*M_PI,100,"w point pointsize 1 pt 7 ");
    g.plot();

    return 0;
}