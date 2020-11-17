#include "../../MathUtils.h"
#include "../../gnuplot.h"
int main() {

    using namespace Calculus;
    typedef Algebric::MultiDimPoint<double> Vec3;
    intergrator<double,double> i = Calculus::SingleVar::Integration<double,double>::SimpsonDefiniteClosed;

    std::function<double(Vec3)> intergrand = [](Vec3 vec)->double { // vec = {K,C,V}
        double value = ((8)/std::sqrt(2*vec[0])) * (1/std::sqrt(1-std::pow(vec[1]*std::sin(vec[2]),2)));
        return value;
    };
    auto l = Calculus::MultiVar::Scalar::Integration::Partial(intergrand,i,2,0.001,(double)0);
    auto k = Calculus::MultiVar::withConstIndex(l,2,(double)M_PI*0.5);
    plt::ThreeDim::FunctionPlot3D(k,{{1,10,30},{0.01,0.5,30}},"","");
    return 0;
}
