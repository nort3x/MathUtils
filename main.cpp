#include <iostream>
#include "library.h"  // for Calculations
#include "DataType.h" // for Types
#include "cmath"
#include "gnuplot.h"
#include "Utils.h"

namespace sig = Calculus::SingleVar;
using namespace Algebric;
namespace sca = Calculus::MultiVar;


int main() {
    typedef MultiDimPoint<double> Vec2;
    const int MolsOfGas = 10;
    const int idealGasR = 8;

    // P = P(T,V)
    auto Pressure = [&MolsOfGas, &idealGasR](const Vec2 &vec) -> double {
        // PV = nRT
        // P = nRT/V
        return MolsOfGas * idealGasR * vec.getValue(0) / vec.getValue(1);
    };


    auto samples = sca::Scalar::Sampler<double, double>::FunctionSampler(Pressure, {{100, +101, 100},
                                                                                    {1,   +2,   100}});
    plt::DataPlot<double, double>(samples, "with l", "set pm3d ; splot", "/root/Desktop/d.dat");

    return 0;
}
