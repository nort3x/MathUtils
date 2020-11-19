#include "../../MathUtils.h"
#include "../../gnuplot.h"


int main() {

    using namespace Calculus::SingleVar;

    typedef Algebric::MultiDimPoint<double> Vec3;
    std::function<Vec3(Vec3,Vec3,double)> NewtonLaw = [](const Vec3& velocity,const Vec3& position,double time){
        //define F(v,r,t)/m
        // for demo i will use gravity of a very large body
        double m0 = 1;
        double G_M_of_large_body = 10; // and also keep in minde large body is fixed at origin
        return -((G_M_of_large_body*m0)/(std::pow(position.Norm(),3)))*position;  // F = (GMm/r^3)r_

    };

    double epsilon = 0.5;  // for tweaking veloctiy (so it gonna be ellipse)
    auto l = ODE<Vec3,double>::EulerMethodSecondOrder(NewtonLaw,0,{0,sqrt(10)+epsilon,0},{1,0,0},0.001);

    std::string p = "\'"+Utils::WriteDataToFile<double>({{0,0,0}})+"\'"; // largebody
    plt::ThreeDim::CurvePlot3D<double,double>(l,0,3,100,
"u 2:3:4 w l,"+p+"w p pt 7 ps 5 lc rgb \'gold\' notitle","set xr [-2:2];set yr [-2:2];set zr [-2:2]",0,10,"/root/Desktop/d.dat");
    return 0;
}
