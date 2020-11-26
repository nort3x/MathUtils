#include <iostream>
#include "MathUtils.h"
#include "transform.h"
#include "gnuplot.h"




int main() {
//    using namespace Algebric;
//    using namespace trns;
//    typedef MultiDimPoint<double> Vec3;
//
//    Vec3 p1 = {-2,-1,0};
//    Vec3 p2 = {-1,1,0};
//    Vec3 p3 = {1,1,0};
//    Vec3 p4 = {1,-1,0};
//
//    std::string adress1 = Utils::WriteDataToFile<double>({p1,p2,p3,p4},"/root/Desktop/d1.dat");
//
//    // no transform these!
//
//    RotationTransform<double> r;
//    r.rotate_y(M_PI/4);
//    //r.rotate_x(M_PI/4);
//    //r.rotate_z(M_PI/4);
//    ProjectionTransform<double> p;
//    p.ProjectOnPlane({0,0,1});
//
//    //std::cout<<r.build();
//    Matrix<double> whole_transform = p.build()*r.build();
//
//    Vec3 q1 = whole_transform*p1;
//    Vec3 q2 = whole_transform*p2;
//    Vec3 q3 = whole_transform*p3;
//    Vec3 q4 = whole_transform*p4;
//
//    std::string adress2 = Utils::WriteDataToFile<double>({q1,q2,q3,q4},"/root/Desktop/d2.dat");
//
//    plt::TwoDim::Gnuplot2D<double,double> neg;
//    //std::cout<<whole_transform;
//    neg.addPlotCommand('\''+adress1+'\''+" w p pt 7 ps 3");
//    neg.addPlotCommand('\''+adress2+'\''+" w p pt 7 ps 3");
//
//    neg.plot();
//
//    std::cout<<Angle_deg<double>(q1-q2,q3-q4)<<"\n";
//    std::cout<<Angle_deg<double>(p1-p2,p3-p4)<<"\n";
//    //std::cout<<Angle_deg<double>({1,1},{1,0});
//


    return 0;
}