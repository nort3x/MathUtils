#include <iostream>
#include "library.h"
#include "cmath"




double h(double x) {
    return (pow(x, 3) - 25);
}

int main() {

    // How to use Differentiation lib
//    std::cout<<std::setprecision(16)<<"sin(0.5): "<<sin(0.5)<<
//    " \ncos(0.5): "<<cos(0.5)<<
//    "\n\nWe Differentiate Sin function:\n3point:"
//    <<Differentiation::ThreePointMidPoint1D(sin,0.001,0.5)<<"\n5point:"<<Differentiation::FivePointMidPoint1D(sin,0.001,0.5)<<std::endl;
//    std::cout<<Differentiation::NthDiff(Differentiation::FivePointMidPoint1D,sin,0.01,0,2);
//
//    std::cout<<"\n\nWe now differentiate sin function 3 times! \nsin'''(0):";
//    std::cout<<Differentiation::NthDiff(Differentiation::Forward1D,sin,0.0001,0,3); //Nth Diff!




    //How to use Interpolation methods


//    double x[] = {-1,0,1};
//    double y[] = {1,0,1};
//    RealFunction1D func = Interpolation::LagrangePolynomialFromDataSetWithNPoint(3,x,y);
//    std::cout<<func(0.5)<<std::endl;


    // sorry for this unpreferable design but u have to allocate Array before using FromFunctionInterpolation
//    Point x[3];
//    Point y[3];
//    auto func1 = Interpolation::LagrangePolynomialFromFunction(sin,3,0.5,+1,x,y);
//    std::cout<<func1(0.1)<<std::endl;


//integration:

//std::cout<<Integration::SimpsonIndefiniteClosed(sin,0,0.01)(M_PI)<<std::endl;

// ArrayList stuff:
//    ArrayList<Point> arrayList;
//    arrayList.add(-1);
//    arrayList.add(1);
//    arrayList.add(2);
//    arrayList.add(3);
//    arrayList.add(4);
//    arrayList.subtract(1,2);
//    arrayList.purge();


    std::vector<MDPoint> listOfPoints;  // a list to hold our Points

//    MDPoint f = MDPoint(2);         // dimension of data (2) means (x,y)    //just for introducing
//    f.setValue(1.5, 0);   // fx = 1.5
//    f.setValue(2, 1);     // fy = 2

    //MDPoint f1 = MDPoint((const Point[]){1,3,4},3);

    listOfPoints.emplace_back((const Point[]){1.5,2},2); // correspond to = { (1.5,2) }
    listOfPoints.push_back(MDPoint(2).setValue(1, 0).setValue(1, 1));   // {(1.5,2) , (1,1)}
    listOfPoints.emplace_back((const Point[]){0.5,0.5},2);   // { (1.5,2), (1,1), (0.5,0.5)} final set

    MDPoint lineParams = Interpolation::LinearRegression2D_parameters(listOfPoints, 0.0000000001,
                                                                      Differentiation::FivePointMidPoint1D,
                                                                      1000); //(a,b) corresponding to y = ax+b
    std::cout << lineParams << std::endl;
    return 0;
}

