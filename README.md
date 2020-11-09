 # MathUtils  
  

### Numerical Analysis Methods in nutshell  
  


* fully templated for DataTypes

* functional style (like you can do Nth Derivative with desired Method)

* support plotting by using GNUPlot

* support inline functions

* support file I/Os

* support MultiVar Calculus and LinearAlgebra

* not fully documented ! (yet)  




```cpp
/*
*    example:
*/

#include "library.h"  // for Calculations
#include "DataType.h" // for Types
#include "cmath"
#include "gnuplot.h"


int main() {

    using namespace Calculus::SingleVar;

    auto function = [](double x){
        if((int)x%2 == 0 ){
            return 1;
        } else{
            return -1;
        }
    };

    auto function_diff = Differentiation<double,double>::FivePointMidPoint1D(function,0.001);

    plt::TwoDim<double,double>::FunctionPlot(function,-10,10,1000,"with l");
    plt::TwoDim<double,double>::FunctionPlot(function_diff,-10,10,1000,"with l");

    auto composit_function = [&function](Algebric::MultiDimPoint<double> vec){
        return function(vec[0])*sin(vec[1])*vec[0];
    };
    auto samples = Calculus::MultiVar::Scalar::Sampler<double,double>::FunctionSampler
    (composit_function,{{-5,5,100},{-5,5,100}});
    plt::DataPlot<double,double>(samples,"with l","set hidden; splot ");
    return 0;
}

```
![Alt text](https://github.com/nort3x/MathUtils/blob/master/box_func.png)


![Alt text](https://github.com/nort3x/MathUtils/blob/master/box_func_diff.png)


![Alt text](https://github.com/nort3x/MathUtils/blob/master/comp_func.png)
