 # MathUtils  
  

### Numerical Analysis Methods in nutshell  
  


* fully templated for DataTypes

* functional style (like you can do Nth Derivative with desired Method)

* support plotting by using GNUPlot

* support inline functions

* support file I/Os

* support MultiVar Calculus and LinearAlgebra

* not fully documented ! (yet)  



#How To Use?

* on Linux:

1 - install gnuplot 

2 - copy headerfiles and cpp files of project to dir of your project

3 Done!





--------
* on Windows:

  1- go to http://www.gnuplot.info/download.html and install gnuplot (i prefer MinGW)
  
  
  2- [add GnuPlot to cmd](https://superuser.com/questions/689333/how-to-add-installed-program-to-command-prompt-in-windows)
  
  
  3- copy headerfiles and cpp files of project to dir of your project
  
  
  4- Done!


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
![alt text](https://github.com/nort3x/MathUtils/blob/master/imgs/box_func.png)


![alt text](https://github.com/nort3x/MathUtils/blob/master/imgs/box_func_diff.png)


![alt text](https://github.com/nort3x/MathUtils/blob/master/imgs/comp_func.png)
