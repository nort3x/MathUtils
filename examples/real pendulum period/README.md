

```cpp
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

```

```
#DATA BEGINS

1	0.01	8.881483103
1	0.02633333333	8.882800656
1	0.04266666667	8.885305223
1	0.059	8.889001327
1	0.07533333333	8.893895666
1	0.09166666667	8.89999715
1	0.108	8.907316952
1	0.1243333333	8.915868571
1	0.1406666667	8.925667907
...
```


![alt text](https://raw.githubusercontent.com/nort3x/MathUtils/master/examples/real%20pendulum%20period/Screenshot_2020-11-15_13-03-50.png)
![alt text](https://raw.githubusercontent.com/nort3x/MathUtils/master/examples/real%20pendulum%20period/Screenshot_2020-11-15_13-15-58.png)


