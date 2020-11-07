#ifndef MATHUTILS_LIBRARY_CPP
#define MATHUTILS_LIBRARY_CPP

#include "DataType.h"
#include "cmath"
#include "library.h"
#include "limits.h"

tempT void Calculus::Misc<T>::StartBeforeEnd(T &start, T &end) {
    if (start > end) {  // just to make sure a is after b :)
        T temp = start;
        start = end;
        end = temp;
    }
}


tempTK K
Calculus::SingleVar::EquationSolver<T, K>::BisectionZeroFinder(RealFunction1D &function1D, K a, K b, T accuracy) {

    if (function1D(a) * function1D(b) < 0) {
        K answer = (a + b) * 0.5;
        while (std::abs((function1D(answer))) > accuracy) {
            if (function1D(answer) * function1D(a) < 0) {
                b = answer;
            } else {
                a = answer;
            }
            answer = (a + b) * 0.5;
        }
        return answer;
    } else {
        throw;
    }

}

tempTK K
Calculus::SingleVar::EquationSolver<T, K>::Newton_Raphson(RealFunction1D &function1D, K init_p, T accuracy, int MaxTry,
                                                          T(*diifer)(const RealFunction1D &function, K s, K p)) {
    RealFunction1D f_prime = [&function1D, &accuracy, &diifer](K x) {
        return diifer(function1D, accuracy, x);
    };
    int n = 0;
    while (std::abs(function1D(init_p)) > accuracy) {
        if (n >= MaxTry) {
            break;
        }
        init_p -= function1D(init_p) / f_prime(init_p);
        n++;
    }
    return init_p;

}


tempTK T Calculus::SingleVar::Differentiation<T, K>::Forward1D(const RealFunction1D &rf, K StepSize, K atPoint) {
    return (rf(atPoint + StepSize) - rf(atPoint)) / StepSize;
}

tempTK std::function<T(K)> Calculus::SingleVar::Differentiation<T, K>::Forward1D(const RealFunction1D &rf, K StepSize) {
    return [&rf, StepSize](K k) {
        return Forward1D(rf, StepSize, k);
    };
}

tempTK T Calculus::SingleVar::Differentiation<T, K>::Backward1D(const RealFunction1D &rf, K StepSize, K atPoint) {
    return Differentiation::Forward1D(rf, -StepSize, atPoint);
}

tempTK std::function<T(K)> Calculus::SingleVar::Differentiation<T, K>::Backward1D(const RealFunction1D &rf, K StepSize) {
    return [&rf, StepSize](K x) {
        return Backward1D(rf, StepSize, x);
    };
}

tempTK T Calculus::SingleVar::Differentiation<T, K>::ThreePointMidPoint1D(const RealFunction1D &rf, K StepSize, K atPoint) {
    return (rf(atPoint + StepSize) - rf(atPoint - StepSize)) / (2 * StepSize);
}

tempTK std::function<T(K)> Calculus::SingleVar::Differentiation<T, K>::ThreePointMidPoint1D(const RealFunction1D &rf, K StepSize) {
    return [&rf, StepSize](K x) {
        return ThreePointMidPoint1D(rf, StepSize, x);
    };
}

tempTK T Calculus::SingleVar::Differentiation<T, K>::FivePointMidPoint1D(const RealFunction1D &rf, K StepSize, K atPoint) {

    return (double(4) * ThreePointMidPoint1D(rf, StepSize, atPoint) - ThreePointMidPoint1D(rf, 2 * StepSize, atPoint)) /
           double(3);
}

tempTK std::function<T(K)> Calculus::SingleVar::Differentiation<T, K>::FivePointMidPoint1D(const RealFunction1D &rf, K StepSize) {
    return [&rf, StepSize](K x) {
        return FivePointMidPoint1D(rf, StepSize, x);
    };
}

tempTK T Calculus::SingleVar::Differentiation<T, K>::NthDiff(T (*differentiator)(const RealFunction1D &rf, K StepSize, K atPoint),
                                        const RealFunction1D &rf, K StepSize, K atPoint, int n) {

    if (n > 0) {

        auto func = [differentiator, StepSize, rf](K x) {
            return differentiator(rf, StepSize, x);
        };
        return NthDiff(differentiator, func, StepSize, atPoint, n - 1);
    } else {
        return rf(atPoint);
    }
}

tempTK std::function<T(K)> Calculus::SingleVar::Differentiation<T, K>::NthDiff(T (*differentiator)(const RealFunction1D &, K, K),
                                                          const RealFunction1D &rf, K StepSize, int n) {
    return [&differentiator, &rf, StepSize, n](K x) {
        return NthDiff(differentiator, rf, StepSize, x, n);
    };
}

tempTK K Calculus::SingleVar::Interpolation<T, K>::ChebyshevNode(int n, K a, K b, int samples) {
    Misc<T>::StartBeforeEnd(b, a);
    return 0.5*(a+b) + 0.5*(b-a)*cos(  (2*n - 1)*M_PI/((K)2*samples) );
}

tempTK std::function<T(K)> Calculus::SingleVar::Interpolation<T, K>::LagrangeCoefficient(std::vector<K> PointSet,int k_point) {

    return [PointSet,k_point](K x) {
        T answer = T(1);
        K k = PointSet.at(k_point);
        for (int i = 0; i < PointSet.size(); i++) {
            if (i != k_point) {
                answer = answer * ((x - PointSet.at(i)) / (PointSet.at(k_point) - PointSet.at(i)));
            }
        }
        return answer;
    };
}

tempTK std::function<T(K)> Calculus::SingleVar::Interpolation<T, K>::LagrangePolynomialFromDataSetWithNPoint (typename Sampler<T,K>::SampledFunction sampledFunction) {
    auto func = [sampledFunction](K x) {
        K answer = 0;
        sampledFunction.fitSize();
        sampledFunction.avoidinf();
        for (int i = 0; i < sampledFunction.x.size(); i++) {
            answer += (sampledFunction.y.at(i)) * (LagrangeCoefficient(sampledFunction.x,i)(x));
        }
        return answer;
    };
    return func;
}

tempTK std::function<T(K)>
Calculus::SingleVar::Interpolation<T, K>::LagrangePolynomialFromFunctionSimpleStep(const RealFunction1D& function1D, int samples,
                                                              K start, K end) {
    Misc<T>::StartBeforeEnd(start, end);
    return LagrangePolynomialFromDataSetWithNPoint(Calculus::SingleVar::Sampler<T,K>::FunctionSampler(function1D,start,end,samples));
}

tempTK std::function<T(K)> Calculus::SingleVar::Interpolation<T, K>::LagrangePolynomialFromFunctionChebyshevNodes(const RealFunction1D &function1D,
                                                                                             int samples, T start,
                                                                                             T end) {
    Misc<T>::StartBeforeEnd(start, end);
    return LagrangePolynomialFromDataSetWithNPoint(Calculus::SingleVar::Sampler<T,K>::FunctionSamplerChebyshevNodes(function1D,start,end,samples));
}

tempTK Algebric::MultiDimPoint<T>
Calculus::SingleVar::Interpolation<T, K>::LinearRegression2D_parameters(std::vector<Algebric::MultiDimPoint<T>> arr, T accuracy,
                                                   T(*diifer)(const RealFunction1D &function, K s, K p), int MaxTry) {
    if (arr.at(0).getDim() == 2) { //its 2d reg
        int NumberOfPoints = arr.size();
        T X = 0;
        T Y = 0;
        T MostHigh = 1;
        T MostLow = 1;
        for (int i = 0; i < NumberOfPoints; i++) {
            X += arr.at(i).getValue(0); // x values
            Y += arr.at(i).getValue(1); // y values
            if (arr.at(i).getValue(1) / arr.at(i).getValue(0) > MostHigh) {
                MostHigh = arr.at(i).getValue(1) / arr.at(i).getValue(0);
            } else if (arr.at(i).getValue(1) / arr.at(i).getValue(0) < MostLow) {
                MostLow = arr.at(i).getValue(1) / arr.at(i).getValue(0);
            }
        }

        RealFunction1D finalFunction_a = [&arr, &NumberOfPoints, &X, &Y](T a) {
            T FirstSigma = 0;
            T SecondSigma = 0;
            T temp;
            for (int i = 0; i < NumberOfPoints; i++) {
                temp = (arr.at(i).getValue(1) -
                        (a * arr.at(i).getValue(0) + (1 / double(NumberOfPoints)) * (Y - a * X)));
                FirstSigma += pow(temp, 2);
                SecondSigma += arr.at(i).getValue(0) * temp;
            }
            return ((a / (1 + pow(a, 2))) * FirstSigma + SecondSigma);
        };
        Algebric::MultiDimPoint<T> answer = Algebric::MultiDimPoint<T>(2);
        answer.setValue(EquationSolver<T, K>::Newton_Raphson(finalFunction_a, MostHigh, accuracy, MaxTry, diifer), 0);
        answer.setValue((1 / (T) (NumberOfPoints)) * (Y - answer.getValue(0) * X), 1);
        return answer;
    }
}


tempTK std::function<T(K)>
Calculus::SingleVar::Integration<T, K>::SimpleIndefinite(const RealFunction1D &function1D, K lowerbound, K stepSize) {
    return [function1D, lowerbound, stepSize](K x) {
        int chunks = std::abs(x - lowerbound) / stepSize;
        K answer = 0;
        for (int i = 0; i < chunks; i++) {
            answer += function1D(lowerbound + (i + 0.5) * stepSize) * (stepSize);
        }
        return answer;
    };
}

tempTK T Calculus::SingleVar::Integration<T, K>::SimpleDefinite(const RealFunction1D &function1D, K LowerBound, K UpperBound, K stepsize) {
    return Integration::SimpleIndefinite(function1D, LowerBound, stepsize)(UpperBound);
}

tempTK std::function<T(K)>
Calculus::SingleVar::Integration<T, K>::SimpsonIndefiniteClosed(const RealFunction1D &function1D, K LowerBound, K stepSize) {
    return [function1D, LowerBound, stepSize](K x) {
        int n = (x - LowerBound) / stepSize;
        K answer = 0;
        for (int i = 0; i < n; i++) {
            answer += (stepSize / (double) 6) * (function1D(LowerBound + i * stepSize)
                                                 + 4 * function1D(LowerBound + (i + 0.5) * stepSize)
                                                 + function1D(LowerBound + (i + 1) * stepSize));
        }
        return answer;
    };
}

tempTK T Calculus::SingleVar::Integration<T, K>::SimpsonDefiniteClosed(const RealFunction1D &function1D, K LowerBound, K UpperBound,
                                                  K stepSize) {
    return SimpsonIndefiniteClosed(function1D, LowerBound, stepSize)(UpperBound);
}


tempTK std::vector<K> Calculus::SingleVar::Sampler<T, K>::linespace(K a, K b, int n) {
    std::vector<K> lin;
    lin.reserve(n);
    for (int i = 0; i <= n; ++i) {
        lin.insert(lin.begin() + i, a + (i * ((b-a)/ n)));
    }
    return lin;
}

tempTK typename Calculus::SingleVar::Sampler<T, K>::SampledFunction
Calculus::SingleVar::Sampler<T, K>::FunctionSampler(const std::function<T(K)> &rf, K a, K b, int n) {
    SampledFunction s;
    s.x = linespace(a, b, n);
    s.y.reserve(n);
    for (int i = 0; i < s.x.size(); ++i) {
        s.y.insert(s.y.begin() + i, rf(s.x.at(i)));
    }
    return s;
}
tempTK void Calculus::SingleVar::Sampler<T,K>::SampledFunction::fitSize() const{
    x.resize(std::min(x.size(),y.size()));
    y.resize(std::min(x.size(),y.size()));
}

tempTK typename Calculus::SingleVar::Sampler<T,K>::SampledFunction Calculus::SingleVar::Sampler<T,K>::FunctionSamplerChebyshevNodes(const std::function<T(K)> &rf,
                                                                                        K a, K b, int n) {
    SampledFunction s;
    K p;
    for (int i = 1; i <= n; ++i) {
        p = Calculus::SingleVar::Interpolation<T,K>::ChebyshevNode(i,a,b,n);
        s.x.push_back(p);
        s.y.push_back(rf(p));
    }
    return s;
}

tempTK void Calculus::SingleVar::Sampler<T,K>::SampledFunction::avoidinf() const {
    for (int i = 0; i < std::min(y.size(),x.size()); ++i) {
        if(std::isinf(x.at(i)) || std::isinf(y.at(i))){
            x.erase(x.begin()+ i);
            y.erase(y.begin()+i);
        }
    }
}
#endif