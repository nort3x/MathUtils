#ifndef MATHUTILS_LIBRARY_H
#define MATHUTILS_LIBRARY_H
#include "DataType.h"
#define tempT template<typename T>
#define tempTK template<typename T,typename K>

namespace Calculus{
    tempT class Misc{
    public:
        static void StartBeforeEnd(T &start,T &end);
    };
    namespace SingleVar{

        tempTK class Sampler{
        public:
            struct SampledFunction{
                mutable std::vector<T> y;
                mutable std::vector<K> x;
                void fitSize() const;
                void avoidinf() const;
            };

            static std::vector<K> linespace(K a,K b,int n);
            static SampledFunction FunctionSampler(const std::function<T(K)> &rf,K a,K b,int n);
            static SampledFunction FunctionSamplerChebyshevNodes(const std::function<T(K)> &rf,K a,K b,int n);
        };

        tempTK class Differentiation{

            typedef std::function<T(K)> RealFunction1D;
        public:
            static RealFunction1D Forward1D(const RealFunction1D& rf,K StepSize);
            static T Forward1D(const RealFunction1D& rf,K StepSize,K atPoint);

            static T Backward1D(const RealFunction1D& rf,K StepSize,K atPoint);
            static RealFunction1D Backward1D(const RealFunction1D& rf,K StepSize);

            static T ThreePointMidPoint1D(const RealFunction1D& rf,K StepSize,K atPoint);
            static RealFunction1D ThreePointMidPoint1D(const RealFunction1D& rf,K StepSize);

            static T FivePointMidPoint1D(const RealFunction1D& rf,K StepSize,K atPoint);
            static RealFunction1D FivePointMidPoint1D(const RealFunction1D& rf,K StepSize);

            static T NthDiff(T (*differentiator)(const RealFunction1D& rf,K StepSize,K atPoint),const RealFunction1D& rf,K StepSize,K atPoint,int n);
            static RealFunction1D NthDiff(T (*differentiator)(const RealFunction1D& rf,K StepSize,K atPoint),const RealFunction1D& rf,K StepSize,int n);
        };

        tempTK class Interpolation{
            typedef std::function<T(K)> RealFunction1D;
        private:
            static RealFunction1D LagrangeCoefficient(std::vector<K> PointSet,int k_point);
        public:
            static K ChebyshevNode(int n,K a ,K b,int samples);
            static RealFunction1D LagrangePolynomialFromDataSetWithNPoint(typename Sampler<T,K>::SampledFunction sampledFunction);
            static RealFunction1D LagrangePolynomialFromFunctionSimpleStep(const RealFunction1D& function1D,int samples,K start,K end);
            static RealFunction1D LagrangePolynomialFromFunctionChebyshevNodes(const RealFunction1D& realFunction1D,int samples,T start,T end);
            static Algebric::MultiDimPoint<T>  LinearRegression2D_parameters(std::vector<Algebric::MultiDimPoint<T>> arr,T accuracy,T(*diifer)(const RealFunction1D& function, K s, K p),int MaxTry);
        };

        tempTK class Integration{
            typedef std::function<T(K)> RealFunction1D;
        public:
            static RealFunction1D SimpleIndefinite(const RealFunction1D& function1D,K LowerBound,K stepSize);
            static T SimpleDefinite(const RealFunction1D& function1D,K LowerBound,K UpperBound,K stepsize);

            static RealFunction1D SimpsonIndefiniteClosed(const RealFunction1D& function1D,K LowerBound,K stepSize);
            static T SimpsonDefiniteClosed(const RealFunction1D& function1D,K LowerBound,K UpperBound,K stepSize);
        };

        tempTK class EquationSolver{
            typedef std::function<T(K)> RealFunction1D;
        public:
            static K BisectionZeroFinder(RealFunction1D& function1D, K a, K b, T accuracy);
            static K Newton_Raphson(RealFunction1D &function1D, K init_p, T accuracy,int maxtry,T(*diifer)(const RealFunction1D& function, K stepsize, K point));

        };


    }
    namespace MultiVar{
        namespace Scalar{
            tempTK class Sampler{
            public:
                struct SampledScalarFunction{
                    mutable std::vector<std::vector<K>> variables;
                    mutable std::vector<T> function;
                    mutable std::vector<int> numberOfsamples;
                    void fitSize() const;
                    void avoidinf() const;
                };
                struct linspace { K a;K b ; int n;};
                static std::vector<std::vector<K>> MlineSpace(const std::initializer_list<const linspace>&);
                static SampledScalarFunction FunctionSampler(const std::function<T(Algebric::MultiDimPoint<K>)> &rf,const std::initializer_list<const linspace>&);
                static SampledScalarFunction FunctionSamplerChebyshevNodes(const std::function<T(Algebric::MultiDimPoint<K>)> &rf,std::initializer_list<const linspace>&);
            };
            namespace Differentiation{
                tempTK std::function<T(Algebric::MultiDimPoint<K>)> Partial(const std::function<T(Algebric::MultiDimPoint<K>)>&,
                        T(*Differ)(const std::function<T(K)>& rf,K StepSize,K atPoint),int Index,K accuracy);
            }

        }
    }

}
#include "library.cpp"
#endif //MATHUTILS_LIBRARY_H
