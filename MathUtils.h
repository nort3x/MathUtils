#ifndef MATHUTILS_MATHUTILS_H
#define MATHUTILS_MATHUTILS_H
#include "DataType.h"
#ifndef tempT
#define tempT template<typename T>
#define tempTK template<typename T,typename K>
#endif

namespace Calculus{
    tempTK using Function1D = std::function<T(K)>;
    tempTK using differ = T(*)(const Function1D<T,K>& rf,K accuracy,K atPoint);
    tempTK using intergrator = T(*)(const Function1D<T,K>& function1D,K LowerBound,K UpperBound,K stepsize);

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
            static SampledFunction FunctionSampler(const Function1D<T,K> &rf,K a,K b,int n);
            static SampledFunction FunctionSamplerChebyshevNodes(const Function1D<T,K> &rf,K a,K b,int n);
        };
        tempTK class Differentiation{


        public:
            static Function1D<T,K> Forward1D(const Function1D<T,K>& rf,K StepSize);
            static T Forward1D(const Function1D<T,K>& rf,K StepSize,K atPoint);

            static T Backward1D(const Function1D<T,K>& rf,K StepSize,K atPoint);
            static Function1D<T,K> Backward1D(const Function1D<T,K>& rf,K StepSize);

            static T ThreePointMidPoint1D(const Function1D<T,K>& rf,K StepSize,K atPoint);
            static Function1D<T,K> ThreePointMidPoint1D(const Function1D<T,K>& rf,K StepSize);

            static T FivePointMidPoint1D(const Function1D<T,K>& rf,K StepSize,K atPoint);
            static Function1D<T,K> FivePointMidPoint1D(const Function1D<T,K>& rf,K StepSize);

            static T NthDiff(T(*differentiator)(const Function1D<T,K>& rf,K StepSize,K atPoint),const Function1D<T,K>& rf,K StepSize,K atPoint,int n);
            static Function1D<T,K> NthDiff(T (*differentiator)(const Function1D<T,K>& rf,K StepSize,K atPoint),const Function1D<T,K>& rf,K StepSize,int n);
        };
        tempTK class Interpolation{

        private:
            static Function1D<T,K> LagrangeCoefficient(std::vector<K> PointSet,int k_point);
        public:
            static K ChebyshevNode(int n,K a ,K b,int samples);
            static Function1D<T,K> LagrangePolynomialFromDataSetWithNPoint(typename Sampler<T,K>::SampledFunction sampledFunction);
            static Function1D<T,K> LagrangePolynomialFromFunctionSimpleStep(const Function1D<T,K>& function1D,int samples,K start,K end);
            static Function1D<T,K> LagrangePolynomialFromFunctionChebyshevNodes(const Function1D<T,K>& Function1D,int samples,T start,T end);
            static Algebric::MultiDimPoint<T>  LinearRegression2D_parameters(std::vector<Algebric::MultiDimPoint<T>> arr,T accuracy,T(*diifer)(const Function1D<T,K>& function, K s, K p),int MaxTry);
        };
        tempTK class Integration{

        public:
            static Function1D<T,K> SimpleIndefinite(const Function1D<T,K>& function1D,K LowerBound,K stepSize);
            static T SimpleDefinite(const Function1D<T,K>& function1D,K LowerBound,K UpperBound,K stepsize);

            static Function1D<T,K> SimpsonIndefiniteClosed(const Function1D<T,K>& function1D,K LowerBound,K stepSize);
            static T SimpsonDefiniteClosed(const Function1D<T,K>& function1D,K LowerBound,K UpperBound,K stepSize);
        };
        tempTK class EquationSolver{

        public:
            static K BisectionZeroFinder(Function1D<T,K>& function1D, K a, K b, T accuracy);
            static K Newton_Raphson(Function1D<T,K> &function1D, K init_p, T accuracy,int maxtry,T(*diifer)(const Function1D<T,K>& function, K stepsize, K point));

        };
        tempTK class ODE{
        public:
            static std::function<T(K)> EulerMethodFirstOrder(const std::function<T(T,K)> &F_y_t ,K t0,T y0, K stepsize);
            static std::function<T(K)> Runge_KuttaFirstOrder(const std::function<T(T,K)> &F_y_t,K t0,T y0,K stepsize);

            static std::function<T(K)> EulerMethodSecondOrder(const std::function<T(T,T,K)> &F_yp_y_t ,K t0,T yp0,T y0, K stepsize);

        };

    }
    namespace MultiVar{
        tempT struct linspace { T a;T b ; int n;};
        tempT std::vector<std::vector<T>> MlineSpace(const std::initializer_list<const linspace<T>>&);
        tempTK std::function<T(Algebric::MultiDimPoint<K>)> withConstIndex(const std::function<T(Algebric::MultiDimPoint<K>)> &rf,int index,K value);
        namespace Scalar{
            tempTK std::function<T(K)> asSingleVar(const std::function<T(Algebric::MultiDimPoint<K>)> &rf,int index,const Algebric::MultiDimPoint<K> &init_vec);
            tempTK class Sampler{
            public:
                struct SampledScalarFunction{
                    mutable std::vector<std::vector<K>> variables;
                    mutable std::vector<T> function;
                    mutable std::vector<int> numberOfsamples;
                    void fitSize() const;
                    void avoidinf() const;
                };
                static SampledScalarFunction FunctionSampler(const std::function<T(Algebric::MultiDimPoint<K>)> &rf,const std::initializer_list<const linspace<K>>&);
                static SampledScalarFunction FunctionSamplerChebyshevNodes(const std::function<T(Algebric::MultiDimPoint<K>)> &rf,std::initializer_list<const linspace<K>>&);
            };
            namespace Differentiation{
                tempTK std::function<T(Algebric::MultiDimPoint<K>)> Partial(
                        const std::function<T(Algebric::MultiDimPoint<K>)> &rf, differ<T,K> d,int Index, K accuracy);
                tempTK std::function<Algebric::MultiDimPoint<T>(Algebric::MultiDimPoint<K>)> Gradiant(const std::function<T(Algebric::MultiDimPoint<K>)>&,
                                                                             differ<T,K> d,K accuracy);
            }
            namespace Integration{
               tempTK std::function<T(Algebric::MultiDimPoint<K>)> Partial(const std::function<T(Algebric::MultiDimPoint<K>)> &rf,const intergrator<T,K> &d,int Index, K accuracy,K lowerBand);
            }

        }
        namespace VectorValued{
            tempTK class Sampler{
            public:
                struct SampledVectorFunction{
                    mutable std::vector<std::vector<K>> variables;
                    mutable std::vector<std::vector<T>> function;
                    mutable std::vector<int> numberOfsamples;
                    void fitSize() const;
                    void avoidinf() const;
                };
                struct Curve{
                    struct range{
                        T min;
                        T max;
                    };
                    mutable std::vector<K> variable;
                    mutable std::vector<std::vector<T>> function;
                    mutable std::vector<range> ranges;
                    mutable std::vector<int> numberOfsamples;
                    void fitSize() const;
                    void avoidinf() const;
                };
                static SampledVectorFunction FunctionSampler(const std::function<Algebric::MultiDimPoint<T>(Algebric::MultiDimPoint<K>)> &rf,const std::initializer_list<const linspace<K>>&);
                static SampledVectorFunction FunctionSamplerChebyshevNodes(const std::function<Algebric::MultiDimPoint<T>(Algebric::MultiDimPoint<K>)> &rf,const std::initializer_list<const linspace<K>>&);

                static Curve CurveSampler(const std::function<Algebric::MultiDimPoint<T>(K)> &rf,const K& a,const K& b,int n);
            };
        }
    }

}
#undef tempT
#undef tempTK
#include "MathUtils.cpp"
#endif //MATHUTILS_MATHUTILS_H
