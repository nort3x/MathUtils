//
// Created by root on 11/7/20.
//

#ifndef MATHUTILS_UTILS_H
#define MATHUTILS_UTILS_H

#include <ostream>
#include "chrono"
#include "iomanip"

namespace Utils{
    class Timer{
       mutable std::chrono::time_point<std::chrono::system_clock, std::chrono::duration<double>> time;
    public:
        Timer(){
            time = std::chrono::system_clock::now();
        }
        std::chrono::duration<double> readNow() const{

            auto end = std::chrono::system_clock::now();
            auto t = end - time;
            reset();
            return t ;
        }
        void reset() const{
            time = std::chrono::system_clock::now();
        }
        friend std::ostream& operator<<(std::ostream& stream,const Timer& t){
            stream<<std::setprecision(10)<<t.readNow().count()<<" s\n";
            return stream;
        }
    };

    tempT
    class BenchMark{
    public:
        static double methodDuration( std::function<void(T)> func,T data){
            Timer t;
            func(data);
            return t.readNow().count();
        }

        static std::function<double(T)> methodDurationFunction(std::function<void(T)> func){
            return [func](T data)->double{
                return methodDuration(func,data);
            };
        }

    };


}
#endif //MATHUTILS_UTILS_H
