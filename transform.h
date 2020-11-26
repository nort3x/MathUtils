//
// Created by root on 11/21/20.
//

#ifndef IMAGE_MESSURE_TRANSFORM_H
#define IMAGE_MESSURE_TRANSFORM_H
#include "DataType.h"
#include "cmath"

#ifndef tempT
#define tempT template<typename T>
#define tempTK template<typename T,typename K>
#endif
#define Vector(T) Algebric::MultiDimPoint<T>

namespace trns{
    tempT T Dot(const Vector(T) &v1,const Vector(T) &v2){
        T ans = 0;
        for (int i = 0; i < std::min(v1.getDim(),v2.getDim()); ++i) {
            T v = v1.getValue(i)*v2.getValue(i);
            ans += v;
        }
        return ans;
    }
    tempT double Angle_deg(const Vector(T) &v1,const Vector(T) &v2){
        double d = Dot(v1,v2);
        double k = (v1.Norm())*(v2.Norm());
        double m =  d/k;
        return (180/M_PI )*std::acos( m);
    }

    template <typename T> class RotationTransform{
        Algebric::Matrix<T> t = Algebric::Matrix<T>(3,3,0);
    public:
        RotationTransform(){
            t(0,0) =1;
            t(1,1) =1;
            t(2,2) =1;
        }
        Algebric::Matrix<T> build(){
            return t;
        }

        void rotate_y(T radian){
            Algebric::Matrix<T> temp(3,3,0);
            temp(0,0) = cos(radian);
            temp(1,1) = 1;
            temp(0,2) = sin(radian);
            temp(2,0) = -temp(0,2);
            temp(2,2) = temp(0,0);
            t = temp*t;
        }



        void rotate_z(T radian){
            Algebric::Matrix<T> temp(3,3,0);
            temp(0,0) = cos(radian);
            temp(1,1) = temp(0,0);
            temp(0,1) = sin(radian);
            temp(1,0) = -temp(0,1);
            temp(2,2) = 1;
            t = (temp.operator*(t));

        }


        void rotate_x(T radian){
            Algebric::Matrix<T> temp(3,3,0);
            temp(0,0) = 1;
            temp(1,1) = cos(radian);
            temp(1,2) = sin(radian);
            temp(2,1) = -temp(1,2);
            temp(2,2) = temp(2,2);
            t = temp*t;
        }
    };
    template <typename T> class ProjectionTransform{
        Algebric::Matrix<T> t = Algebric::Matrix<T>(3,3,0);
    public:
        void ProjectOnVector(const Algebric::MultiDimPoint<T>& vector){
            t = (vector^vector)*(1/std::pow(vector.Norm(),2));
        }
        void ProjectOnPlane(const Algebric::MultiDimPoint<T>& perpendicularVec){
            Algebric::Matrix<T> I(3,3,0);
            I(0,0) = 1;
            I(1,1) = 1;
            I(2,2) = 1;
            t = I - ((perpendicularVec^perpendicularVec)*(1/std::pow(perpendicularVec.Norm(),2)));
        }
        Algebric::Matrix<T> build(){
            return t;
        }
    };

}
#undef Vector
#undef tempT
#undef tempTK
#endif //IMAGE_MESSURE_TRANSFORM_H
