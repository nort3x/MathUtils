//
// Created by root on 11/2/20.
//

#ifndef MATHUTILS_DATATYPE_H
#define MATHUTILS_DATATYPE_H
#include <ostream>
#include <functional>
#include <vector>

namespace Algebric{

    template <typename T> class MultiDimPoint {
        std::vector<T> mdp;
        int dim;
    public:
        MultiDimPoint(const std::initializer_list<T> &arg);
        MultiDimPoint(const int &dim);
        MultiDimPoint();
        MultiDimPoint(const std::vector<T>&);

        MultiDimPoint& setValue(T p,const int& index);
        T getValue(int index) const;
        int getDim() const;
        std::vector<T> getVec()const;
        const T& operator[](const int &index)const{
            return getValue(index);
        };
        T& operator[](const int &index){
            return (mdp.at(index));
        };

        T Norm()const;

        friend std::ostream& operator<<(std::ostream& o ,const MultiDimPoint& mdPoint){
            o<<"(";
            for(int i=0 ; i<mdPoint.getDim();i++){
                if(i!=mdPoint.getDim()-1){
                    o<<mdPoint.getValue(i) <<" ,";
                } else{
                    o<<mdPoint.getValue(i)<<")";
                }
            }
            return o;
        };
    };

    template <typename T> class Matrix  {
    private:
        std::vector<std::vector<T> > mat;
        unsigned rows;
        unsigned cols;

    public:
        Matrix(unsigned _rows, unsigned _cols, const T& _initial);
        Matrix(const Matrix<T>& rhs);
        virtual ~Matrix();

        // Operator overloading, for "standard" mathematical matrix operations
        Matrix<T>& operator=(const Matrix<T>& rhs);

        // Matrix mathematical operations
        Matrix<T> operator+(const Matrix<T>& rhs);
        Matrix<T>& operator+=(const Matrix<T>& rhs);
        Matrix<T> operator-(const Matrix<T>& rhs);
        Matrix<T>& operator-=(const Matrix<T>& rhs);
        Matrix<T> operator*(const Matrix<T>& rhs);
        Matrix<T>& operator*=(const Matrix<T>& rhs);
        Matrix<T> transpose();

        // Matrix/scalar operations
        Matrix<T> operator+(const T& rhs);
        Matrix<T> operator-(const T& rhs);
        Matrix<T> operator*(const T& rhs);
        Matrix<T> operator/(const T& rhs);

        // Matrix/vector operations
        std::vector<T> operator*(const std::vector<T>& rhs);
        std::vector<T> diag_vec();

        // Access the individual elements
        T& operator()(const unsigned& row, const unsigned& col);
        const T& operator()(const unsigned& row, const unsigned& col) const;

        // Access the row and column sizes
        unsigned get_rows() const;
        unsigned get_cols() const;

        // better print!
         friend std::ostream& operator<<(std::ostream& stream,const Matrix& m){
            for(int i=0;i<m.get_rows();i++){
                stream<<"| ";
                for(int j=0;j<m.get_cols();j++){
                    stream<<m(i,j);
                    if(j == m.get_cols()-1){
                        stream<<" |\n";
                    } else{
                        stream<<", ";
                    }
                }
            }
            return stream;
        }

    };
};

#include "DataType.cpp"
#endif //MATHUTILS_DATATYPE_H
