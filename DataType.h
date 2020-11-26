#ifndef IMAGE_MESSURE_DATATYPE_H
#define IMAGE_MESSURE_DATATYPE_H
#include <ostream>
#include <functional>
#include <vector>
#include "iomanip"

namespace Algebric{
    template <typename T> class Matrix;
    template <typename T> class MultiDimPoint {
    private:
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
        const T& operator()(const int &index)const{
            return getValue(index);
        };
        T& operator()(const int &index){
            return (mdp.at(index));
        };
        template<typename K> friend MultiDimPoint<T> operator*(const K& k,const MultiDimPoint<T> &o){
            std::vector<T> t = o.getVec();
            for(T &s:t)
                s = k*s;
            return t;
        }
        MultiDimPoint<T> operator+(const MultiDimPoint<T> &vec2)const{
            if(getDim() == vec2.getDim())
            {
                std::vector<T> res;
                for (int i = 0; i < getDim(); ++i) {
                    res.push_back(getValue(i)+vec2.getValue(i));
                }
                return res;
            }
        }
        MultiDimPoint<T> operator-(const MultiDimPoint<T> &vec2)const {
            return operator+(-1*vec2);
        }
        Matrix<T> operator^(const MultiDimPoint<T>& vec2)const;

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
        MultiDimPoint<T> operator*(const MultiDimPoint<T>& rhs);
        std::vector<T> diag_vec();

        // Access the individual elements
        T& operator()(const unsigned& row, const unsigned& col);
        const T& operator()(const unsigned& row, const unsigned& col) const;

        // Access the row and column sizes
        unsigned get_rows() const;
        unsigned get_cols() const;

        void pretty_row(int which,std::ostream& stream)const{
            using namespace std;
            stream<<setprecision(3)<<"|  ";
            for (int i = 0; i < get_cols(); ++i) {
                if(mat[which][i]>=0){
                    stream<<setw(6)<<left<<mat[which][i]<<"  ";
                    if(i == get_cols() -1)
                        stream<<"\b\b\b";
                }
                else {
                    stream << '\b' << setw(6) << left << mat[which][i] << "   ";
                    if(i == get_cols() -1)
                        stream<<"\b\b\b";
                }

            }
            stream<<"|";
        }
        // better print!
         friend std::ostream& operator<<(std::ostream& stream,const Matrix& m){
            for(int i=0;i<m.get_rows();i++){
                m.pretty_row(i,stream);
                stream<<"\n";
            }
            return stream;
        }

    };
};

#include "DataType.cpp"
#endif //MATHUTILS_DATATYPE_H
