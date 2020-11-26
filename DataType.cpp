#ifndef MATHUTILS_DATATYPE_CPP
#define MATHUTILS_DATATYPE_CPP

#include <cmath>
#include "DataType.h"
#ifndef tempT
#define tempT template<typename T>
#define tempTK template<typename T,typename K>
#endif

tempT Algebric::MultiDimPoint<T>::MultiDimPoint(const int& d) {
    mdp.reserve(d);
     dim = d;
}

tempT Algebric::MultiDimPoint<T>::MultiDimPoint() {
    dim = 0;
}
tempT Algebric::MultiDimPoint<T>::MultiDimPoint(const std::initializer_list<T>& arg) {
    mdp.reserve(arg.size());
    this->dim = arg.size();
    for(auto p: arg){
        mdp.push_back(p);
    }
}
tempT Algebric::MultiDimPoint<T>& Algebric::MultiDimPoint<T>::setValue(T p, const int& index) {
    mdp.insert(mdp.begin()+index,p);
    return *this;
}
tempT T Algebric::MultiDimPoint<T>::getValue(int index) const {
    return mdp.at(index);
}
tempT int Algebric::MultiDimPoint<T>::getDim() const {
    return dim;
}

template<typename T>
Algebric::MultiDimPoint<T>::MultiDimPoint(const std::vector<T>& vec ) {
    mdp = vec;
    dim = vec.size();
}

tempT T Algebric::MultiDimPoint<T>::Norm() const{
    T t =0;
    for(auto i: mdp){
        t += std::pow(i,2);
    }
    return std::sqrt(t);
}
tempT std::vector<T> Algebric::MultiDimPoint<T>::getVec() const {
    return mdp;
}
tempT Algebric::Matrix<T> Algebric::MultiDimPoint<T>::operator^(const MultiDimPoint <T> &vec2)const{
    // v^w = v * wt
    Matrix<T> m(getDim(),vec2.getDim(),0);
    for (int i = 0; i <getDim() ; ++i) {
        for (int j = 0; j < vec2.getDim(); ++j) {
            m(i,j) = getValue(i)*vec2.getValue(j);
        }
    }
    return m;
};

///////////// matrix///////////////


// Parameter Constructor
template<typename T>
Algebric::Matrix<T>::Matrix(unsigned _rows, unsigned _cols, const T& _initial) {
    mat.resize(_rows);
    for (unsigned i=0; i<mat.size(); i++) {
        mat[i].resize(_cols, _initial);
    }
    rows = _rows;
    cols = _cols;
}

// Copy Constructor
template<typename T>
Algebric::Matrix<T>::Matrix(const Algebric::Matrix<T>& rhs) {
    mat = rhs.mat;
    rows = rhs.get_rows();
    cols = rhs.get_cols();
}

// (Virtual) Destructor
template<typename T>
Algebric::Matrix<T>::~Matrix() {}

// Assignment Operator
template<typename T>
Algebric::Matrix<T>& Algebric::Matrix<T>::operator=(const Algebric::Matrix<T>& rhs) {
    if (&rhs == this)
        return *this;

    unsigned new_rows = rhs.get_rows();
    unsigned new_cols = rhs.get_cols();

    mat.resize(new_rows);
    for (unsigned i=0; i<mat.size(); i++) {
        mat[i].resize(new_cols);
    }

    for (unsigned i=0; i<new_rows; i++) {
        for (unsigned j=0; j<new_cols; j++) {
            mat[i][j] = rhs(i, j);
        }
    }
    rows = new_rows;
    cols = new_cols;

    return *this;
}

// Addition of two matrices
template<typename T>
Algebric::Matrix<T> Algebric::Matrix<T>::operator+(const Algebric::Matrix<T>& rhs) {
    Matrix result(rows, cols, 0.0);

    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result(i,j) = this->mat[i][j] + rhs(i,j);
        }
    }

    return result;
}

// Cumulative addition of this matrix and another
template<typename T>
Algebric::Matrix<T>& Algebric::Matrix<T>::operator+=(const Algebric::Matrix<T>& rhs) {
    unsigned rows = rhs.get_rows();
    unsigned cols = rhs.get_cols();

    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            this->mat[i][j] += rhs(i,j);
        }
    }

    return *this;
}

// Subtraction of this matrix and another
template<typename T>
Algebric::Matrix<T> Algebric::Matrix<T>::operator-(const Algebric::Matrix<T>& rhs) {
    unsigned rows = rhs.get_rows();
    unsigned cols = rhs.get_cols();
    Matrix result(rows, cols, 0.0);

    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result(i,j) = this->mat[i][j] - rhs(i,j);
        }
    }

    return result;
}

// Cumulative subtraction of this matrix and another
template<typename T>
Algebric::Matrix<T>& Algebric::Matrix<T>::operator-=(const Algebric::Matrix<T>& rhs) {
    unsigned rows = rhs.get_rows();
    unsigned cols = rhs.get_cols();

    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            this->mat[i][j] -= rhs(i,j);
        }
    }

    return *this;
}

// Left multiplication of this matrix and another
template<typename T>
Algebric::Matrix<T> Algebric::Matrix<T>::operator*(const Algebric::Matrix<T>& rhs) {
    unsigned rows = rhs.get_rows();
    unsigned cols = rhs.get_cols();
    Matrix result(rows, cols, 0.0);

    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            for (unsigned k=0; k<rows; k++) {
                result(i,j) += this->mat[i][k] * rhs(k,j);
            }
        }
    }

    return result;
}

// Cumulative left multiplication of this matrix and another
template<typename T>
Algebric::Matrix<T>& Algebric::Matrix<T>::operator*=(const Algebric::Matrix<T>& rhs) {
   Matrix result = (*this) * rhs;
    (*this) = result;
    return *this;
}

// Calculate a transpose of this matrix
template<typename T>
Algebric::Matrix<T> Algebric::Matrix<T>::transpose() {
    Matrix result(rows, cols, 0.0);

    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result(i,j) = this->mat[j][i];
        }
    }

    return result;
}

// Matrix/scalar addition
template<typename T>
Algebric::Matrix<T> Algebric::Matrix<T>::operator+(const T& rhs) {
    Matrix result(rows, cols, 0.0);

    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result(i,j) = this->mat[i][j] + rhs;
        }
    }

    return result;
}

// Matrix/scalar subtraction
template<typename T>
Algebric::Matrix<T> Algebric::Matrix<T>::operator-(const T& rhs) {
   Matrix result(rows, cols, 0.0);

    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result(i,j) = this->mat[i][j] - rhs;
        }
    }

    return result;
}

// Matrix/scalar multiplication
template<typename T>
Algebric::Matrix<T> Algebric::Matrix<T>::operator*(const T& rhs) {
    Matrix result(rows, cols, 0.0);

    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result(i,j) = this->mat[i][j] * rhs;
        }
    }

    return result;
}

// Matrix/scalar division
template<typename T>
Algebric::Matrix<T> Algebric::Matrix<T>::operator/(const T& rhs) {
    Matrix result(rows, cols, 0.0);

    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result(i,j) = this->mat[i][j] / rhs;
        }
    }

    return result;
}

// Multiply a matrix with a vector
template<typename T>
std::vector<T> Algebric::Matrix<T>::operator*(const std::vector<T>& rhs) {
    std::vector<T> result(rhs.size(), 0.0);

    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result[i] += this->mat[i][j] * rhs[j];
        }
    }

    return result;
}
tempT Algebric::MultiDimPoint<T> Algebric::Matrix<T>::operator*(const Algebric::MultiDimPoint<T>& vec){
    return operator*(vec.getVec());
}

// Obtain a vector of the diagonal elements
template<typename T>
std::vector<T> Algebric::Matrix<T>::diag_vec() {
    std::vector<T> result(rows, 0.0);

    for (unsigned i=0; i<rows; i++) {
        result[i] = this->mat[i][i];
    }

    return result;
}

// Access the individual elements
template<typename T>
T& Algebric::Matrix<T>::operator()(const unsigned& row, const unsigned& col) {
    return this->mat[row][col];
}

// Access the individual elements (const)
template<typename T>
const T& Algebric::Matrix<T>::operator()(const unsigned& row, const unsigned& col) const {
    return this->mat[row][col];
}

// Get the number of rows of the matrix
template<typename T>
unsigned Algebric::Matrix<T>::get_rows() const {
    return this->rows;
}

// Get the number of columns of the matrix
template<typename T>
unsigned Algebric::Matrix<T>::get_cols() const {
    return this->cols;
}

#undef tempT
#undef tempTK
#endif
