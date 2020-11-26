//
// Created by root on 11/7/20.
//

#ifndef MATHUTILS_UTILS_H
#define MATHUTILS_UTILS_H

#include <iostream>
#include <fstream>
#include "chrono"
#include "iomanip"

#ifndef tempT
#define tempT template<typename T>
#define tempTK template<typename T,typename K>
#endif

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
    class ScopeTimer: Timer{
        std::ostream* s;
    public:
        ScopeTimer(std::ostream& strm ){
            this->s = &strm;
        }

        ~ScopeTimer(){
            (*s)<< *this;
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

    tempT std::string WriteDataToFile(const std::initializer_list<Algebric::MultiDimPoint<T>> &args,const std::string &path="",const int& precision=10){
        std::fstream file;
        std::string rpath;

        file<<std::setprecision(precision);
        if (path.empty()){
            rpath = tempnam("1","1234567890");
            file.open(rpath,std::ios::out);
        } else{
            rpath = path;
            file.open(rpath,std::ios::out);
        }
        if(file.is_open()){
            file<<"#DATA BEGINS\n";
            for(auto vec : args){
                for (int i = 0; i < vec.getDim(); ++i) {
                    file<<vec.getValue(i)<<"\t";
                }
                file<<"\n";
            }

            file<<"#DATA ENDS\n\n";
        }
        file.close();
        return rpath;

    }

    template<typename T,size_t N> std::string WriteDataToFile(const std::vector<T> (&data)[N],const std::string &path="",const int& precision=10){
        std::fstream file;
        std::string rpath;

        file<<std::setprecision(precision);
        if (path.empty()){
            rpath = tempnam("1","1234567890");
            file.open(rpath,std::ios::out);
        } else{
            rpath = path;
            file.open(rpath,std::ios::out);
        }
        if(file.is_open()){
            file<<"#DATA BEGINS\n";
            int min_n_data = data[0].size();
            for (int i = 0; i < N; ++i) {
                min_n_data = std::min(min_n_data,(int)data[i].size());
            }

            for(int i=0;i<min_n_data;i++) {
                for (int j = 0; j < N; ++j) {
                    file<<data[j].at(i)<<"\t";
                }
                file<<"\n";
            }
            file<<"#DATA ENDS\n\n";
        }
        file.close();
        return rpath;
    }

    tempTK std::string WriteDataToFile(const typename Calculus::SingleVar::Sampler<T,K>::SampledFunction s,
                                       const std::string &path="",const int& precision=10){

        std::fstream file;
        std::string rpath;
        file<<std::setprecision(10);

        if (path.empty()){
            rpath = tempnam("1","1234567890");
            file.open(rpath,std::ios::out);
        } else{
            rpath = path;
            file.open(rpath,std::ios::out);
        }
        if(file.is_open()){
            file<<"#DATA BEGINS\n";
            s.fitSize();

            for (int i = 0; i < s.y.size(); ++i) {
                file<<s.x.at(i)<<"\t";
                file<<s.y.at(i)<<"\n";
            }
            file<<"#DATA ENDS\n\n";
        }
        file.close();
        return rpath;
    }
    tempTK std::string WriteDataToFile(const typename Calculus::MultiVar::Scalar::Sampler<T,K>::SampledScalarFunction s,
                                       const std::string &path="",const int& precision=10){

        std::fstream file;
        std::string rpath;
        file<<std::setprecision(10);

        if (path.empty()){
            rpath = tempnam("1","1234567890");
            file.open(rpath,std::ios::out);
        } else{
            rpath = path;
            file.open(rpath,std::ios::out);
        }
        if(file.is_open()){
            file<<"#DATA BEGINS\n";
            s.fitSize();

            for (int i = 0; i < s.function.size(); ++i) {
                if( i%(s.numberOfsamples.at(0)+1) == 0 )
                    file<<"\n";
                for(auto cor:s.variables)
                    file<<cor.at(i)<<"\t";
                file<<s.function.at(i)<<"\n";


            }
            file<<"#DATA ENDS\n\n";
        }
        file.close();
        return rpath;
    }
    tempTK std::string WriteDataToFile(const typename Calculus::MultiVar::VectorValued::Sampler<T,K>::SampledVectorFunction s,
                                       const std::string &path="",const int& precision=10){

        std::fstream file;
        std::string rpath;
        file<<std::setprecision(10);

        if (path.empty()){
            rpath = tempnam("1","1234567890");
            file.open(rpath,std::ios::out);
        } else{
            rpath = path;
            file.open(rpath,std::ios::out);
        }
        if(file.is_open()){
            file<<"#DATA BEGINS\n";
            s.fitSize();

            for (int i = 0; i < s.function.size(); ++i) {
                if( i%(s.numberOfsamples.at(0)+1) == 0 )
                    file<<"\n";
                for(auto cor:s.variables)
                    file<<cor.at(i)<<"\t";
                for(auto p:s.function.at(i))
                    file<<p<<"\t";
                file<<"\n";


            }
            file<<"#DATA ENDS\n\n";
        }
        file.close();
        return rpath;
    }
    tempTK std::string WriteDataToFile(const typename Calculus::MultiVar::VectorValued::Sampler<T,K>::Curve s,
                                       const std::string &path="",const int& precision=10){

        std::fstream file;
        std::string rpath;
        file<<std::setprecision(10);

        if (path.empty()){
            rpath = tempnam("1","1234567890");
            file.open(rpath,std::ios::out);
        } else{
            rpath = path;
            file.open(rpath,std::ios::out);
        }
        if(file.is_open()){
            file<<"#DATA BEGINS\n";
            s.fitSize();

            for (int i = 0; i < s.function.size(); ++i) {
                file<<s.variable.at(i)<<"\t";
                for(auto p:s.function.at(i))
                    file<<p<<"\t";
                file<<"\n";
            }
            file<<"#DATA ENDS\n\n";
        }
        file.close();
        return rpath;
    }


    tempT void CopyVectorToItSelf(std::vector<T> &v, int n){
        v.reserve(v.size()*n);
        std::vector<T> a = v;
        for (int i = 0; i < n; ++i) {
            v.insert( v.end(), a.begin(), a.end());
        }
    }
    tempT void growthVector(std::vector<T>& v,int n){
        v.reserve(n*v.size());
        std::vector<T> a = v;
        v.clear();
        for (int i = 0; i < a.size(); ++i) {
            for (int j = 0; j < n; ++j) {
                v.insert(v.begin()+i*n+j,a.at(i));
            }
        }
    }

    tempT void CartiProductsimple ( std::vector<T> v1, std::vector<T> v2){
        std::vector<std::vector<T>> ans;
        int i = v1.size();
        growthVector(v1,v2.size());
        CopyVectorToItSelf(v2,i);
    }

    tempT void CartiProduct( std::vector<std::vector<T>> &vecs,int n=1){
        if(vecs.size()-1 != n-1){
            int s = vecs.at(0).size();
            for (int i = 0; i < n; ++i) {
                growthVector(vecs.at(i),(int)vecs.at(n).size());
            }
            CopyVectorToItSelf(vecs.at(n),s-1);
            int q = n+1;
            CartiProduct(vecs,q);
        }
    }

    tempT std::vector<std::vector<T>> CartiProduct( std::initializer_list<std::vector<T>> vecs,int n=1){
        std::vector<std::vector<T>> temp;
        for(auto vec: vecs )
            temp.push_back(vec);
        CartiProduct(temp);
        return temp;
    }
    tempT std::vector<T> ConcateVector(const std::vector<std::vector<T>> &ves){
        std::vector<T> ans;
        for(auto i: ves){
            ans.insert(ans.end(),i.begin(),i.end());
        }
        return ans;
    }
}
#undef tempT
#undef tempTK
#endif //MATHUTILS_UTILS_H
