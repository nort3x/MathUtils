//
// Created by root on 11/6/20.
//

#ifndef MATHUTILS_GNUPLOT_H
#define MATHUTILS_GNUPLOT_H


#include <string>
#include <vector>
#include <fstream>
namespace plt{

    void gpc(const std::string& command){
        system(("gnuplot -p -e \"" + command +" \"").c_str());
    }
    void gpf(const std::string& filepath){
        system(("gnuplot -p  \"" + filepath +" \"").c_str());
    }
    tempT void DataPlot(const std::string& b_command, std::vector<T> *arr, int cols,const std::string& a_command){
        std::string rpath = tempnam("1","1234567890");
        std::fstream file;
        file.open(rpath,std::ios::out);
        if(file.is_open()){
            file<<"#DATA BEGINS\n";
            int min_n_data = arr[0].size();
            for(int i=0;i<cols;i++)
                min_n_data = std::min(min_n_data,(int)arr[i].size());

            for(int i=0;i<min_n_data;i++) {
                for (int j = 0; j < cols; ++j) {
                    file << arr[j][i] << "\t";
                }
                file << "\n";
            }
            file<<"#DATA ENDS\n\n";
        }
        file.close();
        gpc(b_command+" \'"+rpath+"\' " + a_command);
        std::remove(rpath.c_str());
    };

    tempTK class TwoDim{
    public:
        static void DataPlot2D(const std::string& b_command,const std::vector<K> X,const std::vector<T> Y,const std::string& a_command){
            std::string rpath = tempnam("1","1234567890");
            std::fstream file;
            file.open(rpath,std::ios::out);
            if(file.is_open()){
                file<<"#DATA BEGINS\n";
                int min_n_data = std::min((int)X.size(),(int)Y.size());
                for(int i=0;i<min_n_data;i++) {
                    file<<X.at(i)<<"\t"<<Y.at(i)<<"\n";
                }
                file<<"#DATA ENDS\n\n";
            }
            file.close();
            gpc(b_command+" \'"+rpath+"\' " + a_command);
        }
        static void FunctionPlot(const std::string& pre_q, std::function<T(K)> rf,K a, K b,int n,const std::string& pas_q){
            auto s = Calculus::SingleVar::Sampler<T,K>::FunctionSampler(rf,a,b,n);
            DataPlot2D("plot "+pre_q ,s.x,s.y,pas_q);
        }
    };
}
#endif //MATHUTILS_GNUPLOT_H