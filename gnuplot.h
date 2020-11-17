//
// Created by root on 11/6/20.
//

#ifndef MATHUTILS_GNUPLOT_H
#define MATHUTILS_GNUPLOT_H


#include <string>
#include <vector>
#include <fstream>
#include "Utils.h"
namespace plt{

    void gpc(const std::string& command){
        system(("gnuplot -p -e \"" + command +" ;pause -1\" ").c_str());
    }
    void gpf(const std::string& filepath){
        system(("gnuplot -p  \"" + filepath +" \"").c_str());
    }
    template <typename T,size_t N> void DataPlot(const std::vector<T>(&data)[N],const std::string& a_command="",const std::string& b_command="",const int &i = 10,const std::string path=""){
        std::string rpath = Utils::WriteDataToFile(data,path,i);
        gpc(b_command+" \'"+rpath+"\' " + a_command);
        if(path.empty())
            std::remove(rpath.c_str());
    };
    template <typename T,typename K> void DataPlot(typename Calculus::MultiVar::Scalar::Sampler<T,K>::SampledScalarFunction data,const std::string& a_command="",const std::string& b_command="",const int&i=10,const std::string path=""){
        std::string rpath = Utils::WriteDataToFile<T,K>(data,path,i);
        gpc(b_command+" \'"+rpath+"\' " + a_command);
        if(path.empty())
            std::remove(rpath.c_str());
    };

    namespace TwoDim{
        tempTK void DataPlot2D(const std::vector<K> X,const std::vector<T> Y,const std::string& a_command="",const std::string& b_command="",const int&i=10,const std::string path=""){

            std::string rpath = Utils::WriteDataToFile((std::vector<T>[2]){X,Y},path,i);
            gpc(b_command+" \'"+rpath+"\' " + a_command);
            if(path.empty())
                std::remove(rpath.c_str());
        }
        tempTK void FunctionPlot(std::function<T(K)> rf,K a, K b,int n,const std::string& a_command="",const std::string& b_command=""){
            auto s = Calculus::SingleVar::Sampler<T,K>::FunctionSampler(rf,a,b,n);
            DataPlot2D(s.x,s.y,a_command ,"plot "+b_command+" ");
        }
    };

    namespace ThreeDim{

          tempTK void DataPlot3D(typename Calculus::MultiVar::Scalar::Sampler<T,K>::SampledScalarFunction data,const std::string& a_command="",const std::string& b_command="",const int& mode=0,const int&i=10,const std::string path=""){
             switch (mode) {
                 case 0:
                     DataPlot<T,K>(data,"with l " + a_command,"set hidden; set grid ; set contour base; set cntrparam levels auto 10;set pm3d at b;  "+b_command+" splot",i,path);
                     break;
                 case 1:
                     DataPlot<T,K>(data,a_command,"set hidden; set grid ;  "+b_command+" splot",i,path);
                     break;
                 case 2:
                     DataPlot<T,K>(data,a_command,"set hidden; set grid ; set contour base; set cntrparam levels auto 10;"+b_command+" splot",i,path);
                     break;
                 case 4:
                     DataPlot<T,K>(data,"with pm3d "+a_command,"set hidden; set grid ;set pm3d at b;  "+b_command+" splot",i,path);

             }
         }
          tempTK  void FunctionPlot3D(const std::function<T(Algebric::MultiDimPoint<K>)>& rf,const std::initializer_list<const typename Calculus::MultiVar::linspace<K>> &arg,const std::string& a_command="",const std::string& b_command="",const int& mode=0,const int&i=10,const std::string &path=""){
            auto samples = Calculus::MultiVar::Scalar::Sampler<T,K>::FunctionSampler(rf,arg);
             switch (mode) {
                 case 0:
                     DataPlot<T,K>(samples,"with l " + a_command,"set hidden; set grid ; set contour base; set cntrparam levels auto 10;set pm3d at b;  "+b_command+"; splot",i,path);
                     break;
                 case 1:
                     DataPlot<T,K>(samples,a_command,"set hidden; set grid ;  "+b_command+"; splot",i,path);
                     break;
                 case 2:
                     DataPlot<T,K>(samples,a_command,"set hidden; set grid ; set contour base; set cntrparam levels auto 10;"+b_command+"; splot",i,path);
                     break;
                 case 4:
                     DataPlot<T,K>(samples,"with pm3d "+a_command,"set hidden; set grid ;set pm3d at b;  "+b_command+" splot",i,path);

             }

         }
     };



}
#endif //MATHUTILS_GNUPLOT_H