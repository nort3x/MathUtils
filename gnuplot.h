//
// Created by root on 11/6/20.
//

#ifndef MATHUTILS_GNUPLOT_H
#define MATHUTILS_GNUPLOT_H


#include <string>
#include <vector>
#include <fstream>
#include <cstdio>
#include "Utils.h"
#define tempT template<typename T>
#define tempTK template<typename T,typename K>

namespace plt{
    void gpc(const std::string& command) {
        system(("gnuplot  -e \"" + command +" ;pause -1\" ").c_str());
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
        tempTK class Gnuplot2D{
        private:
            std::string command;
            std::vector<std::string> paths;
            std::string rpath;
        public:
            explicit Gnuplot2D(const std::string& b_command=""){
                if(b_command.empty())
                    command ="plot ";
                else
                    command = b_command + "; plot";
                rpath = "";
            }
            explicit Gnuplot2D(const std::string &path,const std::string& b_command=""){
                if(b_command.empty())
                    command ="plot ";
                else
                    command = b_command + "; plot";
                rpath = path;
            }
            void plot(){
                command = command.substr(0,command.rfind(','));
                gpc(command);
            }
            void clear(){
                command = "";
                plot();
            }
            void dispose(){
                for(const std::string &p: paths)
                    std::remove(p.c_str());
            }
            void addFunctionPlot(std::function<T(K)> rf,K a, K b,int n,const std::string& a_command="",const std::string& b_command="",const std::string &name=""){
                auto s = Calculus::SingleVar::Sampler<T,K>::FunctionSampler(rf,a,b,n);
                addDataPlot(s.x,s.y,a_command ," "+b_command+" ",10,name);
            }
            void addDataPlot(const std::vector<K> X,const std::vector<T> Y,const std::string& a_command="",const std::string& b_command="",const int&i=10,std::string name=""){
                if(name.empty())
                    name =tempnam(rpath.c_str(),"12345");
                else
                    name = rpath+"/"+name;
                std::string wpath = Utils::WriteDataToFile((std::vector<T>[2]){X,Y},name,i);
                paths.push_back(wpath);
                command += b_command+" \'"+wpath+"\' " + a_command+",";
            }
            void addPlotCommand(const std::string& s){
                command += s+" ,";
            }
        };
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
        tempTK void CurvePlot3D_anim(const std::function<Algebric::MultiDimPoint<T>(K)>& curve_func,K a,K b,int n,const std::string& anim_gif_dir,const std::string& a_command="",const std::string& b_command="",const int& mode=0,const int&i=10,const std::string &path=""){
            auto c =Calculus::MultiVar::VectorValued::Sampler<T,K>::CurveSampler(curve_func,a,b,n);
            std::string rpath = Utils::WriteDataToFile<T,K>(c,path,i);
            system(("mkdir "+anim_gif_dir).c_str());
            std::string script = b_command+";set xr ["+std::to_string(c.ranges.at(0).min)+':'+std::to_string(c.ranges.at(0).max)+"]";
            if(c.ranges.size()>1)
                script += ";set yr ["+std::to_string(c.ranges.at(1).min)+':'+std::to_string(c.ranges.at(1).max)+"]";
            if(c.ranges.size()>2)
                script += ";set zr ["+std::to_string(c.ranges.at(2).min)+':'+std::to_string(c.ranges.at(2).max)+"]";

            script+=                    ";set term png ;do for[ii=1:"+std::to_string(n)+":1]{"+
                                           ";fname = sprintf(\'"+anim_gif_dir+"/"+"img%i.png\',ii)"
                                           ";set output fname"
                                           ";splot "+" \'"+rpath+"\' " + a_command + " every ::1::ii"+
                                           ";}";
            gpc(script);
            if(path.empty())
                std::remove(rpath.c_str());
        }
        tempTK void CurvePlot3D(const std::function<Algebric::MultiDimPoint<T>(K)>& curve_func,K a,K b,int n,const std::string& a_command="",const std::string& b_command="",const int& mode=0,const int&i=10,const std::string &path=""){
            auto c =Calculus::MultiVar::VectorValued::Sampler<T,K>::CurveSampler(curve_func,a,b,n);
            std::string rpath = Utils::WriteDataToFile<T,K>(c,path,i);
            std::string script = b_command;
            script+=(";splot  \'"+rpath+"\' " + a_command );
            gpc(script);
            if(path.empty())
                std::remove(rpath.c_str());
        }

    };



}
#undef tempT
#undef tempTK
#endif //MATHUTILS_GNUPLOT_H