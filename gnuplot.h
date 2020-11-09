//
// Created by root on 11/6/20.
//

#ifndef MATHUTILS_GNUPLOT_H
#define MATHUTILS_GNUPLOT_H


#include <string>
#include <vector>
#include <fstream>
#include "Utils.h"

namespace plt {

    void gpc(const std::string &command) {
        system(("gnuplot -p -e \"" + command + " ;pause -1\" ").c_str());
    }

    void gpf(const std::string &filepath) {
        system(("gnuplot -p  \"" + filepath + " \"").c_str());
    }

    template<typename T, size_t N>
    void DataPlot(const std::vector<T>(&data)[N], const std::string &a_command = "", const std::string &b_command = "",
                  const std::string path = "") {
        std::string rpath = Utils::WriteDataToFile(data, path);
        gpc(b_command + " \'" + rpath + "\' " + a_command);
        if (path.empty())
            std::remove(rpath.c_str());
    };

    template<typename T, typename K>
    void DataPlot(typename Calculus::MultiVar::Scalar::Sampler<T, K>::SampledScalarFunction data,
                  const std::string &a_command = "", const std::string &b_command = "", const std::string path = "") {
        std::string rpath = Utils::WriteDataToFile<T, K>(data, path);
        gpc(b_command + " \'" + rpath + "\' " + a_command);
        if (path.empty())
            std::remove(rpath.c_str());
    };

    tempTK
    class TwoDim {
    public:
        static void DataPlot2D(const std::vector<K> X, const std::vector<T> Y, const std::string &a_command = "",
                               const std::string &b_command = "", const std::string path = "") {

            std::string rpath = Utils::WriteDataToFile((std::vector<T>[2]) {X, Y}, path);
            gpc(b_command + " \'" + rpath + "\' " + a_command);
            if (path.empty())
                std::remove(rpath.c_str());
        }

        static void FunctionPlot(std::function<T(K)> rf, K a, K b, int n, const std::string &a_command = "",
                                 const std::string &b_command = "") {
            auto s = Calculus::SingleVar::Sampler<T, K>::FunctionSampler(rf, a, b, n);
            DataPlot2D(s.x, s.y, a_command, "plot " + b_command + " ");
        }
    };


}
#endif //MATHUTILS_GNUPLOT_H