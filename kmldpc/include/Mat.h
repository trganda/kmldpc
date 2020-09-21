#ifndef KMLDPC_MAT_H
#define KMLDPC_MAT_H

#include <matio.h>
#include <string>
#include <vector>
#include <complex>

namespace kmldpc
{
    // Encapsulation of matio https://github.com/tbeu/matio
    class Mat
    {
    public:
        Mat(std::string& filename);
        ~Mat();

    public:
        void open();
        void close();
        void writeInt(const std::string& varname, int32_t data);
        void writeDouble(const std::string& varname, double data);
        void writeComplex(const std::string& varname, std::complex<double> data);
        void writeVector(const std::string& varname, const std::vector<int32_t>& data);
        void writeVector(const std::string& varname, const std::vector<double>& data);
        void writeVector(const std::string& varname, const std::vector<std::complex<double>>& data);

    private:
        std::string _filename;
        mat_t* _matfp;
    };
}

#endif //KMLDPC_MAT_H
