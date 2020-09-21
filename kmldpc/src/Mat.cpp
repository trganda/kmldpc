#include "Mat.h"
#include "Log.h"

namespace kmldpc
{
    Mat::Mat(std::string &filename) :
        _filename(filename), _matfp(nullptr)
    {}

    Mat::~Mat() {
        if (nullptr != _matfp) {
            Mat_Close(_matfp);
        }
    }

    void Mat::open() {
        _matfp = Mat_CreateVer(_filename.c_str(), nullptr,
                MAT_FT_DEFAULT);
        if (nullptr == _matfp) {
            kmldpc::Log::get().setLevel(Error);
            LOG(Error, true) << "Creating file failed, file name is "
                                    << _filename << std::endl;
            exit(-1);
        }
    }

    void Mat::close() {
        if (nullptr != _matfp) {
            Mat_Close(_matfp);
            _matfp = nullptr;
        }
    }

    void Mat::writeInt(const std::string& varname, int32_t data) {
        if (nullptr == _matfp) {
            kmldpc::Log::get().setLevel(Error);
            LOG(Error, true) << "Create of open the file first" << std::endl;
            exit(-1);
        }
        size_t dims[2] = {1, 1};
        int32_t datas[1] = {data};

        matvar_t *matvar = Mat_VarCreate(varname.c_str(), MAT_C_INT32, MAT_T_INT32, 2, dims, datas, 0);
        if (nullptr == _matfp) {
            Log::get().setLevel(Error);
            LOG(Error, true) << "Error creating variable for " << varname << std::endl;
            exit(-1);
        } else {
            Mat_VarWrite(_matfp, matvar, MAT_COMPRESSION_NONE);
            Mat_VarFree(matvar);
        }

        Log::get().setLevel(Info);
        LOG(Info, false) << "Writed " << varname << " to " << _filename << std::endl;
    }

    void Mat::writeDouble(const std::string &varname, double data) {
        if (nullptr == _matfp) {
            Log::get().setLevel(Error);
            LOG(Error, true) << "Create of open the file first" << std::endl;
            exit(-1);
        }
        size_t dims[2] = {1, 1};
        double datas[1] = {data};

        matvar_t *matvar = Mat_VarCreate(
                varname.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, datas, 0);
        if (nullptr == matvar) {
            Log::get().setLevel(Error);
            LOG(Error, true) << "Error creating variable for " << varname << std::endl;
            exit(-1);
        } else {
            Mat_VarWrite(_matfp, matvar, MAT_COMPRESSION_NONE);
            Mat_VarFree(matvar);
        }

        Log::get().setLevel(Info);
        LOG(Info, false) << "Writed " << varname << " to " << _filename << std::endl;
    }

    void Mat::writeComplex(const std::string& varname, std::complex<double> data) {
        if (nullptr == _matfp) {
            Log::get().setLevel(Error);
            LOG(Error, true) << "Create of open the file first" << std::endl;
            exit(-1);
        }
        size_t dims[2] = {1, 1};
        double x[1] = {data.real()},
               y[1] = {data.imag()};

        struct mat_complex_split_t datas = {x, y};

        matvar_t *matvar = Mat_VarCreate(
                varname.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &datas, MAT_F_COMPLEX);
        if (nullptr == matvar) {
            Log::get().setLevel(Error);
            LOG(Error, true) << "Error creating variable for " << varname << std::endl;
            exit(-1);
        } else {
            Mat_VarWrite(_matfp, matvar, MAT_COMPRESSION_NONE);
            Mat_VarFree(matvar);
        }

        Log::get().setLevel(Info);
        LOG(Info, false) << "Writed " << varname << " to " << _filename << std::endl;
    }

    void Mat::writeVector(const std::string &varname, const std::vector<int32_t> &data) {
        if (nullptr == _matfp) {
            Log::get().setLevel(Error);
            LOG(Error, true) << "Create of open the file first" << std::endl;
            exit(-1);
        }
        size_t dims[2] = {data.size(), 1};
        int32_t datas[data.size()] = {0};
        for (int i=0; i<data.size(); i++) {
            datas[i] = data[i];
        }

        matvar_t *matvar = Mat_VarCreate(varname.c_str(), MAT_C_INT32, MAT_T_INT32, 2, dims, datas, 0);
        if (nullptr == matvar) {
            Log::get().setLevel(Error);
            LOG(Error, true) << "Error creating variable for " << varname << std::endl;
            exit(-1);
        } else {
            Mat_VarWrite(_matfp, matvar, MAT_COMPRESSION_NONE);
            Mat_VarFree(matvar);
        }

        Log::get().setLevel(Info);
        LOG(Info, false) << "Writed " << varname << " to " << _filename << std::endl;
    }

    void Mat::writeVector(const std::string &varname, const std::vector<double> &data) {
        if (nullptr == _matfp) {
            Log::get().setLevel(Error);
            LOG(Error, true) << "Create of open the file first" << std::endl;
            exit(-1);
        }
        size_t dims[2] = {data.size(), 1};
        double datas[data.size()] = {0};
        for (int i=0; i<data.size(); i++) {
            datas[i] = data[i];
        }

        matvar_t *matvar = Mat_VarCreate(
                varname.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, datas, 0);
        if (nullptr == matvar) {
            Log::get().setLevel(Error);
            LOG(Error, true) << "Error creating variable for " << varname << std::endl;
            exit(-1);
        } else {
            Mat_VarWrite(_matfp, matvar, MAT_COMPRESSION_NONE);
            Mat_VarFree(matvar);
        }

        Log::get().setLevel(Info);
        LOG(Info, false) << "Writed " << varname << " to " << _filename << std::endl;
    }

    void Mat::writeVector(const std::string &varname, const std::vector<std::complex<double> > &data) {
        if (nullptr == _matfp) {
            kmldpc::Log::get().setLevel(Error);
            LOG(Error, true) << "Create of open the file first" << std::endl;
            exit(-1);
        }
        size_t dims[2] = {data.size(), 1};
        double x[data.size()],
               y[data.size()];
        for (int i=0; i<data.size(); i++) {
            x[i] = data[i].real();
            y[i] = data[i].imag();
        }

        struct mat_complex_split_t datas = {x, y};
        matvar_t* matvar = Mat_VarCreate(varname.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &datas, MAT_F_COMPLEX);
        if (nullptr == matvar) {
            Log::get().setLevel(Error);
            LOG(Error, true) << "Error creating variable for " << varname << std::endl;
            exit(-1);
        } else {
            Mat_VarWrite(_matfp, matvar, MAT_COMPRESSION_NONE);
            Mat_VarFree(matvar);
        }

        Log::get().setLevel(Info);
        LOG(Info, false) << "Writed " << varname << " to " << _filename << std::endl;
    }
}