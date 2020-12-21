#include "modem.h"

namespace lab {
    Modem::Modem(const toml::value& arguments)
            : symbol_num_(0), input_len_(0), output_len_(0),
              input_symbol_(nullptr), output_symbol_(nullptr), symbol_prob_(nullptr),
              energy_(0.0) {
        const auto modem = toml::find(arguments, "modem");
        std::string modem_file = toml::find<std::string>(modem, "modem_file");
        this->init(modem_file);
    }

    Modem::~Modem() {
        for (int i = 0; i < symbol_num_; i++) {
            delete[] output_symbol_[i];
            delete[] input_symbol_[i];
        }
        delete[] output_symbol_;
        delete[] input_symbol_;
        delete[] symbol_prob_;
    }

    void Modem::Malloc(const std::string &modem_file) {
        int i, j;
        char temp_str[255] = {' '};
        int sym;
        int temp;
        double energy;
        FILE *fp;

        if ((fp = fopen(modem_file.c_str(), "r")) == nullptr) {
            LOG(lab::logger::Error, true) << "Cannot Open " << modem_file << std::endl;
            exit(-1);
        }

        fscanf(fp, "%s", temp_str);
        fscanf(fp, "%d", &input_len_);

        fscanf(fp, "%s", temp_str);
        fscanf(fp, "%d", &output_len_);

        fscanf(fp, "%s", temp_str);
        energy_ = 0.0;
        symbol_num_ = 1 << input_len_;
        input_symbol_ = new int *[symbol_num_];
        output_symbol_ = new double *[symbol_num_];
        symbol_prob_ = new double[symbol_num_];

        for (i = 0; i < symbol_num_; i++) {
            fscanf(fp, "%d", &sym);

            input_symbol_[i] = new int[input_len_];
            temp = 0;
            for (j = 0; j < input_len_; j++) {
                fscanf(fp, "%d", &input_symbol_[i][j]);
                //temp = (temp<<1) | input_symbol_[i][j];
                temp = (temp << 1) + input_symbol_[i][j];
            }
            if (sym != temp || sym != i) {
                fprintf(stderr, "\nSym = %d is not the binary expression of %d!\n", sym, temp);
                fprintf(stderr, "\nSome symbols are missed or not in the right order!\n");
                system("pause");
                exit(3);
            }
            output_symbol_[i] = new double[output_len_];
            for (j = 0; j < output_len_; j++) {
                fscanf(fp, "%lf", &output_symbol_[i][j]);
            }

            for (j = 0; j < output_len_; j++) {
                energy = output_symbol_[i][j] * output_symbol_[i][j];
                energy_ += energy;
            }
        }
        fclose(fp);

        energy_ = energy_ / symbol_num_;
        for (i = 0; i < symbol_num_; i++) {
            for (j = 0; j < output_len_; j++) {
                output_symbol_[i][j] /= sqrt(energy_);
            }
        }
    }

    void Modem::Mapping(const int *bin_cc, double *xx, int num_bit_in_blk) const {
        int i, j, r, t;
        int num_point_in_blk;
        int symbol;

        num_point_in_blk = num_bit_in_blk / input_len_;

        r = 0;
        t = 0;
        for (i = 0; i < num_point_in_blk; i++) {
            symbol = 0;
            for (j = 0; j < input_len_; j++) {
                //symbol = (symbol << 1) | bin_cc[r];
                symbol = (symbol << 1) + bin_cc[r];
                r++;
            }
            for (j = 0; j < output_len_; j++) {
                xx[t] = output_symbol_[symbol][j];
                t++;
            }
        }
    }

    std::vector<std::complex<double>> Modem::Mapping(const int* bin_cc, int bin_cc_len) {
        int blk_symbol_num = bin_cc_len / input_len_;
        std::vector<std::complex<double>> ret(blk_symbol_num);
        for (int i=0; i<blk_symbol_num; i++) {
            int symbol_idx = 0;
            for (int j=0; j<input_len_; j++) {
                symbol_idx = (symbol_idx << 1) + bin_cc[j + i*input_len_];
            }
            ret[i] = symbol_out_[symbol_idx];
        }

        return ret;
    }

    void Modem::DeMapping(double *bitLin, double *symRin,
                          double *bitLout, int num_bit_in_blk) const {
        int i, j, q, t;
        int num_point_in_blk;
        double P_0, P_1;
        double sum;

        num_point_in_blk = num_bit_in_blk / input_len_;

        utility::ProbClip(bitLin, num_bit_in_blk);
        utility::ProbClip(symRin, symbol_num_ * num_point_in_blk);

        for (i = 0; i < num_point_in_blk; i++) //���� prob(x=0 | y);
        {
            for (q = 0; q < symbol_num_; q++) {
                symbol_prob_[q] = 1.0;
            }

            //convert bit extrinsic message to symbol message
            t = i * input_len_;
            for (j = 0; j < input_len_; j++) {
                for (q = 0; q < symbol_num_; q++) {
                    if (input_symbol_[q][j] == 0)
                        symbol_prob_[q] *= bitLin[t];
                    else
                        symbol_prob_[q] *= 1.0 - bitLin[t];
                }
                t++;
            }

            //compute symbol full message
            sum = 0.0;
            t = i * symbol_num_;
            for (q = 0; q < symbol_num_; q++) {
                symbol_prob_[q] *= symRin[t];
                sum += symbol_prob_[q];
                t++;
            }

            //normalized
            for (q = 0; q < symbol_num_; q++)
                symbol_prob_[q] /= sum;

            //compute bit extrinsic message for output
            t = i * input_len_;
            for (j = 0; j < input_len_; j++) {
                P_0 = 0;
                P_1 = 0;
                for (q = 0; q < symbol_num_; q++) {
                    if (input_symbol_[q][j] == 0)
                        P_0 += symbol_prob_[q];
                    else
                        P_1 += symbol_prob_[q];
                }
                P_0 /= bitLin[t];
                P_1 /= (1.0 - bitLin[t]);
                bitLout[t] = P_0 / (P_0 + P_1);
                t++;
            }
        } //end of for (i = 0; i < num_point_in_blk; i++)

        utility::ProbClip(bitLout, num_bit_in_blk);
    }

    std::vector<std::complex<double>> Modem::GetConstellations() const {
//        std::vector<std::complex<double>> result(symbol_num_);
//        for (int i = 0; i < symbol_num_; i++) {
//            result[i] = std::complex<double>(output_symbol_[i][0], output_symbol_[i][1]);
//        }
//        return result;
        return symbol_out_;
    }

    int Modem::GetNumSymbol() const {
        return symbol_num_;
    }

    int Modem::GetInputLen() const {
        return input_len_;
    }

    int Modem::GetOutputLen() const {
        return output_len_;
    }

    double **Modem::GetOutputSymbol() const {
        return output_symbol_;
    }

    void Modem::init(const std::string& modem_file) {
        std::ifstream ifs(modem_file, std::ios_base::binary);
        if (!ifs.is_open()) {
            LOG(lab::logger::Error, true) << "Cannot Open " << modem_file << std::endl;
            exit(-1);
        }

        std::string temp;
        ifs >> temp;
        ifs >> input_len_;
        ifs >> temp;
        ifs >> output_len_;
        ifs >> temp;

        symbol_num_ = 1<<input_len_;
        symbol_in_ = std::vector<std::vector<int>>(symbol_num_,
                                                   std::vector<int>(input_len_));
        symbol_out_ = std::vector<std::complex<double>>(symbol_num_);

        int sym_dec = 0;
        double energies = 0;
        for (int i=0; i<symbol_num_; i++) {
            ifs >> sym_dec;
            int temp_dec = 0;
            for (int j=0; j<input_len_; j++) {
                ifs >> symbol_in_[i][j];
                temp_dec = (temp_dec << 1) + symbol_in_[i][j];
            }
            if (sym_dec != temp_dec || sym_dec != i) {
                LOG(lab::logger::Error, true) << sym_dec << " is not the binary expression of "
                                              << temp_dec << std::endl;
                exit(-1);
            }
            double real, imag;
            ifs >> real >> imag;
            symbol_out_[i] = std::complex<double>(real, imag);
            energies += pow(abs(symbol_out_[i]), 2);
        }
        ifs.close();

        energies /= symbol_num_;
        for (int i=0; i<symbol_num_; i++) {
            symbol_out_[i] /= sqrt(energies);
        }
    }
}