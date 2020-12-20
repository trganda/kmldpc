#include "modem.h"

namespace lab {
    Modem::Modem()
            : num_symbol_(0), input_len_(0), output_len_(0),
              input_symbol_(nullptr), output_symbol_(nullptr), symbol_prob_(nullptr),
              energy_(0.0) {}

    Modem::~Modem() {
        for (int i = 0; i < num_symbol_; i++) {
            delete[] output_symbol_[i];
            delete[] input_symbol_[i];
        }
        delete[] output_symbol_;
        delete[] input_symbol_;
        delete[] symbol_prob_;
    }

    void Modem::Malloc(int code_no, const std::string &modem_file) {
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
        num_symbol_ = 1 << input_len_;
        input_symbol_ = new int *[num_symbol_];
        output_symbol_ = new double *[num_symbol_];
        symbol_prob_ = new double[num_symbol_];

        for (i = 0; i < num_symbol_; i++) {
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

        energy_ = energy_ / num_symbol_;
        for (i = 0; i < num_symbol_; i++) {
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

    void Modem::DeMapping(double *bitLin, double *symRin,
                          double *bitLout, int num_bit_in_blk) const {
        int i, j, q, t;
        int num_point_in_blk;
        double P_0, P_1;
        double sum;

        num_point_in_blk = num_bit_in_blk / input_len_;

        utility::ProbClip(bitLin, num_bit_in_blk);
        utility::ProbClip(symRin, num_symbol_ * num_point_in_blk);

        for (i = 0; i < num_point_in_blk; i++) //���� prob(x=0 | y);
        {
            for (q = 0; q < num_symbol_; q++) {
                symbol_prob_[q] = 1.0;
            }

            //convert bit extrinsic message to symbol message
            t = i * input_len_;
            for (j = 0; j < input_len_; j++) {
                for (q = 0; q < num_symbol_; q++) {
                    if (input_symbol_[q][j] == 0)
                        symbol_prob_[q] *= bitLin[t];
                    else
                        symbol_prob_[q] *= 1.0 - bitLin[t];
                }
                t++;
            }

            //compute symbol full message
            sum = 0.0;
            t = i * num_symbol_;
            for (q = 0; q < num_symbol_; q++) {
                symbol_prob_[q] *= symRin[t];
                sum += symbol_prob_[q];
                t++;
            }

            //normalized
            for (q = 0; q < num_symbol_; q++)
                symbol_prob_[q] /= sum;

            //compute bit extrinsic message for output
            t = i * input_len_;
            for (j = 0; j < input_len_; j++) {
                P_0 = 0;
                P_1 = 0;
                for (q = 0; q < num_symbol_; q++) {
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
        std::vector<std::complex<double>> result(num_symbol_);
        for (int i = 0; i < num_symbol_; i++) {
            result[i] = std::complex<double>(output_symbol_[i][0], output_symbol_[i][1]);
        }
        return result;
    }

    int Modem::GetNumSymbol() const {
        return num_symbol_;
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

}