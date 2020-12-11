#ifndef LAB_RAND_NUM_HPP
#define LAB_RAND_NUM_HPP

#include <cstdio>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>

namespace lab {

class RandNum {
    public:
        explicit RandNum() = default;

        virtual ~RandNum() = default;

        virtual void SetSeed(int flag) = 0;

        virtual void PrintState(FILE *fp) = 0;

        virtual double Uniform() = 0;

        virtual void Normal(double *, int) = 0;
};

/*
 * The following generator employs the linear-congruential method,
 * and specifically uses a choice of multiplier that was proposed
 * as a standard by Stephen K. Park et al. in "Technical correspondence,"
 * Communications of the ACM36(1993), number 7, 108-110
 */

class CLCRandNum : public RandNum {
    public:
        explicit CLCRandNum()
                : state(0), A(48271),
                  M(2147483647), Q(M / A), R(M % A) {}
        ~CLCRandNum() override = default;

        void SetSeed(int flag) final {
            if (flag < 0)
                state = 17;
            else if (flag == 0) {
                state = 0;
                while (state == 0) {
                    srand((unsigned) time(nullptr));
                    state = rand();
                }
            } else {
                fprintf(stdout, "\nEnter the initial state: ");
                fscanf(stdin, "%ld", &state);
            }
        }

        void PrintState(FILE *fp) final {
            fprintf(fp, "\n***init_state = %ld***\n", state);
        }

        double Uniform() final {
            double u;

            int tmpState = A * (state % Q) - R * (state / Q);
            if (tmpState >= 0)
                state = tmpState;
            else
                state = tmpState + M;

            u = state / (double) M;

            return u;
        }

        void Normal(double *nn, int len_nn) final {
            double x1, x2, w;
            int t;

            for (t = 0; 2 * t + 1 < len_nn; t++) {
                w = 2.0;
                while (w > 1.0) {
                    x1 = 2.0 * Uniform() - 1.0;
                    x2 = 2.0 * Uniform() - 1.0;

                    w = x1 * x1 + x2 * x2;
                }

                w = sqrt(-2.0 * log(w) / w);

                nn[2 * t] = x1 * w;
                nn[2 * t + 1] = x2 * w;
            }

            if (len_nn % 2 == 1) {
                w = 2.0;
                while (w > 1.0) {
                    x1 = 2.0 * Uniform() - 1.0;
                    x2 = 2.0 * Uniform() - 1.0;

                    w = x1 * x1 + x2 * x2;
                }

                w = sqrt(-2.0 * log(w) / w);

                nn[len_nn - 1] = x1 * w;
            }
        };
        // Singleton pattern
        static CLCRandNum& Get() {
            static CLCRandNum instance;
            return instance;
        }

    private:
        long int state;
        int A;
        long M;
        int Q;
        int R;
};
/*
 *The following generator employs the Wichman-Hill algorithm
 */
class CWHRandNum : public RandNum {
    public:
        explicit CWHRandNum()
                : X(0), Y(0), Z(0) {}
        ~CWHRandNum() override = default;

        void SetSeed(int flag) final {
            if (flag < 0)
            {
                X = 13;
                Y = 37;
                Z = 91;
            }
            else if (flag == 0)
            {
                X = 0;
                Y = 0;
                Z = 0;
                while (X == 0 || Y == 0 || Z == 0)
                {
                    srand((unsigned)time(nullptr));
                    X = rand();
                    Y = rand();
                    Z = rand();
                }
            }
            else
            {
                fprintf(stdout, "\nEnter the initial state (X Y Z): ");
                fscanf(stdin, "%d %d %d", &X, &Y, &Z);
            }
        }

        void PrintState(FILE *fp) final {
            fprintf(fp, "\n***init_state (X Y Z) = %d %d %d***\n", X, Y, Z);
        }

        double Uniform() final {
            double U;

            X = 171 * X % 30269;
            Y = 172 * Y % 30307;
            Z = 170 * Z % 30323;

            U = X / 30269.0 + Y / 30307.0 + Z / 30323.0;
            U = U - int(U);

            return U;
        }

        void Normal(double *nn, int len_nn) final {
            double x1, x2, w;
            int t;

            for (t = 0; 2 * t + 1 < len_nn; t++)
            {
                w = 2.0;
                while (w > 1.0)
                {
                    x1 = 2.0 * Uniform() - 1.0;
                    x2 = 2.0 * Uniform() - 1.0;

                    w = x1 * x1 + x2 * x2;
                }

                w = sqrt(-2.0 * log(w) / w);

                nn[2 * t] = x1 * w;
                nn[2 * t + 1] = x2 * w;
            }

            if (len_nn % 2 == 1)
            {
                w = 2.0;
                while (w > 1.0)
                {
                    x1 = 2.0 * Uniform() - 1.0;
                    x2 = 2.0 * Uniform() - 1.0;

                    w = x1 * x1 + x2 * x2;
                }

                w = sqrt(-2.0 * log(w) / w);

                nn[len_nn - 1] = x1 * w;
            }
        }

        static CWHRandNum& Get() {
            static CWHRandNum instance;
            return instance;
        }

    private:
        int X;
        int Y;
        int Z;
};

}

#endif // RAND_NUM_H
