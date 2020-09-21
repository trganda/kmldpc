#ifndef KMLDPC_BINARY5GLDPCCODEC_H
#define KMLDPC_BINARY5GLDPCCODEC_H

/**********************************************************
类名: EdgeLDPC
功能: Tanner图上的边
***********************************************************/
#include <cstdio>
#include <cstring>

typedef struct EdgeLDPC
{
    int m_row_no;
    int m_col_no;

    double m_alpha[2];
    double m_beta[2];
    double m_v2c[2];
    double m_c2v[2];

    struct EdgeLDPC *left;
    struct EdgeLDPC *right;
    struct EdgeLDPC *up;
    struct EdgeLDPC *down;
} EdgeLDPC;


/*******************************************************************************************
类名：CBinaryLDPCCodec
功能 ：二进制LDPC码编译码
*******************************************************************************************/

class CBinary5GLDPCCodec
{
public:
//code parameters
    int m_codedim;//码的维数
    int m_codelen_no_puncture;//码的总长度
    int m_codelen_puncture;  //码的发送长度，打孔之后
    int m_codechk;//校验位的个数
    double m_coderate;//码率
    int m_encoder_active;//是否编码   0：不编码  1：编码
    char m_file_name_of_H[255];

    int m_lifting_factor;

    int *m_cc_nopuncture;
    double *m_cc_soft_nopuncture;

//parity-check matrix
    int m_num_row;//校验矩阵的行数
    int m_num_col;//校验矩阵的列数
    char **m_decH;//用于译码的校验矩阵
    char **m_encH;//用于编码校验矩阵

    // soft syndrom
    double *m_syndromsoft;

//graph
    EdgeLDPC *m_row_head;
    EdgeLDPC *m_col_head;
//
    int *m_cc_hat;
    int m_max_iter;
    int m_success;
    int m_detection_error;//added by liangchulong @2013-01-23-19-04

    void Free();
    void Malloc(int code_no, char *file_name);

    void Encoder(int *uu, int *cc);

    void Encoder_5G(int *uu, int *cc);

    int Decoder(double *M2V, int *uu_hat);

    int Decoder_5G(double *M2V, int *uu_hat, int iter_cnt);

    int SoftInSoftOut(double *M2V, int *uu_hat, double* cc_hat_soft);

    int Parity_Checking();

    int parityCheck(int *cc);

    void InitMsg();//modified by liang chu long
    void PrintCodeParameter(FILE *fp);
private:
    void SystH();
    void SystH_1();
    void SystH2();
    void SystH_5G();
    void Free_Tanner_Graph();

};

#endif //KMLDPC_BINARY5GLDPCCODEC_H
