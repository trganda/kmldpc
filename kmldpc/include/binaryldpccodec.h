#ifndef BINARY_LDPC_CODEC_H
#define BINARY_LDPC_CODEC_H

#include <cstdio>
#include <cstring>

#include "utility.h"
#include "randnum.h"

/**********************************************************
类名: Edge
功能: Tanner图上的边
***********************************************************/

typedef struct Edge
{
	int m_row_no;
	int m_col_no;
	
	double m_alpha[2];
	double m_beta[2];
	double m_v2c[2];
	double m_c2v[2];

	struct Edge *left;
	struct Edge *right;
	struct Edge *up;
	struct Edge *down;
} Edge;


/*******************************************************************************************
类名：CBinaryLDPCCodec
功能 ：二进制LDPC码编译码
*******************************************************************************************/

class CBinaryLDPCCodec  
{
public:

//code parameters
	int m_codedim;//码的维数
	int m_codelen;//码的长度
	int m_codechk;//校验位的个数
	double m_coderate;//码率
	int m_encoder_active;//是否编码   0：不编码  1：编码
	char m_file_name_of_H[255];

//parity-check matrix
	int m_num_row;//校验矩阵的行数
	int m_num_col;//校验矩阵的列数
	char **m_decH;//用于译码的校验矩阵
	char **m_encH;//用于编码校验矩阵

	// soft syndrom
	double *m_syndromsoft;

//graph
	Edge *m_row_head;
	Edge *m_col_head;
//
	int *m_cc_hat;
	int m_max_iter;
	int m_success;
	int m_detection_error;//added by liangchulong @2013-01-23-19-04

	void Free();
	void Malloc(int code_no, char *file_name);
	
	void Encoder(int *uu, int *cc);
	void OnlineEncoder(int *uu, int *cc);
	
	int Decoder(double *M2V, int *uu_hat, int iter_count);
	int parityCheck(double *M2V);
	int parityCheck(int* rr);
	void SoftInHardOutDetectorDecoder(double *U2V, double *M2V, int *uu_hat);
	void Kite_Encoder(int *uu, int *cc);
	//////////////////////////////////////////////////////////////////////////
	//added by liang chu long
	int SoftInSoftOut(double *M2V, int *uu_hat, double* cc_hat_soft);
	int SoftInSoftOut(double *U2V, double *M2V, double* uu_soft_out, double* cc_soft_out, char *arrow);
	//////////////////////////////////////////////////////////////////////////
	void InitMsg();//modified by liang chu long
	void PrintCodeParameter(FILE *fp);
private:
	void SystH();
	void SystH2();
	void Free_Tanner_Graph();

};

#endif