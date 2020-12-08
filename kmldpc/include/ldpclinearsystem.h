#ifndef LDPC_LINEAR_SYSTEM_H
#define LDPC_LINEAR_SYSTEM_H

#include "sourcesink.h"
#include "xorsegcodec.h"
#include "modemlinearsystem.h"
#include "kmeans.h"
#include "randnum.h"
#include <complex>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>


class LDPC_Linear_System  
{
public:
	LDPC_Linear_System() {}
	virtual ~LDPC_Linear_System(){}

	//read Setup_of_LDPC_Linear_System*.txt
	double m_min_snr;
	double m_max_snr;
	double m_inc_snr;
	int m_max_blk_err;
	int m_max_blk_num;

	int m_total_angle;
		
	int * m_uu;
	int * m_uu_hat;
	int   m_len_uu;
	
	int * m_cc;
	int * m_cc_hat;
	int   m_len_cc;

	double * m_sym_prob;
	
	struct tm * ptr;
	time_t loctime;

	CSourceSink m_source_sink;
	CSourceSink m_source_extrabits;

	XORSegCodec m_codec;
	Modem_Linear_System Modem_Lin_Sym;
				
	void StartSimulator();
	void EndSimulator();
	void Simulator();
};

#endif