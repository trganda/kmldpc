#ifndef SOURCE_SINK_H
#define SOURCE_SINK_H

#include <complex>
#include <vector>

class CSourceSink  
{
public:
	CSourceSink();
	virtual ~CSourceSink();

	double m_num_tot_blk{};
	double m_num_tot_bit{};
	int m_num_err_blk{};
	int m_num_err_bit{};
	int m_temp_err{};

	double m_ber{};
	double m_fer{};
	
	void GetBitStr(int *uu, int len);
	void GetSymStr(int *uu, int qary, int len);
	void ClrCnt();
	void CntErr(int *uu, int *uu_hat, int len, int accumulator);
	void PrintResult(FILE *fp) const;
};

#endif