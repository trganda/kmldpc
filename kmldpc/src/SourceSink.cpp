#include "SourceSink.h"
#include "RandNum.h"

extern CLCRandNum rndGen0;

void CSourceSink::GetBitStr(int *uu, int len)
{
	int t; 

	for (t = 0; t < len; t++)
		uu[t] = (rndGen0.Uniform() < 0.5?0:1);
	    //uu[t] = 0;

}

void CSourceSink::GetSymStr(int *uu, int qary, int len)
{
	int t; 

	for (t = 0; t < len; t++){
		uu[t] = qary;
		while (uu[t] == qary)
			uu[t] = (int)(qary * rndGen0.Uniform());
//	    uu[t] = 0;
	}
}

void CSourceSink::ClrCnt()
{
	m_num_tot_blk = 0;
	m_num_tot_bit = 0;
	m_num_err_blk = 0;
	m_num_err_bit = 0;
}

void CSourceSink::CntErr(int *uu, int *uu_hat, int len, int accumulator)
{
	int t;

	m_temp_err = 0;
	for (t = 0; t < len; t++){
		if (uu_hat[t] != uu[t])
			m_temp_err++;
	}

	if (accumulator == 1){
		if (m_temp_err > 0){
			m_num_err_bit += m_temp_err;
			m_num_err_blk += 1;
		}
		
		m_num_tot_blk += 1.0;
		m_num_tot_bit += len;

		m_ber = m_num_err_bit / m_num_tot_bit;
		m_fer = m_num_err_blk / m_num_tot_blk;
	}
}

void CSourceSink::PrintResult(double snr) const
{
	LOG(kmldpc::Info, true) << std::fixed << std::setprecision(0) << std::setfill('0')
	                           << "SNR = "
	                           << std::setw(3) << std::right << snr << ' '
	                           << "Total blk = "
	                           << std::setw(7) << std::right << m_num_tot_blk << ' '
	                           << "Error blk = "
	                           << std::setw(7) << std::right << m_num_err_blk << ' '
                               << "Error bit = "
                               << std::setw(7) << std::right << m_num_err_bit << ' '
                               << std::fixed << std::setprecision(14)
                               << "BER = " << m_ber << ' '
                               << "FER = " << m_fer
                               << std::endl;
}

void CSourceSink::PrintResult(FILE *fp) const
{
    fprintf(fp, "tot_blk = %d: err_blk = %d: err_bit = %d: ber = %12.10lf: fer = %12.10lf\n",
            (int)m_num_tot_blk, m_num_err_blk, m_num_err_bit, m_ber, m_fer);
}