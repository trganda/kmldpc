#include "linearsystem.h"

extern CLCRandNum rndGen0;
extern CWHRandNum rndGen1;

void Linear_System::Malloc(CModem * modem, int len_xx, int code_no, char * file_name)
{
	char temp_str[80] = {' '};
	FILE *fp;
	m_modem = modem;
	m_len_xx = len_xx;

	if ((fp = fopen(file_name, "r")) == nullptr){
		fprintf(stderr, "\nCannot open %s", file_name);
		exit(3);
	}
		
	noise = new double[m_len_xx];
	m_yy = new double[2];
	m_Nyy = new double[m_len_xx];
	
	m_symbol_prob = new double[m_modem->num_symbol];
	
	num_symbol = m_modem->num_symbol;
}

void Linear_System::Free()
{
	delete []noise;
	delete []m_yy;
	delete []m_Nyy;
	delete []m_symbol_prob;
}

void Linear_System::AWGN_linear_system(double *xx, double *sym_prob)
{
	int i;
	int xxidx, nnidx;
	int temp;
	
	int num_of_symbol_blk = m_len_xx / 2;
	
	rndGen0.Normal(noise, m_len_xx);

	xxidx = 0;
	nnidx = 0;
	for (i = 0; i < num_of_symbol_blk; i++)
	{
		m_yy[0] = xx[xxidx] + (sigma/SQRT2) * noise[nnidx];
		m_yy[1] = xx[xxidx+1] + (sigma/SQRT2) * noise[nnidx+1];
		m_Nyy[xxidx] = m_yy[0];
		m_Nyy[xxidx + 1] = m_yy[1];
		xxidx += 2;
		nnidx += 2;

		temp = i * m_modem->num_symbol;		
		Soft_AWGN_Demodulation(m_yy, (sym_prob + temp));
	} //end of for (i = 0; i < num_sym_in_blk; i++)
}

void Linear_System::Parition_HAWGN_system(double* xx, std::vector<std::complex<double>>& selectH)
{
	int i;
	int xxidx, nnidx;

	int num_of_symbol_blk = m_len_xx / 2;
	int partLen = selectH.size();

	rndGen0.Normal(noise, m_len_xx);

	for (int i = 0; i < partLen; i++) {
		int num_of_symbol_part_blk = num_of_symbol_blk / partLen;
		xxidx = 0 + i * num_of_symbol_part_blk * 2;
		nnidx = 0 + i * num_of_symbol_part_blk * 2;
		for (int j = 0; j < num_of_symbol_part_blk; j++) {
			double treal = xx[xxidx] * selectH[i].real() - xx[xxidx + 1] * selectH[i].imag();
			double timag = xx[xxidx + 1] * selectH[i].real() + xx[xxidx] * selectH[i].imag();
			m_Nyy[xxidx] = treal + (sigma / SQRT2) * noise[nnidx];
			m_Nyy[xxidx + 1] = timag + (sigma / SQRT2) * noise[nnidx + 1];

			xxidx += 2;
			nnidx += 2;
		}
	}
}

void Linear_System::Soft_AWGN_Demodulation(double *yy, double *sym_prob)
{
	int i, q;
	double sqr_norm, sum;
	double temp;
	
	temp = 0.0;
	for (i = 0; i < num_symbol; i++)
	{
		sqr_norm = 0.0;
		
		double x_real = m_modem->output_symbol[i][0];
		double x_imag = m_modem->output_symbol[i][1];

		sqr_norm += ((yy[0] - x_real) * (yy[0] - x_real)
			      + (yy[1] - x_imag) * (yy[1] - x_imag)) / var;

		m_symbol_prob[i] = -sqr_norm;

	}
	
	temp = Seqmax(m_symbol_prob, num_symbol);

	for (q = 0; q < num_symbol; q++)
	{
		m_symbol_prob[q] = exp(m_symbol_prob[q] - temp);
	}

	//normalization
	sum = 0.0;
	for (q = 0; q < num_symbol; q++)
	{
		sum += m_symbol_prob[q];
	}

	for (q = 0; q < num_symbol; q++)
	{
		m_symbol_prob[q] /= sum;
	}

	for (q = 0; q < num_symbol; q++)
	{
		sym_prob[q] = m_symbol_prob[q];
	}

	ProbClip(sym_prob, num_symbol);
}
