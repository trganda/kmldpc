#include "ModemLinearSystem.h"
#include "RandNum.h"

extern CLCRandNum rndGen0;
extern CWHRandNum rndGen1;

void Modem_Linear_System::Malloc(int len_cc, int code_no, char *file_name)
{
	int temp0, temp;
	char temp_str[80] = {' '};
	char mark[80];
	FILE *fp;
	
	m_len_cc = len_cc;

	if ((fp = fopen(file_name, "r")) == NULL){
		fprintf(stderr, "\nCannot open %s", file_name);
		exit(3);
	}
	
	sprintf(mark, "Modem_Linear_System***%d***PARAMETERS", code_no);
	while (strcmp(temp_str, mark) != 0)
		fscanf(fp, "%s", temp_str);
	
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%s", &m_constellation_file);
	fclose(fp);

	if ((fp = fopen("snrber.txt", "a+")) == nullptr){
		fprintf(stderr, "\n Cannot open the file!!!\n");
		exit(1);
	}
	fprintf(fp, "\n The constellation: %s", m_constellation_file);
	fclose(fp);

	if ((fp = fopen("snrresult.txt", "a+")) == nullptr){
		fprintf(stderr, "\n Cannot open the file!!!\n");
		exit(1);
	}
	fprintf(fp, "\n The constellation: %s", m_constellation_file);
	fclose(fp);

	m_modem.Malloc(0, m_constellation_file);	
		
	temp0 = m_len_cc / m_modem.len_input;
	if (m_len_cc % m_modem.len_input != 0)
	{
		fprintf(stderr, "\n(m_len_cc = %d) %% (len_input = %d) != 0 !\n", m_len_cc, m_modem.len_input);
		system("pause");
		exit(3);
	}

	m_len_xx = temp0 * m_modem.len_output;
	m_xx = new double[m_len_xx];

	temp = temp0 * m_modem.num_symbol;
	m_sym_prob = new double[temp];
	
	Lin_Sym.Malloc(&m_modem, m_len_xx, 0, file_name);

}

void Modem_Linear_System::Free()
{
	delete []m_xx;
	delete []m_sym_prob;
	
	m_modem.Free();
	Lin_Sym.Free();

}

void Modem_Linear_System::modem_linear_system(int *cc, double *m_sym_prob)
{	
	m_modem.Mapping(cc, m_xx, m_len_cc);

	Lin_Sym.AWGN_linear_system(m_xx, m_sym_prob);
}

void Modem_Linear_System::modem_linear_system_parition(int* cc, std::vector<std::complex<double>>& selectH)
{
	m_modem.Mapping(cc, m_xx, m_len_cc);

	Lin_Sym.Parition_HAWGN_system(m_xx, selectH);
}

void Modem_Linear_System::Soft_Demodulation(std::vector<std::pair<int, std::complex<double>>>& thetaList)
{
	std::vector<std::complex<double>> um = getRSymbol();

	int symbolPerPart = um.size() / thetaList.size();

	std::vector<std::complex<double>> tum(um.size());
	for (int i = 0; i < thetaList.size(); i++) {
		for (int j = 0; j < symbolPerPart; j++) {
			tum[j + i * symbolPerPart] = um[j + i * symbolPerPart] / thetaList[i].second;
			Lin_Sym.m_yy[0] = tum[j + i * symbolPerPart].real();
			Lin_Sym.m_yy[1] = tum[j + i * symbolPerPart].imag();
			int temp = (j + i * symbolPerPart) * Lin_Sym.m_modem->num_symbol;
			Lin_Sym.Soft_AWGN_Demodulation(Lin_Sym.m_yy, (m_sym_prob + temp));
		}
	}
}

std::vector<std::complex<double>> Modem_Linear_System::getRSymbol()
{
	std::vector<std::complex<double>> um(m_len_xx / 2);
	for (int i = 0; i < um.size(); i++) {
		um[i] = std::complex<double>(Lin_Sym.m_Nyy[i * 2], Lin_Sym.m_Nyy[i * 2 + 1]);
	}
	return um;
}

