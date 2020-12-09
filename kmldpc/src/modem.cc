#include "modem.h"

void CModem::Malloc(int code_no, char *file_name)
{
	int i, j;
	char temp_str[255] = {' '};
	int sym;
	int temp;
	double energy;
	FILE *fp;

	if ((fp = fopen(file_name, "r")) == nullptr)
	{
		fprintf(stderr, "\nCannot open %s", file_name);
		exit(3);
	}

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &len_input);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &len_output);

	fscanf(fp, "%s", temp_str);
	m_Es = 0.0;
	num_symbol = 1 << len_input;
	input_symbol = new int *[num_symbol];
	output_symbol = new double *[num_symbol];
	m_symbol_prob = new double[num_symbol];

	for (i = 0; i < num_symbol; i++)
	{
		fscanf(fp, "%d", &sym);

		input_symbol[i] = new int[len_input];
		temp = 0;
		for (j = 0; j < len_input; j++)
		{
			fscanf(fp, "%d", &input_symbol[i][j]);
			//temp = (temp<<1) | input_symbol[i][j];
			temp = (temp << 1) + input_symbol[i][j];
		}
		if (sym != temp || sym != i)
		{
			fprintf(stderr, "\nSym = %d is not the binary expression of %d!\n", sym, temp);
			fprintf(stderr, "\nSome symbols are missed or not in the right order!\n");
			system("pause");
			exit(3);
		}
		output_symbol[i] = new double[len_output];
		for (j = 0; j < len_output; j++)
		{
			fscanf(fp, "%lf", &output_symbol[i][j]);
		}

		for (j = 0; j < len_output; j++)
		{
			energy = output_symbol[i][j] * output_symbol[i][j];
			m_Es += energy;
		}
	}
	fclose(fp);

	m_Es = m_Es / num_symbol;
	for (i = 0; i < num_symbol; i++)
	{
		for (j = 0; j < len_output; j++)
		{
			output_symbol[i][j] /= sqrt(m_Es);
		}
	}

}

void CModem::Free()
{
	int i;
	for (i = 0; i < num_symbol; i++)
	{
		delete[] output_symbol[i];
		delete[] input_symbol[i];
	}
	delete[] output_symbol;
	delete[] input_symbol;
	delete[] m_symbol_prob;
}

void CModem::PrintCodeParameter(FILE *fp) const
{
	int i, j;
	FILE *fpp;

	fprintf(fp, "\n%%CModem Parameters\n");
	fprintf(fp, "%%%-20s = %d;\n", "bits_per_constellation_point", len_input);
	fprintf(fp, "%%%-20s = %d;\n", "symbols_per_constellation_point", len_output);
	fprintf(fp, "%%%-20s = %lf;\n", "Es", m_Es);
	fprintf(fp, "%%%-20s = %s;\n", "spacial_constell_file", "'spacial_constellation.txt'");

	if ((fpp = fopen("spacial_constellation.txt", "w")) == nullptr)
	{
		fprintf(stderr, "\nCannot open %s", "spacial_constellation.txt");
		exit(3);
	}

	fprintf(fpp, "\"symbol\"____\"binary_expression\"____\"mapping(real_image_..._real_image)\"\n");
	for (i = 0; i < num_symbol; i++)
	{
		fprintf(fpp, "%-20d", i);
		for (j = 0; j < len_input; j++)
		{
			fprintf(fpp, "%d ", input_symbol[i][j]);
		}

		fprintf(fpp, "%-5s", " ");
		for (j = 0; j < len_output; j++)
		{
			fprintf(fpp, "%13.10lf ", output_symbol[i][j]);
		}
		fprintf(fpp, "\n");
	}

	fclose(fpp);
}

void CModem::Mapping(int *bin_cc, double *xx, int num_bit_in_blk)
{
	int i, j, r, t;
	int num_point_in_blk;
	int symbol;

	num_point_in_blk = num_bit_in_blk / len_input;

	r = 0;
	t = 0;
	for (i = 0; i < num_point_in_blk; i++)
	{
		symbol = 0;
		for (j = 0; j < len_input; j++)
		{
			//symbol = (symbol << 1) | bin_cc[r];
			symbol = (symbol << 1) + bin_cc[r];
			r++;
		}
		for (j = 0; j < len_output; j++)
		{
			xx[t] = output_symbol[symbol][j];
			t++;
		}
	}
}

void CModem::Demapping(double *bitLin, double *symRin, double *bitLout, int num_bit_in_blk)
{
	int i, j, q, t;
	int num_point_in_blk;
	double P_0, P_1;
	double sum;

	num_point_in_blk = num_bit_in_blk / len_input;

	ProbClip(bitLin, num_bit_in_blk);
	ProbClip(symRin, num_symbol * num_point_in_blk);

	for (i = 0; i < num_point_in_blk; i++) //���� prob(x=0 | y);
	{
		for (q = 0; q < num_symbol; q++)
		{
			m_symbol_prob[q] = 1.0;
		}

		//convert bit extrinsic message to symbol message
		t = i * len_input;
		for (j = 0; j < len_input; j++)
		{
			for (q = 0; q < num_symbol; q++)
			{
				if (input_symbol[q][j] == 0)
					m_symbol_prob[q] *= bitLin[t];
				else
					m_symbol_prob[q] *= 1.0 - bitLin[t];
			}
			t++;
		}

		//compute symbol full message
		sum = 0.0;
		t = i * num_symbol;
		for (q = 0; q < num_symbol; q++)
		{
			m_symbol_prob[q] *= symRin[t];
			sum += m_symbol_prob[q];
			t++;
		}

		//normalized
		for (q = 0; q < num_symbol; q++)
			m_symbol_prob[q] /= sum;

		//compute bit extrinsic message for output
		t = i * len_input;
		for (j = 0; j < len_input; j++)
		{
			P_0 = 0;
			P_1 = 0;
			for (q = 0; q < num_symbol; q++)
			{
				if (input_symbol[q][j] == 0)
					P_0 += m_symbol_prob[q];
				else
					P_1 += m_symbol_prob[q];
			}
			P_0 /= bitLin[t];
			P_1 /= (1.0 - bitLin[t]);
			bitLout[t] = P_0 / (P_0 + P_1);
			t++;
		}
	} //end of for (i = 0; i < num_point_in_blk; i++)

	ProbClip(bitLout, num_bit_in_blk);
}

std::vector<std::complex<double>> CModem::getConstellation()
{
	std::vector<std::complex<double>> result(num_symbol);
	for (int i = 0; i < num_symbol; i++)
	{
		result[i] = std::complex<double>(output_symbol[i][0], output_symbol[i][1]);
	}
	return result;
}
