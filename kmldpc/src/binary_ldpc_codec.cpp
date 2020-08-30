#include "binary_ldpc_codec.h"
#include "utility.h"
#include "randnum.h"

extern CLCRandNum rndGen0;

CBinaryLDPCCodec::CBinaryLDPCCodec()
= default;

CBinaryLDPCCodec::~CBinaryLDPCCodec()
= default;

void CBinaryLDPCCodec::Malloc(int code_no, char *file_name)
{	
	int i, j;
	int row_no, row_deg, col_no;
	Edge *temp_edge;
	FILE *fq=nullptr;
	
	char temp_str[80] = {' '};
	char mark[80];
	FILE *fp;

	sprintf(mark, "LDPC***%d***PARAMETERS",code_no);

	if ((fp = fopen(file_name, "r")) == nullptr){
		fprintf(stderr, "\nCannot open %s", file_name);
		exit(3);
	}

	while (strcmp(temp_str, mark) != 0)
		fscanf(fp, "%s", temp_str);

//decdoing	
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &m_max_iter);

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &m_encoder_active);

	//file_name_of_H	
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%s", m_file_name_of_H);

	fclose(fp);

//Read H from file temp_str
	if ((fp = fopen(m_file_name_of_H, "r")) == nullptr){
		fprintf(stderr, "\nCannot open %s", m_file_name_of_H);
		exit(0);
	}

	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d %d %d", &m_num_row, &m_num_col, &m_codechk);
	m_codelen = m_num_col;
	m_codedim = m_codelen - m_codechk;
	m_coderate = (double)m_codedim / m_codelen;

	m_row_head = new Edge[m_num_row];
	m_col_head = new Edge[m_num_col];

	for (i = 0; i < m_num_row; i++){
		(m_row_head+i)->m_row_no = i;
		(m_row_head+i)->m_col_no = -1;
		(m_row_head+i)->left = m_row_head+i;
		(m_row_head+i)->right = m_row_head+i;
		(m_row_head+i)->up = m_row_head+i;
		(m_row_head+i)->down = m_row_head+i;
	}

	for (i = 0; i < m_num_col; i++){
		(m_col_head+i)->m_row_no = -1;
		(m_col_head+i)->m_col_no = i;
		(m_col_head+i)->left = m_col_head+i;
		(m_col_head+i)->right = m_col_head+i;
		(m_col_head+i)->up = m_col_head+i;
		(m_col_head+i)->down = m_col_head+i;
	}

	fscanf(fp, "%s", temp_str);
	for (i = 0; i < m_num_row; i++){
		fscanf(fp, "%d %d", &row_no, &row_deg);
		for (j = 0; j < row_deg; j++){
			temp_edge = new Edge;
			temp_edge->m_row_no = row_no;
			fscanf(fp, "%d", &col_no);
			temp_edge->m_col_no = col_no;

			temp_edge->right = (m_row_head+i)->right;
			(m_row_head+i)->right = temp_edge;
			temp_edge->left = m_row_head+i;
			(temp_edge->right)->left = temp_edge;

			temp_edge->down = (m_col_head+col_no)->down;
			(m_col_head+col_no)->down = temp_edge;
			temp_edge->up = m_col_head+col_no;
			(temp_edge->down)->up = temp_edge;
		}
	}

	fclose(fp);
#ifdef _DEBUG
	if ((fq = fopen("a.txt", "a+")) == NULL){
		fprintf(stderr, "\nCannot open %s", "a.txt");
		exit(0);
	}

	for(i = 0; i < m_num_row; i++){
		fprintf(fp,  "%d: %d* %d --- %d* %d --- %d* %d --- %d* %d --- %d* %d\n",i,
		             (m_row_head + i) -> m_row_no,(m_row_head + i)-> m_col_no, 
		             (m_row_head + i) -> left -> m_row_no,(m_row_head + i) -> left -> m_col_no,
					 (m_row_head + i) -> right -> m_row_no,(m_row_head + i) -> right -> m_col_no,
					 (m_row_head + i) -> up -> m_row_no,(m_row_head + i) -> up -> m_col_no,
					 (m_row_head + i) -> down -> m_row_no,(m_row_head + i) -> down -> m_col_no);
	}

    for(i = 0; i < m_num_col; i++){
		fprintf(fp,  "%d: %d* %d --- %d* %d --- %d* %d --- %d* %d --- %d* %d\n",i,
		             (m_col_head + i) -> m_row_no,(m_col_head + i)-> m_col_no, 
		             (m_col_head + i) -> left -> m_row_no,(m_col_head + i) -> left -> m_col_no,
					 (m_col_head + i) -> right -> m_row_no,(m_col_head + i) -> right -> m_col_no,
					 (m_col_head + i) -> up -> m_row_no,(m_col_head + i) -> up -> m_col_no,
					 (m_col_head + i) -> down -> m_row_no,(m_col_head + i) -> down -> m_col_no);
	}
	fclose(fq);
#endif
	if (m_encoder_active == 1)
		SystH();

	m_cc_hat = new int[m_codelen];

	return;
}

void CBinaryLDPCCodec::PrintCodeParameter( FILE *fp )
{
	fprintf(fp, "%%#######CBinaryLDPCCodec Parameters#######\n");
	fprintf(fp, "%%%-20s = [%d, %d, %d, %g];\n", "[n, k, m, r]",
		m_codelen, m_codedim, m_codechk, m_coderate);
	fprintf(fp, "%%%m_max_iter = %d;\n", m_max_iter);
	fprintf(fp, "%%%m_encoder_active = %d;\n", m_encoder_active);
	fprintf(fp, "%%%m_file_name_of_H = '%s';\n", m_file_name_of_H);
}

void CBinaryLDPCCodec::Free_Tanner_Graph()
{
	int i;
	Edge *temp_edge;

	for (i = 0; i < m_num_row; i++)
	{
		while ((m_row_head+i)->right->m_col_no != -1)
		{
			temp_edge = (m_row_head+i)->right;
			(m_row_head+i)->right = temp_edge->right;
			delete temp_edge;
		}
	}

	delete []m_row_head;
	delete []m_col_head;
}

void CBinaryLDPCCodec::Free()
{
	int i=0;
	Free_Tanner_Graph();
	
	if (m_encoder_active != 0)
	{
		for (i = 0; i < m_num_row; i++)
		{
			delete []m_encH[i];
		}
		delete []m_encH;
	}
	
	delete []m_cc_hat;

	return;
}

void CBinaryLDPCCodec::SystH()
{
	int i, j, ii, jj, m, n; 
	int temp;
	int flag;
	int *tempP;
	char **tempH;
	Edge *p_edge;

	m_decH = new char*[m_num_row];
	for (i = 0; i < m_num_row; i++)
		m_decH[i] = new char[m_num_col];

	for (i = 0; i < m_num_row; i++){
		for (j = 0; j < m_num_col; j++)
			m_decH[i][j] = 0;
	}

	for (i = 0; i < m_num_row; i++){
		p_edge = (m_row_head + i)->right;
		while (p_edge->m_col_no != -1){
			m_decH[i][p_edge->m_col_no] = 1;
			p_edge = p_edge->right;
		}
	}

	tempP = new int[m_num_col];
	for (j = 0; j < m_num_col; j++)
		tempP[j] = j;
	tempH = new char*[m_num_row];
	for (i = 0; i < m_num_row; i++)
		tempH[i] = new char[m_num_col];

	for (i = 0; i < m_num_row; i++){
		for (j = 0; j < m_num_col; j++)
			tempH[i][j] = m_decH[i][j];
	}

	m_encH = new char*[m_num_row];
	for (i = 0; i < m_num_row; i++)
	{
		m_encH[i] = new char[m_num_col];
	}
	for (i = 0; i < m_num_row; i++)
	{
		for (j = 0; j < m_num_col; j++)
		{
			m_encH[i][j] = m_decH[i][j];
		}
	}

	m_codechk = 0;

	for (i = 0; i < m_num_row; i++)
	{
		flag = 0;
		for (jj = i; jj < m_num_col; jj++)
		{
			for (ii = i; ii < m_num_row; ii++)
			{
				if (m_encH[ii][jj] != 0)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 1)
			{
				m_codechk++;
				break;
			}
		}

		if (flag == 0)
			break;
		else
		{
//swap i and ii row
			if (ii != i)
			{
				for (n = 0; n < m_num_col; n++)
				{
					temp = m_encH[i][n];
					m_encH[i][n] = m_encH[ii][n];
					m_encH[ii][n] = temp;
				}
			}
//swap i and jj col
			if (jj != i)
			{
				temp = tempP[i];
				tempP[i] = tempP[jj];
				tempP[jj] = temp;

				for (m = 0; m < m_num_row; m++)
				{
					temp = m_encH[m][i];
					m_encH[m][i] = m_encH[m][jj];
					m_encH[m][jj] = temp;
				}
			}
//elimination
			for (m = 0; m < m_num_row; m++)
			{
				if (m != i && m_encH[m][i] == 1)
				{
					for (n = 0; n < m_num_col; n++)
						m_encH[m][n] ^= m_encH[i][n];
				}
			}
		}
	}

	//for (i = 0; i < m_num_col; i++)
	//{
	//	if (tempP[i] != i)
	//	{
	//		fprintf(stderr, "\nWarning: please make sure the information bits on the right!");
	//	}
	//}

	
	for (j = 0; j < m_num_col; j++)
	{
		for (i = 0; i < m_num_row; i++)
		{
			m_decH[i][j] = tempH[i][tempP[j]];
		}
	}

	Free_Tanner_Graph();

	m_row_head = new Edge[m_num_row];
	m_col_head = new Edge[m_num_col];

	for (i = 0; i < m_num_row; i++)
	{
		(m_row_head + i)->m_row_no = i;
		(m_row_head + i)->m_col_no = -1;
		(m_row_head + i)->left = m_row_head+i;
		(m_row_head + i)->right = m_row_head+i;
		(m_row_head + i)->up = m_row_head+i;
		(m_row_head + i)->down = m_row_head+i;
	}

	for (i = 0; i < m_num_col; i++)
	{
		(m_col_head + i)->m_row_no = -1;
		(m_col_head + i)->m_col_no = i;
		(m_col_head + i)->left = m_col_head+i;
		(m_col_head + i)->right = m_col_head+i;
		(m_col_head + i)->up = m_col_head+i;
		(m_col_head + i)->down = m_col_head+i;
	}

	for (i = 0; i < m_num_row; i++)
	{
		for (j = 0; j < m_num_col; j++)
		{
			if (m_decH[i][j] != 0)
			{
				p_edge = new Edge;
				p_edge->m_row_no = i;
				p_edge->m_col_no = j;

				p_edge->right = (m_row_head + i)->right;
				(m_row_head + i)->right = p_edge;
				p_edge->left = m_row_head + i;
				(p_edge->right)->left = p_edge;

				p_edge->down = (m_col_head + j)->down;
				(m_col_head + j)->down = p_edge;
				p_edge->up = m_col_head + j;
				(p_edge->down)->up = p_edge;
			}
		}
	}
	
	m_codelen = m_num_col;
	m_codedim = m_codelen - m_codechk;
	m_coderate = (double)m_codedim / m_codelen;

	delete []tempP;
	for (i = 0; i < m_num_row; i++)
	{
		delete[]tempH[i];
		delete[]m_decH[i];
	}
	delete []tempH;
	delete []m_decH;
}

void CBinaryLDPCCodec::SystH2()
{
	int i, j, ii, jj, m, n; 
	int temp;
	int flag;
	int *tempP;
	char **tempH;
	Edge *p_edge;

	m_decH = new char*[m_num_row];
	for (i = 0; i < m_num_row; i++)
		m_decH[i] = new char[m_num_col];

	for (i = 0; i < m_num_row; i++){
		for (j = 0; j < m_num_col; j++)
			m_decH[i][j] = 0;
	}

	for (i = 0; i < m_num_row; i++){
		p_edge = (m_row_head + i)->right;
		while (p_edge->m_col_no != -1){
			m_decH[i][p_edge->m_col_no] = 1;
			p_edge = p_edge->right;
		}
	}

	tempP = new int[m_num_col];
	for (j = 0; j < m_num_col; j++)
		tempP[j] = j;
	tempH = new char*[m_num_row];
	for (i = 0; i < m_num_row; i++)
		tempH[i] = new char[m_num_col];

	for (i = 0; i < m_num_row; i++){
		for (j = 0; j < m_num_col; j++)
			tempH[i][j] = m_decH[i][j];
	}

	m_encH = new char*[m_num_row];
	for (i = 0; i < m_num_row; i++)
	{
		m_encH[i] = new char[m_num_col];
	}
	for (i = 0; i < m_num_row; i++)
	{
		for (j = 0; j < m_num_col; j++)
		{
			m_encH[i][j] = m_decH[i][j];
		}
	}

	m_codechk = 0;

	for (i = 0; i < m_num_row; i++)
	{
		flag = 0;
		for (jj = i; jj < m_num_col; jj++)
		{
			for (ii = i; ii < m_num_row; ii++)
			{
				if (m_encH[ii][jj] != 0)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 1)
			{
				m_codechk++;
				break;
			}
		}

		if (flag == 0)
			break;
		else
		{
			//swap i and ii row
			if (ii != i)
			{
				for (n = 0; n < m_num_col; n++)
				{
					temp = m_encH[i][n];
					m_encH[i][n] = m_encH[ii][n];
					m_encH[ii][n] = temp;
				}
			}
			//swap i and jj col
			if (jj != i)
			{
				temp = tempP[i];
				tempP[i] = tempP[jj];
				tempP[jj] = temp;

				for (m = 0; m < m_num_row; m++)
				{
					temp = m_encH[m][i];
					m_encH[m][i] = m_encH[m][jj];
					m_encH[m][jj] = temp;
				}
			}
			//elimination
			for (m = 0; m < m_num_row; m++)
			{
				if (m != i && m_encH[m][i] == 1)
				{
					for (n = 0; n < m_num_col; n++)
						m_encH[m][n] ^= m_encH[i][n];
				}
			}
		}
	}

	//for (i = 0; i < m_num_col; i++)
	//{
	//	if (tempP[i] != i)
	//	{
	//		fprintf(stderr, "\nWarning: please make sure the information bits on the right!");
	//	}
	//}


	for (j = 0; j < m_num_col; j++)
	{
		for (i = 0; i < m_num_row; i++)
		{
			m_decH[i][j] = tempH[i][tempP[j]];
		}
	}

	Free_Tanner_Graph();

	m_row_head = new Edge[m_num_row];
	m_col_head = new Edge[m_num_col];

	for (i = 0; i < m_num_row; i++)
	{
		(m_row_head + i)->m_row_no = i;
		(m_row_head + i)->m_col_no = -1;
		(m_row_head + i)->left = m_row_head+i;
		(m_row_head + i)->right = m_row_head+i;
		(m_row_head + i)->up = m_row_head+i;
		(m_row_head + i)->down = m_row_head+i;
	}

	for (i = 0; i < m_num_col; i++)
	{
		(m_col_head + i)->m_row_no = -1;
		(m_col_head + i)->m_col_no = i;
		(m_col_head + i)->left = m_col_head+i;
		(m_col_head + i)->right = m_col_head+i;
		(m_col_head + i)->up = m_col_head+i;
		(m_col_head + i)->down = m_col_head+i;
	}

	for (i = 0; i < m_num_row; i++)
	{
		for (j = 0; j < m_num_col; j++)
		{
			if (m_decH[i][j] != 0)
			{
				p_edge = new Edge;
				p_edge->m_row_no = i;
				p_edge->m_col_no = j;

				p_edge->right = (m_row_head + i)->right;
				(m_row_head + i)->right = p_edge;
				p_edge->left = m_row_head + i;
				(p_edge->right)->left = p_edge;

				p_edge->down = (m_col_head + j)->down;
				(m_col_head + j)->down = p_edge;
				p_edge->up = m_col_head + j;
				(p_edge->down)->up = p_edge;
			}
		}
	}

	m_codelen = m_num_col;
	m_codedim = m_codelen - m_codechk;
	m_coderate = (double)m_codedim / m_codelen;

	delete []tempP;
	for (i = 0; i < m_num_row; i++)
	{
		delete[]tempH[i];
		delete[]m_decH[i];
	}
	delete []tempH;
	delete []m_decH;
}

void CBinaryLDPCCodec::Kite_Encoder(int *uu, int *cc)
{
	int i,j=0,t;

	Edge *p_edge;
//codeword = [parity_check_bits information_bits]
	switch (m_encoder_active){
	case 0:
		for (i = 0; i < m_codedim; i++)
			uu[i] = 0;
		for (i = 0; i < m_codelen; i++)
			cc[i] = 0;
		break;
	case 1:
		for (t = 0; t < m_codedim; t++)
			cc[t] = uu[t];

		for(t = m_codelen - m_codechk ; t < m_codelen ; t ++)
		{
			cc[t] = 0;
			p_edge = this->m_row_head[t - m_codelen + m_codechk ].left;
			while(p_edge->m_col_no != -1)
			{ 
				if(p_edge->m_col_no < t)//encoding care!!
				{
					cc[t] ^= cc[p_edge->m_col_no];
				}
				p_edge = p_edge->left;
			}
		}
		/*for (t = 0; t < m_codechk; t++){
			cc[t] = 0;
			for (j = m_codechk; j < m_codelen; j++)
				cc[t] ^= (cc[j] & m_encH[t][j]);
		}*/
		break;
	default:
		break;
	}

	return;
}

void CBinaryLDPCCodec::Encoder(int *uu, int *cc)
{
	int i, j, t;

//codeword = [parity_check_bits information_bits]
	switch (m_encoder_active){
	case 0:
		for (i = 0; i < m_codedim; i++)
			uu[i] = 0;
		for (i = 0; i < m_codelen; i++)
			cc[i] = 0;
		break;
	case 1:
		for (t = m_codechk; t < m_codelen; t++)
			cc[t] = uu[t- m_codechk];

		for (t = 0; t < m_codechk; t++){
			cc[t] = 0;
			for (j = m_codechk; j < m_codelen; j++)
				cc[t] ^= (cc[j] & m_encH[t][j]);
		}
		break;
	default:
		break;
	}

	return;
}

void CBinaryLDPCCodec::OnlineEncoder(int *uu, int *cc)
{
	int i, row_no;
	Edge *p_edge;

//codeword = [information_bits parity_check_bits]
	for (i = 0; i < m_codedim; i++)
		cc[i] = uu[i];

	for (row_no = 0; row_no < m_codechk; row_no++){
		cc[i] = 0;
		p_edge = (m_row_head+row_no)->left;
		while (p_edge->m_col_no != -1){
			if (p_edge->m_col_no < i)
				cc[i] ^= cc[p_edge->m_col_no];
			p_edge = p_edge->left;
		}
		i++;
	}
	
	return;
}

void CBinaryLDPCCodec::InitMsg()
{
	int i;
	Edge *p_edge;
	
	for (i = 0; i < m_num_col; i++){
		p_edge = (m_col_head+i)->down;
		while (p_edge->m_row_no != -1){
			p_edge->m_c2v[0] = 0.5;
			p_edge->m_c2v[1] = 0.5;
			p_edge = p_edge->down;
		}
	}

	return;
}

int CBinaryLDPCCodec::Decoder(double *M2V, int *uu_hat)
{
	int i;
	int iter;
	int parity_check;
	double temp0, temp1, temp_sum;

	Edge *p_edge;

//initialization
	InitMsg();
//iteration
	for (iter = 0; iter < m_max_iter; iter++){
//from vnode to cnode
		for (i = 0; i < m_num_col; i++){
//forward
			p_edge = (m_col_head+i)->down;
			p_edge->m_alpha[0] = M2V[i];
			p_edge->m_alpha[1] = 1.0 - M2V[i];

			while (p_edge->m_row_no != -1){
				p_edge->down->m_alpha[0] = p_edge->m_alpha[0] * p_edge->m_c2v[0];
				p_edge->down->m_alpha[1] = p_edge->m_alpha[1] * p_edge->m_c2v[1];
				temp_sum = p_edge->down->m_alpha[0] + p_edge->down->m_alpha[1];
				p_edge->down->m_alpha[0] /= temp_sum;
				p_edge->down->m_alpha[1] /= temp_sum;
				p_edge = p_edge->down;
			}
//hard decision
			if ((m_col_head+i)->m_alpha[0] > (m_col_head+i)->m_alpha[1])
				m_cc_hat[i] = 0;
			else
				m_cc_hat[i] = 1;

//backward
			p_edge = (m_col_head+i)->up;
			p_edge->m_beta[0] = 1.0;
			p_edge->m_beta[1] = 1.0;
			while (p_edge->m_row_no != -1){
				temp0 = p_edge->m_alpha[0] * p_edge->m_beta[0];
				temp1 = p_edge->m_alpha[1] * p_edge->m_beta[1];
				temp_sum = temp0 + temp1;

				p_edge->m_v2c[0] = temp0 / temp_sum;
				p_edge->m_v2c[1] = temp1 / temp_sum;

				
				p_edge->up->m_beta[0] = p_edge->m_beta[0] * p_edge->m_c2v[0];
				p_edge->up->m_beta[1] = p_edge->m_beta[1] * p_edge->m_c2v[1];

				temp_sum = p_edge->up->m_beta[0] + p_edge->up->m_beta[1];
				p_edge->up->m_beta[0] /= temp_sum;
				p_edge->up->m_beta[1] /= temp_sum;

				p_edge = p_edge->up;
			}
		}

		for (i = 0; i < m_codedim; i++)
		{
			uu_hat[i] = m_cc_hat[i+m_codechk];
//			uu_hat[i] = m_cc_hat[i];//kiteʹ��
		}
//parity checking
		m_success = 1;
		for (i = 0; i < m_num_row; i++){
			parity_check = 0;
			p_edge = (m_row_head+i)->right;
			while (p_edge->m_col_no != -1){
				parity_check = parity_check ^ m_cc_hat[p_edge->m_col_no];
				p_edge = p_edge->right;
			}
			if (parity_check != 0){
				m_success = 0;
				break;
			}
		}

		if (m_success == 1)
			break;			

//from c node to v node
		for (i = 0; i < m_num_row; i++){
//forward
			p_edge = (m_row_head + i)->right;
			p_edge->m_alpha[0] = 1.0;
			p_edge->m_alpha[1] = 0.0;
			while (p_edge->m_col_no != -1){//over trellis with two states
				p_edge->right->m_alpha[0] = p_edge->m_alpha[0] * p_edge->m_v2c[0] + p_edge->m_alpha[1] * p_edge->m_v2c[1];
				p_edge->right->m_alpha[1] = p_edge->m_alpha[0] * p_edge->m_v2c[1] + p_edge->m_alpha[1] * p_edge->m_v2c[0];
				temp_sum = p_edge->right->m_alpha[0] + p_edge->right->m_alpha[1];
				p_edge->right->m_alpha[0] /= temp_sum;
				p_edge->right->m_alpha[1] /= temp_sum;

				p_edge = p_edge->right;
			}
//backward
			p_edge = (m_row_head + i)->left;
			p_edge->m_beta[0] = 1.0;
			p_edge->m_beta[1] = 0.0;
			while (p_edge->m_col_no != -1){
				temp0 = p_edge->m_alpha[0] * p_edge->m_beta[0] + p_edge->m_alpha[1] * p_edge->m_beta[1];
				temp1 = p_edge->m_alpha[0] * p_edge->m_beta[1] + p_edge->m_alpha[1] * p_edge->m_beta[0];
				temp_sum = temp0 + temp1;
				p_edge->m_c2v[0] = temp0 / temp_sum;
				p_edge->m_c2v[1] = temp1 / temp_sum;
				
				if (p_edge->m_c2v[0] > 1.0 - SMALLPROB)
					p_edge->m_c2v[0] = 1.0 - SMALLPROB;
				if (p_edge->m_c2v[0] < SMALLPROB)
					p_edge->m_c2v[0] = SMALLPROB;

				p_edge->m_c2v[1] = 1.0 - p_edge->m_c2v[0];

				p_edge->left->m_beta[0] = p_edge->m_beta[0] * p_edge->m_v2c[0] + p_edge->m_beta[1] * p_edge->m_v2c[1];
				p_edge->left->m_beta[1] = p_edge->m_beta[0] * p_edge->m_v2c[1] + p_edge->m_beta[1] * p_edge->m_v2c[0];
				
				temp_sum = p_edge->left->m_beta[0] + p_edge->left->m_beta[1];
				p_edge->left->m_beta[0] /= temp_sum;
				p_edge->left->m_beta[1] /= temp_sum;
				
				p_edge = p_edge->left;
			}
		}
	}

	m_detection_error = 1-m_success;//added by liangchulong @2013-01-23-19-04
	return iter + (iter<m_max_iter);
}

/*
 * // ͳ��У�������m_cc_hat�ĳ˻�H*m_cc_hat'��1�ĸ���
 */
int CBinaryLDPCCodec::parityCheck(double * M2V)
{
	int i;
	int iter;
	int parity_check;
	double temp0, temp1, temp_sum;

	Edge *p_edge;

	//initialization
	InitMsg();
	//from vnode to cnode
	for (i = 0; i < m_num_col; i++) {
		//hard decision
		if ((m_col_head + i)->m_alpha[0] > (m_col_head + i)->m_alpha[1])
			m_cc_hat[i] = 0;
		else
			m_cc_hat[i] = 1;
	}

	// m_success = 1;
	int errCnt = 0;
	for (i = 0; i < m_num_row; i++) {
		parity_check = 0;
		p_edge = (m_row_head + i)->right;
		while (p_edge->m_col_no != -1) {
			parity_check = parity_check ^ m_cc_hat[p_edge->m_col_no];
			p_edge = p_edge->right;
		}
		errCnt += parity_check;
	}
	
	return errCnt;
}

int CBinaryLDPCCodec::parityCheck(int* rr)
{
	int parity_check, count, i;
	Edge* p_edge;
	count = 0;
	for (i = 0; i < m_num_row; i++) {
		parity_check = 0;
		p_edge = (m_row_head + i)->right;
		while (p_edge->m_col_no != -1) {
			if (p_edge->m_col_no < m_num_col)
				parity_check = parity_check ^ rr[p_edge->m_col_no];
			p_edge = p_edge->right;
		}
		if (parity_check != 0) {
			count++;
		}
	}
	return count;
}

/***************************************************************************
������SoftInHardOutDecoder(double *U2V, double *M2V, int *uu_hat)
���ܣ�
���������*U2V    ---- �������
          *M2V    ---- ������Ϣ�������м�����غ������
���������*uu_hat ---- ������   
****************************************************************************/

void CBinaryLDPCCodec::SoftInHardOutDetectorDecoder(double *U2V, double *M2V, int *uu_hat)
{
	int i;
	int iter;
	int parity_check;
	double temp0, temp1, temp_sum;

	Edge *p_edge;

//initialization
	InitMsg();

	for (iter = 0; iter < m_max_iter; iter++){

//from vnode to cnode
		for (i = 0; i < m_num_col; i++){
//forward
			p_edge = (m_col_head+i)->down;
			if (i < m_codedim){
				temp0 = U2V[i] * M2V[i];
				temp1 = (1.0 - U2V[i]) * (1.0 - M2V[i]);
				temp_sum = temp0 + temp1;
				p_edge->m_alpha[0] = temp0 / temp_sum;
				p_edge->m_alpha[1] = temp1 / temp_sum;
			}			
			else{
				p_edge->m_alpha[0] = M2V[i];
				p_edge->m_alpha[1] = 1.0 - M2V[i];
			}

			while (p_edge->m_row_no != -1){
				if (p_edge->m_row_no < m_num_row){
					p_edge->down->m_alpha[0] = p_edge->m_alpha[0] * p_edge->m_c2v[0];
					p_edge->down->m_alpha[1] = p_edge->m_alpha[1] * p_edge->m_c2v[1];
					temp_sum = p_edge->down->m_alpha[0] + p_edge->down->m_alpha[1];
					p_edge->down->m_alpha[0] /= temp_sum;
					p_edge->down->m_alpha[1] /= temp_sum;
				}
				else{
					p_edge->down->m_alpha[0] = p_edge->m_alpha[0];
					p_edge->down->m_alpha[1] = p_edge->m_alpha[1];
				}
				p_edge = p_edge->down;
			}
//hard decision
			if ((m_col_head+i)->m_alpha[0] > (m_col_head+i)->m_alpha[1])
				m_cc_hat[i] = 0;
			else
				m_cc_hat[i] = 1;

//backward
			p_edge = (m_col_head+i)->up;
			p_edge->m_beta[0] = 1.0;
			p_edge->m_beta[1] = 1.0;
			while (p_edge->m_row_no != -1){
				if (p_edge->m_row_no < m_num_row){
					temp0 = p_edge->m_alpha[0] * p_edge->m_beta[0];
					temp1 = p_edge->m_alpha[1] * p_edge->m_beta[1];
					temp_sum = temp0 + temp1;
					p_edge->m_v2c[0] = temp0 / temp_sum;
					p_edge->m_v2c[1] = temp1 / temp_sum;
				
					p_edge->up->m_beta[0] = p_edge->m_beta[0] * p_edge->m_c2v[0];
					p_edge->up->m_beta[1] = p_edge->m_beta[1] * p_edge->m_c2v[1];

					temp_sum = p_edge->up->m_beta[0] + p_edge->up->m_beta[1];
					p_edge->up->m_beta[0] /= temp_sum;
					p_edge->up->m_beta[1] /= temp_sum;
				}
				else{
					p_edge->up->m_beta[0] = p_edge->m_beta[0];
					p_edge->up->m_beta[1] = p_edge->m_beta[1];
				}
				p_edge = p_edge->up;
			}
		}

		for (i = 0; i < m_codedim; i++)
			uu_hat[i] = m_cc_hat[m_codechk+i];

//parity checking
		m_success = 1;
		for (i = 0; i < m_num_row; i++){
			parity_check = 0;
			p_edge = (m_row_head+i)->right;
			while (p_edge->m_col_no != -1){
				if (p_edge->m_col_no < m_num_col)
					parity_check = parity_check ^ m_cc_hat[p_edge->m_col_no];
				p_edge = p_edge->right;
			}
			if (parity_check != 0){
				m_success = 0;
				break;
			}
		}

		if (m_success == 1)	
			break;			

//from c node to v node
		for (i = 0; i < m_num_row; i++){
//forward
			p_edge = (m_row_head + i)->right;
			p_edge->m_alpha[0] = 1.0;
			p_edge->m_alpha[1] = 0.0;
			while (p_edge->m_col_no != -1){//over trellis with two states
				p_edge->right->m_alpha[0] = p_edge->m_alpha[0] * p_edge->m_v2c[0] + p_edge->m_alpha[1] * p_edge->m_v2c[1];
				p_edge->right->m_alpha[1] = p_edge->m_alpha[0] * p_edge->m_v2c[1] + p_edge->m_alpha[1] * p_edge->m_v2c[0];
				temp_sum = p_edge->right->m_alpha[0] + p_edge->right->m_alpha[1];
				p_edge->right->m_alpha[0] /= temp_sum;
				p_edge->right->m_alpha[1] /= temp_sum;
	
				p_edge = p_edge->right;
			}
//backward
			p_edge = (m_row_head + i)->left;
			p_edge->m_beta[0] = 1.0;
			p_edge->m_beta[1] = 0.0;
			while (p_edge->m_col_no != -1){
				temp0 = p_edge->m_alpha[0] * p_edge->m_beta[0] + p_edge->m_alpha[1] * p_edge->m_beta[1];
				temp1 = p_edge->m_alpha[0] * p_edge->m_beta[1] + p_edge->m_alpha[1] * p_edge->m_beta[0];
				temp_sum = temp0 + temp1;
				p_edge->m_c2v[0] = temp0 / temp_sum;
				p_edge->m_c2v[1] = temp1 / temp_sum;
				
				if (p_edge->m_c2v[0] > 1.0 - SMALLPROB)
					p_edge->m_c2v[0] = 1.0 - SMALLPROB;
				if (p_edge->m_c2v[0] < SMALLPROB)
					p_edge->m_c2v[0] = SMALLPROB;

				p_edge->m_c2v[1] = 1.0 - p_edge->m_c2v[0];

				p_edge->left->m_beta[0] = p_edge->m_beta[0] * p_edge->m_v2c[0] + p_edge->m_beta[1] * p_edge->m_v2c[1];
				p_edge->left->m_beta[1] = p_edge->m_beta[0] * p_edge->m_v2c[1] + p_edge->m_beta[1] * p_edge->m_v2c[0];
				
				temp_sum = p_edge->left->m_beta[0] + p_edge->left->m_beta[1];
				p_edge->left->m_beta[0] /= temp_sum;
				p_edge->left->m_beta[1] /= temp_sum;
				
				p_edge = p_edge->left;
			}
		}
	}

	m_detection_error = 1-m_success;
}

int CBinaryLDPCCodec::SoftInSoftOut(double *M2V, int *uu_hat,  double* cc_hat_soft)
{
	int i;
	int iter;
	int parity_check;
	double temp0, temp1, temp_sum;

	Edge *p_edge;

	//initialization
	InitMsg();
	//iteration
	for (iter = 0; iter < m_max_iter; iter++){
		//from vnode to cnode �Ⱥŵ��Ӻ�
		for (i = 0; i < m_num_col; i++){
			//forward
			p_edge = (m_col_head+i)->down;
			p_edge->m_alpha[0] = M2V[i];
			p_edge->m_alpha[1] = 1.0 - M2V[i];

			while (p_edge->m_row_no != -1){
				p_edge->down->m_alpha[0] = p_edge->m_alpha[0] * p_edge->m_c2v[0];
				p_edge->down->m_alpha[1] = p_edge->m_alpha[1] * p_edge->m_c2v[1];
				temp_sum = p_edge->down->m_alpha[0] + p_edge->down->m_alpha[1];
				p_edge->down->m_alpha[0] /= temp_sum;
				p_edge->down->m_alpha[1] /= temp_sum;
				p_edge = p_edge->down;
			}
			//hard decision
			if ((m_col_head+i)->m_alpha[0] > (m_col_head+i)->m_alpha[1])
				m_cc_hat[i] = 0;
			else
				m_cc_hat[i] = 1;

			//backward
			p_edge = (m_col_head+i)->up;
			p_edge->m_beta[0] = 1.0;
			p_edge->m_beta[1] = 1.0;
			while (p_edge->m_row_no != -1){
				temp0 = p_edge->m_alpha[0] * p_edge->m_beta[0];
				temp1 = p_edge->m_alpha[1] * p_edge->m_beta[1];
				temp_sum = temp0 + temp1;

				p_edge->m_v2c[0] = temp0 / temp_sum;
				p_edge->m_v2c[1] = temp1 / temp_sum;


				p_edge->up->m_beta[0] = p_edge->m_beta[0] * p_edge->m_c2v[0];
				p_edge->up->m_beta[1] = p_edge->m_beta[1] * p_edge->m_c2v[1];

				temp_sum = p_edge->up->m_beta[0] + p_edge->up->m_beta[1];
				p_edge->up->m_beta[0] /= temp_sum;
				p_edge->up->m_beta[1] /= temp_sum;

				p_edge = p_edge->up;
			}
		}

		for (i = 0; i < m_codedim; i++)
		{
			uu_hat[i] = m_cc_hat[i+m_codechk];
//			uu_hat[i] = m_cc_hat[i];
		}
		//parity checking
		m_success = 1;
		for (i = 0; i < m_num_row; i++){
			parity_check = 0;
			p_edge = (m_row_head+i)->right;
			while (p_edge->m_col_no != -1){
				parity_check = parity_check ^ m_cc_hat[p_edge->m_col_no];
				p_edge = p_edge->right;
			}
			if (parity_check != 0){
				m_success = 0;
				break;
			}
		}

		if (m_success == 1)	
			break;			

		//from c node to v node   �Ӻŵ��Ⱥ�
		for (i = 0; i < m_num_row; i++){
			//forward
			p_edge = (m_row_head + i)->right;
			p_edge->m_alpha[0] = 1.0;
			p_edge->m_alpha[1] = 0.0;
			while (p_edge->m_col_no != -1){//over trellis with two states
				p_edge->right->m_alpha[0] = p_edge->m_alpha[0] * p_edge->m_v2c[0] + p_edge->m_alpha[1] * p_edge->m_v2c[1];
				p_edge->right->m_alpha[1] = p_edge->m_alpha[0] * p_edge->m_v2c[1] + p_edge->m_alpha[1] * p_edge->m_v2c[0];
				temp_sum = p_edge->right->m_alpha[0] + p_edge->right->m_alpha[1];
				p_edge->right->m_alpha[0] /= temp_sum;
				p_edge->right->m_alpha[1] /= temp_sum;

				p_edge = p_edge->right;
			}
			//backward
			p_edge = (m_row_head + i)->left;
			p_edge->m_beta[0] = 1.0;
			p_edge->m_beta[1] = 0.0;
			while (p_edge->m_col_no != -1){
				temp0 = p_edge->m_alpha[0] * p_edge->m_beta[0] + p_edge->m_alpha[1] * p_edge->m_beta[1];
				temp1 = p_edge->m_alpha[0] * p_edge->m_beta[1] + p_edge->m_alpha[1] * p_edge->m_beta[0];
				temp_sum = temp0 + temp1;
				p_edge->m_c2v[0] = temp0 / temp_sum;
				p_edge->m_c2v[1] = temp1 / temp_sum;

				if (p_edge->m_c2v[0] > 1.0 - SMALLPROB)
					p_edge->m_c2v[0] = 1.0 - SMALLPROB;
				if (p_edge->m_c2v[0] < SMALLPROB)
					p_edge->m_c2v[0] = SMALLPROB;

				p_edge->m_c2v[1] = 1.0 - p_edge->m_c2v[0];

				p_edge->left->m_beta[0] = p_edge->m_beta[0] * p_edge->m_v2c[0] + p_edge->m_beta[1] * p_edge->m_v2c[1];
				p_edge->left->m_beta[1] = p_edge->m_beta[0] * p_edge->m_v2c[1] + p_edge->m_beta[1] * p_edge->m_v2c[0];

				temp_sum = p_edge->left->m_beta[0] + p_edge->left->m_beta[1];
				p_edge->left->m_beta[0] /= temp_sum;
				p_edge->left->m_beta[1] /= temp_sum;

				p_edge = p_edge->left;
			}
		}
	}

//compute soft message
	if (cc_hat_soft!=nullptr)
	{
		for (i = 0; i < m_num_col; i++){
			//forward
			p_edge = (m_col_head+i)->down;
			p_edge->m_alpha[0] = 0.5;//����������Ϣ
			p_edge->m_alpha[1] = 0.5;

			while (p_edge->m_row_no != -1){
				p_edge->down->m_alpha[0] = p_edge->m_alpha[0] * p_edge->m_c2v[0];
				p_edge->down->m_alpha[1] = p_edge->m_alpha[1] * p_edge->m_c2v[1];
				temp_sum = p_edge->down->m_alpha[0] + p_edge->down->m_alpha[1];
				p_edge->down->m_alpha[0] /= temp_sum;
				p_edge->down->m_alpha[1] /= temp_sum;
				p_edge = p_edge->down;
			}
			//soft messaage
			cc_hat_soft[i] = (m_col_head+i)->m_alpha[0];
		}
	}

	m_detection_error = 1-m_success;//added by liangchulong @2013-01-23-19-04
	return iter + (iter<m_max_iter);
}

/***************************************************************************
������SoftInHardOutDecoder(double *U2V, double *M2V, int *uu_hat)
���ܣ�
���������*U2V    ---- �������
          *M2V    ---- ������Ϣ�������м�����غ������
		  *uu_soft_out ���������Ϣ
���������*uu_hat ---- ������   
****************************************************************************/
//added by liang chu long
int CBinaryLDPCCodec::SoftInSoftOut(double *U2V, double *M2V, double* uu_soft_out, double* cc_soft_out, char *arraw)
{
	int i, j;
	int iter;
	int parity_check;
	double temp0, temp1, temp_sum;

	Edge *p_edge;

//initialization
	InitMsg();

	for (iter = 0; iter < m_max_iter; iter++){

//from vnode to cnode
		j = 0;
		for (i = 0; i < m_num_col; i++){
//forward
			p_edge = (m_col_head+i)->down;
			if (i < m_codechk){
				p_edge->m_alpha[0] = M2V[i];
				p_edge->m_alpha[1] = (1.0 - M2V[i]);
			}			
			else{
				temp0 = U2V[j] * M2V[i];
				temp1 = (1.0 - U2V[j]) * (1.0 - M2V[i]);
				temp_sum = temp0 + temp1;
				p_edge->m_alpha[0] = temp0 / temp_sum;
				p_edge->m_alpha[1] = temp1 / temp_sum;
				j++;
			}

			while (p_edge->m_row_no != -1){
				if (p_edge->m_row_no < m_num_row){
					p_edge->down->m_alpha[0] = p_edge->m_alpha[0] * p_edge->m_c2v[0];
					p_edge->down->m_alpha[1] = p_edge->m_alpha[1] * p_edge->m_c2v[1];
					temp_sum = p_edge->down->m_alpha[0] + p_edge->down->m_alpha[1];
					p_edge->down->m_alpha[0] /= temp_sum;
					p_edge->down->m_alpha[1] /= temp_sum;
				}
				else{
					p_edge->down->m_alpha[0] = p_edge->m_alpha[0];
					p_edge->down->m_alpha[1] = p_edge->m_alpha[1];
				}
				p_edge = p_edge->down;
			}
//hard decision
			if ((m_col_head+i)->m_alpha[0] > (m_col_head+i)->m_alpha[1])
				m_cc_hat[i] = 0;
			else
				m_cc_hat[i] = 1;

//backward
			p_edge = (m_col_head+i)->up;
			p_edge->m_beta[0] = 1.0;
			p_edge->m_beta[1] = 1.0;
			while (p_edge->m_row_no != -1){
				if (p_edge->m_row_no < m_num_row){
					temp0 = p_edge->m_alpha[0] * p_edge->m_beta[0];
					temp1 = p_edge->m_alpha[1] * p_edge->m_beta[1];
					temp_sum = temp0 + temp1;
					p_edge->m_v2c[0] = temp0 / temp_sum;
					p_edge->m_v2c[1] = temp1 / temp_sum;
				
					p_edge->up->m_beta[0] = p_edge->m_beta[0] * p_edge->m_c2v[0];
					p_edge->up->m_beta[1] = p_edge->m_beta[1] * p_edge->m_c2v[1];

					temp_sum = p_edge->up->m_beta[0] + p_edge->up->m_beta[1];
					p_edge->up->m_beta[0] /= temp_sum;
					p_edge->up->m_beta[1] /= temp_sum;
				}
				else{
					p_edge->up->m_beta[0] = p_edge->m_beta[0];
					p_edge->up->m_beta[1] = p_edge->m_beta[1];
				}
				p_edge = p_edge->up;
			}
		}

//parity checking
		m_success = 1;
		for (i = 0; i < m_num_row; i++){
			parity_check = 0;
			p_edge = (m_row_head+i)->right;
			while (p_edge->m_col_no != -1){
				if (p_edge->m_col_no < m_num_col)
					parity_check = parity_check ^ m_cc_hat[p_edge->m_col_no];
				p_edge = p_edge->right;
			}
			if (parity_check != 0){
				m_success = 0;
				break;
			}
		}

		if (m_success == 1)	
			break;			

//from c node to v node
		for (i = 0; i < m_num_row; i++){
//forward
			p_edge = (m_row_head + i)->right;
			p_edge->m_alpha[0] = 1.0;
			p_edge->m_alpha[1] = 0.0;
			while (p_edge->m_col_no != -1){//over trellis with two states
				p_edge->right->m_alpha[0] = p_edge->m_alpha[0] * p_edge->m_v2c[0] + p_edge->m_alpha[1] * p_edge->m_v2c[1];
				p_edge->right->m_alpha[1] = p_edge->m_alpha[0] * p_edge->m_v2c[1] + p_edge->m_alpha[1] * p_edge->m_v2c[0];
				temp_sum = p_edge->right->m_alpha[0] + p_edge->right->m_alpha[1];
				p_edge->right->m_alpha[0] /= temp_sum;
				p_edge->right->m_alpha[1] /= temp_sum;
	
				p_edge = p_edge->right;
			}
//backward
			p_edge = (m_row_head + i)->left;
			p_edge->m_beta[0] = 1.0;
			p_edge->m_beta[1] = 0.0;
			while (p_edge->m_col_no != -1){
				temp0 = p_edge->m_alpha[0] * p_edge->m_beta[0] + p_edge->m_alpha[1] * p_edge->m_beta[1];
				temp1 = p_edge->m_alpha[0] * p_edge->m_beta[1] + p_edge->m_alpha[1] * p_edge->m_beta[0];
				temp_sum = temp0 + temp1;
				p_edge->m_c2v[0] = temp0 / temp_sum;
				p_edge->m_c2v[1] = temp1 / temp_sum;
				
				if (p_edge->m_c2v[0] > 1.0 - SMALLPROB)
					p_edge->m_c2v[0] = 1.0 - SMALLPROB;
				if (p_edge->m_c2v[0] < SMALLPROB)
					p_edge->m_c2v[0] = SMALLPROB;

				p_edge->m_c2v[1] = 1.0 - p_edge->m_c2v[0];

				p_edge->left->m_beta[0] = p_edge->m_beta[0] * p_edge->m_v2c[0] + p_edge->m_beta[1] * p_edge->m_v2c[1];
				p_edge->left->m_beta[1] = p_edge->m_beta[0] * p_edge->m_v2c[1] + p_edge->m_beta[1] * p_edge->m_v2c[0];
				
				temp_sum = p_edge->left->m_beta[0] + p_edge->left->m_beta[1];
				p_edge->left->m_beta[0] /= temp_sum;
				p_edge->left->m_beta[1] /= temp_sum;
				
				p_edge = p_edge->left;
			}
		}
	}

	if(uu_soft_out!=nullptr)
	{
		//output msg
		//Kite��
		//for (i = 0; i < m_codedim; i++){
		//	//forward
		//	p_edge = (m_col_head+i)->down;
		//	//if (i < m_codedim){
		//	//	temp0 = U2V[i] * M2V[i];
		//	//	temp1 = (1.0 - U2V[i]) * (1.0 - M2V[i]);
		//	//	temp_sum = temp0 + temp1;
		//	//	p_edge->m_alpha[0] = temp0 / temp_sum;
		//	//	p_edge->m_alpha[1] = temp1 / temp_sum;
		//	//}			
		//	//else{
		//		p_edge->m_alpha[0] = M2V[i];//�����������Ϣ
		//		p_edge->m_alpha[1] = 1.0 - M2V[i];
		//	//}

		//	while (p_edge->m_row_no != -1){
		//		if (p_edge->m_row_no < m_num_row){
		//			p_edge->down->m_alpha[0] = p_edge->m_alpha[0] * p_edge->m_c2v[0];
		//			p_edge->down->m_alpha[1] = p_edge->m_alpha[1] * p_edge->m_c2v[1];
		//			temp_sum = p_edge->down->m_alpha[0] + p_edge->down->m_alpha[1];
		//			p_edge->down->m_alpha[0] /= temp_sum;
		//			p_edge->down->m_alpha[1] /= temp_sum;
		//		}
		//		else{
		//			p_edge->down->m_alpha[0] = p_edge->m_alpha[0];
		//			p_edge->down->m_alpha[1] = p_edge->m_alpha[1];
		//		}
		//		p_edge = p_edge->down;
		//	}
		//	//soft msg
		//	uu_soft_out[i] = (m_col_head+i)->m_alpha[0];
		//}
		//Kite��

		for (j = 0, i = m_codechk; i < m_num_col; i++, j++){
			//forward
			p_edge = (m_col_head+i)->down;
			p_edge->m_alpha[0] = M2V[i];//�����������Ϣ
			p_edge->m_alpha[1] = 1.0 - M2V[i];

			while (p_edge->m_row_no != -1){
				if (p_edge->m_row_no < m_num_row){
					p_edge->down->m_alpha[0] = p_edge->m_alpha[0] * p_edge->m_c2v[0];
					p_edge->down->m_alpha[1] = p_edge->m_alpha[1] * p_edge->m_c2v[1];
					temp_sum = p_edge->down->m_alpha[0] + p_edge->down->m_alpha[1];
					p_edge->down->m_alpha[0] /= temp_sum;
					p_edge->down->m_alpha[1] /= temp_sum;
				}
				else{
					p_edge->down->m_alpha[0] = p_edge->m_alpha[0];
					p_edge->down->m_alpha[1] = p_edge->m_alpha[1];
				}
				p_edge = p_edge->down;
			}
			//soft msg
			uu_soft_out[j] = (m_col_head+i)->m_alpha[0];
		}
	}

	if (cc_soft_out!=nullptr)
	{
		//output msg
		j = 0;
		for (i = 0; i < m_num_col; i++){
			//forward
			p_edge = (m_col_head+i)->down;
			if (i < m_codechk){// ֮ǰ��i<m_codedim, 2013����������2014��7��20�����֣�liangchulong����2014��7��20����
				p_edge->m_alpha[0] = 0.5;//�����������Ϣ
				p_edge->m_alpha[1] = 0.5;
			}
			else{
				p_edge->m_alpha[0] = U2V[j];//�����������Ϣ
				p_edge->m_alpha[1] = 1.0 - U2V[j];
				j++;
			}

			while (p_edge->m_row_no != -1){
				if (p_edge->m_row_no < m_num_row){
					p_edge->down->m_alpha[0] = p_edge->m_alpha[0] * p_edge->m_c2v[0];
					p_edge->down->m_alpha[1] = p_edge->m_alpha[1] * p_edge->m_c2v[1];
					temp_sum = p_edge->down->m_alpha[0] + p_edge->down->m_alpha[1];
					p_edge->down->m_alpha[0] /= temp_sum;
					p_edge->down->m_alpha[1] /= temp_sum;
				}
				else{
					p_edge->down->m_alpha[0] = p_edge->m_alpha[0];
					p_edge->down->m_alpha[1] = p_edge->m_alpha[1];
				}
				p_edge = p_edge->down;
			}
			//soft msg
			cc_soft_out[i] = (m_col_head+i)->m_alpha[0];
		}
	}

	m_detection_error = 1-m_success;//added by liangchulong @2013-01-23-19-04
	return iter + (iter<m_max_iter);
}
