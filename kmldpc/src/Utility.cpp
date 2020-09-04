/*******************************************************************************
Programmed by Xiao Ma, Sun Yat-sen University. 
If you have any suggestions, please contact me at maxiao@mail.sysu.edu.cn
The program can only be employed for academic research.
******************************************************************************/

#include "Utility.h"
#include "RandNum.h"

extern CLCRandNum rndGen0;

void MatrixProd(int* uu, int* cc, int** G, int dim, int len) {
	for (int n = 0; n < len; ++n) {
		cc[n] = 0;
	}
	for (int k = 0; k < dim; ++k) {
		for (int n = 0; n < len; ++n) {
			cc[n] ^= (uu[k] & G[k][n]);
		}
	}
}

double Seqmax(double *x, int num)
{
	int i;
	double temp;

	temp = x[0];
	for (i = 1; i < num; i++)
	{
		if (x[i] > temp)
			temp = x[i];
	}
	return temp;
}

double Seqmin(double *x, int num)
{
	int i;
	double temp;

	temp = x[0];
	for (i = 1; i < num; i++)
	{
		if (x[i] < temp)
			temp = x[i];
	}
	return temp;
}

double FunctionQ(double x)
{
	double t,dq;
	double q;
	if(x >= 8.5)
	{
		q = -log(sqrt( 2 * m_PI) * x) - x * x / 2;
	}
	else
	{
		if(x < -7)
			return(0.0);

		t = 1.0 / (1.0 + 0.2316419 * fabs(x));

		dq = ((((1.330274429 * t - 1.821255978) * t + 1.78147793) * t - 0.356563782) * t
			+ 0.319381530) * t * exp(-x*x*0.5) * 0.39894228;
		
		if(x > 0.0)
			q = log(dq);
		else
			q = log(1.0-dq);

	}

	return(q);
}

int sgn(double x)
{
	if (x > 0)
		return 1;
	else if (x < 0)
		return -1;
	else
		return 0;
}

double erf_inv(double x)
{
	double y, a;

	a = 0.147;

	y = sgn(x) * sqrt(sqrt(pow((2/(a*m_PI) + 0.5*log(1-x*x)), 2) - (log(1-x*x)/a)) 
		             - (2./(a*m_PI) + 0.5*log(1-x*x)));

	return y;
}

void ProbClip(double *xx, int len_xx)
{
	int i;

	for (i = 0; i < len_xx; i++){
		if (xx[i] < SMALLPROB)
			xx[i] = SMALLPROB;
		else if (xx[i] > 1.0 - SMALLPROB)
			xx[i] = 1.0 - SMALLPROB;
	}

	return;
}


/**************************************************************************************
         binary (b[0] b[1] b[2] ......), where b[0] is the MSB, b[end-1] is the LSB
         in other words, b is the binary expression of dec number d;
                     d     --> b = (b[0] b[1] b[2])  len_b = 3   
         for example d = 0 --> b = (0 0 0)
                     d = 1 --> b = (0 0 1)
                     d = 2 --> b = (0 1 0)   
                     d = 3 --> b = (0 1 1)
                     d = 4 --> b = (1 0 0)   
                     d = 5 --> b = (1 0 1)
                     d = 6 --> b = (1 1 0)   
                     d = 7 --> b = (1 1 1)
                     d = 8 --> b = (0 0 0)   
**************************************************************************************/

void Dec2Bin(int d, int *b, int len_b)
{
	int i;

	for (i = 0; i < len_b; i++)
		b[len_b-i-1] = (d >> i) % 2;

	return;
}

void SeqDec2Bin(int *bin, int *dec, int len_dec, int len_symbol)
{
	int i, j, t;

	t = 0;
	for (i = 0; i < len_dec; i++){
		for (j = 0; j < len_symbol; j++){
			bin[t] = (dec[i] >> (len_symbol-1-j)) & 1; 
			t++;
		}
	}

	return;
}

void SeqBin2Dec(int *bin, int *dec, int len_dec, int len_symbol)
{
	int i, j, t;

	t = 0;
	for (i = 0; i < len_dec; i++){
		dec[i] = 0;
		for (j = 0; j < len_symbol; j++){
			dec[i] = (dec[i] << 1) ^ bin[t];
			t++;
		}
	}

	return;
}

int BitDotProd(int a, int b, int len)
{
	int i; 
	int temp;
	int prod;

	temp = a & b;

	prod = 0;
	for (i = 0; i < len; i++)
		prod += (temp >> i) % 2;
	
	return (prod % 2);
}

int Systemizer(int num_row, int num_col, int **H, int **sysH, int *pai)
{
	int i, j, ii, jj, m, n;
	int flag;
	int temp;
	int rank;

	for (j = 0; j < num_col; j++)
		pai[j] = j;//pai[0] means that the first column (0th column) of sysH comes from the pai[0]-th of H

//copy H to sysH;
	for (i = 0; i < num_row; i++){
		for (j = 0; j < num_col; j++)
			sysH[i][j] = H[i][j];
	}

//initialize the rank
	rank = 0;
//find the most left_up non_zero element
	for (i = 0; i < num_row; i++){
		flag = 0;
		for (jj = i; jj < num_col; jj++){
			for (ii = i; ii < num_row; ii++){
				if (sysH[ii][jj] != 0){
					flag = 1;
					break;
				}
			}
			if (flag == 1){
				rank++;
				break;
			}
		}

		if (flag == 0)
			break;
		else{//swap to ensure the (1,1) elment is non_zero
//swap i and ii row
			if (ii != i){
				for (n = 0; n < num_col; n++){
					temp = sysH[i][n];
					sysH[i][n] = sysH[ii][n];
					sysH[ii][n] = temp;
				}
			}
//swap i and jj col
			if (jj != i){
				temp = pai[i];
				pai[i] = pai[jj];
				pai[jj] = temp;

				for (m = 0; m < num_row; m++){
					temp = sysH[m][i];
					sysH[m][i] = sysH[m][jj];
					sysH[m][jj] = temp;
				}
			}
//elimination
			for (m = 0; m < num_row; m++){
				if (m != i && sysH[m][i] == 1){
					for (n = 0; n < num_col; n++)
						sysH[m][n] ^= sysH[i][n];
				}
			}
		}
	}

	return rank;
}

int min(int x, int y)
{
	return x < y? x : y;
}

int max(int x, int y)
{
	return x > y? x : y;
}

//interleavers
void RandomIntl(int *pai, int period)
{
	int t;
	double random_num;
	int temp;
	int position;

	for (t = 0; t < period; t++)
		pai[t] = t;

	for (t = period - 1; t > 0; t--){
		random_num = rndGen0.Uniform();
		position = (int) (random_num * t);
		temp = pai[position];
		pai[position] = pai[t];
		pai[t] = temp;
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void N1N2RandomIntl(int *pai, int period, int N1, int N2)
{
	int i, j, k;
	int position;
	int count;
	int *Zn;
	int *tail;
	int *temp_tail;
	int *temp_pai;
	int body_len;
	int tail_len;
	int valid_len;
	double rand_number;
	int success;
  
	success = 0;
  
	Zn = new int[period];
	for (i = 0; i < period; i++)
		Zn[i] = i;
  
	valid_len = period;//number of element that can be selected for the current position
	for (i = 0; i < period; i++){
		rand_number = 1.0;
		while (rand_number >= 1.0){
			rand_number = rndGen0.Uniform();
		}
		position = (int)(valid_len * rand_number);
		count = -1;
		for (j = 0; j < period; j++){
			if (Zn[j] >= 0)
				count++;
			if (count == position)
				break;
		}
		pai[i] = j;//the i-th position of pai is filled by j
		Zn[j] = -1;//the j-th position is set to be -1 to mean that j cannot be selected again 
		valid_len--;
    
//in order to generate pai[i+1], we need tell which number is OK
//if Zn[j] = -1, j is not OK
//integers in interval [pai[i]-N1+1, pai[i]+N1-1] are not OK, set to be -2 
//because |(i+1)-i| < N2 requires |pai[i+1]-pai[i]|>=N1
//integers in interval [pai[i-1]-N1+1, pai[i-1]+N1-1] are not OK, set to be -3
//integers in interval [pai[i-2]-N1+1, pai[i-2]+N1-1] are not OK, set to be -4
//...
//integers in interval [pai[i+2-N2]-N1+1, pai[i+2-N2]+N1-1] are not OK, set to be -N2
//other integers which is not -1 are OK
		
		for (j = 0; j < period; j++){
			if (Zn[j] == - (N2 + 1)){
				Zn[j] = j;
				valid_len++;
			}
			else if(Zn[j] < -1)
				Zn[j]--;
		}

		for (j = max(0, pai[i]-N1+1); j <= min(pai[i]+N1-1, period-1); j++){
			if (Zn[j] >= 0){
				Zn[j] = -2;
				valid_len--;
			}
		}
		if (valid_len <= 0)
			break;
	}

	tail_len = 0;
	for (i = 0; i < period; i++){
		if (Zn[i] < -1)
			tail_len++;
	}
	body_len = period - tail_len;

	if (tail_len == 0){
		fprintf(stdout, "\nS_RandomIntl is built successfully\n");
		success = 1;
	}
	else{
		temp_tail = new int[tail_len];
		temp_pai = new int[tail_len];
		tail = new int[tail_len];
		RandomIntl(temp_pai, tail_len);
		
		j = 0;
		for (i = 0; i < period; i++){
			if (Zn[i] < -1){
				temp_tail[j] = i;
				j++;
			}
		}
		for (i = 0; i < tail_len; i++)
			tail[i] = temp_tail[temp_pai[i]];
    
//integers which have not been used are written in tail

		for (j = 0; j < tail_len; j++){
			//try to insert tail[j] between pai[i-1] and pai[i]
			//check pai[k] corresponding to the intervals [i-N2+1, i-1]
			//[i, i+N2-2] would change to [i+1, i+N2-1]
			//|i-k| < N2 requires that |tail[j]-pai[k]| >= N1
			for (i = 0; i < body_len; i++){
				count = 0;
				for (k = max(0, i-N2+1); k <= min(body_len, i+N2-2); k++){
					if (abs(tail[j] - pai[k]) < N1){
						count++;
						break;
					}
				}
				if (count == 0){
					for (k = period - 1; k > i; k--)
						pai[k] = pai[k - 1];
						pai[i] = tail[j];
						body_len++;
						break;
				}
			}
		}

		if (body_len == period){
			fprintf(stdout, "\nS_RandomIntl is built successfully by inserting tail\n");
			success = 1;
		}
		delete []temp_tail;
		delete []temp_pai;
		delete []tail;
	}

	delete []Zn;
  
	if (success == 0)
		N1N2RandomIntl(pai, period, N1, N2);

	return;
}

void ShiftIntl( int *pai, int shifted_pos, int period )
{
	for (int i=0; i<period; i++)
	{
		pai[i] = (i+shifted_pos)%period;
	}
}

void LLRClip( double *xx, int len_xx )
{
	int i;

	for (i = 0; i < len_xx; i++){
		if (xx[i] < SMALL_LLR)
			xx[i] = SMALL_LLR;
		else if (xx[i] > LARGE_LLR)
			xx[i] = LARGE_LLR;
	}

	return;
}

double maxstar( double x, double y )
{
	if (x>y)
		return x + log( 1.0 + exp(y-x) );
	else
		return y + log( 1.0 + exp(x-y) );
}

void DFT_gen(int size_t, double** F_re, double** F_im)
{
	int i, j;
	for (i = 0; i < size_t; i++)
		for (j = 0; j < size_t; j++)
		{
			F_re[i][j] = 1 / sqrt(size_t*1.0) * cos(2. * m_PI * i * j / size_t);
			F_im[i][j] = -1 / sqrt(size_t*1.0) * sin(2. * m_PI * i * j / size_t);
		}

	return;
}