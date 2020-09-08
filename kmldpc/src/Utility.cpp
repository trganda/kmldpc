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

}

void ShiftIntl( int *pai, int shifted_pos, int period )
{
	for (int i=0; i<period; i++)
	{
		pai[i] = (i+shifted_pos)%period;
	}
}

std::ostream & operator<<( std::ostream & os,const std::complex<double> & c) {
    os << std::fixed << std::setprecision(14)
       << std::resetiosflags(std::ios::showpos)
       << '('
       << std::showpos << c.real() << std::showpos << c.imag() << 'i'
       << ')'
       << std::resetiosflags(std::ios::showpos);
    return os;
}