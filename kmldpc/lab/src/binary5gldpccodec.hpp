#ifndef LAB_BINARY_5GLDPC_CODEC_HPP
#define LAB_BINARY_5GLDPC_CODEC_HPP

#include <cstdio>
#include <cstring>

#include "utility.hpp"

namespace lab {

/************************************************************
类名：CBinaryLDPCCodec
功能 ：二进制LDPC码编译码
*************************************************************/

class Binary5GLDPCCodec
{
public:
    void Free() {
        Free_Tanner_Graph();

        if (m_encoder_active != 0) {
            for (int i = 0; i < m_num_row; i++) {
                delete []m_encH[i];
            }
            delete []m_encH;
        }

        delete []m_cc_hat;
        delete []m_syndromsoft;
    }

    void Malloc(int code_no, char *file_name) {
        int i, j;
        int row_no, row_deg, col_no;
        Edge *temp_edge;
        FILE *fq = nullptr;

        char temp_str[80] = {' '};
        char mark[80];
        FILE *fp;

        sprintf(mark, "LDPC***%d***PARAMETERS",code_no);

        if ((fp = fopen(file_name, "r")) == nullptr){
            fprintf(stderr, "\nCannot Open %s", file_name);
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
            fprintf(stderr, "\nCannot Open %s", m_file_name_of_H);
            exit(0);
        }

        fscanf(fp, "%s", temp_str);
        fscanf(fp, "%d %d %d %d", &m_num_row, &m_num_col, &m_codechk, &m_lifting_factor);
        m_codelen_no_puncture = m_num_col;
        m_codelen_puncture = m_num_col - m_lifting_factor * 2;
        m_codedim = m_codelen_no_puncture - m_codechk;
        m_coderate = (double)m_codedim / m_codelen_puncture;

        m_cc_nopuncture = new int[m_codelen_no_puncture];

        m_cc_soft_nopuncture = new double[m_codelen_no_puncture];

        m_row_head = new Edge[m_num_row];
        m_col_head = new Edge[m_num_col];

        m_syndromsoft = new double[m_num_row];

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
            SystH_5G();

        m_cc_hat = new int[m_codelen_no_puncture];
    }

    void Encoder(int *uu, int *cc) {
        int i, j, t;

//codeword = [parity_check_bits information_bits]
        switch (m_encoder_active){
            case 0:
                for (i = 0; i < m_codedim; i++)
                    uu[i] = 0;
                for (i = 0; i < m_codelen_no_puncture; i++)
                    cc[i] = 0;
                break;
            case 1:
                for (t = m_codechk; t < m_codelen_no_puncture; t++)
                    cc[t] = uu[t- m_codechk];

                for (t = 0; t < m_codechk; t++){
                    cc[t] = 0;
                    for (j = m_codechk; j < m_codelen_no_puncture; j++)
                        cc[t] ^= (cc[j] & m_encH[t][j]);
                }
                break;
            default:
                break;
        }
    }

    void Encoder_5G(int *uu, int *cc) {
        int i, j, t;

        for (i = 0; i < m_codelen_puncture; i++)
            m_cc_nopuncture[i] = 0;

        //codeword = [parity_check_bits information_bits]
        switch (m_encoder_active) {
            case 0:
                for (i = 0; i < m_codedim; i++)
                    uu[i] = 0;
                for (i = 0; i < m_codelen_puncture; i++)
                    cc[i] = 0;
                break;
            case 1:
                for (t = 0; t < m_codedim; t++)
                    m_cc_nopuncture[t] = uu[t];
                for (t = 0; t < m_codechk; t++) {
                    m_cc_nopuncture[m_codedim + t] = 0;
                    for (j = 0; j < m_codedim; j++)
                        m_cc_nopuncture[m_codedim + t] ^= (m_cc_nopuncture[j] & m_encH[t][j]);
                }
                for (t = 0; t < m_codelen_puncture; t++) {
                    cc[t] = m_cc_nopuncture[t + m_lifting_factor * 2];
                }
                break;
            default:
                break;
        }
    }

    int Decoder(double *M2V, int *uu_hat) {
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
//			uu_hat[i] = cc_hat_[i];//kite使用
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

                    if (p_edge->m_c2v[0] > 1.0 - kSmallestProb)
                        p_edge->m_c2v[0] = 1.0 - kSmallestProb;
                    if (p_edge->m_c2v[0] < kSmallestProb)
                        p_edge->m_c2v[0] = kSmallestProb;

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

        return iter + (iter<m_max_iter);
    }

    int Decoder_5G(double *M2V, int *uu_hat, int iter_cnt) {
        int i;
        int iter;
        int parity_check;
        double temp0, temp1, temp_sum;

        Edge *p_edge;

        //initialization
        InitMsg();
        //iteration
        for (iter = 0; iter < iter_cnt; iter++) {
            //from vnode to cnode
            for (i = 0; i < m_num_col; i++) {
                //forward
                if (i < m_lifting_factor * 2) {
                    p_edge = (m_col_head + i)->down;
                    p_edge->m_alpha[0] = 0.5;
                    p_edge->m_alpha[1] = 1.0 - 0.5;
                }
                else {
                    p_edge = (m_col_head + i)->down;
                    p_edge->m_alpha[0] = M2V[i - m_lifting_factor * 2];
                    p_edge->m_alpha[1] = 1.0 - M2V[i - m_lifting_factor * 2];
                }

                while (p_edge->m_row_no != -1) {
                    p_edge->down->m_alpha[0] = p_edge->m_alpha[0] * p_edge->m_c2v[0];
                    p_edge->down->m_alpha[1] = p_edge->m_alpha[1] * p_edge->m_c2v[1];
                    temp_sum = p_edge->down->m_alpha[0] + p_edge->down->m_alpha[1];
                    p_edge->down->m_alpha[0] /= temp_sum;
                    p_edge->down->m_alpha[1] /= temp_sum;
                    p_edge = p_edge->down;
                }
                //hard decision
                if ((m_col_head + i)->m_alpha[0] > (m_col_head + i)->m_alpha[1])
                    m_cc_hat[i] = 0;
                else
                    m_cc_hat[i] = 1;

                //backward
                p_edge = (m_col_head + i)->up;
                p_edge->m_beta[0] = 1.0;
                p_edge->m_beta[1] = 1.0;
                while (p_edge->m_row_no != -1) {
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
                //			uu_hat[i] = cc_hat_[i+m_codechk];
                uu_hat[i] = m_cc_hat[i];//kite使用
            }
            //parity checking
            m_success = 1;
            for (i = 0; i < m_num_row; i++) {
                parity_check = 0;
                p_edge = (m_row_head + i)->right;
                while (p_edge->m_col_no != -1) {
                    parity_check = parity_check ^ m_cc_hat[p_edge->m_col_no];
                    p_edge = p_edge->right;
                }
                if (parity_check != 0) {
                    m_success = 0;
                    break;
                }
            }

            if (m_success == 1)
                break;

            //from c node to v node
            for (i = 0; i < m_num_row; i++) {
                //forward
                p_edge = (m_row_head + i)->right;
                p_edge->m_alpha[0] = 1.0;
                p_edge->m_alpha[1] = 0.0;
                while (p_edge->m_col_no != -1) {//over trellis with two states
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
                while (p_edge->m_col_no != -1) {
                    temp0 = p_edge->m_alpha[0] * p_edge->m_beta[0] + p_edge->m_alpha[1] * p_edge->m_beta[1];
                    temp1 = p_edge->m_alpha[0] * p_edge->m_beta[1] + p_edge->m_alpha[1] * p_edge->m_beta[0];
                    temp_sum = temp0 + temp1;
                    p_edge->m_c2v[0] = temp0 / temp_sum;
                    p_edge->m_c2v[1] = temp1 / temp_sum;

                    if (p_edge->m_c2v[0] > 1.0 - kSmallestProb)
                        p_edge->m_c2v[0] = 1.0 - kSmallestProb;
                    if (p_edge->m_c2v[0] < kSmallestProb)
                        p_edge->m_c2v[0] = kSmallestProb;

                    p_edge->m_c2v[1] = 1.0 - p_edge->m_c2v[0];

                    p_edge->left->m_beta[0] = p_edge->m_beta[0] * p_edge->m_v2c[0] + p_edge->m_beta[1] * p_edge->m_v2c[1];
                    p_edge->left->m_beta[1] = p_edge->m_beta[0] * p_edge->m_v2c[1] + p_edge->m_beta[1] * p_edge->m_v2c[0];

                    temp_sum = p_edge->left->m_beta[0] + p_edge->left->m_beta[1];
                    p_edge->left->m_beta[0] /= temp_sum;
                    p_edge->left->m_beta[1] /= temp_sum;

                    p_edge = p_edge->left;
                }

                m_syndromsoft[i] = (m_row_head + i)->m_alpha[0];
            }
        }

        return iter + (iter<m_max_iter);
    }

    int Parity_Checking() {
        int i;
        int iter;
        int parity_check;
        double temp0, temp1, temp_sum;

        Edge *p_edge;

        //parity checking
        m_success = 1;
        for (i = 0; i < m_num_row; i++) {
            parity_check = 0;
            p_edge = (m_row_head + i)->right;
            while (p_edge->m_col_no != -1) {
                parity_check = parity_check ^ m_cc_hat[p_edge->m_col_no];
                p_edge = p_edge->right;
            }
            if (parity_check != 0) {
                m_success = 0;
                break;
            }
        }

        return m_success;
    }

    int parityCheck(int *cc) {
        int i;
        int parity_check;
        int count = 0;

        Edge *p_edge;

        //parity checking
        m_success = 1;
        for (i = 0; i < m_num_row; i++) {
            parity_check = 0;
            p_edge = (m_row_head + i)->right;
            while (p_edge->m_col_no != -1) {
                parity_check = parity_check ^ cc[p_edge->m_col_no];
                p_edge = p_edge->right;
            }
            if (parity_check != 0) {
                count++;
            }
        }

        return count;
    }

    void InitMsg() {
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
    }

private:
    void SystH() {
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

/* By Zhaosc
-----由于在Gauss消元过程中有行列变换
-----需要删除原有的Tanner Graph，重新建立.
*/
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

        m_codelen_no_puncture = m_num_col;
        m_codedim = m_codelen_no_puncture - m_codechk;
        m_coderate = (double)m_codedim / m_codelen_no_puncture;

        delete []tempP;
        for (i = 0; i < m_num_row; i++)
        {
            delete[]tempH[i];
            delete[]m_decH[i];
        }
        delete []tempH;
        delete []m_decH;
    }
    /***************************************************************************
    函数：SystH_1()
    功能：把5G LDPC 校验矩阵化为系统形式 单位阵在后面 系统位在前面
    ****************************************************************************/
    void SystH_1() {
        int i, j, ii, jj, m, n;
        int temp;
        int flag;
        int *tempP;
        char **tempH;
        Edge *p_edge;

        m_decH = new char*[m_num_row];
        for (i = 0; i < m_num_row; i++)
            m_decH[i] = new char[m_num_col];

        for (i = 0; i < m_num_row; i++) {
            for (j = 0; j < m_num_col; j++)
                m_decH[i][j] = 0;
        }

        for (i = 0; i < m_num_row; i++) {
            p_edge = (m_row_head + i)->right;
            while (p_edge->m_col_no != -1) {
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

        for (i = 0; i < m_num_row; i++) {
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

        for (i = m_num_row - 1; i >= 0; --i)
        {
            flag = 0;
            for (jj = i + m_num_col - m_num_row; jj >= 0; --jj)
            {
                for (ii = i; ii >= 0; --ii)
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
                if (jj != i + m_num_col - m_num_row)
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
                for (m = m_num_row - 1; m >= 0; --m)
                {
                    if (m != i && m_encH[m][i + m_num_col - m_num_row] == 1)
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

        /* By Zhaosc
        -----由于在Gauss消元过程中有行列变换
        -----需要删除原有的Tanner Graph，重新建立.
        */
        Free_Tanner_Graph();

        m_row_head = new Edge[m_num_row];
        m_col_head = new Edge[m_num_col];

        for (i = 0; i < m_num_row; i++)
        {
            (m_row_head + i)->m_row_no = i;
            (m_row_head + i)->m_col_no = -1;
            (m_row_head + i)->left = m_row_head + i;
            (m_row_head + i)->right = m_row_head + i;
            (m_row_head + i)->up = m_row_head + i;
            (m_row_head + i)->down = m_row_head + i;
        }

        for (i = 0; i < m_num_col; i++)
        {
            (m_col_head + i)->m_row_no = -1;
            (m_col_head + i)->m_col_no = i;
            (m_col_head + i)->left = m_col_head + i;
            (m_col_head + i)->right = m_col_head + i;
            (m_col_head + i)->up = m_col_head + i;
            (m_col_head + i)->down = m_col_head + i;
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

        m_codelen_no_puncture = m_num_col;
        m_codedim = m_codelen_no_puncture - m_codechk;
        m_coderate = (double)m_codedim / m_codelen_puncture;

        delete[]tempP;
        for (i = 0; i < m_num_row; i++)
        {
            delete[]tempH[i];
            delete[]m_decH[i];
        }
        delete[]tempH;
        delete[]m_decH;
    }
    void SystH2() {
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

        /* By Zhaosc
        -----由于在Gauss消元过程中有行列变换
        -----需要删除原有的Tanner Graph，重新建立.
        */
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

        m_codelen_no_puncture = m_num_col;
        m_codedim = m_codelen_no_puncture - m_codechk;
        m_coderate = (double)m_codedim / m_codelen_no_puncture;

        delete []tempP;
        for (i = 0; i < m_num_row; i++)
        {
            delete[]tempH[i];
            delete[]m_decH[i];
        }
        delete []tempH;
        delete []m_decH;
    }
    void SystH_5G() {
        int i, j, ii, jj, m, n;
        int temp;
        int flag;
        int *tempP;
        char **tempH;
        Edge *p_edge;

        m_decH = new char*[m_num_row];
        for (i = 0; i < m_num_row; i++)
            m_decH[i] = new char[m_num_col];

        for (i = 0; i < m_num_row; i++) {
            for (j = 0; j < m_num_col; j++)
                m_decH[i][j] = 0;
        }

        for (i = 0; i < m_num_row; i++) {
            p_edge = (m_row_head + i)->right;
            while (p_edge->m_col_no != -1) {
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

        for (i = 0; i < m_num_row; i++) {
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

        for (i = m_num_row - 1; i >= 0; --i)
        {
            flag = 0;
            for (jj = i + m_num_col - m_num_row; jj >= 0; --jj)
            {
                for (ii = i; ii >= 0; --ii)
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
                //swap (i + m_num_col - m_num_row) and jj col
                if (jj != i + m_num_col - m_num_row)
                {
                    temp = tempP[i + m_num_col - m_num_row];
                    tempP[i + m_num_col - m_num_row] = tempP[jj];
                    tempP[jj] = temp;

                    for (m = 0; m < m_num_row; m++)
                    {
                        temp = m_encH[m][i + m_num_col - m_num_row];
                        m_encH[m][i + m_num_col - m_num_row] = m_encH[m][jj];
                        m_encH[m][jj] = temp;
                    }
                }
                //elimination
                for (m = m_num_row - 1; m >= 0; --m)
                {
                    if (m != i && m_encH[m][i + m_num_col - m_num_row] == 1)
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

        /* By Zhaosc
        -----由于在Gauss消元过程中有行列变换
        -----需要删除原有的Tanner Graph，重新建立.
        */
        Free_Tanner_Graph();

        m_row_head = new Edge[m_num_row];
        m_col_head = new Edge[m_num_col];

        for (i = 0; i < m_num_row; i++)
        {
            (m_row_head + i)->m_row_no = i;
            (m_row_head + i)->m_col_no = -1;
            (m_row_head + i)->left = m_row_head + i;
            (m_row_head + i)->right = m_row_head + i;
            (m_row_head + i)->up = m_row_head + i;
            (m_row_head + i)->down = m_row_head + i;
        }

        for (i = 0; i < m_num_col; i++)
        {
            (m_col_head + i)->m_row_no = -1;
            (m_col_head + i)->m_col_no = i;
            (m_col_head + i)->left = m_col_head + i;
            (m_col_head + i)->right = m_col_head + i;
            (m_col_head + i)->up = m_col_head + i;
            (m_col_head + i)->down = m_col_head + i;
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

        m_codelen_no_puncture = m_num_col;
        m_codedim = m_codelen_no_puncture - m_codechk;
        m_coderate = (double)m_codedim / m_codelen_puncture;

        delete[]tempP;
        for (i = 0; i < m_num_row; i++)
        {
            delete[]tempH[i];
            delete[]m_decH[i];
        }
        delete[]tempH;
        delete[]m_decH;
    }
    void Free_Tanner_Graph() {
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

public:
    //code parameters
    int m_codedim;//码的维数
    int m_codelen_no_puncture;//码的总长度
    int m_codelen_puncture;  //码的发送长度，打孔之后
    int m_codechk;//校验位的个数
    double m_coderate;//码率
    int m_encoder_active;//是否编码   0：不编码  1：编码
    char m_file_name_of_H[255];

    int m_lifting_factor;

    int *m_cc_nopuncture;
    double *m_cc_soft_nopuncture;
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

    int *m_cc_hat;
    int m_max_iter;
    int m_success;
};

}

#endif //KMLDPC_BINARY5GLDPCCODEC_H
