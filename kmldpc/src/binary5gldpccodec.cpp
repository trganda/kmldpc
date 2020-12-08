#include "binary5gldpccodec.h"
#include "utility.h"
#include "randnum.h"

extern CLCRandNum rndGen0;


/***************************************************************************
函数：ReadHMalloc(int code_no, char *file_name)
功能：初始化，读入基本参数，读入校验矩阵，构建Tanner图
输入参数：code_no     ---- LDPC码编号
          *file_name  ---- 校验矩阵文件名
****************************************************************************/
void CBinary5GLDPCCodec::Malloc(int code_no, char *file_name)
{
    int i, j;
    int row_no, row_deg, col_no;
    EdgeLDPC *temp_edge;
    FILE *fq=NULL;

    char temp_str[80] = {' '};
    char mark[80];
    FILE *fp;

    sprintf(mark, "LDPC***%d***PARAMETERS",code_no);

    if ((fp = fopen(file_name, "r")) == NULL){
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
    if ((fp = fopen(m_file_name_of_H, "r")) == NULL){
        fprintf(stderr, "\nCannot open %s", m_file_name_of_H);
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

    m_row_head = new EdgeLDPC[m_num_row];
    m_col_head = new EdgeLDPC[m_num_col];

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
            temp_edge = new EdgeLDPC;
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

    return;
}

void CBinary5GLDPCCodec::PrintCodeParameter( FILE *fp )
{
    fprintf(fp, "%%#######CBinary5GLDPCCodec Parameters#######\n");
    fprintf(fp, "%%%-20s = [%d, %d, %d, %g];\n", "[n, k, m, r]",
            m_codelen_no_puncture, m_codedim, m_codechk, m_coderate);
    fprintf(fp, "%%%m_max_iter = %d;\n", m_max_iter);
    fprintf(fp, "%%%m_encoder_active = %d;\n", m_encoder_active);
    fprintf(fp, "%%%m_file_name_of_H = '%s';\n", m_file_name_of_H);
}

/***************************************************************************
By Zhaosc
函数：Free_Tanner_Graph()
功能：释放译码过程中用到的Tanner图的空间
****************************************************************************/
void CBinary5GLDPCCodec::Free_Tanner_Graph()
{
    int i;
    EdgeLDPC *temp_edge;

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


/***************************************************************************
函数：Free()
功能：释放空间
****************************************************************************/

void CBinary5GLDPCCodec::Free()
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
    delete []m_syndromsoft;

    return;
}



/***************************************************************************
函数：SystH()
功能：把校验矩阵化为系统形式
****************************************************************************/

void CBinary5GLDPCCodec::SystH()
{
    int i, j, ii, jj, m, n;
    int temp;
    int flag;
    int *tempP;
    char **tempH;
    EdgeLDPC *p_edge;

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

    m_row_head = new EdgeLDPC[m_num_row];
    m_col_head = new EdgeLDPC[m_num_col];

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
                p_edge = new EdgeLDPC;
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

    return;
}


/***************************************************************************
函数：SystH_1()
功能：把5G LDPC 校验矩阵化为系统形式 单位阵在后面 系统位在前面
****************************************************************************/
void CBinary5GLDPCCodec::SystH_1()
{
    int i, j, ii, jj, m, n;
    int temp;
    int flag;
    int *tempP;
    char **tempH;
    EdgeLDPC *p_edge;

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

    m_row_head = new EdgeLDPC[m_num_row];
    m_col_head = new EdgeLDPC[m_num_col];

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
                p_edge = new EdgeLDPC;
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

    return;
}

void CBinary5GLDPCCodec::SystH_5G()
{
    int i, j, ii, jj, m, n;
    int temp;
    int flag;
    int *tempP;
    char **tempH;
    EdgeLDPC *p_edge;

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

    m_row_head = new EdgeLDPC[m_num_row];
    m_col_head = new EdgeLDPC[m_num_col];

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
                p_edge = new EdgeLDPC;
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

    return;
}

void CBinary5GLDPCCodec::SystH2()
{
    int i, j, ii, jj, m, n;
    int temp;
    int flag;
    int *tempP;
    char **tempH;
    EdgeLDPC *p_edge;

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

    m_row_head = new EdgeLDPC[m_num_row];
    m_col_head = new EdgeLDPC[m_num_col];

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
                p_edge = new EdgeLDPC;
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

    return;
}




/***************************************************************************
函数：Encoder(int *uu, int *cc)
功能：编码
输入参数：*uu   ----  信息
          *cc   ----  码字
****************************************************************************/

void CBinary5GLDPCCodec::Encoder(int *uu, int *cc)
{
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

    return;
}

/***************************************************************************
函数：Encoder_5G(int *uu, int *cc)  信息位在前面
功能：编码
输入参数：*uu   ----  信息
*cc   ----  码字
****************************************************************************/
void CBinary5GLDPCCodec::Encoder_5G(int *uu, int *cc)
{
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

    return;
}

/***************************************************************************
函数：InitMsg()
功能：迭代起始，将校验节点传给信息节点的信息初始化
****************************************************************************/
void CBinary5GLDPCCodec::InitMsg()
{
    int i;
    EdgeLDPC *p_edge;

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

/***************************************************************************
函数：Decoder(double *M2V, int *uu_hat)
功能：SPA算法
输入参数：*M2V    ---- 根据信息接收序列计算比特后验概率Message to Variable Nodes[M2V]
输出参数：*uu_hat ---- 译码结果
****************************************************************************/

int CBinary5GLDPCCodec::Decoder(double *M2V, int *uu_hat)
{
    int i;
    int iter;
    int parity_check;
    double temp0, temp1, temp_sum;

    EdgeLDPC *p_edge;

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
//			uu_hat[i] = m_cc_hat[i];//kite使用
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

/***************************************************************************
函数：Decoder_5G(double *M2V, int *uu_hat)
功能：SPA算法
输入参数：*M2V    ---- 根据信息接收序列计算比特后验概率Message to Variable Nodes[M2V]
输出参数：*uu_hat ---- 译码结果
****************************************************************************/

int CBinary5GLDPCCodec::Decoder_5G(double *M2V, int *uu_hat, int iter_cnt)
{
    int i;
    int iter;
    int parity_check;
    double temp0, temp1, temp_sum;

    EdgeLDPC *p_edge;

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
            //			uu_hat[i] = m_cc_hat[i+m_codechk];
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

            m_syndromsoft[i] = (m_row_head + i)->m_alpha[0];
        }
    }

    m_detection_error = 1 - m_success;//added by liangchulong @2013-01-23-19-04
    return iter + (iter<m_max_iter);
}

/***************************************************************************
函数：SoftOutHardOutDecoder(double *M2V, double* cc_hat_soft, int *uu_hat)
功能：软输出的SPA算法
输入参数：*M2V    ---- 根据信息接收序列计算比特后验概率Message to Variable Nodes[M2V]
  *cc_hat_soft    ---- 码字中每个比特的后验概率
输出参数：*uu_hat ---- 译码结果
****************************************************************************/

int CBinary5GLDPCCodec::SoftInSoftOut(double *M2V, int *uu_hat,  double* cc_hat_soft)
{
    int i;
    int iter;
    int parity_check;
    double temp0, temp1, temp_sum;

    EdgeLDPC *p_edge;

    //initialization
    //InitMsg();
    //iteration
    for (iter = 0; iter < m_max_iter; iter++){
        //from vnode to cnode 等号到加号
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

        //from c node to v node   加号到等号
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
    if (cc_hat_soft!=NULL)
    {
        for (i = 0; i < m_num_col; i++){
            //forward
            p_edge = (m_col_head+i)->down;
            p_edge->m_alpha[0] = 0.5;//不含输入信息
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


int CBinary5GLDPCCodec::Parity_Checking()
{
    int i;
    int iter;
    int parity_check;
    double temp0, temp1, temp_sum;

    EdgeLDPC *p_edge;

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


int CBinary5GLDPCCodec::parityCheck(int *cc)
{
    int i;
    int parity_check;
    int count = 0;

    EdgeLDPC *p_edge;

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
