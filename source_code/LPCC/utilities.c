#define _CRT_SECURE_NO_DEPRECATE
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "fcgimisc.h"
#include "LPCC.h"
#include "utilities.h"
extern int leftnodes;
//reading data
int readarray (FILE *in,
			   int *num_p,
			   double **data_p)
{
	int  status = 0;
	int  max, num;
	char ch;

	num = 0;
	max = 10;

	*data_p = (double*)malloc(max * sizeof(double));
	if ( *data_p == NULL ) {
		status = NO_MEMORY;
		goto TERMINATE;
	}

	for (;;) {
		fscanf (in, "%c", &ch);
		if ( ch == '\t' ||
			ch == '\r' ||
			ch == ' '  ||
			ch == '\n'   ) continue;
		if ( ch == '[' ) break;
		status = -1;
		goto TERMINATE;
	}

	for(;;) {
		int read;
		read = fscanf (in, "%lg", (*data_p)+num);
		if ( read == 0 ) {
			status = -1;
			goto TERMINATE;
		}
		num++;
		if ( num >= max ) {
			max *= 2;
			*data_p = (double*)realloc(*data_p, max * sizeof(double));
			if ( *data_p == NULL ) {
				status = NO_MEMORY;
				goto TERMINATE;
			}
		}
		do {
			fscanf (in, "%c", &ch);
		} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
		if ( ch == ']' ) break;
		else if ( ch != ',' ) {
			status = -1;
			goto TERMINATE;
		}
	}

	*num_p = num;

TERMINATE:

	return (status);
}

int readarray_int (FILE *in,
			   int *num_p,
			   int **data_p)
{
	int  status = 0;
	int  max, num;
	char ch;

	num = 0;
	max = 10;

	*data_p = (int*)malloc(max * sizeof(int));
	if ( *data_p == NULL ) {
		status = NO_MEMORY;
		goto TERMINATE;
	}

	for (;;) {
		fscanf (in, "%c", &ch);
		if ( ch == '\t' ||
			ch == '\r' ||
			ch == ' '  ||
			ch == '\n'   ) continue;
		if ( ch == '[' ) break;
		status = -1;
		goto TERMINATE;
	}

	for(;;) {
		int read;
		read = fscanf (in, "%d", (*data_p)+num);
		if ( read == 0 ) {
			status = -1;
			goto TERMINATE;
		}
		num++;
		if ( num >= max ) {
			max *= 2;
			*data_p = (int*)realloc(*data_p, max * sizeof(int));
			if ( *data_p == NULL ) {
				status = NO_MEMORY;
				goto TERMINATE;
			}
		}
		do {
			fscanf (in, "%c", &ch);
		} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
		if ( ch == ']' ) break;
		else if ( ch != ',' ) {
			status = -1;
			goto TERMINATE;
		}
	}

	*num_p = num;

TERMINATE:

	return (status);
}

int readarray_row (FILE *in,
			   int *num_p,
			   double **data_p)
{
	int  status = 0;
	int  max, num;
	char ch;

	num = 0;
	max = 10;

	*data_p = (double*)malloc(max * sizeof(double));
	if ( *data_p == NULL ) {
		status = NO_MEMORY;
		goto TERMINATE;
	}

	for(;;) {
		int read;
		do {
			fscanf (in, "%c", &ch);
		} while (ch == ' ' || ch == '\t' || ch=='\n' || ch=='\r');
		if(feof(in)) break;
		fseek(in,-1L,SEEK_CUR);
		
		read = fscanf (in, "%lg", (*data_p)+num);
		//fprintf(stdout,"%f ",(double)(*data_p)[num]);
		if ( read == 0 ) {
			status = -1;
			goto TERMINATE;
		}
		num++;
		if ( num >= max ) {
			max *= 2;
			*data_p = (double*)realloc(*data_p, max * sizeof(double));
			if ( *data_p == NULL ) {
				status = NO_MEMORY;
				goto TERMINATE;
			}
		}
		do {
			fscanf (in, "%c", &ch);
		} while ((ch == ' ' || ch == '\t') && !feof(in));
		if ( ch == '\n' || ch=='\r' || feof(in)) break;
		fseek(in,-1L,SEEK_CUR);
	}
	//fprintf(stdout,"\n");
	*num_p = num;

TERMINATE:
	return (status);
}

int readarray_col(FILE *in,
				   int *num_p,
				   double **data_p)
{
	int  status = 0;
	int  max, num;
	char ch;

	num = 0;
	max = 10;

	*data_p = (double*)malloc(max * sizeof(double));
	if ( *data_p == NULL ) {
		status = NO_MEMORY;
		goto TERMINATE;
	}

	for(;;) {
		int read;
		read = fscanf (in, "%lg", (*data_p)+num);
		if ( read == 0 ) {
			status = -1;
			goto TERMINATE;
		}
		//fprintf(stdout,"%f\n",(double)(*data_p)[num]);
		num++;
		if ( num >= max ) {
			max *= 2;
			*data_p = (double*)realloc(*data_p, max * sizeof(double));
			if ( *data_p == NULL ) {
				status = NO_MEMORY;
				goto TERMINATE;
			}
		}
		do {
			fscanf (in, "%c", &ch);
		} while ((ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r') && !feof(in));
		if (feof(in)) break;
		fseek(in,-1L,SEEK_CUR);
	}

	*num_p = num;

TERMINATE:

	return (status);
}

void usage (char *progname)
{
	fprintf (stderr,"Usage: %s <method> <datafile>\n", progname);
	fprintf (stderr," Exiting...\n");
} 
int get_matrix(const double **matrix_A, 
			   const int n_row,
			   const int n_col, 
			   int param, 
			   MATRIX *matrix_A_p) 
{
	int status = 0;
    int i, j;
    int nnz = 0;
    int count = 0;
    int max = 100;
    int *matbeg = NULL;
    int *matcnt = NULL;
    int *matind = NULL;
    double *matval = NULL;
	if ( param == 1) 
	{
		matbeg = (int *)malloc(n_row * sizeof(int));
		matcnt = (int *)malloc(n_row * sizeof(int)); 
		if ( matbeg == NULL || matcnt == NULL ) 
		{   
			status = NO_MEMORY;  
			fprintf(stderr,"get_matrix(): could not allocate memory\n");
			goto TERMINATE;
		}
		nnz = 0;
		matind = (int *)malloc(max * sizeof(int));
		matval = (double *)malloc(max * sizeof(double));
		if ( matind== NULL || matval == NULL ) 
		{   
			status = NO_MEMORY;   
			fprintf(stderr,"get_matrix(): could not allocate memory\n");
			goto TERMINATE;
		}

		for (i = 0; i < n_row; i++) 
		{
			count = 0; 
			for (j = 0; j < n_col; j++) 
			{   
				if (fabs(*(matrix_A[i]+j)) > ZERO_TOLERANCE) 
				{ 
					if ( nnz >= max ) 
					{
						max *= 2;
						matind = (int*)realloc(matind, max * sizeof(int));
						matval = (double *)realloc(matval, max * sizeof(double));
						if ( matind == NULL || matval == NULL) 
						{
							status = NO_MEMORY;
							goto TERMINATE;
						}
					}
					matind[nnz] = j;
				    matval[nnz] = *(matrix_A[i]+j);
					nnz++;
					count++;
				}
			}
			matbeg[i] = nnz - count;
			matcnt[i] = count;
		}
	}else if ( param == 0) 
	{
		matbeg = (int *)malloc(n_col * sizeof(int));
		matcnt = (int *)malloc(n_col * sizeof(int)); 
		if ( matbeg == NULL || matcnt == NULL ) 
		{   
			status = NO_MEMORY;  
			fprintf(stderr,"get_matrix(): could not allocate memory\n");
			goto TERMINATE;
		}

		nnz = 0;
		matind = (int *)malloc(max * sizeof(int));
		matval = (double *)malloc(max * sizeof(double));
		if ( matind== NULL || matval == NULL ) 
		{   
			status = NO_MEMORY;   
			fprintf(stderr,"get_matrix(): could not allocate memory\n");
			goto TERMINATE;
		}

		for (i = 0; i < n_col; i++) 
		{
			count = 0; 
			for (j = 0; j < n_row; j++) 
			{   
				if (fabs(*(matrix_A[j]+i)) > ZERO_TOLERANCE)
				{ 
					if ( nnz >= max ) 
					{
						max *= 2;
						matind = (int*)realloc(matind, max * sizeof(int));
						matval = (double *)realloc(matval, max * sizeof(double));
						if ( matind == NULL || matval == NULL) 
						{
							status = NO_MEMORY;
							goto TERMINATE;
						}
					}
					matind[nnz] = j;
				    matval[nnz] = *(matrix_A[j]+i);
					nnz++;
					count++;
				}
			}
			matbeg[i] = nnz - count;
			matcnt[i] = count;
		}
	}else 
	{
		status = -1;
		fprintf(stderr,"get_martix(): incorrect matrix_A.param\n");
		goto TERMINATE;
	}
	free_matrix(matrix_A_p);
	matrix_A_p->matbeg	= matbeg;
	matrix_A_p->matcnt	= matcnt;
	matrix_A_p->matind	= matind;
	matrix_A_p->matval	= matval;
	matrix_A_p->nnz		= nnz;
	matrix_A_p->param	= param;
	matrix_A_p->n_row	= n_row;
	matrix_A_p->n_col	= n_col;
TERMINATE:
    return (status);
}
int combine_matrix(MATRIX matrix_A, 
				   MATRIX matrix_B,
				   MATRIX *matrix_AB, 
				   int method )
{
	int status=0;
	int i,j;
	int count;
	int mat_AB_col;
	int mat_AB_row;
	int mat_AB_nnz;
	int mat_AB_param;
	int *matbeg=NULL;
	int *matind=NULL;
	int *matcnt=NULL;
	double *matval=NULL;
	switch(method)
	{
	case 0:
		//[A B]
		//check whether A.param == B.param
		if (matrix_A.param!=matrix_B.param)
		{
			status=-1;
			fprintf(stderr,"combine_matrix_r(): matrixA.param!=matrixB.param\n");
			goto TERMINATE;
		}else
		{
			mat_AB_param=matrix_A.param;
		}	
		//check whether A.row==B.row
		if (matrix_A.n_row!=matrix_B.n_row)
		{
			status=-1;
			fprintf(stderr,"combine_matrix_r(): matrixA.row!=matrixB.row\n");
			goto TERMINATE;
		}else
		{
			mat_AB_row=matrix_A.n_row;
		}
		mat_AB_nnz=matrix_A.nnz+matrix_B.nnz;
		mat_AB_col=matrix_A.n_col+matrix_B.n_col;
		matind=(int*) malloc(mat_AB_nnz*sizeof(int));
		matval=(double*) malloc(mat_AB_nnz*sizeof(double));
		if (matval==NULL || matind==NULL)
		{
			status=-1;
			fprintf(stderr,"combine_matrix_r(): unable to allocate memory\n");
			goto TERMINATE;
		}
		switch(mat_AB_param)
		{
		case 0:
			//group by column
			matbeg=(int*) malloc(mat_AB_col*sizeof(int));
			matcnt=(int*) malloc(mat_AB_col*sizeof(int));
			if (matbeg==NULL || matcnt==NULL)
			{
				status=NO_MEMORY;
				fprintf(stderr,"combine_matrix_r(): unable to allocate memory\n");
				goto TERMINATE;
			}
			for (i=0;i<mat_AB_nnz;i++)
			{
				if (i<matrix_A.nnz)
				{
					matind[i]=matrix_A.matind[i];
					matval[i]=matrix_A.matval[i];
				}else
				{
					matind[i]=matrix_B.matind[i-matrix_A.nnz];
					matval[i]=matrix_B.matval[i-matrix_A.nnz];
				}
			}
			for (i=0;i<mat_AB_col;i++)
			{
				if (i<matrix_A.n_col)
				{
					matbeg[i]=matrix_A.matbeg[i];
					matcnt[i]=matrix_A.matcnt[i];
				}else
				{
					matbeg[i]=matrix_A.nnz+matrix_B.matbeg[i-matrix_A.n_col];
					matcnt[i]=matrix_B.matcnt[i-matrix_A.n_col];
				}
			}
			break;
		case 1:
			//group by row
			matbeg=(int*) malloc(mat_AB_row*sizeof(int));
			matcnt=(int*) malloc(mat_AB_row*sizeof(int));
			if (matbeg==NULL || matcnt==NULL)
			{
				status=NO_MEMORY;
				fprintf(stderr,"combine_matrix_r(): unable to allocate memory\n");
				goto TERMINATE;
			}
			matbeg[0]=0;
			count=0;
			for (i=0;i<mat_AB_row;i++)
			{
				for (j=0;j<matrix_A.matcnt[i];j++)
				{
					matind[count]=matrix_A.matind[matrix_A.matbeg[i]+j];
					matval[count]=matrix_A.matval[matrix_A.matbeg[i]+j];
					count++;
				}
				for (j=0;j<matrix_B.matcnt[i];j++)
				{
					matind[count]=matrix_A.n_col+matrix_B.matind[matrix_B.matbeg[i]+j];
					matval[count]=matrix_B.matval[matrix_B.matbeg[i]+j];
					count++;
				}
				matcnt[i]=matrix_A.matcnt[i]+matrix_B.matcnt[i];
				if(i<(mat_AB_row-1))matbeg[i+1]=matbeg[i]+matcnt[i];
			}
			break;
		default:
			status=-1;
			fprintf(stderr," combine_matrix_r(): unexpected matrix.param: %d\n",mat_AB_param);
			goto TERMINATE;
		}
		free_matrix(matrix_AB);
		matrix_AB->n_col=mat_AB_col;
		matrix_AB->n_row=mat_AB_row;
		matrix_AB->nnz=mat_AB_nnz;
		matrix_AB->param=mat_AB_param;
		matrix_AB->matbeg=matbeg;
		matrix_AB->matcnt=matcnt;
		matrix_AB->matind=matind;
		matrix_AB->matval=matval;
		break;
	case 1:
		//[A
		// B]
		//check whether A.param == B.param
		if (matrix_A.param!=matrix_B.param)
		{
			status=-1;
			fprintf(stderr,"combine_matrix_r(): matrixA.param!=matrixB.param\n");
			goto TERMINATE;
		}else
		{
			mat_AB_param=matrix_A.param;
		}
		//check whether A.col==B.col
		if (matrix_A.n_col!=matrix_B.n_col)
		{
			status=-1;
			fprintf(stderr,"combine_matrix_r(): matrixA.col!=matrixB.col\n");
			goto TERMINATE;
		}else
		{
			mat_AB_col=matrix_A.n_col;
		}
		mat_AB_nnz=matrix_A.nnz+matrix_B.nnz;
		mat_AB_row=matrix_A.n_row+matrix_B.n_row;
		matind=(int*) malloc(mat_AB_nnz*sizeof(int));
		matval=(double*) malloc(mat_AB_nnz*sizeof(double));
		if (matval==NULL || matind==NULL)
		{
			status=-1;
			fprintf(stderr,"combine_matrix_r(): unable to allocate memory\n");
			goto TERMINATE;
		}
		switch(mat_AB_param)
		{
		case 0:
			//group by column
			matbeg=(int*) malloc(mat_AB_col*sizeof(int));
			matcnt=(int*) malloc(mat_AB_col*sizeof(int));
			if (matbeg==NULL || matcnt==NULL)
			{
				status=NO_MEMORY;
				fprintf(stderr,"combine_matrix_r(): unable to allocate memory\n");
				goto TERMINATE;
			}
			matbeg[0]=0;
			count=0;
			for (i=0;i<mat_AB_col;i++)
			{
				for (j=0;j<matrix_A.matcnt[i];j++)
				{
					matind[count]=matrix_A.matind[matrix_A.matbeg[i]+j];
					matval[count]=matrix_A.matval[matrix_A.matbeg[i]+j];
					count++;
				}
				for (j=0;j<matrix_B.matcnt[i];j++)
				{
					matind[count]=matrix_A.n_row+matrix_B.matind[matrix_B.matbeg[i]+j];
					matval[count]=matrix_B.matval[matrix_B.matbeg[i]+j];
					count++;
				}
				matcnt[i]=matrix_A.matcnt[i]+matrix_B.matcnt[i];
				if(i<(mat_AB_col-1))matbeg[i+1]=matbeg[i]+matcnt[i];
			}
			break;
		case 1:
			//group by row
			matbeg=(int*) malloc(mat_AB_row*sizeof(int));
			matcnt=(int*) malloc(mat_AB_row*sizeof(int));
			if (matbeg==NULL || matcnt==NULL)
			{
				status=NO_MEMORY;
				fprintf(stderr,"combine_matrix_r(): unable to allocate memory\n");
				goto TERMINATE;
			}
			for (i=0;i<mat_AB_nnz;i++)
			{
				if (i<matrix_A.nnz)
				{
					matind[i]=matrix_A.matind[i];
					matval[i]=matrix_A.matval[i];
				}else
				{
					matind[i]=matrix_B.matind[i-matrix_A.nnz];
					matval[i]=matrix_B.matval[i-matrix_A.nnz];
				}
			}
			for (i=0;i<mat_AB_row;i++)
			{
				if (i<matrix_A.n_row)
				{
					matbeg[i]=matrix_A.matbeg[i];
					matcnt[i]=matrix_A.matcnt[i];
				}else
				{
					matbeg[i]=matrix_A.nnz+matrix_B.matbeg[i-matrix_A.n_row];
					matcnt[i]=matrix_B.matcnt[i-matrix_A.n_row];
				}
			}
			break;
		default:
			status=-1;
			fprintf(stderr," combine_matrix_r(): unexpected matrix.param: %d\n",mat_AB_param);
			goto TERMINATE;
		}
		free_matrix(matrix_AB);
		matrix_AB->n_col=mat_AB_col;
		matrix_AB->n_row=mat_AB_row;
		matrix_AB->nnz=mat_AB_nnz;
		matrix_AB->param=mat_AB_param;
		matrix_AB->matbeg=matbeg;
		matrix_AB->matcnt=matcnt;
		matrix_AB->matind=matind;
		matrix_AB->matval=matval;
		break;
	default:
		status=-1;
		fprintf(stderr," combine_matrix(): unexpected method:%d\n",method);
		goto TERMINATE;
	}
TERMINATE:
	return(status);
}
int print_matrix ( MATRIX *matrix_A_p ) 
{
	int status = 0;
	int i,j;
	if ( matrix_A_p->param ==0) 
	{
		fprintf(stdout,"row: %d, col: %d, nonzero: %d\n",matrix_A_p->n_row,matrix_A_p->n_col,matrix_A_p->nnz);
		fprintf(stdout,"group by column\n");
		for ( i = 0; i < matrix_A_p->n_col; i++)
		{
			for ( j = matrix_A_p->matbeg[i]; j < (matrix_A_p->matbeg[i]+matrix_A_p->matcnt[i]); j++) 
			{
				fprintf(stdout,"(%d,%d):%f\n",matrix_A_p->matind[j],i,matrix_A_p->matval[j]);
			}
		}
	}else if (matrix_A_p->param == 1) 
	{
		fprintf(stdout,"row: %d, col: %d, nonzero: %d\n",matrix_A_p->n_row,matrix_A_p->n_col,matrix_A_p->nnz);
		fprintf(stdout,"group by row\n");
		for ( i = 0; i < matrix_A_p->n_row; i++)
		{
			for ( j = matrix_A_p->matbeg[i]; j < (matrix_A_p->matbeg[i]+matrix_A_p->matcnt[i]); j++) 
			{
				fprintf(stdout,"(%d,%d):%f\n",i,matrix_A_p->matind[j],matrix_A_p->matval[j]);
			}
		}
	}else 
	{
		status = -1;
		fprintf(stderr,"print_martix(): incorrect matrix_A.param\n");
		goto TERMINATE;
	}
TERMINATE:
	return (status);
}
int print_full_matrix(MATRIX *matrix_A_p)
{
	int status=0;
	int i,j;
	const int row=matrix_A_p->n_row;
	const int col=matrix_A_p->n_col;
	double **matrix_p=NULL;
	matrix_p=(double**)malloc(row*sizeof(double*));
	if (matrix_p==NULL)
	{
		status=NO_MEMORY;
		fprintf(stderr,"printf_full_matrix(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<row;i++)
	{
		matrix_p[i]=(double*)calloc(col,sizeof(double));
		if (matrix_p[i]==NULL)
		{
			status=NO_MEMORY;
			fprintf(stderr,"printf_full_matrix(): unable to allocate memory\n");
			goto TERMINATE;
		}
	}
	if ( matrix_A_p->param ==0) 
	{
		fprintf(stdout,"row: %d, col: %d, nonzero: %d\n",matrix_A_p->n_row,matrix_A_p->n_col,matrix_A_p->nnz);
		fprintf(stdout,"group by column\n");
		for ( i = 0; i < col; i++)
		{
			for ( j = matrix_A_p->matbeg[i]; j < (matrix_A_p->matbeg[i]+matrix_A_p->matcnt[i]); j++) 
			{
				matrix_p[matrix_A_p->matind[j]][i]=matrix_A_p->matval[j];
			}
		}
		for (i=0;i<row;i++)
		{
			for (j=0;j<col;j++)
			{
				fprintf(stdout,"%f,",matrix_p[i][j]);
			}
			fprintf(stdout,"\n");
		}
	}else if (matrix_A_p->param == 1) 
	{
		fprintf(stdout,"row: %d, col: %d, nonzero: %d\n",matrix_A_p->n_row,matrix_A_p->n_col,matrix_A_p->nnz);
		fprintf(stdout,"group by row\n");
		for ( i = 0; i < row; i++)
		{
			for ( j = matrix_A_p->matbeg[i]; j < (matrix_A_p->matbeg[i]+matrix_A_p->matcnt[i]); j++) 
			{
				matrix_p[i][matrix_A_p->matind[j]]=matrix_A_p->matval[j];
			}
		}
		for (i=0;i<row;i++)
		{
			for (j=0;j<col;j++)
			{
				fprintf(stdout,"%f,",matrix_p[i][j]);
			}
			fprintf(stdout,"\n");
		}
	}else 
	{
		status = -1;
		fprintf(stderr,"print_martix(): incorrect matrix_A.param\n");
		goto TERMINATE;
	}
TERMINATE:
	if ( matrix_p != NULL ) {
		for (i = 0; i <row ; ++i) {
			free_and_null ((char **) &(matrix_p[i]));
		}
	}
	free_and_null((char**)&matrix_p);
	return (status);
}
int transpose_matrix (MATRIX *matrix_A_p)
{
	int status = 0;
	int i,j;
	int n_row_t = matrix_A_p->n_col;
	int n_col_t = matrix_A_p->n_row;
	int nnz_t = matrix_A_p->nnz;
	int param_t = matrix_A_p->param;
	int *matbeg_t = NULL;
    int *matcnt_t = NULL;
    int *matind_t = NULL;
    double *matval_t = NULL;
	int *count=NULL;
	if ( matrix_A_p->param == 0 )
	{
		matbeg_t = (int *)malloc( n_col_t* sizeof(int));
		matcnt_t = (int *)calloc( n_col_t , sizeof(int)); 
		matind_t = (int *)malloc( nnz_t * sizeof(int));
		matval_t = (double *)malloc( nnz_t * sizeof(double));
		count = (int *) calloc ( n_col_t , sizeof(int));
		if ( matbeg_t == NULL || matcnt_t == NULL || (matind_t == NULL && nnz_t!=0) || (matval_t == NULL && nnz_t!=0) || count == NULL) 
		{   
			status = NO_MEMORY;  
			fprintf(stderr,"transpose_matrix(): could not allocate memory\n");
			goto TERMINATE;
		}
		for ( i = 0; i < nnz_t; i++) 
		{
			matcnt_t[matrix_A_p->matind[i]]++;
		}
		matbeg_t[0] = 0;
		for ( i = 1; i < n_col_t; i++) 
		{
			matbeg_t[i] = matbeg_t[i-1]+matcnt_t[i-1];
		}
		for ( i = 0; i < matrix_A_p->n_col; i++) 
		{
			for ( j = matrix_A_p->matbeg[i]; j < (matrix_A_p->matbeg[i]+matrix_A_p->matcnt[i]); j++) 
			{				
				matval_t[matbeg_t[matrix_A_p->matind[j]]+count[matrix_A_p->matind[j]]] = matrix_A_p->matval[j];
				matind_t[matbeg_t[matrix_A_p->matind[j]]+count[matrix_A_p->matind[j]]] = i;
				count[matrix_A_p->matind[j]]++;
			}
		}		
	}else if ( matrix_A_p->param == 1 ) 
	{
		matbeg_t = (int *)malloc( n_row_t* sizeof(int));
		matcnt_t = (int *)calloc( n_row_t , sizeof(int)); 
		matind_t = (int *)malloc( nnz_t * sizeof(int));
		matval_t = (double *)malloc( nnz_t * sizeof(double));
		count = (int *) calloc ( n_row_t , sizeof(int));
		if ( matbeg_t == NULL || matcnt_t == NULL ||  (matind_t == NULL && nnz_t!=0) || (matval_t == NULL && nnz_t!=0) || count == NULL) 
		{   
			status = NO_MEMORY;  
			fprintf(stderr,"transpose_matrix(): could not allocate memory\n");
			goto TERMINATE;
		}
		for ( i = 0; i < nnz_t; i++) 
		{
			matcnt_t[matrix_A_p->matind[i]]++;
		}
		matbeg_t[0] = 0;
		for ( i = 1; i < n_row_t; i++) 
		{
			matbeg_t[i] = matbeg_t[i-1]+matcnt_t[i-1];
		}
		for ( i = 0; i < matrix_A_p->n_row; i++) 
		{
			for ( j = matrix_A_p->matbeg[i]; j < (matrix_A_p->matbeg[i]+matrix_A_p->matcnt[i]); j++) 
			{				
				matval_t[matbeg_t[matrix_A_p->matind[j]]+count[matrix_A_p->matind[j]]] = matrix_A_p->matval[j];
				matind_t[matbeg_t[matrix_A_p->matind[j]]+count[matrix_A_p->matind[j]]] = i;
				count[matrix_A_p->matind[j]]++;
			}
		}	
	}else 
	{
		status = -1;
		fprintf(stderr,"transpose_matrix(): incorrect matrix_A.param\n");
		goto TERMINATE;
	}
	matrix_A_p->n_col = n_col_t;
	matrix_A_p->n_row = n_row_t;
	free_and_null ((char **)&(matrix_A_p->matbeg));
	free_and_null ((char **)&(matrix_A_p->matcnt));
	free_and_null ((char **)&(matrix_A_p->matind));
	free_and_null ((char **)&(matrix_A_p->matval));
	matrix_A_p->matbeg = matbeg_t;
	matrix_A_p->matcnt = matcnt_t;
	matrix_A_p->matind = matind_t;
	matrix_A_p->matval = matval_t;

TERMINATE:
	free_and_null ((char **)&count);
	return (status);
}



int get_matrix_element(MATRIX *matrix_A_p,
					   const int index_i,
					   const int index_j,
					   double *val_p) 
{
	int status = 0;
	int i;
	int check_flag = 1;
	if (( index_i < 0 ) || ( index_i >= matrix_A_p->n_row)) 
	{
		status = -1;
		fprintf(stderr,"get_matrix_element(): index i is not in the range of matrix\n");
		goto TERMINATE;
	}
	if (( index_j < 0 ) || ( index_j >= matrix_A_p->n_col)) 
	{
		status = -1;
		fprintf(stderr,"get_matrix_element(): index j is not in the range of matrix\n");
		goto TERMINATE;
	}
	if ( matrix_A_p->matbeg == NULL || matrix_A_p->matcnt == NULL || matrix_A_p->matind == NULL || matrix_A_p->matval == NULL) 
	{
		status = -1;
		fprintf(stderr,"get_matrix_element(): matrix_A is empty\n");
		goto TERMINATE;
	}
	if ( matrix_A_p->param == 0) 
	{	
		for ( i = 0; i < matrix_A_p->matcnt[index_j]; i++) 
		{
			if (index_i == matrix_A_p->matind[matrix_A_p->matbeg[index_j]+i]) 
			{
				*val_p = matrix_A_p->matval[matrix_A_p->matbeg[index_j]+i];
				check_flag = 0;
				break;
			}
		}
		if (check_flag) *val_p = 0;
	}else if (matrix_A_p->param == 1) 
	{
		for ( i = 0; i < matrix_A_p->matcnt[index_i]; i++) 
		{
			if (index_j == matrix_A_p->matind[matrix_A_p->matbeg[index_i]+i]) 
			{
				*val_p = matrix_A_p->matval[matrix_A_p->matbeg[index_i]+i];
				check_flag = 0;
				break;
			}
		}
		if (check_flag) *val_p = 0;
	}else
	{
		status = -1;
		fprintf(stderr,"get_martix_element(): incorrect matrix_A.param\n");
		goto TERMINATE;
	}
TERMINATE:
	return (status);
}

void init_matrix ( MATRIX *ptr) //done
{
	//initialize the matrix as gourp by column
	if (ptr != NULL)
	{
		ptr->n_row = 0;
		ptr->n_col = 0;
		ptr->nnz = 0;
		ptr->param = 0;
		ptr->matbeg=NULL;
		ptr->matcnt=NULL;
		ptr->matind=NULL;
		ptr->matval=NULL;
	}
}
void free_matrix ( MATRIX *ptr) //done
{
	if (ptr != NULL)
	{
		ptr->n_row = 0;
		ptr->n_col = 0;
		ptr->nnz = 0;
		ptr->param = 0;
		free_and_null ((char **) &(ptr->matbeg));
		free_and_null ((char **) &(ptr->matcnt));
		free_and_null ((char **) &(ptr->matind));
		free_and_null ((char **) &(ptr->matval));
	}
}
//void free_and_null_matrix ( MATRIX **ptr) //done
//{
//	if (*ptr != NULL)
//	{
//		free_and_null ((char **) &((*ptr)->matbeg));
//		free_and_null ((char **) &((*ptr)->matcnt));
//		free_and_null ((char **) &((*ptr)->matind));
//		free_and_null ((char **) &((*ptr)->matval));
//		free(*ptr);
//		*ptr=NULL;
//	}
//}

int identity_matrix (const int param,
					 const int matrix_size,
					 const double val,
					 MATRIX *matrix_A)
{
	int status = 0;
	int i;
	const int nnz=matrix_size;
	const int n_col=matrix_size;
	const int n_row=matrix_size;
	int *matbeg = NULL;
	int *matcnt = NULL;
	int *matind = NULL;
	double *matval = NULL;
	matbeg = (int *)malloc( matrix_size* sizeof(int));
	matcnt = (int *)malloc( matrix_size* sizeof(int)); 
	matind = (int *)malloc( matrix_size * sizeof(int));
	matval = (double *)malloc( matrix_size * sizeof(double));
	if ( matbeg == NULL || matcnt == NULL || matind == NULL || matval == NULL )
	{   
		status = NO_MEMORY;  
		fprintf(stderr," identity_matrix(): could not allocate memory\n");
		goto TERMINATE;
	}
	switch(param)
	{
	case 0:
		for (i=0;i<matrix_size;i++)
		{
			matbeg[i]=i;
			matcnt[i]=1;
		}
		for (i=0;i<nnz;i++)
		{
			matind[i]=i;
			matval[i]=val;
		}
		break;
	case 1:
		for (i=0;i<matrix_size;i++)
		{
			matbeg[i]=i;
			matcnt[i]=1;
		}
		for (i=0;i<nnz;i++)
		{
			matind[i]=i;
			matval[i]=val;
		}
		break;
	default:
		status=-1;
		fprintf(stderr," identity_matrix(): matrix_param is not 0 or 1\n");
		goto TERMINATE;
	}
	free_matrix(matrix_A);
	matrix_A->param = param;
	matrix_A->n_row = n_row;
	matrix_A->n_col = n_col;
	matrix_A->nnz = nnz;
	matrix_A->matbeg = matbeg;
	matrix_A->matcnt = matcnt;
	matrix_A->matind = matind;
	matrix_A->matval = matval;
TERMINATE:
	return (status);
}

int diag_matrix (const int param,
				 const int matrix_size,
				 const double *val,
				 MATRIX *matrix_A)
{
	int status = 0;
	int i;
	const int nnz=matrix_size;
	const int n_col=matrix_size;
	const int n_row=matrix_size;
	int *matbeg = NULL;
	int *matcnt = NULL;
	int *matind = NULL;
	double *matval = NULL;
	matbeg = (int *)malloc( matrix_size* sizeof(int));
	matcnt = (int *)malloc( matrix_size* sizeof(int)); 
	matind = (int *)malloc( matrix_size * sizeof(int));
	matval = (double *)malloc( matrix_size * sizeof(double));
	if ( matbeg == NULL || matcnt == NULL || matind == NULL || matval == NULL )
	{   
		status = NO_MEMORY;  
		fprintf(stderr," diag_matrix(): could not allocate memory\n");
		goto TERMINATE;
	}
	switch(param)
	{
	case 0:
		for (i=0;i<matrix_size;i++)
		{
			matbeg[i]=i;
			matcnt[i]=1;
		}
		for (i=0;i<nnz;i++)
		{
			matind[i]=i;
			matval[i]=val[i];
		}
		break;
	case 1:
		for (i=0;i<matrix_size;i++)
		{
			matbeg[i]=i;
			matcnt[i]=1;
		}
		for (i=0;i<nnz;i++)
		{
			matind[i]=i;
			matval[i]=val[i];
		}
		break;
	default:
		status=-1;
		fprintf(stderr," diag_matrix(): matrix_param is not 0 or 1\n");
		goto TERMINATE;
	}
	free_matrix(matrix_A);
	matrix_A->param = param;
	matrix_A->n_row = n_row;
	matrix_A->n_col = n_col;
	matrix_A->nnz = nnz;
	matrix_A->matbeg = matbeg;
	matrix_A->matcnt = matcnt;
	matrix_A->matind = matind;
	matrix_A->matval = matval;
TERMINATE:
	return (status);
}

int zero_matrix (const int param,
				 const int n_row,
				 const int n_col,
				 MATRIX *matrix_A)
{
	int status = 0;
	int i;
	const int nnz=0;
	int *matbeg = NULL;
	int *matcnt = NULL;
	switch(param)
	{
	case 0:
		matbeg = (int *)malloc( n_col* sizeof(int));
		matcnt = (int *)malloc( n_col* sizeof(int)); 
		if ( matbeg == NULL || matcnt == NULL)
		{   
			status = NO_MEMORY;  
			fprintf(stderr," zero_matrix(): could not allocate memory\n");
			goto TERMINATE;
		}
		for (i=0;i<n_col;i++)
		{
			matbeg[i]=0;
			matcnt[i]=0;
		}
		break;
	case 1:
		matbeg = (int *)malloc( n_row* sizeof(int));
		matcnt = (int *)malloc( n_row* sizeof(int)); 
		if ( matbeg == NULL || matcnt == NULL)
		{   
			status = NO_MEMORY;  
			fprintf(stderr," zero_matrix(): could not allocate memory\n");
			goto TERMINATE;
		}
		for (i=0;i<n_row;i++)
		{
			matbeg[i]=0;
			matcnt[i]=0;
		}
		break;
	default:
		status=-1;
		fprintf(stderr," zero_matrix(): matrix_param is not 0 or 1\n");
		goto TERMINATE;
	}
	free_matrix(matrix_A);
	matrix_A->param = param;
	matrix_A->n_row = n_row;
	matrix_A->n_col = n_col;
	matrix_A->nnz = nnz;
	matrix_A->matbeg = matbeg;
	matrix_A->matcnt = matcnt;
	matrix_A->matind = NULL;
	matrix_A->matval = NULL;
TERMINATE:
	return (status);
}

int copy_matrix ( MATRIX *matrix_A_p, MATRIX *matrix_A_cpy_p)//done 
{
	int status = 0;
	int i;
	int *matbeg = NULL;
    int *matcnt = NULL;
    int *matind = NULL;
    double *matval = NULL;
	if ( matrix_A_p->param == 0) 
	{
		matbeg = (int *)malloc( matrix_A_p->n_col* sizeof(int));
		matcnt = (int *)malloc( matrix_A_p->n_col * sizeof(int)); 
		matind = (int *)malloc( matrix_A_p->nnz * sizeof(int));
		matval = (double *)malloc( matrix_A_p->nnz * sizeof(double));
		if ( matbeg == NULL || matcnt == NULL || matind == NULL || matval == NULL )
		{   
			status = NO_MEMORY;  
			fprintf(stderr,"copy_matrix(): could not allocate memory\n");
			goto TERMINATE;
		}
		for ( i = 0; i< matrix_A_p->n_col; i++) 
		{
			matbeg[i] = matrix_A_p->matbeg[i];
			matcnt[i] = matrix_A_p->matcnt[i];
		}
		for ( i = 0; i< matrix_A_p->nnz; i++)
		{
			matind[i] = matrix_A_p->matind[i];
			matval[i] = matrix_A_p->matval[i];
		}
	}else if ( matrix_A_p->param == 1) 
	{
		matbeg = (int *)malloc( matrix_A_p->n_row* sizeof(int));
		matcnt = (int *)malloc( matrix_A_p->n_row * sizeof(int)); 
		matind = (int *)malloc( matrix_A_p->nnz * sizeof(int));
		matval = (double *)malloc( matrix_A_p->nnz * sizeof(double));
		if ( matbeg == NULL || matcnt == NULL || matind == NULL || matval == NULL ) 
		{   
			status = NO_MEMORY;  
			fprintf(stderr,"copy_martix(): could not allocate memory\n");
			goto TERMINATE;
		}
		for ( i = 0; i< matrix_A_p->n_row; i++)
		{
			matbeg[i] = matrix_A_p->matbeg[i];
			matcnt[i] = matrix_A_p->matcnt[i];
		}
		for ( i = 0; i< matrix_A_p->nnz; i++) 
		{
			matind[i] = matrix_A_p->matind[i];
			matval[i] = matrix_A_p->matval[i];
		}

	}else {
		status = -1;
		fprintf(stderr,"copy_martix(): incorrect matrix_A.param\n");
		goto TERMINATE;
	}
	free_matrix(matrix_A_cpy_p);
	matrix_A_cpy_p->param = matrix_A_p->param;
	matrix_A_cpy_p->n_row = matrix_A_p->n_row;
	matrix_A_cpy_p->n_col = matrix_A_p->n_col;
	matrix_A_cpy_p->nnz   = matrix_A_p->nnz;
	matrix_A_cpy_p->matbeg = matbeg;
	matrix_A_cpy_p->matcnt = matcnt;
	matrix_A_cpy_p->matind = matind;
	matrix_A_cpy_p->matval = matval;
TERMINATE:
	return (status);
}
void init_startinfo ( STARTINFO *ptr) //done
{
	if (ptr != NULL)
	{
		ptr->ccnt = 0;
		ptr->rcnt = 0;
		ptr->cindices=NULL;
		ptr->rindices=NULL;
		ptr->cstat=NULL;
		ptr->rstat=NULL;
	}
}
void free_startinfo ( STARTINFO *ptr) //done 
{
	if (ptr != NULL)
	{
		ptr->ccnt=0;
		ptr->rcnt=0;
		free_and_null ((char **) &(ptr->cindices));
		free_and_null ((char **) &(ptr->cstat));
		free_and_null ((char **) &(ptr->rindices));
		free_and_null ((char **) &(ptr->rstat));
	}
}
//void free_and_null_startinfo ( STARTINFO **ptr) //done 
//{
//	if (*ptr != NULL)
//	{
//		free_and_null ((char **) &((*ptr)->cindices));
//		free_and_null ((char **) &((*ptr)->cstat));
//		free_and_null ((char **) &((*ptr)->rindices));
//		free_and_null ((char **) &((*ptr)->rstat));
//		free(*ptr);
//		*ptr=NULL;
//	}
//}

int copy_startinfo (STARTINFO *start_p,
					STARTINFO *start_cpy_p) //done
{
	int status = 0;
	//int *cindices = NULL;
	int *cstat = NULL;
	//int *rindices = NULL;
	int *rstat = NULL;
	int i;
	if (start_p->ccnt==0)
	{
		//cindices = NULL;
		cstat = NULL;
	}else 
	{
		//cindices = (int *)malloc( start_p->ccnt* sizeof(int));
		cstat = (int *)malloc( start_p->ccnt * sizeof(int)); 
		if ( cstat == NULL ) //cindices == NULL || 
		{   
			status = NO_MEMORY;  
			fprintf(stderr,"copy_startinfo(): could not allocate memory\n");
			goto TERMINATE;
		}
	}
	if (start_p->rcnt == 0)
	{
		//rindices = NULL;
		rstat = NULL;
	}else 
	{
		//rindices = (int *)malloc( start_p->rcnt * sizeof(int));
		rstat = (int *)malloc( start_p->rcnt * sizeof(int));
		if (  rstat == NULL )//rindices == NULL ||
		{   
			status = NO_MEMORY;  
			fprintf(stderr,"copy_startinfo(): could not allocate memory\n");
			goto TERMINATE;
		}
	}
	for ( i = 0; i < start_p->ccnt; i++ ) 
	{
		//cindices[i] = start_p->cindices[i];
		cstat[i] = start_p->cstat[i];
	}
	for ( i = 0; i < start_p->rcnt; i++) 
	{
		//rindices[i] = start_p->rindices[i];
		rstat[i] = start_p->rstat[i];
	}
	free_startinfo(start_cpy_p);
	start_cpy_p->ccnt = start_p->ccnt;
	start_cpy_p->rcnt = start_p->rcnt;
	//start_cpy_p->cindices = cindices;
	start_cpy_p->cstat = cstat;
	//start_cpy_p->rindices = rindices;
	start_cpy_p->rstat = rstat;
TERMINATE:
	return (status);
}


void init_constraint_set ( CONSTRAINT_SET *ptr) 
{
	if (ptr != NULL) 
	{
		init_matrix(&ptr->cons_matrix);
		ptr->cons_matrix.nnz = 0;
		ptr->cons_matrix.n_col =0;
		ptr->cons_matrix.n_row =0;
		ptr->cons_matrix.param = 1;	
		ptr->rhs=NULL;
		ptr->sen=NULL;
	}
}

void free_constraint_set ( CONSTRAINT_SET *ptr) 
{
	if (ptr != NULL)
	{
		free_and_null((char **) &(ptr->sen));
		free_and_null((char **) &(ptr->rhs));
		free_matrix(&ptr->cons_matrix);
	}
}
//void free_and_null_constraint_set ( CONSTRAINT_SET **ptr) 
//{
//	if (*ptr != NULL)
//	{
//		free_and_null((char **) &((*ptr)->sen));
//		free_and_null((char **) &((*ptr)->rhs));
//		free_matrix(&((*ptr)->cons_matrix));
//		free(*ptr);
//		*ptr=NULL;
//	}
//}

int copy_constraint_set (CONSTRAINT_SET *cons_p, 
						 CONSTRAINT_SET *cons_cpy_p) 
{
	int status = 0;
	int i;
	double *rhs;
	char *sen;
	free_constraint_set(cons_cpy_p);
	init_constraint_set(cons_cpy_p);
	if (cons_p->cons_matrix.matbeg == NULL || cons_p->cons_matrix.matcnt == NULL || cons_p->cons_matrix.matind == NULL || cons_p->cons_matrix.matval == NULL || cons_p->rhs == NULL || cons_p->sen == NULL ) 
	{
		//status = -1;
		//fprintf(stderr,"copy_constraint_set(): cons is empty\n");
		goto TERMINATE;
	}
	
	if (cons_p->cons_matrix.param !=1)
	{
		status = -1;
		fprintf(stderr,"copy_constraint_set(): cons_matrix.param !=1\n");
		goto TERMINATE;
	}

	status = copy_matrix(&(cons_p->cons_matrix),&(cons_cpy_p->cons_matrix));
	if (status) 
	{
		fprintf(stderr,"copy_constraint_set(): unable to copy cons_matrix\n");
		goto TERMINATE;
	}
	rhs = (double *)malloc(cons_p->cons_matrix.n_row *sizeof(double));
	sen = (char *)malloc(cons_p->cons_matrix.n_row *sizeof(char));
	if ( rhs == NULL || sen == NULL) 
	{
		status = NO_MEMORY;  
		fprintf(stderr,"copy_constraint_set(): could not allocate memory\n");
		goto TERMINATE;
	}
	for ( i = 0; i < cons_p->cons_matrix.n_row; i++ )
	{
		rhs[i] = cons_p->rhs[i];
		sen[i] = cons_p->sen[i];
	}
	cons_cpy_p->rhs = rhs;
	cons_cpy_p->sen = sen;
TERMINATE:
	return (status);
}

int get_constraint(CONSTRAINT_SET *cons_p,
				   const int index, 
				   CONSTRAINT_SET **cons_i_p)
{
	int status = 0;
	int i;
	int n_row;
	int n_col;
	int nnz;
	double *rhs = NULL;
	char *sen = NULL;
	int *matbeg = NULL;
	int *matcnt = NULL;
	int *matind = NULL;
	double *matval = NULL;
	if (cons_p->cons_matrix.matbeg == NULL || cons_p->cons_matrix.matcnt == NULL || cons_p->cons_matrix.matind == NULL || cons_p->cons_matrix.matval == NULL || cons_p->rhs == NULL || cons_p->sen == NULL) 
	{
		status = -1;
		fprintf(stderr,"get_constraint(): cons is empty\n");
		goto TERMINATE;
	}
	if (cons_p->cons_matrix.param !=1) 
	{
		status = -1;
		fprintf(stderr," get_constraint(): cons.cons_matrix.param !=1\n");
		goto TERMINATE;
	}
	if ( index >= cons_p->cons_matrix.n_row )
	{
		status = -1;
		fprintf(stderr,"get_constraint(): index > n_row\n");
		goto TERMINATE;
	}
	nnz = cons_p->cons_matrix.matcnt[index];
	n_col = cons_p->cons_matrix.n_col;
	n_row = 1;
	rhs = (double *) malloc(n_row*sizeof(double));
	sen = (char *) malloc(n_row*sizeof(char));
	matbeg = (int *) malloc(n_row*sizeof(int));
	matcnt = (int *) malloc(n_row*sizeof(int));
	matind = (int *) malloc(nnz*sizeof(int));
	matval = (double *) malloc(nnz*sizeof(double));
	if ( rhs == NULL || sen == NULL || matbeg == NULL || matcnt == NULL || matind == NULL || matval == NULL)
	{
		status = NO_MEMORY;  
		fprintf(stderr,"get_constraint_set(): could not allocate memory\n");
		goto TERMINATE;
	}
	rhs[0] = cons_p->rhs[index];
	sen[0] = cons_p->sen[index];
	matbeg[0] = 0;
	matcnt[0] = nnz;
	for( i=0; i< nnz; i++)
	{
		matind[i] = cons_p->cons_matrix.matind[cons_p->cons_matrix.matbeg[index]+i];
		matval[i] = cons_p->cons_matrix.matval[cons_p->cons_matrix.matbeg[index]+i];
	}
	free_constraint_set(*cons_i_p);
	*cons_i_p = (CONSTRAINT_SET*) malloc(sizeof(CONSTRAINT_SET));
	(*cons_i_p)->rhs = rhs;
	(*cons_i_p)->sen = sen;
	(*cons_i_p)->cons_matrix.nnz = nnz;
	(*cons_i_p)->cons_matrix.param = 1;
	(*cons_i_p)->cons_matrix.n_row = n_row;
	(*cons_i_p)->cons_matrix.n_col = n_col;
	(*cons_i_p)->cons_matrix.matbeg = matbeg;
	(*cons_i_p)->cons_matrix.matcnt = matcnt;
	(*cons_i_p)->cons_matrix.matind = matind;
	(*cons_i_p)->cons_matrix.matval = matval;
TERMINATE:
	return (status);
}
int add_constraint(CONSTRAINT_SET *ptr,
				   int nnz, 
				   int *ind, 
				   double *val, 
				   double rhs) //done
{
	int status = 0;
	int init_n_row = ptr->cons_matrix.n_row;
	int init_nnz = ptr->cons_matrix.nnz;
	int i;
	//if (ptr->cons_matrix.param !=1) {
	//	status = -1;
	//	fprintf(stderr,"add_constraint(): cons_matrix.param !=1\n");
	//	goto TERMINATE;	
	//}
	if (nnz < 0) 
	{
		status = -1;
		fprintf(stderr,"add_constraint(): nnz < 0\n");
		goto TERMINATE;
	}
	if ( ind ==NULL || val == NULL ) 
	{
		status = -1;
		fprintf(stderr,"add_constraint(): ind == NULL or val == NULL\n");
		goto TERMINATE;
	}
	ptr->cons_matrix.n_row++;
	ptr->cons_matrix.nnz += nnz;
	ptr->cons_matrix.matbeg = (int*) realloc(ptr->cons_matrix.matbeg,ptr->cons_matrix.n_row*sizeof(int));
	ptr->cons_matrix.matcnt = (int*) realloc(ptr->cons_matrix.matcnt,ptr->cons_matrix.n_row*sizeof(int));
	ptr->cons_matrix.matind = (int*) realloc(ptr->cons_matrix.matind,ptr->cons_matrix.nnz*sizeof(int));
	ptr->cons_matrix.matval = (double*) realloc(ptr->cons_matrix.matval,ptr->cons_matrix.nnz*sizeof(double));
	ptr->rhs = (double*) realloc(ptr->rhs,ptr->cons_matrix.n_row*sizeof(double));
	ptr->sen = (char*) realloc(ptr->sen,ptr->cons_matrix.n_row*sizeof(char));
	if (ptr->cons_matrix.matbeg == NULL || ptr->cons_matrix.matcnt == NULL || ptr->cons_matrix.matind == NULL || ptr->cons_matrix.matval == NULL || ptr->rhs == NULL || ptr->sen == NULL) 
	{
		status = NO_MEMORY;  
		fprintf(stderr,"add_constraint(): could not allocate memory\n");
		goto TERMINATE;
	}
	ptr->sen[init_n_row] = 'G';
	ptr->rhs[init_n_row] = rhs;
	if(init_n_row==0)
	{
		ptr->cons_matrix.matbeg[init_n_row] = 0;
	}else
	{
		ptr->cons_matrix.matbeg[init_n_row] = (ptr->cons_matrix.matbeg[init_n_row-1]) +( ptr->cons_matrix.matcnt[init_n_row-1]);
	}
	
	ptr->cons_matrix.matcnt[init_n_row] = nnz;
	for (i = 0; i< nnz; i++)
	{
		ptr->cons_matrix.matind[init_nnz+i] = ind[i];
		ptr->cons_matrix.matval[init_nnz+i] = val[i];
	}
TERMINATE:
	return (status);
}
int add_cuts(CONSTRAINT_SET *ptr,
			 int new_cols, 
			 int nnz, 
			 int *ind, 
			 double *val, 
			 double rhs) 
{
	int status = 0;
	int init_n_row = ptr->cons_matrix.n_row;
	int init_nnz = ptr->cons_matrix.nnz;
	int i;
	if (ptr->cons_matrix.n_col<new_cols)
	{
		ptr->cons_matrix.n_col=new_cols;
	}
	if (nnz < 0) 
	{
		status = -1;
		fprintf(stderr,"add_cuts(): nnz < 0\n");
		goto TERMINATE;
	}
	if ( ind ==NULL || val == NULL ) 
	{
		status = -1;
		fprintf(stderr,"add_cuts(): ind == NULL or val == NULL\n");
		goto TERMINATE;
	}
	ptr->cons_matrix.n_row++;
	ptr->cons_matrix.nnz += nnz;
	ptr->cons_matrix.matbeg = (int*) realloc(ptr->cons_matrix.matbeg,ptr->cons_matrix.n_row*sizeof(int));
	ptr->cons_matrix.matcnt = (int*) realloc(ptr->cons_matrix.matcnt,ptr->cons_matrix.n_row*sizeof(int));
	ptr->cons_matrix.matind = (int*) realloc(ptr->cons_matrix.matind,ptr->cons_matrix.nnz*sizeof(int));
	ptr->cons_matrix.matval = (double*) realloc(ptr->cons_matrix.matval,ptr->cons_matrix.nnz*sizeof(double));
	ptr->rhs = (double*) realloc(ptr->rhs,ptr->cons_matrix.n_row*sizeof(double));
	ptr->sen = (char*) realloc(ptr->sen,ptr->cons_matrix.n_row*sizeof(char));
	if (ptr->cons_matrix.matbeg == NULL || ptr->cons_matrix.matcnt == NULL || ptr->cons_matrix.matind == NULL || ptr->cons_matrix.matval == NULL || ptr->rhs == NULL || ptr->sen == NULL) 
	{
		status = NO_MEMORY;  
		fprintf(stderr,"add_constraint(): could not allocate memory\n");
		goto TERMINATE;
	}
	ptr->sen[init_n_row] = 'G';
	ptr->rhs[init_n_row] = rhs;
	if(init_n_row==0)
	{
		ptr->cons_matrix.matbeg[init_n_row] = 0;
	}else
	{
		ptr->cons_matrix.matbeg[init_n_row] = (ptr->cons_matrix.matbeg[init_n_row-1]) +( ptr->cons_matrix.matcnt[init_n_row-1]);
	}

	ptr->cons_matrix.matcnt[init_n_row] = nnz;
	for (i = 0; i< nnz; i++)
	{
		ptr->cons_matrix.matind[init_nnz+i] = ind[i];
		ptr->cons_matrix.matval[init_nnz+i] = val[i];
	}
TERMINATE:
	return (status);
}
int add_constraint_a (CONSTRAINT_SET *ptr,
					  int cnt, 
					  double *original_array,
					  double rhs) 
{
	int status =0;
	int nnz;
	int *ind = NULL;
	double *val = NULL;
	status = convertArray(cnt,original_array,&nnz,&ind,&val);
	if (status) goto TERMINATE;
	status = add_constraint(ptr,nnz,ind,val,rhs);
	if (status) goto TERMINATE;
	ptr->cons_matrix.n_col = cnt;
TERMINATE:
	free_and_null ((char **) &ind);
	free_and_null ((char **) &val);
	return (status);
}
int convertArray (int cnt, 
				  double *original_array, 
				  int *nnz_p, 
				  int **ind_p, 
				  double ** val_p)
{
	int status =0;
	int i;
	int nnz = 0;
	int max =10;
	*ind_p = (int *) malloc( max*sizeof(int));
	*val_p = (double *) malloc ( max*sizeof(double));
	if (*ind_p == NULL || *val_p == NULL) 
	{
		status = NO_MEMORY;
		fprintf( stderr,"convertArrary(): could not allocate memory\n");
		goto TERMINATE;
	}
	for ( i=0; i< cnt; i++) 
	{
		if (fabs(original_array[i]) > ZERO_TOLERANCE) 
		{
			if (nnz>=max)
			{
				max*=2;
				*ind_p = (int *) realloc( *ind_p,max*sizeof(int));
				*val_p = (double *) realloc ( *val_p,max*sizeof(double));
				if (*ind_p == NULL || *val_p == NULL) 
				{
					status = NO_MEMORY;
					fprintf( stderr,"convertArrary(): could not allocate memory\n");
					goto TERMINATE;
				}
			}
			(*ind_p)[nnz] = i;
			(*val_p)[nnz] = original_array[i];
			nnz++;
		}
	}
	*nnz_p = nnz;
TERMINATE:
	return (status);
}

int checkComplementary(int index,
					   double *x)
{
	int status = 0;
	if (fabs(x[index])>ZERO_TOLERANCE) status = 1;
	return (status);
}

int checkallComplementary(const int param_n,
						  const int param_m, 
						  const int param_k, 
						  double *x) 
{
	int status = 0;
	int i;
	for ( i = 0; i < param_m; i++) 
	{
		if ((checkComplementary(param_n+i,x)==1) && (checkComplementary(param_n+param_m+i,x)==1))
		{
			status = 1;
			break;
		}
	}
	return (status);
}

int countallComplementary(const int param_n,
						  const int param_m, 
						  const int param_k,
						  double *x) 
{
	int cnt = 0;
	int i;
	for ( i = 0; i < param_m; i++) 
	{
		if ((checkComplementary(param_n+i,x)==1) && (checkComplementary(param_n+param_m+i,x)==1)) 
		{
			cnt++;

		}
	}
	return (cnt);
}

double weaken_value(const double value)
{
	double val;
	if (value>ZERO_TOLERANCE)
	{
		val = value*(1-WEAKEN_SCALE);
	}else if ( value < -ZERO_TOLERANCE) 
	{
		val = value*(1+WEAKEN_SCALE);
	}else {
		val = -ZERO_TOLERANCE;
	}
	return (val);
}

double strengthen_value(const double value)
{
	double val;
	if (value>ZERO_TOLERANCE)
	{
		val = value*(1+STRENGTHEN_SCALE);
	}else if ( value < -ZERO_TOLERANCE) 
	{
		val = value*(1-STRENGTHEN_SCALE);
	}else 
	{
		val = ZERO_TOLERANCE;
	}
	return (val);
}

//double max_value(double a, double b)
//{
//    return (a < b ? b : a);
//}
//
//double min_value(double a, double b)
//{
//    return (a < b ? a : b);
//}

void free_and_null (char **ptr)
{	
	//fprintf(stdout,"free debug1\n");
   if ( *ptr != NULL ) 
   {
	   //fprintf(stdout,"free debug2\n");
      free (*ptr);
      *ptr = NULL;
   }
   //fprintf(stdout,"free debug3\n");
} /* END free_and_null */  


int pushnode_ordered(NODE **head_ptr,
					 NODE newnode)//done
{ 
	int status = 0;
	NODE *ptr= NULL;
	NODE *pre_ptr=NULL;
	NODE *newnode_ptr = NULL;
	leftnodes++;
	newnode_ptr =(NODE*) malloc(sizeof(NODE));
	init_node(newnode_ptr);
	status = copy_node(&newnode,newnode_ptr);
	if (status) goto TERMINATE;

	//addNodeBranchInfo(branch_hist,&newnode);	
	if(*head_ptr==NULL)
	{
		*head_ptr = newnode_ptr;
	}else
	{
		ptr=*head_ptr;
		while(newnode_ptr->nodelb>ptr->nodelb)
		{
			pre_ptr=ptr;
			ptr=ptr->link;
			if (ptr==NULL) break;
		}
		if(pre_ptr==NULL)
		{
			*head_ptr=newnode_ptr;
			newnode_ptr->link=ptr;
		}else
		{
			pre_ptr->link=newnode_ptr;
			newnode_ptr->link=ptr;
		}

	}
TERMINATE:
	return(status);
}
int pushnode_unordered( NODE **head_ptr,
					    NODE newnode)//done
{ 
	int status = 0;
	NODE *ptr= NULL;
	NODE *pre_ptr=NULL;
	NODE *newnode_ptr = NULL;
	leftnodes++;
	newnode_ptr =(NODE*) malloc(sizeof(NODE));
	init_node(newnode_ptr);
	status = copy_node(&newnode,newnode_ptr);
	if (status) goto TERMINATE;
	//addNodeBranchInfo(branch_hist,&newnode);
	if(*head_ptr==NULL)
	{
		*head_ptr = newnode_ptr;
	}else
	{
		ptr=*head_ptr;
		*head_ptr=newnode_ptr;
		newnode_ptr->link=ptr;	
	}
TERMINATE:
	return(status);
}
void popnode(NODE **head_ptr, 
			 NODE **topnode_ptr)//done
{
	if (*head_ptr==NULL)
	{
		*topnode_ptr=NULL;
	}else
	{
		*topnode_ptr=*head_ptr;
		*head_ptr=(*head_ptr)->link;
		(*topnode_ptr)->link = NULL;
		//removeNodeBranchInfo(branch_hist,*topnode_ptr);
	}
	leftnodes--;
}

double getLastnodevalue(NODE **head_ptr)
{
	double value;
	NODE *tmp_node_ptr=NULL;
	if(*head_ptr==NULL)
	{
		value=INF_BOUND;
	}else
	{
		tmp_node_ptr=*head_ptr;
		while(1)
		{
			value=tmp_node_ptr->nodelb;
			tmp_node_ptr=tmp_node_ptr->link;
			if(tmp_node_ptr==NULL)break;
		}	
	}
	return(value);
}
double getaveragenodevalue(NODE **head_ptr)
{
	double value=0;
	double sum = 0;
	int cnt=0;
	NODE *tmp_node_ptr=NULL;
	if(*head_ptr==NULL)
	{
		value=INF_BOUND;
	}else
	{
		tmp_node_ptr=*head_ptr;
		while(1)
		{
			sum+=tmp_node_ptr->nodelb;
			cnt++;
			tmp_node_ptr=tmp_node_ptr->link;
			if(tmp_node_ptr==NULL)break;
		}
		value=sum/cnt;

	}
	return(value);
}
double getmediannodevalue(NODE **head_ptr)
{
	double value=0;
	double sum = 0;
	int cnt=0;
	int mid_index=(int) (leftnodes/2);
	NODE *tmp_node_ptr=NULL;
	if(*head_ptr==NULL)
	{
		value=INF_BOUND;
	}else
	{
		tmp_node_ptr=*head_ptr;
		while(1)
		{
			value=tmp_node_ptr->nodelb;
			if(cnt>mid_index)break;
			cnt++;
			tmp_node_ptr=tmp_node_ptr->link;
			if(tmp_node_ptr==NULL)break;
		}
	}
	return(value);
}

int copy_node(NODE *o_node_p,
			 NODE *cpynode_ptr)//done
{
	int i;
	int status =0;
	free_node(cpynode_ptr);
	cpynode_ptr->link = NULL;
	cpynode_ptr->n = o_node_p->n;
	cpynode_ptr->m=o_node_p->m;
	cpynode_ptr->k = o_node_p->k;
	cpynode_ptr->nodevalue = o_node_p->nodevalue;
	cpynode_ptr->nodelb = o_node_p->nodelb;
	cpynode_ptr->branchDirectionofParentNode=o_node_p->branchDirectionofParentNode;
	cpynode_ptr->branchIndexofParentNode=o_node_p->branchIndexofParentNode;
	cpynode_ptr->parentNodeValue=o_node_p->parentNodeValue;
	cpynode_ptr->parentNodeBranchCoefNorm=o_node_p->parentNodeBranchCoefNorm;
	cpynode_ptr->parentNodeBranchReduceCostNorm=o_node_p->parentNodeBranchReduceCostNorm;
	cpynode_ptr->parentNodeBranchViolation=o_node_p->parentNodeBranchViolation;
	status=copy_startinfo(&(o_node_p->nodestartinfo),&(cpynode_ptr->nodestartinfo));
	if (status) goto TERMINATE;
	status=copy_constraint_set(&(o_node_p->nodecuts),&(cpynode_ptr->nodecuts));
	if (status) goto TERMINATE;
	cpynode_ptr->nodebranchinfo = (int*) malloc(o_node_p->m*sizeof(int));
	cpynode_ptr->nodex = (double*) malloc ((o_node_p->m+o_node_p->n+o_node_p->k+o_node_p->m)*sizeof(double));
	cpynode_ptr->x_start_index=o_node_p->x_start_index;
	cpynode_ptr->y_bar_start_index=o_node_p->y_bar_start_index;
	cpynode_ptr->mcc_index_col=o_node_p->mcc_index_col;
	cpynode_ptr->mcc_index_row=o_node_p->mcc_index_row;
	cpynode_ptr->x_lb=(double*) malloc(o_node_p->n*sizeof(double));
	cpynode_ptr->x_ub=(double*) malloc(o_node_p->n*sizeof(double));
	cpynode_ptr->y_bar_lb=(double*) malloc(o_node_p->n*sizeof(double));
	cpynode_ptr->y_bar_ub=(double*) malloc(o_node_p->n*sizeof(double));
	if(cpynode_ptr->nodebranchinfo==NULL ||
	   cpynode_ptr->nodex==NULL ||
	   cpynode_ptr->x_lb==NULL ||
	   cpynode_ptr->x_ub==NULL ||
	   cpynode_ptr->y_bar_lb==NULL ||
	   cpynode_ptr->y_bar_ub==NULL)
	{
		status =NO_MEMORY;
		fprintf(stderr," copynode(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<o_node_p->n;i++)
	{
		cpynode_ptr->x_lb[i]=o_node_p->x_lb[i];
		cpynode_ptr->x_ub[i]=o_node_p->x_ub[i];
		cpynode_ptr->y_bar_lb[i]=o_node_p->y_bar_lb[i];
		cpynode_ptr->y_bar_ub[i]=o_node_p->y_bar_ub[i];
	}
	for (i=0; i<o_node_p->m; i++)
	{
		cpynode_ptr->nodebranchinfo[i]=o_node_p->nodebranchinfo[i];
		//cpynode_ptr->branchscore[i]=o_node_p->branchscore[i];
	}

	for (i=0; i<(o_node_p->m+o_node_p->n+o_node_p->k+o_node_p->m); i++)
	{
		cpynode_ptr->nodex[i] = o_node_p->nodex[i];
	}
TERMINATE:
	return(status);
}

int init_branch_nodes(NODE activenode,
					  NODE **branchnode_1_ptr, 
					  NODE **branchnode_2_ptr)
{
	int status = 0;
	status=init_branch_node(activenode,branchnode_1_ptr);
	if(status) goto TERMINATE;
	status=init_branch_node(activenode,branchnode_2_ptr);
	if(status) goto TERMINATE;
TERMINATE:
	return(status);
}
int init_branch_node(NODE activenode,
					 NODE **branchnode_ptr)
{
	int i;
	int status = 0;
	*branchnode_ptr =(NODE*) malloc(sizeof(NODE));
	init_node(*branchnode_ptr);
	(*branchnode_ptr)->nodevalue=activenode.nodevalue;
	(*branchnode_ptr)->nodelb = activenode.nodelb;
	(*branchnode_ptr)->link = NULL;
	(*branchnode_ptr)->n = activenode.n;
	(*branchnode_ptr)->m = activenode.m;
	(*branchnode_ptr)->k = activenode.k;
	(*branchnode_ptr)->parentNodeValue = activenode.nodevalue;
	(*branchnode_ptr)->x_start_index=activenode.x_start_index;
	(*branchnode_ptr)->y_bar_start_index=activenode.y_bar_start_index;
	(*branchnode_ptr)->mcc_index_col=activenode.mcc_index_col;
	(*branchnode_ptr)->mcc_index_row=activenode.mcc_index_row;
	(*branchnode_ptr)->nodebranchinfo = (int*) malloc(activenode.m*sizeof(int));
	(*branchnode_ptr)->nodex = (double*) malloc ((activenode.n+activenode.m+activenode.m+activenode.k)*sizeof(double));
	(*branchnode_ptr)->x_lb =(double*) malloc(activenode.n*sizeof(double));
	(*branchnode_ptr)->x_ub =(double*) malloc(activenode.n*sizeof(double));
	(*branchnode_ptr)->y_bar_lb =(double*) malloc(activenode.n*sizeof(double));
	(*branchnode_ptr)->y_bar_ub =(double*) malloc(activenode.n*sizeof(double));
	if((*branchnode_ptr)->nodebranchinfo==NULL ||
	   (*branchnode_ptr)->nodex == NULL ||
	   (*branchnode_ptr)->x_lb==NULL ||
	   (*branchnode_ptr)->x_ub==NULL ||
	   (*branchnode_ptr)->y_bar_lb==NULL ||
	   (*branchnode_ptr)->y_bar_ub==NULL)
	{
		status=NO_MEMORY;
		fprintf(stderr," init_branch_node(): unable to allocate memory\n");
		goto TERMINATE;
	}
	status= copy_constraint_set(&(activenode.nodecuts),&(*branchnode_ptr)->nodecuts);
	if(status) goto TERMINATE;
	for (i=0;i<activenode.n;i++)
	{
		(*branchnode_ptr)->x_lb[i]=activenode.x_lb[i];
		(*branchnode_ptr)->x_ub[i]=activenode.x_ub[i];
		(*branchnode_ptr)->y_bar_lb[i]=activenode.y_bar_lb[i];
		(*branchnode_ptr)->y_bar_ub[i]=activenode.y_bar_ub[i];
	}
	for (i=0; i<activenode.m; i++)
	{
		(*branchnode_ptr)->nodebranchinfo[i] = activenode.nodebranchinfo[i];
	}
TERMINATE:
	return(status);
}
void deletenode(NODE **node_ptr)//done
{
	
	NODE *tempnode = NULL;
	while(*node_ptr!=NULL)
	{
		popnode(node_ptr,&tempnode);
		free_node(tempnode);
		free(tempnode);
		tempnode=NULL;
	}
}

void filternode(NODE **head_ptr, 
				double val)
{ 
	
	NODE *ptr= NULL;
	NODE *pre_ptr=NULL;
	if(*head_ptr!=NULL)
	{
		ptr=*head_ptr;
		while(val>strengthen_value(ptr->nodelb))
		{
			pre_ptr=ptr;
			ptr=ptr->link;
			if (ptr==NULL) break;
		}
		if(pre_ptr==NULL)
		{
			deletenode(head_ptr);
		}else
		{
			pre_ptr->link=NULL;
			deletenode(&ptr);
		}
	}
}

int fixedcnt_node(NODE *node_ptr)
{
	int fixed_cnt=0;
	int i;
	for (i=0;i<node_ptr->m;i++)
	{
		if (node_ptr->nodebranchinfo[i]!=0)
		{
			fixed_cnt++;
		}
	}
	return(fixed_cnt);
}
void init_node(NODE *ptr)//done
{ 
	ptr->link=NULL;
	ptr->n=0;
	ptr->m=0;
	ptr->k = 0;
	ptr->nodebranchinfo = NULL;
	init_startinfo(&(ptr->nodestartinfo));
	init_constraint_set(&(ptr->nodecuts));
	ptr->nodevalue = -INF_BOUND;
	ptr->nodelb = -INF_BOUND;
	ptr->nodex =NULL;
	ptr->x_start_index=-1;
	ptr->y_bar_start_index=-1;
	ptr->mcc_index_col=-1;
	ptr->mcc_index_row=-1;
	ptr->x_lb=NULL;
	ptr->x_ub=NULL;
	ptr->y_bar_lb=NULL;
	ptr->y_bar_ub=NULL;
	ptr->branchDirectionofParentNode=-1;
	ptr->branchIndexofParentNode=-1;
	ptr->parentNodeValue=-INF_BOUND;
	ptr->parentNodeBranchCoefNorm=-1;
	ptr->parentNodeBranchReduceCostNorm=-1;
	ptr->parentNodeBranchViolation=-1;
}

void free_node( NODE *node_ptr)//done
{
	if (node_ptr!=NULL)
	{
		node_ptr->link = NULL;
		free_startinfo(&(node_ptr->nodestartinfo));
		free_constraint_set(&(node_ptr->nodecuts));
		free_and_null((char**) &(node_ptr->nodebranchinfo));
		free_and_null((char**) &(node_ptr->nodex));
		node_ptr->x_start_index=-1;
		node_ptr->y_bar_start_index=-1;
		node_ptr->mcc_index_col=-1;
		node_ptr->mcc_index_row=-1;
		free_and_null((char**) &(node_ptr->x_lb));
		free_and_null((char**) &(node_ptr->x_ub));
		free_and_null((char**) &(node_ptr->y_bar_lb));
		free_and_null((char**) &(node_ptr->y_bar_ub));
		node_ptr->n=0;
		node_ptr->m=0;
		node_ptr->k=0;
		node_ptr->branchDirectionofParentNode=-1;
		node_ptr->branchIndexofParentNode=-1;
		node_ptr->parentNodeValue=-INF_BOUND;
		node_ptr->nodevalue=-INF_BOUND;
		node_ptr->nodelb=-INF_BOUND;
		node_ptr->parentNodeBranchCoefNorm=-1;
		node_ptr->parentNodeBranchReduceCostNorm=-1;
		node_ptr->parentNodeBranchViolation=-1;
	}
}

//void free_and_null_node( NODE **node_ptr)//done
//{
//	if (*node_ptr!=NULL)
//	{
//		free_node(*node_ptr);
//		//(*node_ptr)->link = NULL;
//		//free_startinfo(&((*node_ptr)->nodestartinfo));
//		//free_constraint_set(&((*node_ptr)->nodecuts));
//		//free_and_null((char**) &((*node_ptr)->nodebranchinfo));
//		//free_and_null((char**) &((*node_ptr)->nodex));
//		free(*node_ptr);
//		*node_ptr=NULL;
//	}
//}
void displaynode( NODE *head_ptr)//done
{ 
	NODE *ptr = NULL;
	ptr= head_ptr;
	while(ptr!=NULL)
	{
		fprintf(stdout,"%f\n",ptr->nodevalue);
		ptr=ptr->link;
	}
}
int get_sort_index(const int cnt,
				   double *sort_val,
				   int *sorted_index,
				   const int sort_order)
{
	int status=0;
	int i,j;
	int tmp_index;
	double *val=NULL;
	double tmp_val;
	val=(double*)malloc(cnt*sizeof(double));
	if (val==NULL)
	{
		status=-1;
		fprintf(stderr," get_sort_index(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for(i=0;i<cnt;i++)
	{
		val[i]=sort_val[i];
		sorted_index[i]=i;
	}
	switch(sort_order)
	{
	case 1:// increase order
		for (i=1;i<cnt;i++)
		{
			tmp_val=val[i];
			tmp_index=sorted_index[i];
			j=i-1;
			while(val[j]>tmp_val && j>=0)
			{
				val[j+1]=val[j];
				sorted_index[j+1]=sorted_index[j];
				j--;
			}
			val[j+1]=tmp_val;
			sorted_index[j+1]=tmp_index;
		}
		break;
	case -1://decrease order
		for (i=1;i<cnt;i++)
		{
			tmp_val=val[i];
			tmp_index=sorted_index[i];
			j=i-1;
			while(val[j]<tmp_val && j>=0)
			{
				val[j+1]=val[j];
				sorted_index[j+1]=sorted_index[j];
				j--;
			}
			val[j+1]=tmp_val;
			sorted_index[j+1]=tmp_index;
		}
		break;
	default:
		status=-1;
		fprintf(stderr," get_sort_index(): unexpected sort_order:%d",sort_order);
		goto TERMINATE;
	}

	
	//fprintf(stdout,"initial:\n");
	//for (i=0;i<cnt;i++)
	//{
	//	fprintf(stdout,"%f\n",sort_val[i]);
	//}
	//fprintf(stdout,"sorted:\n");
	//for (i=0;i<cnt;i++)
	//{
	//	fprintf(stdout,"%f\n",val[i]);
	//	fprintf(stdout,"index: %d\n",sorted_index[i]);
	//}
TERMINATE:
	free_and_null((char**)&val);
	return(status);
}

double sum_max_t(const int cnt,
				 double *val,
				 const int t)
{
	int i,j;
	double tmp_val;
	double sum;
	//fprintf(stdout,"initial:\n");
	//for (i=0;i<cnt;i++)
	//{
	//	fprintf(stdout,"%f\n",val[i]);
	//}
	for (i=1;i<cnt;i++)
	{
		tmp_val=val[i];
		j=i-1;
		while(val[j]<tmp_val && j>=0)
		{
			val[j+1]=val[j];
			j--;
		}
		val[j+1]=tmp_val;
	}
	sum=0;
	for (i=0;i<t;i++)
	{
		sum+=val[i];
	}
	//fprintf(stdout,"sorted:\n");
	//for (i=0;i<cnt;i++)
	//{
	//	fprintf(stdout,"%f\n",val[i]);
	//}
	//fprintf(stdout,"sum: %f,t: %d\n",sum,t);
	return(sum);
}

int symmetrizing_matrix(MATRIX *matrix_A_p, 
						MATRIX *sym_matrix_A_p, 
						int *sym_status)
{
	int status=0;
	int i;
	MATRIX matrix_A_t;
	if (matrix_A_p->n_col!=matrix_A_p->n_row)
	{
		status=-1;
		fprintf(stderr," symmetrizing_matrix(): matrix is not square matrix\n");
		goto TERMINATE;
	}
	init_matrix(&matrix_A_t);
	status=copy_matrix(matrix_A_p,&matrix_A_t);
	if(status) goto TERMINATE;
	status=transpose_matrix(&matrix_A_t);
	if(status) goto TERMINATE;
	//status=print_matrix(matrix_A_p);
	//if(status) goto TERMINATE;
	//status=print_matrix(&matrix_A_t);
	//if(status) goto TERMINATE;
	status=matrix_addition(*matrix_A_p,matrix_A_t,sym_matrix_A_p);
	if(status) goto TERMINATE;
	for (i=0;i<sym_matrix_A_p->nnz;i++)
	{
		sym_matrix_A_p->matval[i]=(sym_matrix_A_p->matval[i])/2;
	}
	if (sym_matrix_A_p->nnz==0)
	{
		*sym_status=0;
	}else
	{
		*sym_status=1;
	}
TERMINATE:
	free_matrix(&matrix_A_t);
	return(status);
}

int matrix_addition(MATRIX matrix_A, 
					MATRIX matrix_B,
					MATRIX *matrix_sum)
{
	int status=0;
	int i;
	int cnt;
	int A_cnt;
	int B_cnt;
	int matrix_sum_param;
	int matrix_sum_nnz;
	int matrix_sum_ncol;
	int matrix_sum_nrow;
	int *matrix_sum_matbeg=NULL;
	int *matrix_sum_matcnt=NULL;
	int *matrix_sum_matind=NULL;
	double *matrix_sum_matval=NULL;
	if (matrix_A.param!=matrix_B.param)
	{
		status=-1;
		fprintf(stderr," matrix_addition(): matrix's param is not equal \n");
		goto TERMINATE;
	}
	if (matrix_A.n_col!=matrix_B.n_col || matrix_A.n_row!=matrix_B.n_row)
	{
		status=-1;
		fprintf(stderr," matrix_addition(): two matrix's row or col is not equal\n");
		goto TERMINATE;
	}
	matrix_sum_nrow=matrix_A.n_row;
	matrix_sum_ncol=matrix_A.n_col;
	matrix_sum_param=matrix_A.param;
	if (matrix_A.nnz==0)
	{
		status=copy_matrix(&matrix_B,matrix_sum);
		if(status) goto TERMINATE;
	}else if (matrix_B.nnz==0)
	{
		status=copy_matrix(&matrix_A,matrix_sum);
		if(status) goto TERMINATE;
	}else
	{
		if (matrix_sum_param==0)//group by col
		{
			matrix_sum_matind=(int*)malloc((matrix_A.nnz+matrix_B.nnz)*sizeof(int));
			matrix_sum_matval=(double*)malloc((matrix_A.nnz+matrix_B.nnz)*sizeof(double));
			if (matrix_sum_matind==NULL || matrix_sum_matval==NULL)
			{
				status=NO_MEMORY;
				fprintf(stderr," matrix_addition(): unable to allocate memory\n");
				goto TERMINATE;
			}
			matrix_sum_matbeg=(int*) malloc(matrix_sum_ncol*sizeof(int));
			matrix_sum_matcnt=(int*) malloc(matrix_sum_ncol*sizeof(int));
			if (matrix_sum_matbeg==NULL || matrix_sum_matcnt==NULL)
			{
				status=NO_MEMORY;
				fprintf(stderr," matrix_addition(): unable to allocate memory\n");
				goto TERMINATE;
			}
			matrix_sum_nnz=0;
			for (i=0;i<matrix_sum_ncol;i++)
			{
				cnt=0;
				A_cnt=0;
				B_cnt=0;
				while(A_cnt<matrix_A.matcnt[i] || B_cnt<matrix_B.matcnt[i])
				{
					if (A_cnt<matrix_A.matcnt[i] && B_cnt<matrix_B.matcnt[i])
					{
						if (matrix_A.matind[matrix_A.matbeg[i]+A_cnt]<matrix_B.matind[matrix_B.matbeg[i]+B_cnt])
						{
							matrix_sum_matind[matrix_sum_nnz]=matrix_A.matind[matrix_A.matbeg[i]+A_cnt];
							matrix_sum_matval[matrix_sum_nnz]=matrix_A.matval[matrix_A.matbeg[i]+A_cnt];
							matrix_sum_nnz++;
							cnt++;
							A_cnt++;
						}else if (matrix_A.matind[matrix_A.matbeg[i]+A_cnt]>matrix_B.matind[matrix_B.matbeg[i]+B_cnt])
						{
							matrix_sum_matind[matrix_sum_nnz]=matrix_B.matind[matrix_B.matbeg[i]+B_cnt];
							matrix_sum_matval[matrix_sum_nnz]=matrix_B.matval[matrix_B.matbeg[i]+B_cnt];
							matrix_sum_nnz++;
							cnt++;
							B_cnt++;
						}else
						{	if (fabs(matrix_A.matval[matrix_A.matbeg[i]+A_cnt]+matrix_B.matval[matrix_B.matbeg[i]+B_cnt])>ZERO_TOLERANCE)
							{
								matrix_sum_matind[matrix_sum_nnz]=matrix_A.matind[matrix_A.matbeg[i]+A_cnt];
								matrix_sum_matval[matrix_sum_nnz]=matrix_A.matval[matrix_A.matbeg[i]+A_cnt]+matrix_B.matval[matrix_B.matbeg[i]+B_cnt];
								matrix_sum_nnz++;
								cnt++;
							}
							A_cnt++;
							B_cnt++;
						}
					}else if (A_cnt<matrix_A.matcnt[i])
					{
						matrix_sum_matind[matrix_sum_nnz]=matrix_A.matind[matrix_A.matbeg[i]+A_cnt];
						matrix_sum_matval[matrix_sum_nnz]=matrix_A.matval[matrix_A.matbeg[i]+A_cnt];
						matrix_sum_nnz++;
						cnt++;
						A_cnt++;
					}else
					{
						matrix_sum_matind[matrix_sum_nnz]=matrix_B.matind[matrix_B.matbeg[i]+B_cnt];
						matrix_sum_matval[matrix_sum_nnz]=matrix_B.matval[matrix_B.matbeg[i]+B_cnt];
						matrix_sum_nnz++;
						cnt++;
						B_cnt++;
					}
				}
				if(i==0)matrix_sum_matbeg[i]=0;
				else matrix_sum_matbeg[i]=matrix_sum_matbeg[i-1]+matrix_sum_matcnt[i-1];
				matrix_sum_matcnt[i]=cnt;
			}
		}else if(matrix_sum_param==1)// group by row
		{
			matrix_sum_matind=(int*)malloc((matrix_A.nnz+matrix_B.nnz)*sizeof(int));
			matrix_sum_matval=(double*)malloc((matrix_A.nnz+matrix_B.nnz)*sizeof(double));
			if (matrix_sum_matind==NULL || matrix_sum_matval==NULL)
			{
				status=NO_MEMORY;
				fprintf(stderr," matrix_addition(): unable to allocate memory\n");
				goto TERMINATE;
			}
			matrix_sum_matbeg=(int*) malloc(matrix_sum_nrow*sizeof(int));
			matrix_sum_matcnt=(int*) malloc(matrix_sum_nrow*sizeof(int));
			if (matrix_sum_matbeg==NULL || matrix_sum_matcnt==NULL)
			{
				status=NO_MEMORY;
				fprintf(stderr," matrix_addition(): unable to allocate memory\n");
				goto TERMINATE;
			}
			matrix_sum_nnz=0;
			for (i=0;i<matrix_sum_nrow;i++)
			{
				cnt=0;
				A_cnt=0;
				B_cnt=0;
				while(A_cnt<matrix_A.matcnt[i] || B_cnt<matrix_B.matcnt[i])
				{
					if (A_cnt<matrix_A.matcnt[i] && B_cnt<matrix_B.matcnt[i])
					{
						if (matrix_A.matind[matrix_A.matbeg[i]+A_cnt]<matrix_B.matind[matrix_B.matbeg[i]+B_cnt])
						{
							matrix_sum_matind[matrix_sum_nnz]=matrix_A.matind[matrix_A.matbeg[i]+A_cnt];
							matrix_sum_matval[matrix_sum_nnz]=matrix_A.matval[matrix_A.matbeg[i]+A_cnt];
							matrix_sum_nnz++;
							cnt++;
							A_cnt++;
						}else if (matrix_A.matind[matrix_A.matbeg[i]+A_cnt]>matrix_B.matind[matrix_B.matbeg[i]+B_cnt])
						{
							matrix_sum_matind[matrix_sum_nnz]=matrix_B.matind[matrix_B.matbeg[i]+B_cnt];
							matrix_sum_matval[matrix_sum_nnz]=matrix_B.matval[matrix_B.matbeg[i]+B_cnt];
							matrix_sum_nnz++;
							cnt++;
							B_cnt++;
						}else
						{
							if (fabs(matrix_A.matval[matrix_A.matbeg[i]+A_cnt]+matrix_B.matval[matrix_B.matbeg[i]+B_cnt])>ZERO_TOLERANCE)
							{
								matrix_sum_matind[matrix_sum_nnz]=matrix_A.matind[matrix_A.matbeg[i]+A_cnt];
								matrix_sum_matval[matrix_sum_nnz]=matrix_A.matval[matrix_A.matbeg[i]+A_cnt]+matrix_B.matval[matrix_B.matbeg[i]+B_cnt];
								matrix_sum_nnz++;
								cnt++;
							}
							A_cnt++;
							B_cnt++;
						}
					}else if (A_cnt<matrix_A.matcnt[i])
					{
						matrix_sum_matind[matrix_sum_nnz]=matrix_A.matind[matrix_A.matbeg[i]+A_cnt];
						matrix_sum_matval[matrix_sum_nnz]=matrix_A.matval[matrix_A.matbeg[i]+A_cnt];
						matrix_sum_nnz++;
						cnt++;
						A_cnt++;
					}else
					{
						matrix_sum_matind[matrix_sum_nnz]=matrix_B.matind[matrix_B.matbeg[i]+B_cnt];
						matrix_sum_matval[matrix_sum_nnz]=matrix_B.matval[matrix_B.matbeg[i]+B_cnt];
						matrix_sum_nnz++;
						cnt++;
						B_cnt++;
					}
				}
				if(i==0)matrix_sum_matbeg[i]=0;
				else matrix_sum_matbeg[i]=matrix_sum_matbeg[i-1]+matrix_sum_matcnt[i-1];
				matrix_sum_matcnt[i]=cnt;
			}
		}else
		{
			status=-1;
			fprintf(stderr," matrix_addition(): unexpected matrix_param:%d\n",matrix_sum_param);
			goto TERMINATE;
		}
		//set matrix_sum
		free_matrix(matrix_sum);
		matrix_sum_matind=(int*) realloc(matrix_sum_matind,matrix_sum_nnz*sizeof(int));
		matrix_sum_matval=(double*) realloc(matrix_sum_matval,matrix_sum_nnz*sizeof(double));
		matrix_sum->matbeg=matrix_sum_matbeg;
		matrix_sum->matcnt=matrix_sum_matcnt;
		matrix_sum->matind=matrix_sum_matind;
		matrix_sum->matval=matrix_sum_matval;
		matrix_sum->n_col=matrix_sum_nrow;
		matrix_sum->n_row=matrix_sum_ncol;
		matrix_sum->nnz=matrix_sum_nnz;
		matrix_sum->param=matrix_sum_param;
	}

TERMINATE:
	return(status);
}

int test_matrix_addition_symmetrize()
{
	int status=0;
	int i,j;
	int sym_status;
	double **matrix_A=NULL;
	double **matrix_B=NULL;
	MATRIX m_A;
	MATRIX m_B;
	MATRIX m_AB_sum;

	matrix_A=(double**)malloc(4*sizeof(double*));
	matrix_B=(double**)malloc(4*sizeof(double*));
	for (i=0;i<4;i++)
	{
		matrix_A[i]=(double*)malloc(4*sizeof(double));
		matrix_B[i]=(double*)malloc(4*sizeof(double));
	}
	matrix_A[0][0]=0;
	matrix_A[0][1]=43;
	matrix_A[0][2]=4;
	matrix_A[0][3]=3;
	matrix_A[1][0]=0;
	matrix_A[1][1]=0;
	matrix_A[1][2]=0;
	matrix_A[1][3]=0;
	matrix_A[2][0]=4;
	matrix_A[2][1]=0;
	matrix_A[2][2]=3;
	matrix_A[2][3]=0;
	matrix_A[3][0]=0;
	matrix_A[3][1]=0;
	matrix_A[3][2]=0;
	matrix_A[3][3]=1;


	matrix_B[0][0]=1;
	matrix_B[0][1]=3;
	matrix_B[0][2]=0;
	matrix_B[0][3]=2;
	matrix_B[1][0]=3;
	matrix_B[1][1]=0;
	matrix_B[1][2]=0;
	matrix_B[1][3]=0;
	matrix_B[2][0]=0;
	matrix_B[2][1]=0;
	matrix_B[2][2]=0;
	matrix_B[2][3]=3;
	matrix_B[3][0]=0;
	matrix_B[3][1]=4;
	matrix_B[3][2]=0;
	matrix_B[3][3]=1;
	fprintf(stdout,"matrix_A:\n");
	for (i=0;i<4;i++)
	{
		for (j=0;j<4;j++)
		{
			fprintf(stdout,"%f ",matrix_A[i][j]);
		}
		fprintf(stdout,"\n");
	}
	fprintf(stdout,"matrix_B:\n");
	for (i=0;i<4;i++)
	{
		for (j=0;j<4;j++)
		{
			fprintf(stdout,"%f ",matrix_B[i][j]);
		}
		fprintf(stdout,"\n");
	}
	init_matrix(&m_A);
	init_matrix(&m_B);
	init_matrix(&m_AB_sum);
	status=get_matrix(matrix_A,4,4,0,&m_A);
	if(status) goto TERMINATE;
	status=get_matrix(matrix_B,4,4,0,&m_B);
	if(status) goto TERMINATE;
	status=print_full_matrix(&m_A);
	if(status) goto TERMINATE;
	status=print_full_matrix(&m_B);
	if(status) goto TERMINATE;
	status=matrix_addition(m_A,m_B,&m_AB_sum);
	if(status) goto TERMINATE;
	status=print_full_matrix(&m_AB_sum);
	if(status) goto TERMINATE;
	status=matrix_addition(m_A,m_AB_sum,&m_AB_sum);
	if(status) goto TERMINATE;
	status=print_full_matrix(&m_AB_sum);
	if(status) goto TERMINATE;
	status=symmetrizing_matrix(&m_AB_sum,&m_AB_sum,&sym_status);
	if(status) goto TERMINATE;
	status=print_full_matrix(&m_AB_sum);
	if(status) goto TERMINATE;
TERMINATE:
	free_matrix(&m_A);
	free_matrix(&m_B);
	free_matrix(&m_AB_sum);
	for (i=0;i<4;i++)
	{
		free_and_null((char**)&(matrix_A[i]));
		free_and_null((char**)&(matrix_B[i]));
	}
	free_and_null((char**)&matrix_A);
	free_and_null((char**)&matrix_B);
	return(status);
}

void vector_sum(double *vector_1,
	double *vector_2,
	double cnt,
	double *sum_vector)
{
	int i;
	double coef_1;
	double coef_2;
	for (i=0;i<cnt;i++)
	{
		if (vector_1!=NULL) coef_1=vector_1[i];
		else coef_1=0;
		if (vector_2!=NULL) coef_2=vector_2[i];
		else coef_2=0;
		sum_vector[i]=coef_1+coef_2;
	}
}

void init_branchhist ( BRANCHHIST *ptr) //done
{
	if (ptr != NULL)
	{
		ptr->branch_var_cnt=0;
		ptr->average_obj_improve_1=NULL;
		ptr->average_obj_improve_2=NULL;
		ptr->average_violation_1_cnt=NULL;
		ptr->average_violation_2_cnt=NULL;
		ptr->branch_1_cnt=NULL;
		ptr->branch_2_cnt=NULL;
		//ptr->branch_problem_cnt_1=NULL;
		//ptr->branch_problem_cnt_2=NULL;
	}
}

void free_branchhist ( BRANCHHIST *ptr) //done 
{
	if (ptr != NULL)
	{
		ptr->branch_var_cnt=0;
		free_and_null ((char **) &(ptr->average_obj_improve_1));
		free_and_null ((char **) &(ptr->average_obj_improve_2));
		free_and_null ((char **) &(ptr->average_violation_1_cnt));
		free_and_null ((char **) &(ptr->average_violation_2_cnt));
		free_and_null ((char **) &(ptr->branch_1_cnt));
		free_and_null ((char **) &(ptr->branch_2_cnt));
		//free_and_null ((char **) &(ptr->branch_problem_cnt_1));
		//free_and_null ((char **) &(ptr->branch_problem_cnt_2));
	}
}

int init_branchhist_info(BRANCHHIST *ptr, int branch_var_cnt)
{
	int status=0;
	if (ptr != NULL)
	{
		ptr->branch_var_cnt=branch_var_cnt;
		ptr->average_obj_improve_1=(double*)calloc(branch_var_cnt,sizeof(double));
		ptr->average_obj_improve_2=(double*)calloc(branch_var_cnt,sizeof(double));
		ptr->average_violation_1_cnt=(double*)calloc(branch_var_cnt,sizeof(double));
		ptr->average_violation_2_cnt=(double*)calloc(branch_var_cnt,sizeof(double));
		ptr->branch_1_cnt=(long int*)calloc(branch_var_cnt,sizeof(long int));
		ptr->branch_2_cnt=(long int*)calloc(branch_var_cnt,sizeof(long int));
		//ptr->branch_problem_cnt_1=(long int*)calloc(branch_var_cnt,sizeof(long int));
		//ptr->branch_problem_cnt_2=(long int*)calloc(branch_var_cnt,sizeof(long int));
		if (ptr->average_obj_improve_1==NULL 
			|| ptr->average_obj_improve_2==NULL
			|| ptr->average_violation_1_cnt==NULL
			|| ptr->average_violation_2_cnt==NULL
			|| ptr->branch_1_cnt==NULL
			|| ptr->branch_2_cnt==NULL)
			//|| ptr->branch_problem_cnt_1==NULL
			//|| ptr->branch_problem_cnt_2==NULL)
		{
			status=NO_MEMORY;
			fprintf(stderr,"init_branchhist_info(): unable to allocate memory\n");
			goto TERMINATE;
		}
	}
TERMINATE:
	return(status);
}

int get_branchhist_info(BRANCHHIST *ptr,int index, int direction, double *violation_cost, double *improve_cost)
{
	int status=0;
	int i;
	int cnt;
	double average_violation_cost;
	double average_improve_cost;
	if (ptr !=NULL)
	{
		if (index>=ptr->branch_var_cnt)
		{
			status=-1;
			goto TERMINATE;
		}
		switch (direction)
		{
		case 1 : 
			if(ptr->branch_1_cnt[index]<=8)
			{
				average_violation_cost=0;
				average_improve_cost=0;
				cnt=0;
				for(i=0;i<ptr->branch_var_cnt;i++)
				{
					if (ptr->branch_1_cnt[i]!=0)
					{
						average_violation_cost+=ptr->average_violation_1_cnt[i]/ptr->branch_1_cnt[i];
						average_improve_cost+=ptr->average_obj_improve_1[i]/ptr->branch_1_cnt[i];
						cnt++;
					}
				}
				if (cnt==0)
				{
					*violation_cost=1;
					*improve_cost=1;
				}else
				{
					average_improve_cost/=cnt;
					average_violation_cost/=cnt;
					*violation_cost=average_violation_cost;
					*improve_cost=average_improve_cost;
				}
			}else
			{
				*violation_cost=(ptr->average_violation_1_cnt[index])/ptr->branch_1_cnt[index];		
				*improve_cost=(ptr->average_obj_improve_1[index])/ptr->branch_1_cnt[index];
			}
			break;
		case 2 : 
			if(ptr->branch_2_cnt[index]<=8)
			{
				average_violation_cost=0;
				average_improve_cost=0;
				cnt=0;
				for(i=0;i<ptr->branch_var_cnt;i++)
				{
					if (ptr->branch_2_cnt[i]!=0)
					{
						average_violation_cost+=ptr->average_violation_2_cnt[i]/ptr->branch_2_cnt[i];
						average_improve_cost+=ptr->average_obj_improve_2[i]/ptr->branch_2_cnt[i];
						cnt++;
					}
				}
				if (cnt==0)
				{
					*violation_cost=1;
					*improve_cost=1;
				}else
				{
					average_improve_cost/=cnt;
					average_violation_cost/=cnt;
					*violation_cost=average_violation_cost;
					*improve_cost=average_improve_cost;
				}
			}else
			{
				*violation_cost=ptr->average_violation_2_cnt[index]/ptr->branch_2_cnt[index];		
				*improve_cost=ptr->average_obj_improve_2[index]/ptr->branch_2_cnt[index];	
			}
			break;
		default:
			status=-1;
			goto TERMINATE;
		}
	}else
	{
		status=-1;
		goto TERMINATE;
	}
	*improve_cost=max(*improve_cost,ZERO_TOLERANCE);
TERMINATE:
	return(status);
}
int update_branchhist_info(BRANCHHIST *ptr, int index, int direction, double average_violation, double obj_improvement)
{
	int status=0;
	int period=5;
	double ratio=2/(period+1);
	//fprintf(stdout,"update hist: %f, %f\n",average_violation,obj_improvement);
	if (ptr!=NULL)
	{
		if (index>=ptr->branch_var_cnt)
		{
			status=-1;
			goto TERMINATE;
		}
		switch (direction)
		{
		case 1: 
			ptr->branch_1_cnt[index]++;
			ptr->average_violation_1_cnt[index]+=average_violation;
			ptr->average_obj_improve_1[index]+=obj_improvement;
			break;
		case 2: 
			ptr->branch_2_cnt[index]++;
			ptr->average_violation_2_cnt[index]+=average_violation;
			ptr->average_obj_improve_2[index]+=obj_improvement;
			break;
		default:
			status=-1;
			goto TERMINATE;
		}
	}else
	{
		status=-1;
		goto TERMINATE;
	}
TERMINATE:
	return(status);
}

//void addNodeBranchInfo(BRANCHHIST *ptr,NODE *node_ptr)
//{
//	int i;
//	for(i=0;i<node_ptr->m;i++)
//	{
//		if (node_ptr->nodebranchinfo[i]==1)
//		{
//			ptr->branch_problem_cnt_1[i]++;
//		}else if(node_ptr->nodebranchinfo[i]==2)
//		{
//			ptr->branch_problem_cnt_2[i]++;
//		}
//	}
//}
//
//void removeNodeBranchInfo(BRANCHHIST *ptr,NODE *node_ptr)
//{
//	int i;
//	for(i=0;i<node_ptr->m;i++)
//	{
//		if (node_ptr->nodebranchinfo[i]==1)
//		{
//			ptr->branch_problem_cnt_1[i]--;
//		}else if(node_ptr->nodebranchinfo[i]==2)
//		{
//			ptr->branch_problem_cnt_2[i]--;
//		}
//	}
//}

double vectorproduct(double **Vector_1, 
	double **Vector_2, 
	const int dimension, 
	const int Vector_1_index, 
	const int Vector_2_index)
{
	int i;
	double vector_prod=0;
	for (i=0;i<dimension;i++)
	{
		vector_prod+=Vector_1[Vector_1_index][i]*Vector_2[Vector_2_index][i];
	}
	return(vector_prod);
}

int partition_range(double *x_lb, 
					double *x_ub, 
					double *y_bar_lb,
					double *y_bar_ub,
					double ***partition_x_lb,
					double ***partition_x_ub,
					double ***partition_y_bar_lb,
					double ***partition_y_bar_ub,
					const int var_num,
					int *partition_x_cnt,
					int *partition_y_bar_cnt,
					int *total_partition_cnt)
{
	int status=0;
	int i,j;
	int cnt;
	int n;
	int total_cnt=1;
	double *x_range=NULL;
	double *y_bar_range=NULL;
	x_range=(double*) malloc(var_num*sizeof(double));
	y_bar_range=(double*) malloc(var_num*sizeof(double));
	for (i=0;i<var_num;i++)
	{
		x_range[i]=fabs(x_ub[i]-x_lb[i]);
		y_bar_range[i]=fabs(y_bar_ub[i]-y_bar_lb[i]);
		if (x_range[i]<ZERO_TOLERANCE) partition_x_cnt[i]=1;
		if (y_bar_range[i]<ZERO_TOLERANCE) partition_y_bar_cnt[i]=1;
		x_range[i]/=partition_x_cnt[i];
		y_bar_range[i]/=partition_y_bar_cnt[i];
		total_cnt*=partition_x_cnt[i]*partition_y_bar_cnt[i];
	}
	*total_partition_cnt=total_cnt;
	if (((*partition_x_lb)=(double**) malloc(total_cnt*sizeof(double*)))==NULL ||
		((*partition_x_ub)=(double**) malloc(total_cnt*sizeof(double*)))==NULL ||
		((*partition_y_bar_lb)=(double**) malloc(total_cnt*sizeof(double*)))==NULL ||
		((*partition_y_bar_ub)=(double**) malloc(total_cnt*sizeof(double*)))==NULL)
	{
		status=NO_MEMORY;
		fprintf(stderr,"unable to allocate memory\n");
		goto TERMIANTE;
	}
	cnt=0;
	for (i=0;i<total_cnt;i++)
	{
		if (((*partition_x_lb)[i]=(double*) malloc(var_num*sizeof(double)))==NULL ||
			((*partition_x_ub)[i]=(double*) malloc(var_num*sizeof(double)))==NULL ||
			((*partition_y_bar_lb)[i]=(double*) malloc(var_num*sizeof(double)))==NULL ||
			((*partition_y_bar_ub)[i]=(double*) malloc(var_num*sizeof(double)))==NULL)
		{
			status=NO_MEMORY;
			fprintf(stderr,"unable to allocate memory\n");
			goto TERMIANTE;
		}
		n=i;
		for (j=0;j<var_num;j++)
		{
			(*partition_x_lb)[i][j]=x_lb[j]+fmod(n,partition_x_cnt[j])*x_range[j];
			(*partition_x_ub)[i][j]=(*partition_x_lb)[i][j]+x_range[j];
			n/=partition_x_cnt[j];
			(*partition_y_bar_lb)[i][j]=y_bar_lb[j]+fmod(n,partition_y_bar_cnt[j])*y_bar_range[j];
			(*partition_y_bar_ub)[i][j]=(*partition_y_bar_lb)[i][j]+y_bar_range[j];
			n/=partition_y_bar_cnt[j];
		}
	}

TERMIANTE:
	free_and_null((char**) &x_range);
	free_and_null((char**) &y_bar_range);
	return(status);
}

int readProfile(PARAM *parameter)
{
	int status=0;
	char buffer[1024];
	char *pch=NULL;
	FILE *in=NULL;
	fprintf(stdout,"Reading profile...\n");
	in=fopen(parameter->profilefile,"r");
	if (in ==NULL)
	{
		status = -1;
		fprintf(stderr," unable to open the profile: %s\n",parameter->profilefile);
		goto TERMINATE;
	}else
	{
		while(!feof(in))
		{
			if (fgets(buffer,1024,in)!=NULL)
			{
				pch=strtok(buffer," =#/");
				if (strcmp(pch,"FEASIBILITY_RECOVERY_BREATH")==0)
				{
					//fprintf(stdout,"FEASIBILITY_RECOVERY_BREATH");
					//fprintf(stdout,"=");
					pch=strtok(NULL," =#/");
					//fprintf(stdout,"%s\n",pch);
					parameter->feasibility_recovery_breath=atof(pch);
					//fprintf(stdout,"%f\n",parameter.feasibility_recovery_breath);
				}
				if (strcmp(pch,"FEASIBILITY_RECOVERY_DEPTH")==0)
				{
					//fprintf(stdout,"FEASIBILITY_RECOVERY_DEPTH");
					//fprintf(stdout,"=");
					pch=strtok(NULL," =#/");
					//fprintf(stdout,"%s\n",pch);
					parameter->feasibility_recovery_depth=atoi(pch);
					//fprintf(stdout,"%d\n",parameter.feasibility_recovery_depth);
				}
				if (strcmp(pch,"FEASIBILITY_RECOVERY_REFINEMENT")==0)
				{
					//fprintf(stdout,"FEASIBILITY_RECOVERY_REFINEMENT");
					//fprintf(stdout,"=");
					pch=strtok(NULL," =#/");
					//fprintf(stdout,"%s\n",pch);
					parameter->feasibility_recovery_refinement=atoi(pch);
					//fprintf(stdout,"%d\n",parameter.feasibility_recovery_refinement);
				}
				if (strcmp(pch,"FEASIBILITY_RECOVERY_MAXLOOP")==0)
				{
					//fprintf(stdout,"FEASIBILITY_RECOVERY_MAXLOOP");
					//fprintf(stdout,"=");
					pch=strtok(NULL," =#/");
					//fprintf(stdout,"%s\n",pch);
					parameter->feasibility_recovery_maxloop=atoi(pch);
					//fprintf(stdout,"%d\n",parameter.feasibility_recovery_maxloop);
				}
				if (strcmp(pch,"BOUNDCUT_MAXLOOP")==0)
				{
					//fprintf(stdout,"BOUNDCUT_MAXLOOP");
					//fprintf(stdout,"=");
					pch=strtok(NULL," =#/");
					//fprintf(stdout,"%s\n",pch);
					parameter->boundcut_loop=atoi(pch);
					//fprintf(stdout,"%d\n",parameter.boundcut_loop);
				}
				if (strcmp(pch,"BOUNDCUT_P")==0)
				{
					//fprintf(stdout,"BOUNDCUT_P");
					//fprintf(stdout,"=");
					pch=strtok(NULL," =#/");
					//fprintf(stdout,"%s\n",pch);
					parameter->boundcut_p=atoi(pch);
					//fprintf(stdout,"%d\n",parameter.boundcut_p);
				}
				if (strcmp(pch,"DISJUNCTIVECUT_MAXLOOP")==0)
				{
					//fprintf(stdout,"DISJUNCTIVECUT_MAXLOOP");
					//fprintf(stdout,"=");
					pch=strtok(NULL," =#/");
					//fprintf(stdout,"%s\n",pch);
					parameter->disjunctivecut_loop_param=atof(pch);
					//fprintf(stdout,"%f\n",parameter.disjunctivecut_loop_param);
				}
				if (strcmp(pch,"SIMPLECUT_MAXLOOP")==0)
				{
					//fprintf(stdout,"SIMPLECUT_MAXLOOP");
					//fprintf(stdout,"=");
					pch=strtok(NULL," =#/");
					//fprintf(stdout,"%s\n",pch);
					parameter->simplecut_loop_param=atof(pch);
					//fprintf(stdout,"%f\n",parameter.simplecut_loop_param);
				}
				if (strcmp(pch,"MCCORMICK_REFINEMENT_CNT")==0)
				{
					//fprintf(stdout,"MCCORMICK_REFINEMENT_CNT");
					//fprintf(stdout,"=");
					pch=strtok(NULL," =#/");
					//fprintf(stdout,"%s\n",pch);
					parameter->MccormickRefineCnt=atoi(pch);
					//fprintf(stdout,"%d\n",parameter.MccormickRefineCnt);
				}
				if (strcmp(pch,"TIME_LIMIT")==0)
				{
					//fprintf(stdout,"TIME_LIMIT");
					//fprintf(stdout,"=");
					pch=strtok(NULL," =#/");
					//fprintf(stdout,"%s\n",pch);
					parameter->time_limit=atof(pch);
					//fprintf(stdout,"%f\n",parameter.time_limit);
				}
				if (strcmp(pch,"OPTIMAL_TOLERANCE")==0)
				{
					//fprintf(stdout,"OPTIMAL_TOLERANCE");
					//fprintf(stdout,"=");
					pch=strtok(NULL," =#/");
					//fprintf(stdout,"%s\n",pch);
					parameter->opt_tolerance=atof(pch);
					//fprintf(stdout,"%12.12f\n",parameter.opt_tolerance);
				}
				if (strcmp(pch,"PARTITION_X")==0)
				{
					//fprintf(stdout,"PARTITION_X");
					//fprintf(stdout,"=");
					pch=strtok(NULL," =#/");
					//fprintf(stdout,"%s\n",pch);
					parameter->partition_x=atoi(pch);
					//fprintf(stdout,"%d\n",parameter.partition_x);
				}
				if (strcmp(pch,"PARTITION_Y")==0)
				{
					//fprintf(stdout,"PARTITION_Y");
					//fprintf(stdout,"=");
					pch=strtok(NULL," =#/");
					//fprintf(stdout,"%s\n",pch);
					parameter->partition_y=atoi(pch);
					//fprintf(stdout,"%d\n",parameter.partition_y);
				}
			}
			else
			{
				status=-1;
				fprintf(stderr,"unable to read the profile\n");
				goto TERMINATE;
			}
		}
		fclose(in);
	}
TERMINATE:
	return(status);
}

//
void parse_command_line(int argc, char **argv, PARAM *parameter)
{
	int i;
	// default values
	parameter->boundcut_loop=DEFAULT_BOUNDCUT_MAXLOOP;
	parameter->boundcut_p=DEFAULT_BOUNDCUT_P;
	parameter->disjunctivecut_loop_param=DEFAULT_DISJUNCTIVECUT_MAXLOOP;
	parameter->feasibility_recovery_breath=DEFAULT_FEASIBILITY_RECOVERY_BREATH;
	parameter->feasibility_recovery_depth=DEFAULT_FEASIBILITY_RECOVERY_DEPTH;
	parameter->feasibility_recovery_maxloop=DEFAULT_FEASIBILITY_RECOVERY_MAXLOOP;
	parameter->feasibility_recovery_refinement=DEFAULT_FEASIBILITY_RECOVERY_REFINEMENT;
	parameter->time_limit=DEFAULT_TIME_LIMIT;
	parameter->MccormickRefineCnt=DEFAULT_MCCORMICK_REFINEMENT_CNT;
	parameter->opt_tolerance=DEFAULT_OPTIMAL_TOLERANCE;
	parameter->simplecut_loop_param=DEFAULT_SIMPLECUT_MAXLOOP;
	parameter->partition_x=DEFAULT_PARTITION_X;
	parameter->partition_y=DEFAULT_PARTITION_Y;
	parameter->format_mode=1;
	parameter->solve_mode='c';
	// parse options
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break;
		if(++i>=argc)
			exit_with_help();
		switch(argv[i-1][1])
		{
		case 'm':
			if (strcmp("c",argv[i])==0 || strcmp("b",argv[i])==0)
			{
				parameter->solve_mode=argv[i][0];
				break;
			}else
			{
				fprintf(stdout,"solve mode must be 'b' or 'c' \n");
				exit_with_help();
			}
		case 'f':
			parameter->format_mode=atoi(argv[i]);
			if (parameter->format_mode<1|| parameter->format_mode>4)
			{
				fprintf(stdout,"file format must be 1-4' \n");
				exit_with_help();
			}		
			break;
		case 'p':
			strcpy(parameter->profilefile,argv[i]);
			readProfile(parameter);
			break;
		default:
			fprintf(stdout,"Unknown option: -%c\n", argv[i-1][1]);
			exit_with_help();
		}
	}
	if(i>=argc)
		exit_with_help();

	strcpy(parameter->inputfile, argv[i]);

	if(i<argc-1)
		strcpy(parameter->solnfile,argv[i+1]);
	else
	{
		sprintf(parameter->solnfile,"%s.out",parameter->inputfile);
	}
}

void exit_with_help()
{
	printf(
		"Usage:lpcc [options] inputfile [outputfile]\n"
		"options:\n"
		"-m solve mode:\n"
		"	c -- use cplex to solve MIP formulation of LPCC\n"
		"	b -- use b&c to solve LPCC (Default)\n"
		"-f inputfile_format :\n"
		"	1 -- full matrix format in separate files\n"
		"	2 -- full matrix format in single file\n"
		"	3 -- compact matrix format in single FILE (Default)\n"
		"	4 -- .mps file from Binary Integer Program\n"	
		"-p profile_name\n"
		);
	exit(1);
}

void printParam(PARAM parameter)
{
	fprintf(stdout,"print parameter content:\n");
	fprintf(stdout,"%f\n",parameter.feasibility_recovery_breath);
	fprintf(stdout,"%d\n",parameter.feasibility_recovery_depth);
	fprintf(stdout,"%d\n",parameter.feasibility_recovery_refinement);
	fprintf(stdout,"%d\n",parameter.feasibility_recovery_maxloop);
	fprintf(stdout,"%d\n",parameter.boundcut_loop);
	fprintf(stdout,"%d\n",parameter.boundcut_p);
	fprintf(stdout,"%f\n",parameter.disjunctivecut_loop_param);
	fprintf(stdout,"%f\n",parameter.simplecut_loop_param);
	fprintf(stdout,"%d\n",parameter.MccormickRefineCnt);
	fprintf(stdout,"%c\n",parameter.solve_mode);
	fprintf(stdout,"%d\n",parameter.format_mode);
	fprintf(stdout,"%s\n",parameter.profilefile);
	fprintf(stdout,"%s\n",parameter.inputfile);
	fprintf(stdout,"%s\n",parameter.solnfile);
	fprintf(stdout,"%f\n",parameter.time_limit);
	fprintf(stdout,"%12.12f\n",parameter.opt_tolerance);
	fprintf(stdout,"%d\n",parameter.partition_x);
	fprintf(stdout,"%d\n",parameter.partition_y);
}

