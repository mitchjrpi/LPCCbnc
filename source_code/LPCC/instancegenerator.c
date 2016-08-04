#define _CRT_SECURE_NO_DEPRECATE
#include <stdlib.h>
// #include <direct.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <LPCC.h>
#include <instancegenerator.h>
#include <utilities.h>
#include <rngs.h>
#include <rvgs.h>
int generateLPCCinstance(int param_n,	
						 int param_m,	
						 int param_k, 	
						 int rankM,
						 double spars,
						 double **c_coef_p,		
						 double **d_coef_p,		
						 double **b_coef_p, 		
						 double **q_coef_p,  		
						 double ***matrix_A_p,			
						 double ***matrix_B_p,		
					     double ***matrix_N_p,
						 double ***matrix_M_p)
{
	int status = 0;
	int i,j,k;
	double Mt_ij;
	double** matrix_L=NULL;
	double** Mdelta=NULL;
	double *x_bar=NULL;
	double *y_bar=NULL;
	matrix_L=(double**) calloc(param_m,sizeof(double*));
	Mdelta=(double**)calloc(param_m,sizeof(double*));
	x_bar=(double*)calloc(param_n,sizeof(double));
	y_bar=(double*)calloc(param_m,sizeof(double));
	*c_coef_p=(double*)calloc(param_n,sizeof(double));
	*d_coef_p=(double*)calloc(param_m,sizeof(double));
	*b_coef_p=(double*)calloc(param_k,sizeof(double));
	*q_coef_p=(double*)calloc(param_m,sizeof(double));
	*matrix_A_p=(double**)calloc(param_k,sizeof(double*));
	*matrix_B_p=(double**)calloc(param_k,sizeof(double*));
	*matrix_N_p=(double**)calloc(param_m,sizeof(double*));
	*matrix_M_p=(double**)calloc(param_m,sizeof(double*));
	if (matrix_L==NULL || 
		Mdelta==NULL || 
		x_bar==NULL || 
		y_bar==NULL || 
		matrix_L==NULL || 
		*c_coef_p==NULL || 
		*d_coef_p==NULL || 
		*b_coef_p==NULL || 
		*q_coef_p==NULL || 
		*matrix_A_p==NULL || 
		*matrix_B_p==NULL || 
		*matrix_N_p==NULL || 
		*matrix_M_p==NULL)
	{
		status = NO_MEMORY;
		fprintf(stderr," Unable to allocate enough memory for the input data...\n");
		goto TERMINATE;
	}
	for (i=0;i<param_k;i++)
	{
		(*matrix_A_p)[i]=(double*)calloc(param_n,sizeof(double));
		(*matrix_B_p)[i]=(double*)calloc(param_m,sizeof(double));
		if ((*matrix_A_p)[i]==NULL || (*matrix_B_p)[i]==NULL)
		{
			status = NO_MEMORY;
			fprintf(stderr," Unable to allocate enough memory for the input data...\n");
			goto TERMINATE;
		}
	}
	for (i=0;i<param_m;i++)
	{
		(*matrix_N_p)[i]=(double*)calloc(param_n,sizeof(double));
		(*matrix_M_p)[i]=(double*)calloc(param_m,sizeof(double));
		matrix_L[i]=(double*)calloc(rankM,sizeof(double));
		Mdelta[i]=(double*)calloc(param_m,sizeof(double));
		if ((*matrix_N_p)[i]==NULL || (*matrix_M_p)[i]==NULL || matrix_L[i]==NULL || Mdelta[i]==NULL )
		{
			status = NO_MEMORY;
			fprintf(stderr," Unable to allocate enough memory for the input data...\n");
			goto TERMINATE;
		}
	}
	//x_bar
	for (i=0;i<param_n;i++)
	{
		x_bar[i]=floor(Uniform(0,10));
	}
	//y_bar
	for (i=0;i<param_m;i++)
	{
		if (i<param_m/3)
		{
			y_bar[i]=floor(Uniform(0,10));
		}else
		{
			y_bar[i]=0;
		}
	}
	/*c_coef*/
	for (i=0;i<param_n;i++)
	{
		(*c_coef_p)[i]=floor(Uniform(0,10));
	}
	/*d_coef*/
	for (i=0;i<param_m;i++)
	{
		(*d_coef_p)[i]=floor(Uniform(0,10));
	}
	/*matrix_A*/
	for (i=0;i<param_k;i++)
	{
		for (j=0;j<param_n;j++)
		{
			if (Uniform(0,1)<spars)
			{
				(*matrix_A_p)[i][j]=floor(Uniform(-5,6));
			}else
			{
				(*matrix_A_p)[i][j]=0;
			}			
		}
	}
	/*matrix_B*/
	for (i=0;i<param_k;i++)
	{
		for (j=0;j<param_m;j++)
		{
			
			if (Uniform(0,1)<spars)
			{
				(*matrix_B_p)[i][j]=floor(Uniform(-5,6));
			} 
			else
			{
				(*matrix_B_p)[i][j]=0;
			}
		}
	}
	/*matrix_N*/
	for (i=0;i<param_m;i++)
	{
		for (j=0;j<param_n;j++)
		{
			if (Uniform(0,1)<spars)
			{
				(*matrix_N_p)[i][j]=floor(Uniform(-5,6));
			} 
			else
			{
				(*matrix_N_p)[i][j]=0;
			}
		}
	}
	/*matrix_L*/
	for (i=0;i<param_m;i++)
	{
		for (j=0;j<rankM;j++)
		{
			if (Uniform(0,1)<spars)
			{
				matrix_L[i][j]=floor(Uniform(-5,6));
			} 
			else
			{
				matrix_L[i][j]=0;
			}
		}
	}
	/*matrix_M*/
	for (i=0;i<param_m;i++)
	{
		for(j=0;j<param_m;j++)
		{
			Mt_ij=0;
			Mdelta[i][j]=floor(Uniform(-2,3));
			for (k=0;k<rankM;k++)
			{
				Mt_ij+=matrix_L[i][k]*matrix_L[j][k];
			}
			if (fabs(Mt_ij)<1e-6)
			{
				(*matrix_M_p)[i][j]=0;
			}else if (i==j)
			{
				(*matrix_M_p)[i][j]=Mt_ij;
			} 
			else if (i<j)
			{
				(*matrix_M_p)[i][j]=Mt_ij+Mdelta[i][j];
			} 
			else
			{
				(*matrix_M_p)[i][j]=Mt_ij-Mdelta[j][i];
			}
		}
	}
	/*b_coef*/
	for (i=0;i<param_k;i++)
	{
		(*b_coef_p)[i]=0;
		for (j=0;j<param_n;j++)
		{
			(*b_coef_p)[i]+=(*matrix_A_p)[i][j]*x_bar[j];
		}
		for (j=0;j<param_m;j++)
		{
			(*b_coef_p)[i]+=(*matrix_B_p)[i][j]*y_bar[j];
		}
		(*b_coef_p)[i]-=floor(Uniform(1,11));
	}
	/*q_coef*/
	for (i=0;i<param_m;i++)
	{
		(*q_coef_p)[i]=0;
		for (j=0;j<param_n;j++)
		{
			(*q_coef_p)[i]-=(*matrix_N_p)[i][j]*x_bar[j];
		}
		for (j=0;j<param_m;j++)
		{
			(*q_coef_p)[i]-=(*matrix_M_p)[i][j]*y_bar[j];
		}
		if (i>=(2*param_m)/3)
		{
			(*q_coef_p)[i]+=floor(Uniform(1,11));
		}
	}

TERMINATE:
	if (matrix_L!=NULL)
	{
		for (i=0;i<param_m;i++)
		{
			free_and_null((char**)&(matrix_L[i]));
		}
	}
	if (Mdelta!=NULL)
	{
		for (i=0;i<param_m;i++)
		{
			free_and_null((char**)&(Mdelta[i]));
		}
	}
	free_and_null((char**) &matrix_L);
	free_and_null((char**) &Mdelta);
	free_and_null((char**) &x_bar);
	free_and_null((char**) &y_bar);
	return (status);
}
//write down random generated LPCC instance in the following format
//min c*x+d*y
//st. A*x+B*y>=b
//     N*x+M*y+q=w
//	   y>=0;w>=0;
//	   x>=0;
//	   y comp w
//c: c.txt
//d: d.txt
//A: A.txt
//B: B.txt
//N: N.txt
//M: M.txt
//q: q.txt
//parameter setting: readme.txt
int writeLPCCinstance_format_1(int random_seed,
							   int param_n,	
							   int param_m,	
							   int param_k, 	
							   int rankM,
							   double spars,
							   double *c_coef_p,		
							   double *d_coef_p,		
							   double *b_coef_p, 		
							   double *q_coef_p,  		
							   double **matrix_A_p,			
							   double **matrix_B_p,		
							   double **matrix_N_p,
							   double **matrix_M_p)
{
	int status=0;
	int i,j;
	int int_spars=(int)(spars*100);
	FILE *pFile=NULL;
	char filedir[100];
	char buffer[100];
	strcpy(filedir,".\\");
	strcat(filedir,_itoa(random_seed,buffer,10));
	strcat(filedir,"_");
	strcat(filedir,_itoa(param_n,buffer,10));
	strcat(filedir,"_");
	strcat(filedir,_itoa(param_m,buffer,10));
	strcat(filedir,"_");
	strcat(filedir,_itoa(param_k,buffer,10));
	strcat(filedir,"_");
	strcat(filedir,_itoa(rankM,buffer,10));
	strcat(filedir,"_");
	strcat(filedir,_itoa(int_spars,buffer,10));
	_mkdir(filedir);
	strcat(filedir,"\\");
	fprintf(stdout,"writing LPCC instances:\n random seed: %d, n: %d, m: %d, k:%d, rankM: %d, spars: %1.2f\n",random_seed,param_n,param_m,param_k,rankM,spars);
	
	strcpy(buffer,filedir);
	strcat(buffer,"readme.txt");
	pFile=fopen(buffer,"w");
	fprintf(pFile,"random_seed=%d\n",random_seed);
	fprintf(pFile,"n=%d\n",param_n);
	fprintf(pFile,"m=%d\n",param_m);
	fprintf(pFile,"k=%d\n",param_k);
	fprintf(pFile,"rankM=%d\n",rankM);
	fprintf(pFile,"spars=%f\n",spars);
	fclose(pFile);

	strcpy(buffer,filedir);
	strcat(buffer,"A.txt");
	pFile=fopen(buffer,"w");
	for (i=0;i<param_k;i++)
	{
		for (j=0;j<param_n;j++)
		{
			fprintf(pFile,"%12.12f  ",matrix_A_p[i][j]);
		}
		if(i!=param_k-1)fprintf(pFile,"\n");
	}
	fclose(pFile);

	strcpy(buffer,filedir);
	strcat(buffer,"B.txt");
	pFile=fopen(buffer,"w");
	for (i=0;i<param_k;i++)
	{
		for (j=0;j<param_m;j++)
		{
			fprintf(pFile,"%12.12f  ",matrix_B_p[i][j]);
		}
		if(i!=param_k-1)fprintf(pFile,"\n");
	}
	fclose(pFile);

	strcpy(buffer,filedir);
	strcat(buffer,"c.txt");
	pFile=fopen(buffer,"w");
	for (i=0;i<param_n;i++)
	{
		fprintf(pFile,"%12.12f  ",c_coef_p[i]);
	}
	fclose(pFile);

	strcpy(buffer,filedir);
	strcat(buffer,"d.txt");
	pFile=fopen(buffer,"w");
	for (i=0;i<param_m;i++)
	{
		fprintf(pFile,"%12.12f  ",d_coef_p[i]);
	}
	fclose(pFile);

	strcpy(buffer,filedir);
	strcat(buffer,"M.txt");
	pFile=fopen(buffer,"w");
	for (i=0;i<param_m;i++)
	{
		for (j=0;j<param_m;j++)
		{
			fprintf(pFile,"%12.12f  ",matrix_M_p[i][j]);
		}
		if(i!=param_m-1)fprintf(pFile,"\n");
	}
	fclose(pFile);

	strcpy(buffer,filedir);
	strcat(buffer,"N.txt");
	pFile=fopen(buffer,"w");
	for (i=0;i<param_m;i++)
	{
		for (j=0;j<param_n;j++)
		{
			fprintf(pFile,"%12.12f  ",matrix_N_p[i][j]);
		}
		if(i!=param_m-1)fprintf(pFile,"\n");
	}
	fclose(pFile);

	strcpy(buffer,filedir);
	strcat(buffer,"q.txt");
	pFile=fopen(buffer,"w");
	for (i=0;i<param_m;i++)
	{
		if(i!=param_m-1)fprintf(pFile,"%12.12f\n",q_coef_p[i]);
		else fprintf(pFile,"%12.12f",q_coef_p[i]);
	}
	fclose(pFile);

	strcpy(buffer,filedir);
	strcat(buffer,"f.txt");
	pFile=fopen(buffer,"w");
	for (i=0;i<param_k;i++)
	{
		if (i!=param_k-1)fprintf(pFile,"%12.12f\n",b_coef_p[i]);
		else fprintf(pFile,"%12.12f",b_coef_p[i]);
	}
	fclose(pFile);
	return(status);
}

//write down random generated LPCC instance in the following format:
//
//in full matrix format

int writeLPCCinstance_format_2(int random_seed,
							   int param_n,	
							   int param_m,	
							   int param_k, 	
							   int rankM,
							   double spars,
							   double *c_coef_p,		
							   double *d_coef_p,		
							   double *b_coef_p, 		
							   double *q_coef_p,  		
							   double **matrix_A_p,			
							   double **matrix_B_p,		
							   double **matrix_N_p,
							   double **matrix_M_p)
{
	int status=0;
	int i,j;
	int int_spars=(int)(spars*100);
	FILE *pFile=NULL;
	char filename[100];
	char buffer[100];
	strcpy(filename,"input_full_");
	strcat(filename,_itoa(random_seed,buffer,10));
	strcat(filename,"_");
	strcat(filename,_itoa(param_n,buffer,10));
	strcat(filename,"_");
	strcat(filename,_itoa(param_m,buffer,10));
	strcat(filename,"_");
	strcat(filename,_itoa(param_k,buffer,10));
	strcat(filename,"_");
	strcat(filename,_itoa(rankM,buffer,10));
	strcat(filename,"_");
	strcat(filename,_itoa(int_spars,buffer,10));
	strcat(filename,".dat");
	fprintf(stdout,"writing LPCC instances:\n random seed: %d, n: %d, m: %d, k:%d, rankM: %d, spars: %1.2f\n",random_seed,param_n,param_m,param_k,rankM,spars);

	pFile=fopen(filename,"w");
	//parameter
	fprintf(pFile,"[%d,%d,%d]\n",param_n,param_m,param_k);
	//c
	fprintf(pFile,"[");
	for (i=0;i<param_n;i++)
	{
		if(i==param_n-1)fprintf(pFile,"%12.12f",c_coef_p[i]);
		else fprintf(pFile,"%12.12f,",c_coef_p[i]);
	}
	fprintf(pFile,"]\n");
	//d
	fprintf(pFile,"[");
	for (i=0;i<param_m;i++)
	{
		if(i==param_m-1)fprintf(pFile,"%12.12f",d_coef_p[i]);
		else fprintf(pFile,"%12.12f,",d_coef_p[i]);
	}
	fprintf(pFile,"]\n");	
	//b
	fprintf(pFile,"[");
	for (i=0;i<param_k;i++)
	{
		if (i==param_k-1) fprintf(pFile,"%12.12f",b_coef_p[i]);
		else fprintf(pFile,"%12.12f,",b_coef_p[i]);	
	}
	fprintf(pFile,"]\n");	
	//q
	fprintf(pFile,"[");
	for (i=0;i<param_m;i++)
	{
		if(i==param_m-1)fprintf(pFile,"%12.12f",q_coef_p[i]);
		else fprintf(pFile,"%12.12f,",q_coef_p[i]);
	}
	fprintf(pFile,"]\n");	
	//A
	fprintf(pFile,"[");
	for (i=0;i<param_k;i++)
	{
		fprintf(pFile,"[");
		for (j=0;j<param_n;j++)
		{
			if(j==param_n-1)fprintf(pFile,"%12.12f  ",matrix_A_p[i][j]);
			else fprintf(pFile,"%12.12f,",matrix_A_p[i][j]);
		}
		if(i==param_k-1)fprintf(pFile,"]\n");
		else fprintf(pFile,"],\n");
	}
	fprintf(pFile,"]\n");	
	//B
	fprintf(pFile,"[");
	for (i=0;i<param_k;i++)
	{
		fprintf(pFile,"[");
		for (j=0;j<param_m;j++)
		{
			if(j==param_m-1)fprintf(pFile,"%12.12f  ",matrix_B_p[i][j]);
			else fprintf(pFile,"%12.12f,",matrix_B_p[i][j]);
		}
		if(i==param_k-1)fprintf(pFile,"]\n");
		else fprintf(pFile,"],\n");
	}
	fprintf(pFile,"]\n");	
	//N
	fprintf(pFile,"[");
	for (i=0;i<param_m;i++)
	{
		fprintf(pFile,"[");
		for (j=0;j<param_n;j++)
		{
			if(j==param_n-1)fprintf(pFile,"%12.12f  ",matrix_N_p[i][j]);
			else fprintf(pFile,"%12.12f,",matrix_N_p[i][j]);
		}
		if(i==param_m-1)fprintf(pFile,"]\n");
		else fprintf(pFile,"],\n");
	}
	fprintf(pFile,"]\n");	
	//M
	fprintf(pFile,"[");
	for (i=0;i<param_m;i++)
	{
		fprintf(pFile,"[");
		for (j=0;j<param_m;j++)
		{
			if(j==param_m-1)fprintf(pFile,"%12.12f  ",matrix_M_p[i][j]);
			else fprintf(pFile,"%12.12f,",matrix_M_p[i][j]);
		}
		if(i==param_m-1)fprintf(pFile,"]\n");
		else fprintf(pFile,"],\n");
	}
	fprintf(pFile,"]\n");	
	fclose(pFile);

	return(status);
}

//write down random generated LPCC instance in format 2
//in compact matrix format
int writeLPCCinstance_format_3(int random_seed,
							   int param_n,	
							   int param_m,	
							   int param_k, 	
							   int rankM,
							   double spars,
							   double *c_coef_p,		
							   double *d_coef_p,		
							   double *b_coef_p, 		
							   double *q_coef_p,  		
							   double **matrix_A_p,			
							   double **matrix_B_p,		
							   double **matrix_N_p,
							   double **matrix_M_p)
{
	//note: compact matrix here is always grouped by row
	int status=0;
	int i;
	int int_spars=(int)(spars*100);
	FILE *pFile=NULL;
	char filename[100];
	char buffer[100];
	MATRIX _matrix_A;
	MATRIX _matrix_B;
	MATRIX _matrix_N;
	MATRIX _matrix_M;
	init_matrix(&_matrix_A);
	init_matrix(&_matrix_B);
	init_matrix(&_matrix_N);
	init_matrix(&_matrix_M);
	//get compact matrix
	status=get_matrix(matrix_A_p,param_k,param_n,1,&_matrix_A);
	if (status) goto TERMINATE;
	status=get_matrix(matrix_B_p,param_k,param_m,1,&_matrix_B);
	if (status) goto TERMINATE;
	status=get_matrix(matrix_N_p,param_m,param_n,1,&_matrix_N);
	if (status) goto TERMINATE;
	status=get_matrix(matrix_M_p,param_m,param_m,1,&_matrix_M);
	if (status) goto TERMINATE;
	strcpy(filename,"input_compact_");
	strcat(filename,_itoa(random_seed,buffer,10));
	strcat(filename,"_");
	strcat(filename,_itoa(param_n,buffer,10));
	strcat(filename,"_");
	strcat(filename,_itoa(param_m,buffer,10));
	strcat(filename,"_");
	strcat(filename,_itoa(param_k,buffer,10));
	strcat(filename,"_");
	strcat(filename,_itoa(rankM,buffer,10));
	strcat(filename,"_");
	strcat(filename,_itoa(int_spars,buffer,10));
	strcat(filename,".dat");
	fprintf(stdout,"writing LPCC instances:\n random seed: %d, n: %d, m: %d, k:%d, rankM: %d, spars: %1.2f\n",random_seed,param_n,param_m,param_k,rankM,spars);

	pFile=fopen(filename,"w");
	//parameter
	fprintf(pFile,"[%d,%d,%d]\n",param_n,param_m,param_k);
	//c
	fprintf(pFile,"[");
	for (i=0;i<param_n;i++)
	{
		if(i==param_n-1)fprintf(pFile,"%12.12f",c_coef_p[i]);
		else fprintf(pFile,"%12.12f,",c_coef_p[i]);
	}
	fprintf(pFile,"]\n");
	//d
	fprintf(pFile,"[");
	for (i=0;i<param_m;i++)
	{
		if(i==param_m-1)fprintf(pFile,"%12.12f",d_coef_p[i]);
		else fprintf(pFile,"%12.12f,",d_coef_p[i]);
	}
	fprintf(pFile,"]\n");	
	//b
	fprintf(pFile,"[");
	for (i=0;i<param_k;i++)
	{
		if (i==param_k-1) fprintf(pFile,"%12.12f",b_coef_p[i]);
		else fprintf(pFile,"%12.12f,",b_coef_p[i]);	
	}
	fprintf(pFile,"]\n");	
	//q
	fprintf(pFile,"[");
	for (i=0;i<param_m;i++)
	{
		if(i==param_m-1)fprintf(pFile,"%12.12f",q_coef_p[i]);
		else fprintf(pFile,"%12.12f,",q_coef_p[i]);
	}
	fprintf(pFile,"]\n");	
	//A
	fprintf(pFile,"[");
	//parameter n_row n_col n_nz
	fprintf(pFile,"[%d,%d,%d],\n",_matrix_A.n_row,_matrix_A.n_col,_matrix_A.nnz);
	//matbeg n_row
	fprintf(pFile,"[");
	for (i=0;i<_matrix_A.n_row;i++)
	{
		if (i==_matrix_A.n_row-1) fprintf(pFile,"%d ",_matrix_A.matbeg[i]);
		else fprintf(pFile,"%d,",_matrix_A.matbeg[i]);
	}
	fprintf(pFile,"],\n");
	//matcnt n_row
	fprintf(pFile,"[");
	for (i=0;i<_matrix_A.n_row;i++)
	{
		if(i==_matrix_A.n_row-1) fprintf(pFile,"%d ",_matrix_A.matcnt[i]);
		else fprintf(pFile,"%d,",_matrix_A.matcnt[i]);
	}
	fprintf(pFile,"],\n");
	//matind n_nz
	fprintf(pFile,"[");
	for (i=0;i<_matrix_A.nnz;i++)
	{
		if(i==_matrix_A.nnz-1) fprintf(pFile,"%d ",_matrix_A.matind[i]);
		else fprintf(pFile,"%d,",_matrix_A.matind[i]);
	}
	fprintf(pFile,"],\n");
	//matval n_nz
	fprintf(pFile,"[");
	for (i=0;i<_matrix_A.nnz;i++)
	{
		if(i==_matrix_A.nnz-1) fprintf(pFile,"%12.12f ",_matrix_A.matval[i]);
		else fprintf(pFile,"%12.12f,",_matrix_A.matval[i]);
	}
	fprintf(pFile,"]\n");
	fprintf(pFile,"]\n");	

	//B
	fprintf(pFile,"[");
	//parameter n_row n_col n_nz
	fprintf(pFile,"[%d,%d,%d],\n",_matrix_B.n_row,_matrix_B.n_col,_matrix_B.nnz);
	//matbeg n_row
	fprintf(pFile,"[");
	for (i=0;i<_matrix_B.n_row;i++)
	{
		if (i==_matrix_B.n_row-1) fprintf(pFile,"%d ",_matrix_B.matbeg[i]);
		else fprintf(pFile,"%d,",_matrix_B.matbeg[i]);
	}
	fprintf(pFile,"],\n");
	//matcnt n_row
	fprintf(pFile,"[");
	for (i=0;i<_matrix_B.n_row;i++)
	{
		if(i==_matrix_B.n_row-1) fprintf(pFile,"%d ",_matrix_B.matcnt[i]);
		else fprintf(pFile,"%d,",_matrix_B.matcnt[i]);
	}
	fprintf(pFile,"],\n");
	//matind n_nz
	fprintf(pFile,"[");
	for (i=0;i<_matrix_B.nnz;i++)
	{
		if(i==_matrix_B.nnz-1) fprintf(pFile,"%d ",_matrix_B.matind[i]);
		else fprintf(pFile,"%d,",_matrix_B.matind[i]);
	}
	fprintf(pFile,"],\n");
	//matval n_nz
	fprintf(pFile,"[");
	for (i=0;i<_matrix_B.nnz;i++)
	{
		if(i==_matrix_B.nnz-1) fprintf(pFile,"%12.12f ",_matrix_B.matval[i]);
		else fprintf(pFile,"%12.12f,",_matrix_B.matval[i]);
	}
	fprintf(pFile,"]\n");
	fprintf(pFile,"]\n");	
	//N
	fprintf(pFile,"[");
	//parameter n_row n_col n_nz
	fprintf(pFile,"[%d,%d,%d],\n",_matrix_N.n_row,_matrix_N.n_col,_matrix_N.nnz);
	//matbeg n_row
	fprintf(pFile,"[");
	for (i=0;i<_matrix_N.n_row;i++)
	{
		if(i==_matrix_N.n_row-1) fprintf(pFile,"%d ",_matrix_N.matbeg[i]);
		else fprintf(pFile,"%d,",_matrix_N.matbeg[i]);
	}
	fprintf(pFile,"],\n");
	//matcnt n_row
	fprintf(pFile,"[");
	for (i=0;i<_matrix_N.n_row;i++)
	{
		if(i==_matrix_N.n_row-1) fprintf(pFile,"%d ",_matrix_N.matcnt[i]);
		else fprintf(pFile,"%d,",_matrix_N.matcnt[i]);
	}
	fprintf(pFile,"],\n");
	//matind n_nz
	fprintf(pFile,"[");
	for (i=0;i<_matrix_N.nnz;i++)
	{
		if(i==_matrix_N.nnz-1) fprintf(pFile,"%d ",_matrix_N.matind[i]);
		else fprintf(pFile,"%d,",_matrix_N.matind[i]);
	}
	fprintf(pFile,"],\n");
	//matval n_nz
	fprintf(pFile,"[");
	for (i=0;i<_matrix_N.nnz;i++)
	{
		if(i==_matrix_N.nnz-1) fprintf(pFile,"%12.12f ",_matrix_N.matval[i]);
		else  fprintf(pFile,"%12.12f,",_matrix_N.matval[i]);
	}
	fprintf(pFile,"]\n");
	fprintf(pFile,"]\n");	
	//M
	fprintf(pFile,"[");
	//parameter n_row n_col n_nz
	fprintf(pFile,"[%d,%d,%d],\n",_matrix_M.n_row,_matrix_M.n_col,_matrix_M.nnz);
	//matbeg n_row
	fprintf(pFile,"[");
	for (i=0;i<_matrix_M.n_row;i++)
	{
		if(i==_matrix_M.n_row-1) fprintf(pFile,"%d ",_matrix_M.matbeg[i]);
		else fprintf(pFile,"%d,",_matrix_M.matbeg[i]);
	}
	fprintf(pFile,"],\n");
	//matcnt n_row
	fprintf(pFile,"[");
	for (i=0;i<_matrix_M.n_row;i++)
	{
		if(i==_matrix_M.n_row-1) fprintf(pFile,"%d ",_matrix_M.matcnt[i]);
		else fprintf(pFile,"%d,",_matrix_M.matcnt[i]);
	}
	fprintf(pFile,"],\n");
	//matind n_nz
	fprintf(pFile,"[");
	for (i=0;i<_matrix_M.nnz;i++)
	{
		if(i==_matrix_M.nnz-1) fprintf(pFile,"%d ",_matrix_M.matind[i]);
		else fprintf(pFile,"%d,",_matrix_M.matind[i]);
	}
	fprintf(pFile,"],\n");
	//matval n_nz
	fprintf(pFile,"[");
	for (i=0;i<_matrix_M.nnz;i++)
	{
		if (i==_matrix_M.nnz-1) fprintf(pFile,"%12.12f ",_matrix_M.matval[i]);
		else fprintf(pFile,"%12.12f,",_matrix_M.matval[i]);
	}
	fprintf(pFile,"]\n");
	fprintf(pFile,"]\n");	
	fclose(pFile);	

TERMINATE:
	free_matrix(&_matrix_A);
	free_matrix(&_matrix_B);
	free_matrix(&_matrix_N);
	free_matrix(&_matrix_M);
	return(status);
}
int  writeLPCCinstance(int random_seed,
					   int param_k,
					   int param_m,
					   int param_n,
					   int rankM,
					   double spars,
					   int mode)
{
	int status = 0;
	int printflag = 0;
	int    i;
	double *c_coef	=NULL;
	double *d_coef=NULL;
	double *b_coef=NULL;
	double *q_coef=NULL;
	double **matrix_A= NULL;
	double **matrix_B=NULL;
	double **matrix_N=NULL;
	double **matrix_M=NULL;
	status= generateLPCCinstance(param_n,param_m,param_k,rankM,spars,&c_coef,&d_coef,&b_coef,&q_coef,&matrix_A,&matrix_B,&matrix_N,&matrix_M);
	if(status) goto TERMINATE;
	switch (mode)
	{
	case 0:
		status=writeLPCCinstance_format_1(random_seed,param_n,param_m,param_k,rankM,spars,c_coef,d_coef,b_coef,q_coef,matrix_A,matrix_B,matrix_N,matrix_M);
		if (status) goto TERMINATE;
		break;
	case 1:
		status=writeLPCCinstance_format_2(random_seed,param_n,param_m,param_k,rankM,spars,c_coef,d_coef,b_coef,q_coef,matrix_A,matrix_B,matrix_N,matrix_M);
		if (status) goto TERMINATE;
		break;
	case 2:
		status=writeLPCCinstance_format_3(random_seed,param_n,param_m,param_k,rankM,spars,c_coef,d_coef,b_coef,q_coef,matrix_A,matrix_B,matrix_N,matrix_M);
		if (status) goto TERMINATE;
		break;
	default:
		status=-1;
		fprintf(stderr,"Unknown mode for writing inputfile\n");
		goto TERMINATE;
	}
TERMINATE:
	if ( matrix_A != NULL ) {
		for (i = 0; i <param_k ; ++i) {
			free_and_null ((char **) &(matrix_A[i]));
		}
	}
	if ( matrix_B != NULL ) {
		for (i = 0; i < param_k; ++i) {
			free_and_null ((char **) &(matrix_B[i]));
		}
	}
	if ( matrix_N != NULL ) {
		for (i = 0; i < param_m; ++i) {
			free_and_null ((char **) &(matrix_N[i]));
		}
	}
	if ( matrix_M != NULL ) {
		for (i = 0; i < param_m ; ++i) {
			free_and_null ((char **) &(matrix_M[i]));
		}
	}
	free_and_null ((char **) &matrix_A);
	free_and_null ((char **) &matrix_B);
	free_and_null ((char **) &matrix_N);
	free_and_null ((char **) &matrix_M);
	free_and_null ((char **) &c_coef);
	free_and_null ((char **) &d_coef);
	free_and_null ((char **) &b_coef);
	free_and_null ((char **) &q_coef);

	return (status);
}  

