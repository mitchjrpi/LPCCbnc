#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "LPCC.h"
#include "CPLEXcontrol.h"
#include "utilities.h"
int add_ub_constraint(CPXENVptr env,
					  CPXLPptr lp,
					  double ub)
{
	int status=0;
	const int num_col=CPXgetnumcols(env,lp);
	const int num_row=CPXgetnumrows(env,lp);
	const int int_zero=0;
	int i;
	int *matind=NULL;
	double *obj_coef=NULL;
	matind=(int*) malloc(num_col*sizeof(int));
	obj_coef=(double*) malloc(num_col*sizeof(double));
	if (matind==NULL || obj_coef==NULL)
	{
		status=NO_MEMORY;
		fprintf(stderr,"add_ub_constraint(): unable to allocate memory\n");
		goto TERMINATE;
	}
	status=CPXgetobj(env,lp,obj_coef,0,num_col-1);
	if (status) goto TERMINATE;
	for (i=0;i<num_col;i++)
	{
		matind[i]=i;
	}
	status=CPXaddrows(env,lp,1,1,num_col,&ub,NULL,&int_zero,matind,obj_coef,NULL,NULL);
	if(status) goto TERMINATE;
	status=CPXchgcoef(env,lp,num_row,num_col,1.0);
	if (status) goto TERMINATE;
TERMINATE:
	free_and_null((char**)&matind);
	free_and_null((char**)&obj_coef);
	return(status);
}

int add_lb_constraint(CPXENVptr env,
					  CPXLPptr lp,
					  double lb)
{
	int status=0;
	const int num_col=CPXgetnumcols(env,lp);
	const int num_row=CPXgetnumrows(env,lp);
	const int int_zero=0;
	int i;
	int *matind=NULL;
	double *obj_coef=NULL;
	matind=(int*) malloc(num_col*sizeof(int));
	obj_coef=(double*) malloc(num_col*sizeof(double));
	if (matind==NULL || obj_coef==NULL)
	{
		status=NO_MEMORY;
		fprintf(stderr,"add_lb_constraint(): unable to allocate memory\n");
		goto TERMINATE;
	}
	status=CPXgetobj(env,lp,obj_coef,0,num_col-1);
	if (status) goto TERMINATE;
	for (i=0;i<num_col;i++)
	{
		matind[i]=i;
	}
	status=CPXaddrows(env,lp,1,1,num_col,&lb,NULL,&int_zero,matind,obj_coef,NULL,NULL);
	if(status) goto TERMINATE;
	status=CPXchgcoef(env,lp,num_row,num_col,-1.0);
	if (status) goto TERMINATE;
TERMINATE:
	free_and_null((char**)&matind);
	free_and_null((char**)&obj_coef);
	return(status);

}

int lp_ub_control(CPXENVptr env,
				  CPXLPptr lp,
				  double ub,
				  int reset_flag)
{
	int status=0;
	static int flag=-1;
	static int row_index=-1;
	int num_col;
	int num_row;
	const int int_zero=0;
	int i;
	int *matind=NULL;
	double *obj_coef=NULL;
	double old_ub;
	if (reset_flag)
	{
		flag=-1;
		row_index=-1;
		goto TERMINATE;
	}
	num_col=CPXgetnumcols(env,lp);
	num_row=CPXgetnumrows(env,lp);
	if (flag==-1)
	{
		matind=(int*) malloc(num_col*sizeof(int));
		obj_coef=(double*) malloc(num_col*sizeof(double));
		if (matind==NULL || obj_coef==NULL)
		{
			status=NO_MEMORY;
			fprintf(stderr,"lp_ub_control(): unable to allocate memory\n");
			goto TERMINATE;
		}
		status=CPXgetobj(env,lp,obj_coef,0,num_col-1);
		if (status) goto TERMINATE;
		for (i=0;i<num_col;i++)
		{
			matind[i]=i;
		}
		status=CPXaddrows(env,lp,1,1,num_col,&ub,NULL,&int_zero,matind,obj_coef,NULL,NULL);
		if(status) goto TERMINATE;
		status=CPXchgcoef(env,lp,num_row,num_col,1.0);
		if (status) goto TERMINATE;
		row_index=num_row;
		flag=1;
		//fprintf(stdout,"add ub constraint at row: %d, with bound: %f\n",row_index,ub);
	}else
	{
		status=CPXgetrhs(env,lp,&old_ub,row_index,row_index);
		if(status) goto TERMINATE;
		if (old_ub<=ub)
		{
			//fprintf(stdout,"can not change ub from %f to %f at row: %d\n",old_ub,ub,row_index);
		}else
		{
			//fprintf(stdout,"change ub from %f to %f at row: %d\n",old_ub,ub,row_index);
			status=CPXchgrhs(env,lp,1,&row_index,&ub);
			if(status)goto TERMINATE;
		}
	}
TERMINATE:
	free_and_null((char**)&matind);
	free_and_null((char**)&obj_coef);
	return(status);
}

int zeroobjectivefunc(CPXENVptr env,
					  CPXLPptr lp)
{
	int status=0;
	int i;
	int *indices=NULL;
	double *values=NULL;
	indices=(int*) malloc(CPXgetnumcols(env,lp)*sizeof(int));
	values=(double*) malloc(CPXgetnumcols(env,lp)*sizeof(double));
	if (indices==NULL || values==NULL)
	{
		status=NO_MEMORY;
		fprintf(stderr,"zeroobjectivefunc(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<CPXgetnumcols(env,lp);i++)
	{
		values[i]=0;
		indices[i]=i;
	}
	status=CPXchgobj(env,lp,CPXgetnumcols(env,lp),indices,values);
	if (status) goto TERMINATE;
TERMINATE:
	free_and_null((char**)&indices);
	free_and_null((char**)&values);
	return(status);
}

int setup_initial_CGLP(CPXENVptr env, 
					   CPXLPptr lp,
					   const int param_k,
					   const int param_m,
					   const int param_n,
					   double *b_coef,  
					   double *q_coef, 
					   MATRIX matrix_A, 
					   MATRIX matrix_B, 
					   MATRIX matrix_N, 
					   MATRIX matrix_M)
{
	//matrix need to be grouped by row
	int status=0;
	int i;
	int tmp_int;
	int index;
	int cur_lp_ncol;
	int cur_lp_nrow;
	char lu='L';
	const double negtiveINF=-CPX_INFBOUND; 
	const double one=1;
	//double tmp_double;
	MATRIX matrix_At;
	MATRIX matrix_Bt;
	MATRIX matrix_Mt;
	MATRIX matrix_Nt;
	MATRIX matrix_lp;
	MATRIX matrix_tmp;
	MATRIX matrix_zero;
	MATRIX matrix_identity;
	init_matrix(&matrix_At);
	init_matrix(&matrix_Bt);
	init_matrix(&matrix_Mt);
	init_matrix(&matrix_Nt);
	init_matrix(&matrix_lp);
	init_matrix(&matrix_tmp);
	init_matrix(&matrix_zero);
	init_matrix(&matrix_identity);
	if (matrix_A.param!=1 || matrix_B.param!=1 || matrix_M.param!=1 || matrix_N.param!=1)
	{
		status=-1;
		fprintf(stderr," setup_initial_CGLP():input matrix is not grouped by row\n");
		goto TERMINATE;

	}
	status=copy_matrix(&matrix_A,&matrix_At);
	if(status) goto TERMINATE;
	status=copy_matrix(&matrix_B,&matrix_Bt);
	if(status) goto TERMINATE;
	status=copy_matrix(&matrix_M,&matrix_Mt);
	if(status) goto TERMINATE;
	status=copy_matrix(&matrix_N,&matrix_Nt);
	if(status) goto TERMINATE;
	//transpose matrix
	matrix_At.param=0;
	tmp_int=matrix_At.n_col;
	matrix_At.n_col=matrix_At.n_row;
	matrix_At.n_row=tmp_int;

	matrix_Bt.param=0;
	tmp_int=matrix_Bt.n_col;
	matrix_Bt.n_col=matrix_Bt.n_row;
	matrix_Bt.n_row=tmp_int;

	matrix_Mt.param=0;
	tmp_int=matrix_Mt.n_col;
	matrix_Mt.n_col=matrix_Mt.n_row;
	matrix_Mt.n_row=tmp_int;

	matrix_Nt.param=0;
	tmp_int=matrix_Nt.n_col;
	matrix_Nt.n_col=matrix_Nt.n_row;
	matrix_Nt.n_row=tmp_int;
	//construct initial CGLP
	//order of var:beta, alpha_x, alpha_y,u_11,u_21,u_31,u_12,u_22,u_32
	status = CPXnewrows(env,lp,(2*param_n+2*param_m),NULL,NULL,NULL,NULL);
	if (status) goto TERMINATE;
	//add beta column
	status=CPXnewcols(env,lp,1,NULL,&negtiveINF,NULL,NULL,NULL);
	if(status) goto TERMINATE;
	//alpha_x
	
	status=identity_matrix(0,param_n,-1,&matrix_identity);
	if(status) goto TERMINATE;
	status=zero_matrix(0,param_m,param_n,&matrix_zero);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_identity,matrix_zero,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_tmp,matrix_tmp,&matrix_lp,1);
	if(status) goto TERMINATE;
	//alpha_y
	status=identity_matrix(0,param_m,-1,&matrix_identity);
	if(status) goto TERMINATE;
	status=zero_matrix(0,param_n,param_m,&matrix_zero);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_zero,matrix_identity,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_tmp,matrix_tmp,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_lp,matrix_tmp,&matrix_lp,0);
	if(status) goto TERMINATE;
		
	//u_11
	status=combine_matrix(matrix_At,matrix_Bt,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=zero_matrix(0,(param_n+param_m),param_k,&matrix_zero);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_tmp,matrix_zero,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_lp,matrix_tmp,&matrix_lp,0);
	if(status) goto TERMINATE;
	//u_21
	status=zero_matrix(0,param_n,param_m,&matrix_zero);
	if(status) goto TERMINATE;
	status=identity_matrix(0,param_m,1,&matrix_identity);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_zero,matrix_identity,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=zero_matrix(0,(param_n+param_m),param_m,&matrix_zero);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_tmp,matrix_zero,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_lp,matrix_tmp,&matrix_lp,0);
	if(status) goto TERMINATE;
	//u_31
	status=combine_matrix(matrix_Nt,matrix_Mt,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=zero_matrix(0,(param_n+param_m),param_m,&matrix_zero);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_tmp,matrix_zero,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_lp,matrix_tmp,&matrix_lp,0);
	if(status) goto TERMINATE;
	//u_12
	status=combine_matrix(matrix_At,matrix_Bt,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=zero_matrix(0,(param_n+param_m),param_k,&matrix_zero);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_zero,matrix_tmp,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_lp,matrix_tmp,&matrix_lp,0);
	if(status) goto TERMINATE;
	//u_22	
	status=zero_matrix(0,param_n,param_m,&matrix_zero);
	if(status) goto TERMINATE;
	status=identity_matrix(0,param_m,1,&matrix_identity);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_zero,matrix_identity,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=zero_matrix(0,(param_n+param_m),param_m,&matrix_zero);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_zero,matrix_tmp,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_lp,matrix_tmp,&matrix_lp,0);
	if(status) goto TERMINATE;
	//u_32
	status=combine_matrix(matrix_Nt,matrix_Mt,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=zero_matrix(0,(param_n+param_m),param_m,&matrix_zero);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_zero,matrix_tmp,&matrix_tmp,1);
	if(status) goto TERMINATE;
	status=combine_matrix(matrix_lp,matrix_tmp,&matrix_lp,0);
	if(status) goto TERMINATE;
	//add constraints into lp
	status=CPXaddcols(env,lp,matrix_lp.n_col,matrix_lp.nnz,NULL,matrix_lp.matbeg,matrix_lp.matind,matrix_lp.matval,NULL,NULL,NULL);
	if(status) goto TERMINATE;
	// setting the lb for alpha_x alpha_y as -cpx_infbound
	for (i=0;i<param_n;i++)
	{
		index=i+1;
		lu='L';
		status=CPXchgbds(env,lp,1,&index,&lu,&negtiveINF);
		if(status) goto TERMINATE;
	}
	for (i=0;i<param_m;i++)
	{
		index=param_n+i+1;
		lu='L';
		status=CPXchgbds(env,lp,1,&index,&lu,&negtiveINF);
		if(status) goto TERMINATE;
	}
	//add constraint for beta
	cur_lp_nrow=CPXgetnumrows(env,lp);
	cur_lp_ncol=CPXgetnumcols(env,lp);
	status=CPXnewrows(env,lp,2,NULL,NULL,NULL,NULL);
	if(status) goto TERMINATE;
	status=CPXnewcols(env,lp,2,NULL,NULL,NULL,NULL,NULL);
	if(status) goto TERMINATE;
	status=CPXchgcoef(env,lp,cur_lp_nrow,0,-1);
	if(status) goto TERMINATE;
	status=CPXchgcoef(env,lp,cur_lp_nrow+1,0,-1);
	if(status) goto TERMINATE;
	status=CPXchgcoef(env,lp,cur_lp_nrow,cur_lp_ncol,-1);
	if(status) goto TERMINATE;
	status=CPXchgcoef(env,lp,cur_lp_nrow+1,cur_lp_ncol+1,-1);
	for (i=0;i<param_k;i++)
	{
		index=1+param_n+param_m+i;
		status=CPXchgcoef(env,lp,cur_lp_nrow,index,b_coef[i]);
		if(status) goto TERMINATE;
		index=1+param_n+param_m+param_k+param_m+param_m+i;
		status=CPXchgcoef(env,lp,cur_lp_nrow+1,index,b_coef[i]);
		if (status) goto TERMINATE;
	}
	for (i=0;i<param_m;i++)
	{
		index=1+param_n+param_m+param_k+param_m+i;
		status=CPXchgcoef(env,lp,cur_lp_nrow,index,-q_coef[i]);
		if(status) goto TERMINATE;
		index=1+param_n+param_m+param_k+param_m+param_m+param_k+param_m+i;
		status=CPXchgcoef(env,lp,cur_lp_nrow+1,index,-q_coef[i]);
		if (status) goto TERMINATE;
	}
	//add the last normalize constraint
	cur_lp_nrow=CPXgetnumrows(env,lp);
	cur_lp_ncol=CPXgetnumcols(env,lp);
	//lu='L';
	status=CPXnewrows(env,lp,1,NULL,NULL,NULL,NULL);
	if(status) goto TERMINATE;
	status=CPXchgrhs(env,lp,1,&cur_lp_nrow,&one);
	if(status) goto TERMINATE;
	for (i=(1+param_n+param_m);i<(1+param_n+param_m+2*param_k+2*param_m+2*param_m);i++)
	{
		status=CPXchgcoef(env,lp,cur_lp_nrow,i,1);
		if(status) goto TERMINATE;
	}
TERMINATE:
	free_matrix(&matrix_At);
	free_matrix(&matrix_Bt);
	free_matrix(&matrix_Mt);
	free_matrix(&matrix_Nt);
	free_matrix(&matrix_lp);
	free_matrix(&matrix_tmp);
	free_matrix(&matrix_zero);
	free_matrix(&matrix_identity);
	return(status);
}


int setupGeneralLPCC_lp(CPXENVptr env, 
						CPXLPptr lp,
						double *c_coef, 
						double *d_coef, 
						double *b_coef,  
						double *q_coef, 
						MATRIX  matrix_A, 
						MATRIX matrix_B, 
						MATRIX matrix_N, 
						MATRIX matrix_M,
						int format_mode)
{
	//matrix need to be grouped by row
	//min c*x+d*y
	//st. A*x+B*y		-s=b
	//     N*x+M*y-w      =-q
	//	   y>=0;w>=0;s>=0
	//	   x>=0;
	// and empty bound cuts
	// var: x,y,w,s, slack var for bound cuts

	int status=0;
	int i;
	int x_num;
	int y_num;
	int s_num;
	int w_num;
	double *neg_q_coef=NULL;
	double *x_lb=NULL;
	double *x_ub=NULL;
	MATRIX matrix_AB;
	MATRIX matrix_NM;
	init_matrix(&matrix_AB);
	init_matrix(&matrix_NM);
	//check the matrix dimension and initial the variable dimension
	if (matrix_A.n_col==matrix_N.n_col)
	{
		x_num=matrix_A.n_col;
	}else
	{
		status=-1;
		fprintf(stderr,"setupGeneralLPCC_lp(): A_col!= N_col\n");
		goto TERMINATE;
	}
	if (matrix_B.n_col==matrix_M.n_col)
	{
		y_num=matrix_B.n_col;
	}else
	{
		status=-1;
		fprintf(stderr,"setupGeneralLPCC_lp(): B_col!= M_col\n");
		goto TERMINATE;
	}
	if (matrix_A.n_row==matrix_B.n_row)
	{
		s_num=matrix_A.n_row;
	}else
	{
		status=-1;
		fprintf(stderr,"setupGeneralLPCC_lp(): A_row!= B_row\n");
		goto TERMINATE;
	}
	if (matrix_N.n_row==matrix_M.n_row)
	{
		w_num=matrix_M.n_row;
	}else
	{
		status=-1;
		fprintf(stderr,"setupGeneralLPCC_lp(): N_row!= M_row\n");
		goto TERMINATE;
	}
	//check whether matrix is grouped by row
	if (matrix_A.param!=1 || matrix_B.param!=1 || matrix_M.param!=1 || matrix_N.param!=1)
	{
		status=-1;
		fprintf(stderr," setupGeneralLPCC_lp(): matrix is not grouped by row/n");
		goto TERMINATE;
	}
	neg_q_coef=(double*) malloc(y_num*sizeof(double));
	if (neg_q_coef==NULL)
	{
		status=NO_MEMORY;
		fprintf(stderr," setupGeneralLPCC_lp(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<y_num;i++)
	{
		neg_q_coef[i]=-q_coef[i];
	}
	x_lb=(double*) malloc(x_num*sizeof(double));
	x_ub=(double*) malloc(x_num*sizeof(double));
	if ( x_lb==NULL || x_ub==NULL )	{
		status=NO_MEMORY;
		fprintf(stderr," setupGeneralLPCC_lp(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<x_num;i++)
	{
		if(format_mode==4)x_lb[i]=-CPX_INFBOUND;
		else x_lb[i]=0.0;
		x_ub[i]=CPX_INFBOUND;
	}
	fprintf(stdout,"Setup LPCC problem...\n");
	//print_full_matrix(&matrix_A);
	//print_full_matrix(&matrix_B);
	//var x
	status=CPXnewcols(env,lp,x_num,c_coef,x_lb,x_ub,NULL,NULL);
	if (status) goto TERMINATE; 
	//fprintf(stdout,"Setup LPCC problem...1\n");
	//var y
	status=CPXnewcols(env,lp,y_num,d_coef,NULL,NULL,NULL,NULL);
	if (status) goto TERMINATE;
	//fprintf(stdout,"Setup LPCC problem...2\n");
	// add constraint
	status=combine_matrix(matrix_A,matrix_B,&matrix_AB,0);
	//print_full_matrix(&matrix_AB);
	if (status) goto TERMINATE;
	//fprintf(stdout,"Setup LPCC problem...3\n");
	status=CPXaddrows(env,lp,0,matrix_AB.n_row,matrix_AB.nnz,b_coef,NULL,matrix_AB.matbeg,matrix_AB.matind,matrix_AB.matval,NULL,NULL);
	if (status) goto TERMINATE;
	//fprintf(stdout,"Setup LPCC problem...4\n");
	status=combine_matrix(matrix_N,matrix_M,&matrix_NM,0);
	if (status) goto TERMINATE;
	//fprintf(stdout,"Setup LPCC problem...5\n");
	//print_full_matrix(&matrix_N);
	//print_full_matrix(&matrix_M);
	//print_full_matrix(&matrix_NM);
	status=CPXaddrows(env,lp,0,matrix_NM.n_row,matrix_NM.nnz,neg_q_coef,NULL,matrix_NM.matbeg,matrix_NM.matind,matrix_NM.matval,NULL,NULL);
	if (status) goto TERMINATE;
	//add var w and slack var s
	//var w
	//fprintf(stdout,"Setup LPCC problem...6\n");
	status=CPXnewcols(env,lp,w_num,NULL,NULL,NULL,NULL,NULL);
	if (status) goto TERMINATE;
	//modify the coefficient of w
	//fprintf(stdout,"Setup LPCC problem...7\n");
	for (i=0; i<w_num;i++)
	{
		status = CPXchgcoef (env, lp, s_num+i, x_num+y_num+i, -1.0);
		if (status) goto TERMINATE;
	}
	//slack var
	//fprintf(stdout,"Setup LPCC problem...8\n");
	status=CPXnewcols(env,lp,s_num,NULL,NULL,NULL,NULL,NULL);
	if (status) goto TERMINATE;
	//modify the coefficient of s
	for (i=0; i<s_num;i++)
	{
		status = CPXchgcoef (env, lp, i, x_num+y_num+w_num+i, -1.0);
		if (status) goto TERMINATE;
	}
	//empty bound cuts
	//fprintf(stdout,"Setup LPCC problem...9\n");
	status=CPXnewcols(env,lp,y_num,NULL,NULL,NULL,NULL,NULL);
	if (status) goto TERMINATE;
	//fprintf(stdout,"Setup LPCC problem...10\n");
	status=CPXnewrows(env,lp,y_num,NULL,NULL,NULL,NULL);
	if (status) goto TERMINATE;
	//fprintf(stdout,"Setup LPCC problem...11\n");
	status=lp_ub_control(env,lp,CPX_INFBOUND,1);
	if (status) goto TERMINATE;
	//fprintf(stdout,"Setup LPCC problem...12\n");
TERMINATE:
	if (status) fprintf(stdout,"crashed in Setup LPCC problem...\n");
	free_and_null((char**)&neg_q_coef);
	free_and_null((char**) &x_lb);
	free_and_null((char**) &x_ub);
	free_matrix(&matrix_AB);
	free_matrix(&matrix_NM);
	return(status);
}


int setupLP_GE (CPXENVptr env,
				CPXLPptr lp, 
				double *c_coef, 
				MATRIX matrix_A_c, 
				double *b_coef)
{
	//matrix need to be grouped by column
	//min c*x
	//st. A*x>=b
	int		status		= 0;
	int		i;
	int		n_col_A		= matrix_A_c.n_col;
	int		n_row_A		= matrix_A_c.n_row;
	int		*ind = NULL;
	double	*lb			= NULL; 
	double  *ub			= NULL;
	char	*sense		= NULL;
	if (matrix_A_c.param!=0)
	{
		status=-1;
		fprintf(stderr," setupLP_GE(): matrix is not grouped by column/n");
		goto TERMINATE;
	}
	ind=(int*) malloc(n_col_A*sizeof(int));
	lb = (double*) malloc(n_col_A * sizeof(double));
	ub = (double*) malloc(n_col_A * sizeof(double));
	sense = (char*) malloc(n_row_A * sizeof(char));
	if ( ind == NULL || lb == NULL || sense == NULL || ub == NULL)
	{
		status = NO_MEMORY;
		fprintf (stderr, "setupLP_GE(): could not allocate memory.\n");
		goto TERMINATE;
	}
	for (i = 0; i < n_col_A; i++) 
	{     
		lb[i] = -CPX_INFBOUND;
		ub[i] = CPX_INFBOUND;
		ind[i] =i;
	}
	for (i = 0; i < n_row_A; i++) 
	{
		sense[i] = 'E';
	}
	status = CPXcopylp(env,
					   lp,
					   n_col_A,
					   n_row_A,
					   CPX_MIN,
					   c_coef, 
					   b_coef,
		               sense, 
					   matrix_A_c.matbeg, 
					   matrix_A_c.matcnt,
					   matrix_A_c.matind, 
					   matrix_A_c.matval,
					   lb, 
					   ub,
					   NULL);
	if(status) goto TERMINATE;
	//add slack variables
	status = CPXnewcols(env, lp, n_row_A, NULL, NULL, NULL, NULL, NULL);
	if (status) goto TERMINATE;
	//modify the coefficient of slack variables
	for (i=0; i<n_row_A;i++)
	{
		status = CPXchgcoef (env, lp, i, n_col_A+i, -1.0);
		if (status) goto TERMINATE;
	}
TERMINATE:
	free_and_null((char **) &ind);
	free_and_null ((char **)&sense);
	free_and_null ((char **)&lb);
	free_and_null ((char **)&ub);
	return (status);
}
int setupLP_E (CPXENVptr env,
			   CPXLPptr lp, 
			   double *c_coef, 
			   MATRIX matrix_A_c, 
			   double *b_coef)
{
	//matrix need to be grouped by column
	//min c*x
	//st. A*x=b
	int		status		= 0;
	int		i;
	int		n_col_A		= matrix_A_c.n_col;
	int		n_row_A		= matrix_A_c.n_row;
	int *ind = NULL;
	double	*lb			= NULL; 
	double  *ub			= NULL;
	char	*sense		= NULL;
	if (matrix_A_c.param!=0)
	{
		status=-1;
		fprintf(stderr," setupLP_E(): matrix is not grouped by column/n");
		goto TERMINATE;
	}
	ind=(int*) malloc(n_col_A*sizeof(int));
	lb = (double*) malloc(n_col_A * sizeof(double));
	ub = (double*) malloc(n_col_A * sizeof(double));
	sense = (char*) malloc(n_row_A * sizeof(char));
	if ( ind == NULL || lb == NULL || sense == NULL || ub == NULL)
	{
		status = NO_MEMORY;
		fprintf (stderr, "setupLP_E(): could not allocate memory.\n");
		goto TERMINATE;
	}
	for (i = 0; i < n_col_A; i++) 
	{     
		lb[i] = -CPX_INFBOUND;
		ub[i] = CPX_INFBOUND;
		ind[i] =i;
	}
	for (i = 0; i < n_row_A; i++) 
	{
		sense[i] = 'E';
	}
	status = CPXcopylp(env,
					   lp,
					   n_col_A,
					   n_row_A,
					   CPX_MIN,
					   c_coef,
					   b_coef,
					   sense,
					   matrix_A_c.matbeg,
					   matrix_A_c.matcnt,
					   matrix_A_c.matind,
					   matrix_A_c.matval,
					   lb, 
					   ub,
					   NULL);
	if(status) goto TERMINATE;
TERMINATE:
	free_and_null((char **) &ind);
	free_and_null ((char **)&sense);
	free_and_null ((char **)&lb);
	free_and_null ((char **)&ub);
	return (status);
}

// this routine comibines several cplex routine togethe
// if the cplex solution status shows the problem is infeasible or unbounded, then this routine will turn off the primal reduce and dual reduce preprocess ( is equivalent to turn off presolve ) and resolve the problem.
// this routine will check the cplex solution status, and only consider "infeasible"," unbounded", " optimal" as valid solution status.
int cplexSolveLP(CPXENVptr env, 
				 CPXLPptr lp, 
				 STARTINFO start, 
				 const int num_var,
				 double *x_soln_p, 
				 double *objval_p, 
				 int *solnstat_p, 
				 STARTINFO *start_p,
				 int param_1,
				 int param_2,
				 int pre_rx_lp_condition)
{
	int status = 0;
	int solnstat;
	int solnmethod;
	int solntype;
	int lp_col = CPXgetnumcols (env, lp);
    int lp_row = CPXgetnumrows (env, lp);
	int *cstat = NULL;
	int *rstat = NULL;
	int *cindices = NULL;
	int *rindices = NULL;
	int i; 
	int itcnt =0;
	double *x_soln= NULL;	
	x_soln=(double*) malloc(lp_col*sizeof(double));
	if (x_soln==NULL)
	{
		status = NO_MEMORY;
		fprintf (stderr, "cplexSolveLP(): could not allocate memory.\n");
		goto TERMINATE;
	}
	//cindices = (int *) malloc (lp_col*sizeof(int));
	//rindices = (int *) malloc (lp_row*sizeof(int));
	cindices = (int *) malloc (start.ccnt*sizeof(int));
	rindices = (int *) malloc (start.rcnt*sizeof(int));
			cstat = (int *) malloc (lp_col*sizeof(int));
			rstat = (int *) malloc (lp_row*sizeof(int));
	if (cindices==NULL || rindices==NULL)
	{
		status = NO_MEMORY;
		fprintf (stderr, "cplexSolveLP(): could not allocate memory.\n");
		goto TERMINATE;
	}
	for ( i = 0; i < start.ccnt; i++)
	{
        cstat[i]=i;  //cindices[i]=i;
	}

	for ( i = start.ccnt; i < lp_col; i++)
	{
		cstat[i]=0;
	}

	for ( i = 0; i < start.rcnt; i++)
	{
        rstat[i]=i;   //rindices[i]= i;
	}

	for ( i = start.rcnt; i < lp_row; i++)
	{
		rstat[i]= 1;
	}

	if (start.cstat!=NULL && start.rstat!=NULL && lp_col==start.ccnt)
	{
		/*
            fprintf (stderr, "cplexSolveLP(): start.ccnt:%d\n",start.ccnt);
			fprintf (stderr, "cplexSolveLP(): start.rcnt:%d\n",start.rcnt);
                        fprintf (stderr, "cplexSolveLP(): lp_col:%d\n",lp_col);
                        fprintf (stderr, "cplexSolveLP(): lp_row:%d\n",lp_row);
        fprintf (stderr, "cplexSolveLP(): start.cstat:%d\n",* start.cstat);
        fprintf (stderr, "cplexSolveLP(): start.rstat:%d\n",* start.rstat);
         */

		// status = CPXcopypartialbase(env, lp, start.ccnt, cindices, start.cstat, start.rcnt, rindices, start.rstat);
        // fprintf (stderr, "start.ccnt: %d, start.rcnt: %d, lp_col: %d, lp_row: %d\n",start.ccnt, start.rcnt, lp_col, lp_row);
        //  fprintf (stderr, "lp_col - start.ccnt: %d\n",lp_col - start.ccnt);
		status = CPXcopybase(env, lp, start.cstat, start.rstat);
        // status = CPXcopybase(env, lp, lp_col, lp_row);
		if ( status )
		{
			fprintf (stderr, "cplexSolveLP(): failed to copy basis into LP:%d\n",status);
                        char  errmsg[CPXMESSAGEBUFSIZE];
                        CPXgeterrorstring (env, status, errmsg);
                        fprintf (stderr, "%s", errmsg);
			// goto TERMINATE;
		}
        /*  else
        {
            fprintf (stderr, "cplexSolveLP(): successful copy basis into LP:%d\n",status);
        }   */
	}
    /*  else
    {
        fprintf (stderr, "one of start.cstat or start.rstat is null\n");
    }   */
    status = CPXsetintparam (env, CPX_PARAM_LPMETHOD, param_1);
	if ( status ) 
	{
		fprintf (stderr, "cplexSolveLP(): failure to setup CPX_PARAM_LPMETHOD\n");
		goto TERMINATE;
	}
	status = CPXsetintparam (env, CPX_PARAM_ADVIND, param_2);  
	if ( status ) 
	{
		fprintf (stderr, "cplexSolveLP(): failure to setup CPX_PARAM_ADVIND\n");
		goto TERMINATE;
	}
	status = CPXlpopt (env, lp);
	if ( status ) 
	{
		fprintf (stderr, "cplexSolveLP(): failed to optimize LP.\n");
        goto TERMINATE;
	}
	itcnt = CPXgetitcnt (env, lp);
	solnstat = CPXgetstat (env, lp);
	switch(pre_rx_lp_condition)
	{
	case LP_PRECONDITION_BOUNDED:
		status = CPXsolninfo (env, lp, &solnmethod, &solntype, NULL, NULL);
		if ( status ) 
		{
			fprintf (stderr, "cplexSolveLP(): failed to obtain solution info.\n");
			goto TERMINATE;
		}
		if ( solnstat == CPX_STAT_INFEASIBLE || solnstat==CPX_STAT_INForUNBD)
		{
			*solnstat_p = CPX_STAT_INFEASIBLE;
		}else if ( solnstat == CPX_STAT_OPTIMAL || solnstat == CPX_STAT_OPTIMAL_INFEAS)
		{
			*solnstat_p = CPX_STAT_OPTIMAL;
			//*solnstat_p = solnstat; 
			if ( solntype == CPX_NO_SOLN )
			{
				fprintf (stderr, "cplexSolveLP(): solution not available.\n");
				status = -1;
				goto TERMINATE;
			}
			status = CPXgetobjval (env, lp, objval_p);
			if ( status ) 
			{
				fprintf (stderr, "cplexSolveLP(): failed to obtain objective value.\n");
				goto TERMINATE;
			}
			status = CPXgetx (env, lp, x_soln, 0, lp_col-1);
			if ( status ) 
			{
				fprintf (stderr, "cplexSolveLP(): failed to obtain primal solution.\n");
				goto TERMINATE;
			}
		}else 
		{
			status = -1;
			fprintf (stderr,"cplexSolveLP(): unknow problem happened during CPLEX solving, solnstat == %d\n",solnstat);
			goto TERMINATE;
		}
		if ( solntype == CPX_BASIC_SOLN ) 
		{
			cstat = (int *) malloc (lp_col*sizeof(int));
			rstat = (int *) malloc (lp_row*sizeof(int));
			if ( cstat == NULL || rstat == NULL) 
			{
				status = NO_MEMORY;
				fprintf (stderr, "cplexSolveLP(): could not allocate memory.\n");
				goto TERMINATE; 
			} 
			status = CPXgetbase (env, lp, cstat, rstat);   
			if ( status ) 
			{      
				fprintf (stderr, "cplexSolveLP(): failed to get basis; error %d.\n", status);       
				goto TERMINATE;    
			}
			free_startinfo(start_p);
			start_p->ccnt = lp_col;
			start_p->rcnt = lp_row;
			start_p->cstat = cstat;
			start_p->rstat = rstat;
		}else 
		{ 
			status=copy_startinfo(&start,start_p);
			if(status) goto TERMINATE;
		}
		for (i=0;i< num_var;i++)
		{
			x_soln_p[i]=x_soln[i];
		}
		break;
	case LP_PRECONDITION_UNBOUNDED:
	case LP_PRECONDITION_UNKNOWN:
		if ( solnstat == CPX_STAT_INForUNBD )
		{
			//turn off primal and dual reduction
			status = CPXsetintparam (env, CPX_PARAM_REDUCE, CPX_PREREDUCE_NOPRIMALORDUAL);
			if ( status ) 
			{
				fprintf (stderr,     
					"cplexSolveLP(): failure to turn off primal and dual reduction, error %d.\n", status);
				goto TERMINATE;
			}
			status = CPXsetintparam (env, CPX_PARAM_LPMETHOD, CPX_ALG_PRIMAL);  
			if ( status )
			{
				fprintf (stderr, "cplexSolveLP(): failure to setup CPX_PARAM_LPMETHOD\n");
				goto TERMINATE;
			}
			status = CPXlpopt (env, lp);
			if ( status ) 
			{
				fprintf (stderr, "cplexSolveLP(): failed to optimize LP.\n");
				goto TERMINATE;
			}
			itcnt = CPXgetitcnt (env, lp)+itcnt;
			solnstat = CPXgetstat (env, lp);
			status = CPXsetintparam (env, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALANDDUAL);	
			if ( status ) 
			{
				fprintf (stderr,    
					"cplexSolveLP(): failure to turn on primal and dual reduction, error %d.\n", status);
				goto TERMINATE;
			}  
		}else if (solnstat == CPX_STAT_UNBOUNDED)
		{
			status = CPXsetintparam (env, CPX_PARAM_LPMETHOD, CPX_ALG_PRIMAL);  
			if ( status ) 
			{
				fprintf (stderr, "cplexSolveLP(): failure to setup CPX_PARAM_LPMETHOD\n");
				goto TERMINATE;
			}
			status = CPXlpopt (env, lp);
			if ( status ) 
			{
				fprintf (stderr, "cplexSolveLP(): failed to optimize LP.\n");
				goto TERMINATE;
			}
			solnstat = CPXgetstat (env, lp);
		}
		status = CPXsolninfo (env, lp, &solnmethod, &solntype, NULL, NULL);
		if ( status ) 
		{
			fprintf (stderr, "cplexSolveLP(): failed to obtain solution info.\n");
			goto TERMINATE;
		}
		if ( solnstat == CPX_STAT_UNBOUNDED ) 
		{
			*solnstat_p = CPX_STAT_UNBOUNDED;
			status = CPXgetray(env,lp,x_soln);
			if ( status ) 
			{
				fprintf (stderr, "cplexSolveLP(): failure to get unbounded direction\n");
			}
		}else if ( solnstat == CPX_STAT_INFEASIBLE )
		{
			*solnstat_p = CPX_STAT_INFEASIBLE;
		}else if ( solnstat == CPX_STAT_OPTIMAL || solnstat == CPX_STAT_OPTIMAL_INFEAS)
		{
			*solnstat_p = CPX_STAT_OPTIMAL;
			//*solnstat_p = solnstat; 
			if ( solntype == CPX_NO_SOLN )
			{
				fprintf (stderr, "cplexSolveLP(): solution not available.\n");
				status = -1;
				goto TERMINATE;
			}
			status = CPXgetobjval (env, lp, objval_p);
			if ( status ) 
			{
				fprintf (stderr, "cplexSolveLP(): failed to obtain objective value.\n");
				goto TERMINATE;
			}
			status = CPXgetx (env, lp, x_soln, 0, lp_col-1);
			if ( status ) 
			{
				fprintf (stderr, "cplexSolveLP(): failed to obtain primal solution.\n");
				goto TERMINATE;
			}
		}else 
		{
			status = -1;
			fprintf (stderr,"cplexSolveLP(): unknow problem happened during CPLEX solving, solnstat == %d\n",solnstat);
			goto TERMINATE;
		}
		if ( solntype == CPX_BASIC_SOLN ) 
		{
			cstat = (int *) malloc (lp_col*sizeof(int));
			rstat = (int *) malloc (lp_row*sizeof(int));
			if ( cstat == NULL || rstat == NULL) 
			{
				status = NO_MEMORY;
				fprintf (stderr, "cplexSolveLP(): could not allocate memory.\n");
				goto TERMINATE; 
			} 
			status = CPXgetbase (env, lp, cstat, rstat);   
			if ( status ) 
			{      
				fprintf (stderr, "cplexSolveLP(): failed to get basis; error %d.\n", status);       
				goto TERMINATE;    
			}
			free_startinfo(start_p);
			start_p->ccnt = lp_col;
			start_p->rcnt = lp_row;
			start_p->cstat = cstat;
			start_p->rstat = rstat;
		}else 
		{ 
			status=copy_startinfo(&start,start_p);
			if(status) goto TERMINATE;
		}
		for (i=0;i< num_var;i++)
		{
			x_soln_p[i]=x_soln[i];
		}
		break;
	default:
		status=-1;
		fprintf(stderr,"cplexSolveLP(): unexpected pre_rx_lp_condition\n");
		goto TERMINATE;
	}

TERMINATE:
	free_and_null((char**)&cindices);
	free_and_null((char**)&rindices);
	free_and_null((char**)&x_soln);
	return (status);
}
int relaxLPSolver(CPXENVptr env, 
				  CPXLPptr rx_lp, 
				  CONSTRAINT_SET *cuts_p,
				  const int num_var, 
				  int fix_count, 
				  int* fixed_index_p,
				  STARTINFO *start, 
				  double *x_soln_p, 
				  double *objval_p, 
				  int *solnstat_p, 
				  STARTINFO *start_p,
				  int pre_rx_lp_condition)
{
	int status =0;
	int i;
	const int init_rx_lp_numrows = CPXgetnumrows (env, rx_lp);
	const int init_rx_lp_numcols = CPXgetnumcols (env, rx_lp);
	int cur_rx_lp_numrows;
	int cur_rx_lp_numcols;
	int addcutsflag = 0;
	int new_col;
	char lu;
	double bd;
	double *original_bd=NULL;
	original_bd=(double*)malloc(fix_count*sizeof(double));
	if (original_bd==NULL)
	{
		status=-1;
		fprintf(stderr,"unable to allocate memory\n");
		goto TERMINATE;
	}
	//add cuts
	if (cuts_p!= NULL&& cuts_p->cons_matrix.n_row>0)
	{
		addcutsflag = 1;
		new_col=cuts_p->cons_matrix.n_col-init_rx_lp_numcols;
		status = CPXaddrows(env,rx_lp,new_col,cuts_p->cons_matrix.n_row,cuts_p->cons_matrix.nnz,cuts_p->rhs,NULL,cuts_p->cons_matrix.matbeg,cuts_p->cons_matrix.matind,cuts_p->cons_matrix.matval,NULL,NULL);
		if (status) 
		{
			fprintf(stderr," relaxLPSolver(): add row failed\n");
			goto TERMINATE;
		}
		status = CPXnewcols(env, rx_lp, cuts_p->cons_matrix.n_row, NULL, NULL, NULL, NULL, NULL);
		if (status) goto TERMINATE;
		for (i=0; i<cuts_p->cons_matrix.n_row;i++)
		{
			status = CPXchgcoef (env, rx_lp, init_rx_lp_numrows+i, cuts_p->cons_matrix.n_col+i, -1.0);
			if (status) goto TERMINATE;
		}
	}
	//fix complementary
	for (i=0;i<fix_count;i++)
	{
		lu='U';
		status=CPXgetub(env,rx_lp,&(bd),fixed_index_p[i],fixed_index_p[i]);
		if (status) goto TERMINATE;
		original_bd[i]=bd;
		bd=0.0;
		status = CPXtightenbds(env, rx_lp, 1, &fixed_index_p[i], &lu, &bd);
		if (status) 
		{
			fprintf(stderr," relaxLPSolver(): unable to fix complementary by tightening the upper bound  \n");
			goto TERMINATE;
		}
	}
	cur_rx_lp_numcols = CPXgetnumcols (env, rx_lp);
	status = cplexSolveLP(env,rx_lp,*start,num_var,x_soln_p,objval_p,solnstat_p,start_p,CPX_ALG_DUAL,1,pre_rx_lp_condition);//m
	if (status) goto TERMINATE;
	//recover startinfo
	start_p->rcnt = init_rx_lp_numrows;
	// recover rx_lp
	for (i=0;i<fix_count;i++)
	{
		lu='U';
		bd=original_bd[i];
		status = CPXtightenbds(env, rx_lp, 1, &fixed_index_p[i], &lu, &bd);
		if (status) 
		{
			fprintf(stderr," relaxLPSolver(): unable to recover fixed complementary \n");
			goto TERMINATE;
		}
	}
	if (addcutsflag==1)
	{
		cur_rx_lp_numrows = CPXgetnumrows (env, rx_lp);
		cur_rx_lp_numcols = CPXgetnumcols (env, rx_lp);
		status = CPXdelrows(env, rx_lp,init_rx_lp_numrows,cur_rx_lp_numrows-1);
		if (status)
		{
			fprintf(stderr," relaxLPSolver(): delete row failed\n");
			goto TERMINATE;
		}
		status = CPXdelcols(env, rx_lp,init_rx_lp_numcols,cur_rx_lp_numcols-1);
		if (status)
		{
			fprintf(stderr," relaxLPSolver(): delete col failed\n");
			goto TERMINATE;
		}
	}
TERMINATE:
	free_and_null((char**)&original_bd);
	return (status);
}

int setupnode_rx_lp(CPXENVptr env, 
					CPXLPptr rx_lp, 
					const int param_n,
					const int param_m,
					const int param_k,				
					NODE activenode)
{
	int status =0;
	int i;
	int branch_index_1;
	int branch_index_2;
	int row_indice[4];
	int col_indice[4];
	double rhs[4];
	char lu;
	double bd;
	int *nodebranchinfo=activenode.nodebranchinfo;
	double *x_lb=activenode.x_lb;
	double *x_ub=activenode.x_ub;
	double *y_bar_lb=activenode.y_bar_lb;
	double *y_bar_ub=activenode.y_bar_ub;
	//handle the McCormick constraints
	if (activenode.mcc_index_col!=-1 && activenode.mcc_index_row!=-1)
	{
		for (i=0;i<param_n;i++)
		{
			if (activenode.x_start_index!=-1)
			{
				if (x_lb[i]>-CPX_INFBOUND)
				{
					lu='L';
					branch_index_1=activenode.x_start_index+i;
					bd=x_lb[i];
					status=CPXtightenbds(env,rx_lp,1,&branch_index_1,&lu,&bd);
					if(status) goto TERMINATE;
				}
				if (x_ub[i]<CPX_INFBOUND)
				{
					lu='U';
					branch_index_1=activenode.x_start_index+i;
					bd=x_ub[i];
					status=CPXtightenbds(env,rx_lp,1,&branch_index_1,&lu,&bd);
					if(status) goto TERMINATE;
				}
			}
			if (activenode.y_bar_start_index!=-1)
			{
				if (y_bar_lb[i]>-CPX_INFBOUND)
				{
					lu='L';
					branch_index_2=activenode.y_bar_start_index+i;
					bd=y_bar_lb[i];
					status=CPXtightenbds(env,rx_lp,1,&branch_index_2,&lu,&bd);
					if(status) goto TERMINATE;
				}
				if (y_bar_ub[i]<CPX_INFBOUND)
				{
					lu='U';
					branch_index_2=activenode.y_bar_start_index+i;
					bd=y_bar_ub[i];
					status=CPXtightenbds(env,rx_lp,1,&branch_index_2,&lu,&bd);
					if(status) goto TERMINATE;
				}
			}
			row_indice[0]=activenode.mcc_index_row+4*i;
			row_indice[1]=activenode.mcc_index_row+4*i+1;
			row_indice[2]=activenode.mcc_index_row+4*i+2;
			row_indice[3]=activenode.mcc_index_row+4*i+3;
			col_indice[0]=activenode.mcc_index_col+4*i;
			col_indice[1]=activenode.mcc_index_col+4*i+1;
			col_indice[2]=activenode.mcc_index_col+4*i+2;
			col_indice[3]=activenode.mcc_index_col+4*i+3;
			if (x_lb[i]>-CPX_INFBOUND && y_bar_lb[i]>-CPX_INFBOUND)
			{
				rhs[0]=x_lb[i]*y_bar_lb[i];
				status=CPXchgcoef(env,rx_lp,row_indice[0],i,y_bar_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,rx_lp,row_indice[0],activenode.y_bar_start_index+i,x_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,rx_lp,row_indice[0],col_indice[0],1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,rx_lp,row_indice[0],activenode.y_bar_start_index+param_n+i,-1);
				if (status) goto TERMINATE;
			}else
			{
				fprintf(stderr,"unable to set bound as CPLEX_INFBOUND\n");
				status=-1;
				goto TERMINATE;
			}
			if (x_ub[i]<CPX_INFBOUND && y_bar_ub[i]<CPX_INFBOUND)
			{
				rhs[1]=x_ub[i]*y_bar_ub[i];
				status=CPXchgcoef(env,rx_lp,row_indice[1],i,y_bar_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,rx_lp,row_indice[1],activenode.y_bar_start_index+i,x_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,rx_lp,row_indice[1],col_indice[1],1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,rx_lp,row_indice[1],activenode.y_bar_start_index+param_n+i,-1);
				if (status) goto TERMINATE;
			}else
			{
				fprintf(stderr,"unable to set bound as CPLEX_INFBOUND\n");
				status=-1;
				goto TERMINATE;
			}
			if (x_lb[i]>-CPX_INFBOUND && y_bar_ub[i]<CPX_INFBOUND)
			{
				rhs[2]=x_lb[i]*y_bar_ub[i];
				status=CPXchgcoef(env,rx_lp,row_indice[2],i,y_bar_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,rx_lp,row_indice[2],activenode.y_bar_start_index+i,x_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,rx_lp,row_indice[2],col_indice[2],-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,rx_lp,row_indice[2],activenode.y_bar_start_index+param_n+i,-1);
				if (status) goto TERMINATE;
			}else
			{
				fprintf(stderr,"unable to set bound as CPLEX_INFBOUND\n");
				status=-1;
				goto TERMINATE;
			}

			if (x_ub[i]<CPX_INFBOUND && y_bar_lb[i]>-CPX_INFBOUND)
			{
				rhs[3]=x_ub[i]*y_bar_lb[i];
				status=CPXchgcoef(env,rx_lp,row_indice[3],i,y_bar_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,rx_lp,row_indice[3],activenode.y_bar_start_index+i,x_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,rx_lp,row_indice[3],col_indice[3],-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,rx_lp,row_indice[3],activenode.y_bar_start_index+param_n+i,-1);
				if (status) goto TERMINATE;
			}else
			{
				fprintf(stderr,"unable to set bound as CPLEX_INFBOUND\n");
				status=-1;
				goto TERMINATE;
			}
			status=CPXchgrhs(env,rx_lp,4,row_indice,rhs);
			if(status) goto TERMINATE;
		}
	}
	//handle the complementarity constraints
	for (i= 0; i<param_m; i++)
	{
		if(nodebranchinfo[i]==0)
		{
			branch_index_1=param_n+i;
			branch_index_2=param_n+param_m+i;
			lu='U';
			bd=CPX_INFBOUND;
			status = CPXtightenbds(env, rx_lp, 1, &branch_index_1, &lu, &bd);
			if(status) goto TERMINATE;
			lu='U';
			bd=CPX_INFBOUND;
			status = CPXtightenbds(env, rx_lp, 1, &branch_index_2, &lu, &bd);
			if(status) goto TERMINATE;
		}else if(nodebranchinfo[i]==1)
		{
			branch_index_1=param_n+i;
			branch_index_2=param_n+param_m+i;
			lu='U';
			bd= 0.0;
			status = CPXtightenbds(env, rx_lp, 1, &branch_index_1, &lu, &bd);
			if(status) goto TERMINATE;
			lu='U';
			bd=CPX_INFBOUND;
			status = CPXtightenbds(env, rx_lp, 1, &branch_index_2, &lu, &bd);
			if(status) goto TERMINATE;
		}else if(nodebranchinfo[i]==2)
		{
			branch_index_1=param_n+i;
			branch_index_2=param_n+param_m+i;
			lu='U';
			bd=CPX_INFBOUND;
			status = CPXtightenbds(env, rx_lp, 1, &branch_index_1, &lu, &bd);
			if(status) goto TERMINATE;
			lu='U';
			bd=0.0;
			status = CPXtightenbds(env, rx_lp, 1, &branch_index_2, &lu, &bd);
			if(status) goto TERMINATE;
		}else if (nodebranchinfo[i]==3)
		{
			branch_index_1=param_n+i;
			branch_index_2=param_n+param_m+i;
			lu='U';
			bd=0.0;
			status = CPXtightenbds(env, rx_lp, 1, &branch_index_1, &lu, &bd);
			if(status) goto TERMINATE;
			lu='U';
			bd=0.0;
			status = CPXtightenbds(env, rx_lp, 1, &branch_index_2, &lu, &bd);
			if(status) goto TERMINATE;
		}else
		{
			status = -1;
			fprintf(stderr,"setupnode_rx_lp(): nodebranchinfo[%d] = %d",i,nodebranchinfo[i]);
		}
	}
TERMINATE:
	return(status);
}
int CPX_LPCCSolver(CPXENVptr env,
				   CPXLPptr lp, 
				   const int param_n, 
				   const int param_m,
				   const int param_k,
				   double *start_point)
{
	int status = 0;
	int i;
	const int intzero = 0;
	const int var_cnt=param_n+param_m+param_m+param_k;
	int matind;
	int *indices=NULL;
	const double dblzero = 0.0;
	const double dblone  = 1.0;
	const int lp_col = CPXgetnumcols (env, lp);
	const char     binary = 'B';
	const char sense='E';
   /* Add fixed charge variables */
	for (i = 0; i < param_m; i++) 
	{
		status = CPXaddcols (env, lp, 1, 0, NULL, &intzero, NULL, NULL,&dblzero, &dblone, NULL);
		if ( status )
		{
			 fprintf (stderr, "Failed to add fixed charge variable.\n");
			goto TERMINATE;
		}

		matind = lp_col+i;
		status = CPXchgctype (env, lp, 1, &matind, &binary);
		if ( status ) 
		{
			fprintf (stderr, "Failed to change variable type.\n");
			goto TERMINATE;
		}
	}
    /* Add indicator constraints*/
	for (i = 0; i < param_m; i++) 
	{
		matind=param_n+i;
		status = CPXaddindconstr (env, lp, lp_col+i, 1, 1, 0.0,sense, &matind, &dblone, NULL);
		if ( status ) 
		{
			fprintf (stderr, "Failed to add indicator constraint.");
			goto TERMINATE;
		}
		matind=param_n+param_m+i;
		status = CPXaddindconstr (env, lp, lp_col+i, 0, 1, 0.0, sense, &matind, &dblone, NULL);
		if ( status ) 
		{
			fprintf (stderr, "Failed to add indicator constraint.");
			goto TERMINATE;
		}
	}

	status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
	if ( status )
	{
		fprintf (stderr, 
               "Failure to turn on screen indicator, error %d.\n", status);
		goto TERMINATE;
	}
	status = CPXsetdblparam(env, CPX_PARAM_TILIM,7200);
	if(status) goto TERMINATE;
	indices=(int*)malloc(var_cnt*sizeof(int));
	if (indices==NULL)
	{
		status = NO_MEMORY;	
		fprintf (stderr, "Could not allocate memory.\n");
		goto TERMINATE;
	}
	for (i=0;i<var_cnt;i++)
	{
		indices[i]=i;
	}
	if(start_point!=NULL)
	{
		fprintf(stdout,"copy start feasible solution to CPLEX\n");
		status=CPXcopymipstart(env,lp,var_cnt,indices,start_point);
		if(status) goto TERMINATE;
	}
	status = CPXmipopt (env, lp);
	if ( status ) 
	{
		fprintf (stderr, "Failed to optimize MIP.\n");
		goto TERMINATE;
	}
TERMINATE:
	free_and_null((char**) &indices);
	return(status);
}


int getVARbound(CPXENVptr env, 
				CPXLPptr _lp, 
				STARTINFO  *lp_start, 
				STARTINFO *update_lp_start_p,
				const int var_index,
				const int maxormin,
				double *bound,
				int *condition)
{
	int status=0;
	double one=1;
	double *lp_soln=NULL;
	int solnstat;
	CPXLPptr lp = NULL;
	*condition=1;
	lp = CPXcloneprob(env,_lp,&status);
	if (status) goto TERMINATE;
	lp_soln=(double*) malloc(CPXgetnumcols(env,lp)*sizeof(double));
	if(lp_soln==NULL)
	{
		status=-1;
		fprintf(stderr,"getVARbound(): unable to allocate memory\n");
		goto TERMINATE;
	}
	status=zeroobjectivefunc(env,lp);
	if(status) goto TERMINATE;	
	status=CPXchgobj(env,lp,1,&var_index,&one);
	if (status) goto TERMINATE;
	CPXchgobjsen(env,lp,maxormin);
	status=cplexSolveLP(env,lp,*lp_start,CPXgetnumcols(env,lp),lp_soln,bound,&solnstat,update_lp_start_p,CPX_ALG_AUTOMATIC,1,2);
	switch(solnstat)
	{
	case CPX_STAT_OPTIMAL:
		break;
	case CPX_STAT_OPTIMAL_INFEAS:
		if(maxormin==CPX_MAX) *bound=CPX_INFBOUND;
		else *bound=-CPX_INFBOUND;
		break;
	case CPX_STAT_NUM_BEST:
		if(maxormin==CPX_MAX) *bound=CPX_INFBOUND;
		else *bound=-CPX_INFBOUND;
		*condition=0;
		break;
	case CPX_STAT_INFEASIBLE:
		if(maxormin==CPX_MAX) *bound=-CPX_INFBOUND;
		else *bound=CPX_INFBOUND;
		*condition=-1;
		break;
	case CPX_STAT_UNBOUNDED:
		if(maxormin==CPX_MAX) *bound=CPX_INFBOUND;
		else *bound=-CPX_INFBOUND;
		break;
	default:
		status=-1;
		fprintf(stderr,"getVARbound(): unexpected solnstat: %d\n",solnstat);
		goto TERMINATE;
	}
TERMINATE:
	if ( lp != NULL ) 
	{
		status = CPXfreeprob (env, &lp);
		if ( status ) fprintf (stderr, "CPXfreeprob failed, error code %d.\n",status);
	}
	free_and_null((char**)&lp_soln);
	return(status);
}

int getVARbound_baropt(CPXENVptr env, 
					   CPXLPptr _lp, 
					   const int var_index,
					   const int maxormin,
					   double *bound,
					   int *condition)
{
	int status=0;
	double one=1;
	int lp_solnstat;
	CPXLPptr lp = NULL;
	*condition=1;
	lp = CPXcloneprob(env,_lp,&status);
	if (status) goto TERMINATE;
	status=zeroobjectivefunc(env,lp);
	if(status) goto TERMINATE;	
	status=CPXchgobj(env,lp,1,&var_index,&one);
	if (status) goto TERMINATE;
	CPXchgobjsen(env,lp,maxormin);
	status=CPXbaropt(env,lp);
	status = CPXsolution (env, lp, &lp_solnstat,bound, NULL, NULL, NULL, NULL);
	if ( status ) goto TERMINATE;
	switch(lp_solnstat)
	{
	case CPX_STAT_OPTIMAL:
		break;
	case CPX_STAT_OPTIMAL_INFEAS:
		if(maxormin==CPX_MAX) *bound=CPX_INFBOUND;
		else *bound=-CPX_INFBOUND;
		break;
	case CPX_STAT_NUM_BEST:
		if(maxormin==CPX_MAX) *bound=CPX_INFBOUND;
		else *bound=-CPX_INFBOUND;
		*condition=0;
		break;
	case CPX_STAT_INFEASIBLE:
		if(maxormin==CPX_MAX) *bound=-CPX_INFBOUND;
		else *bound=CPX_INFBOUND;
		*condition=-1;
		break;
	case CPX_STAT_UNBOUNDED:
		if(maxormin==CPX_MAX) *bound=CPX_INFBOUND;
		else *bound=-CPX_INFBOUND;
		break;
	case CPX_STAT_ABORT_PRIM_OBJ_LIM:
		if(maxormin==CPX_MAX) *bound=CPX_INFBOUND;
		else *bound=-CPX_INFBOUND;
		break;
	default:
		status=-1;
		fprintf(stderr,"getVARbound_baropt(): unexpected solnstat: %d\n",lp_solnstat);
		goto TERMINATE;
	}
TERMINATE:
	if ( lp != NULL ) 
	{
		status = CPXfreeprob (env, &lp);
		if ( status ) fprintf (stderr, "CPXfreeprob failed, error code %d.\n",status);
	}
	return(status);
}
