#define _CRT_SECURE_NO_DEPRECATE
#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "fcgimisc.h"
#include "LPCC.h"
#include "cutsgenerator.h"
#include "utilities.h"
#include "CPLEXcontrol.h"
#include <time.h>
extern int nodecnt;
extern int leftnodes;

int CGLP_cutgenerator(CPXENVptr env, 
					  CPXLPptr cg_lp, 
					  STARTINFO *start_p, 
					  double *q_coef,
					  MATRIX matrix_N,
					  MATRIX matrix_M,
					  const int violation_index, 
					  const int param_k,
					  const int param_m,
					  const int param_n,
					  double *cutoff_point, 
					  double *alpha, 
					  double *beta, 
					  int *solnstat_p)
{
	int status=0;
	int i;
	const int cg_lp_n_col=CPXgetnumcols(env,cg_lp);
	const int cg_lp_n_row=CPXgetnumrows(env,cg_lp);
	int lp_solnstat;
	double lp_obj;
	double *lp_soln=NULL;
	STARTINFO updated_start;
	lp_soln=(double*)malloc((param_n+param_m+1)*sizeof(double));
	if (lp_soln==NULL)
	{
		status=-1;
		fprintf(stderr," CGLP_cutgenerator(): unable to allocate memory\n");
		goto TERMINATE;
	}
	init_startinfo(&updated_start);
	//change objective coef
	status=CPXchgcoef(env,cg_lp,-1,0,-1);
	if(status) goto TERMINATE;
	for (i=0;i<(param_n+param_m);i++)
	{
		status=CPXchgcoef(env,cg_lp,-1,(i+1),cutoff_point[i]);
		if(status) goto TERMINATE;
	}
	//add new col related with violation complementarity
	status=CPXnewcols(env,cg_lp,2,NULL,NULL,NULL,NULL,NULL);
	if(status) goto TERMINATE;
	//v1
	status=CPXchgcoef(env,cg_lp,(param_n+violation_index),cg_lp_n_col,-1);
	if(status) goto TERMINATE;
	status=CPXchgcoef(env,cg_lp,(2*param_n+2*param_m+2),cg_lp_n_col,1);
	if(status) goto TERMINATE;
	//v2
	for (i=0;i<matrix_N.matcnt[violation_index];i++)
	{
		status=CPXchgcoef(env,cg_lp,(param_n+param_m+matrix_N.matind[matrix_N.matbeg[violation_index]+i]),cg_lp_n_col+1,-matrix_N.matval[matrix_N.matbeg[violation_index]+i]);
		if(status) goto TERMINATE;
	}
	for (i=0;i<matrix_M.matcnt[violation_index];i++)
	{
		status=CPXchgcoef(env,cg_lp,(param_n+param_m+param_n+matrix_M.matind[matrix_M.matbeg[violation_index]+i]),cg_lp_n_col+1,-matrix_M.matval[matrix_M.matbeg[violation_index]+i]);
		if(status) goto TERMINATE;
	}
	status=CPXchgcoef(env,cg_lp,(param_n+param_m+param_n+param_m+1),cg_lp_n_col+1,q_coef[violation_index]);
	if(status) goto TERMINATE;
	status=CPXchgcoef(env,cg_lp,(2*param_n+2*param_m+2),cg_lp_n_col+1,1);
	if(status) goto TERMINATE;
	//
	status=cplexSolveLP(env,cg_lp,*start_p,(param_n+param_m+1),lp_soln,&lp_obj,&lp_solnstat,&updated_start,CPX_ALG_AUTOMATIC,1,2);
	if(status) goto TERMINATE;
	if (lp_solnstat == CPX_STAT_OPTIMAL) 
	{
		if (lp_obj < -ZERO_TOLERANCE )
		{
			*solnstat_p = 1;
			for (i=0; i< (param_n+param_m); i++) 
			{
				alpha[i] = lp_soln[1+i];
			}
			//weaken the cuts
			*beta = weaken_value(lp_soln[0]);
		}else
		{
			*solnstat_p = 0;
		}
	}else
	{
		status = -1;
		fprintf(stderr," CGLP_cutgenerator(): unknown problem happened during CPLEX solving, solnstat == %d",lp_solnstat);
		goto TERMINATE;
	}
	status=copy_startinfo(&updated_start,start_p);
	if(status) goto TERMINATE;
	//recover cg_lp
	status=CPXdelcols(env,cg_lp,cg_lp_n_col,cg_lp_n_col+1);
	if(status) goto TERMINATE;

TERMINATE:
	free_and_null((char**)&lp_soln);
	free_startinfo(&updated_start);
	return(status);
}

int compute_UB(CPXENVptr env, 
			   CPXLPptr rx_lp, 	
			   STARTINFO rx_lp_start,								 
			   const int param_n,								 
			   const int param_m, 								 
			   const int param_k,	
			   const double lpcc_ub,
			   const int num_var,
			   const char y_w,
			   const int i_index,
			   int param_p,
			   double *UB)
{
	int status =0;
	int index_1;
	int index_2;
	int i,j;
	int lp_solnstat;
	int tmp_index;
	int violation_index;
	int divided_value;
	int* select_indice=NULL;
	char lu;
	double* lp_soln=NULL;
	double* select_violation=NULL;
	double* initial_bd=NULL;
	double tmp_violation;
	double violation;
	double bd;
	double lp_obj;
	double max_UB;
	const double one=1;
	const double zero=0;
	const double negative_one=-1;
	CPXLPptr copy_rx_lp = NULL;
	STARTINFO lp_start_1;
	STARTINFO lp_start_2;
	lp_soln=(double*) malloc(num_var*sizeof(double));
	if (lp_soln==NULL)
	{
		status=-1;
		fprintf(stderr," compute_yw_UB(): unable to allocate memory\n");
		goto TERMINATE;
	}
	if (param_p>=1)
	{
		select_indice=(int*) calloc(param_p,sizeof(int));
		select_violation=(double*) malloc(param_p*sizeof(double));
		initial_bd=(double*) malloc(param_p*sizeof(double));
		if (select_indice==NULL || select_violation==NULL || initial_bd==NULL)
		{
			status=-1;
			fprintf(stderr," compute_UB(): unable to allocate memory\n");
			goto TERMINATE;
		}
		for (i=0;i<param_p;i++)
		{
			select_violation[i]=-CPX_INFBOUND;
		}
	}
	init_startinfo(&lp_start_1);
	init_startinfo(&lp_start_2);
	status=copy_startinfo(&rx_lp_start,&lp_start_1);
	if(status) goto TERMINATE;
	copy_rx_lp = CPXcloneprob(env,rx_lp,&status);
	if (status) {
		fprintf(stderr," compute_yw_UB(): unable to copy rx_lp \n");
		goto TERMINATE;
	}
	if (lpcc_ub<CPX_INFBOUND)
	{
		status=add_ub_constraint(env,copy_rx_lp,lpcc_ub);
		if(status) goto TERMINATE;
	}
	//CPXchgobjsen(env,copy_rx_lp,CPX_MAX);
	status=zeroobjectivefunc(env,copy_rx_lp);
	if(status) goto TERMINATE;

	switch(y_w)
	{
	case 'y':
		index_1=param_n+i_index;
		index_2=param_n+param_m+i_index;
		break;
	case 'w':
		index_1=param_n+param_m+i_index;
		index_2=param_n+i_index;
		break;
	default:
		status=-1;
		fprintf(stderr," compute_UB(): unexpected y_w: %c",y_w);
		goto TERMINATE;
	}
	status=CPXchgobj(env,copy_rx_lp,1,&index_1,&negative_one);
	if(status) goto TERMINATE;
	lu='U';
	bd=0.0;
	status=CPXtightenbds(env,copy_rx_lp,1,&index_2,&lu,&bd);
	if(status) goto TERMINATE;
	status=cplexSolveLP(env,copy_rx_lp,lp_start_1,num_var,lp_soln,&lp_obj,&lp_solnstat,&lp_start_2,CPX_ALG_AUTOMATIC,1,2);
	if(status) goto TERMINATE;
	lp_obj=-lp_obj;
	if (lp_solnstat==CPX_STAT_OPTIMAL)
	{
		*UB=lp_obj;
		if (param_p>=1)
		{
			//select the most param_p violation
			for (i=0;i<param_m;i++)
			{
				if (i!=i_index)
				{
					violation=lp_soln[param_n+i]*lp_soln[param_n+param_m+i];
					violation_index=i;
					//fprintf(stdout,"%d: %f\n",i,violation);
					for (j=0;j<param_p;j++)
					{
						if(violation>select_violation[j])
						{
							tmp_violation=select_violation[j];
							select_violation[j]=violation;
							violation=tmp_violation;
							tmp_index=select_indice[j];
							select_indice[j]=violation_index;
							violation_index=tmp_index;
						}
					}
				}
			}
			//fprintf(stdout," biggest %d violation:\n",param_p);
			//for (i=0;i<param_p;i++)
			//{
			//	fprintf(stdout," %d: %f\n",select_indice[i],select_violation[i]);
			//}
			max_UB=-CPX_INFBOUND;
			for (i=0;i<pow(2,param_p);i++)
			{
				divided_value=i;
				//fprintf(stdout,"%d: ",i);
				for (j=0;j<param_p;j++)
				{
					if ((divided_value%2)==0)
					{
						lu='U';
						bd=0.0;
						tmp_index=param_n+select_indice[j];
						status=CPXgetub(env,copy_rx_lp,&(initial_bd[j]),tmp_index,tmp_index);
						if(status) goto TERMINATE;
						status=CPXtightenbds(env,copy_rx_lp,1,&tmp_index,&lu,&bd);
						if(status) goto TERMINATE;
					}else
					{
						lu='U';
						bd=0.0;
						tmp_index=param_n+param_m+select_indice[j];
						status=CPXgetub(env,copy_rx_lp,&(initial_bd[j]),tmp_index,tmp_index);
						if(status) goto TERMINATE;
						status=CPXtightenbds(env,copy_rx_lp,1,&tmp_index,&lu,&bd);
						if(status) goto TERMINATE;
					}
					//fprintf(stdout," %d",divided_value%2);
					divided_value/=2;
				}
				//fprintf(stdout,"\n");
				status=cplexSolveLP(env,copy_rx_lp,lp_start_1,num_var,lp_soln,&lp_obj,&lp_solnstat,&lp_start_2,CPX_ALG_AUTOMATIC,1,2);
				if(status) goto TERMINATE;
				lp_obj=-lp_obj;
				if (lp_solnstat==CPX_STAT_OPTIMAL)
				{
					status=copy_startinfo(&lp_start_2,&lp_start_1);
					if(status) goto TERMINATE;
					if (max_UB<lp_obj)
					{
						max_UB=lp_obj;
					}
				}else if (lp_solnstat==CPX_STAT_INFEASIBLE)
				{

				}else
				{
					status=-1;
					fprintf(stderr," compute_yw_UB(): unexpected solnstat: %d",lp_solnstat);
					goto TERMINATE;
				}
				//recover lp
				divided_value=i;
				//fprintf(stdout,"%d: ",i);
				for (j=0;j<param_p;j++)
				{
					if ((divided_value%2)==0)
					{
						lu='U';
						bd=initial_bd[j];
						tmp_index=param_n+select_indice[j];
						status=CPXtightenbds(env,copy_rx_lp,1,&tmp_index,&lu,&bd);
					}else
					{
						lu='U';
						bd=initial_bd[j];
						tmp_index=param_n+param_m+select_indice[j];
						status=CPXtightenbds(env,copy_rx_lp,1,&tmp_index,&lu,&bd);
					}
					//fprintf(stdout," %d",divided_value%2);
					divided_value/=2;
				}

			}
			if (max_UB>-CPX_INFBOUND)
			{
				//fprintf(stdout," %f, imroved bd: %f\n",*UB,max_UB);
				*UB=max_UB;
			}else
			{
				*UB=0;
				//fprintf(stdout," %f\n",*UB);
			}
		}
	}else if (lp_solnstat==CPX_STAT_UNBOUNDED)
	{
		*UB=CPX_INFBOUND;
		//fprintf(stdout," %f\n",*UB);
	}else if (lp_solnstat==CPX_STAT_INFEASIBLE)
	{
		*UB=0;
		//fprintf(stdout," %f\n",*UB);
	}else
	{
		status=-1;
		fprintf(stderr," compute_yw_UB(): unexpected solnstat: %d",lp_solnstat);
		goto TERMINATE;
	}
TERMINATE:
	free_and_null((char**)&select_indice);
	free_and_null((char**)&lp_soln);
	free_and_null((char**)&select_violation);
	free_and_null((char**)&initial_bd);
	free_startinfo(&lp_start_1);
	free_startinfo(&lp_start_2);
	if ( copy_rx_lp != NULL ) {
		status = CPXfreeprob (env, &copy_rx_lp);
		if ( status ) {  
			fprintf (stderr, " compute_yw_UB(): CPXfreeprob failed\n");
		}
	}
	return (status);
}
int generateNodeSimpleCuts(CPXENVptr env,
						   CPXLPptr rx_lp, 
						   const int param_n,
						   const int param_m, 
						   const int param_k,
						   const int violate_complementary_index, 
						   const double violation_complementary_1, 
						   const double violation_complementary_2,
						   NODE **branchnode_ptr)
{
	int status=0;
	const int num_cols=CPXgetnumcols(env,rx_lp);
	const int num_rows=CPXgetnumrows(env,rx_lp);
	int i;
	int row_id;
	int *ind = NULL;
	int intzero=0;
	int nnz=0;
	double one=1.0;
	double *val = NULL;
	double *tableau_1=NULL;
	double *tableau_2=NULL;
	double coef_1;
	double coef_2;
	tableau_1 = (double*) malloc(num_cols*sizeof(double));
	tableau_2 = (double*) malloc(num_cols*sizeof(double));
	ind = (int*) malloc(num_cols*sizeof(int));
	val = (double*) malloc(num_cols*sizeof(double));
	if (tableau_1==NULL || tableau_2==NULL || ind ==NULL || val== NULL)
	{
		status=NO_MEMORY;
		fprintf(stderr," generateSimpleCuts():could not allocate memory.\n");
		goto TERMINATE;
	}
	status = CPXgetijrow(env,rx_lp,-1,param_n+violate_complementary_index,&row_id);
	if (status) goto TERMINATE;
	status = CPXbinvarow(env,rx_lp,row_id,tableau_1);
	if (status) goto TERMINATE;
	status = CPXgetijrow(env,rx_lp,-1,param_n+param_m+violate_complementary_index,&row_id);
	if (status) goto TERMINATE;
	status = CPXbinvarow(env,rx_lp,row_id,tableau_2);
	for (i=0;i<num_cols;i++)
	{
		if (i==param_n+violate_complementary_index || i==param_n+param_m+violate_complementary_index)
		{
			continue;
		}
		//if (tableau_1[i]<EP && tableau_2[i]<EP)
		//{
		//	continue;
		//}
		coef_1= tableau_1[i]/violation_complementary_1;
		coef_2= tableau_2[i]/violation_complementary_2;
		ind[nnz]=i;
		val[nnz]=coef_1>coef_2 ? coef_1:coef_2;
		nnz++;
	}
	status=add_cuts(&((*branchnode_ptr)->nodecuts),num_cols,nnz,ind,val,1);
	if(status) goto TERMINATE;
	status = CPXaddrows(env,rx_lp,0,1,nnz,&one,NULL,&intzero,ind,val,NULL,NULL);
	if (status) goto TERMINATE;
	status = CPXnewcols(env,rx_lp,1,NULL,NULL,NULL,NULL,NULL);
	if (status) goto TERMINATE;
	status = CPXchgcoef(env,rx_lp,CPXgetnumrows(env,rx_lp)-1,CPXgetnumcols(env,rx_lp)-1,-1);
	if (status) goto TERMINATE;
TERMINATE:
	free_and_null((char**)&tableau_1);
	free_and_null((char**)&tableau_2);
	free_and_null((char**)&ind);
	free_and_null((char**)&val);
	return(status);
}

int simplecuts_generator(CPXENVptr env, 
						 CPXLPptr rx_lp, 
						 const int param_n,
						 const int param_m,
						 const int param_k,
						 const int violate_complementary_index,
						 const double violation_complementary_1, 
						 const double violation_complementary_2, 
						 CONSTRAINT_SET *cut_ptr)
{
	int status=0;
	const int num_cols=CPXgetnumcols(env,rx_lp);
	const int num_rows=CPXgetnumrows(env,rx_lp);
	int i;
	int row_id;
	int *ind = NULL;
	int intzero=0;
	int nnz=0;
	double one=1.0;
	double *val = NULL;
	double *tableau_1=NULL;
	double *tableau_2=NULL;
	double coef_1;
	double coef_2;
	double var_lb;
	tableau_1 = (double*) malloc(num_cols*sizeof(double));
	tableau_2 = (double*) malloc(num_cols*sizeof(double));
	ind = (int*) malloc(num_cols*sizeof(int));
	val = (double*) malloc(num_cols*sizeof(double));
	if (tableau_1==NULL || tableau_2==NULL || ind ==NULL || val== NULL)
	{
		status=NO_MEMORY;
		fprintf(stderr," SimpleCuts_generator():could not allocate memory.\n");
		goto TERMINATE;
	}
	status = CPXgetijrow(env,rx_lp,-1,param_n+violate_complementary_index,&row_id);
	if (status) goto TERMINATE;
	status = CPXbinvarow(env,rx_lp,row_id,tableau_1);
	if (status) goto TERMINATE;
	status = CPXgetijrow(env,rx_lp,-1,param_n+param_m+violate_complementary_index,&row_id);
	if (status) goto TERMINATE;
	status = CPXbinvarow(env,rx_lp,row_id,tableau_2);
	for (i=0;i<num_cols;i++)
	{
		if (i==param_n+violate_complementary_index || i==param_n+param_m+violate_complementary_index) continue;
		//if (tableau_1[i]<EP && tableau_2[i]<EP)
		//{
		//	continue;
		//}
		status=CPXgetlb(env,rx_lp,&var_lb,i,i);
		if (var_lb<-ZERO_TOLERANCE && (fabs(tableau_1[i])>ZERO_TOLERANCE || fabs(tableau_2[i])>ZERO_TOLERANCE))
		{
			fprintf(stdout,"nonbasic variable is free variable and its coef is nonzero, unable to generate simple cuts\n");
			goto TERMINATE;
		}
		coef_1= tableau_1[i]/violation_complementary_1;
		coef_2= tableau_2[i]/violation_complementary_2;
		//if (fabs(coef_1>coef_2 ? coef_1:coef_2)<=1e-6) 
		//{
		//	//fprintf(stdout,"%f\n",coef_1>coef_2 ? coef_1:coef_2);
		//	continue;
		//}
		ind[nnz]=i;
		val[nnz]=coef_1>coef_2 ? coef_1:coef_2;
		nnz++;	
	}
	status=add_cuts(cut_ptr,num_cols,nnz,ind,val,1-0.001);
	if(status) goto TERMINATE;
TERMINATE:
	free_and_null((char**)&tableau_1);
	free_and_null((char**)&tableau_2);
	free_and_null((char**)&ind);
	free_and_null((char**)&val);
	return(status);
}


int add_disjunctivecuts(CPXENVptr env, 
					   CPXLPptr rx_lp, 	
					   CPXLPptr cg_lp,
					   double *q_coef,
					   MATRIX matrix_N,
					   MATRIX matrix_M,
					   STARTINFO *lp_start,								 
					   const int param_n,								 
					   const int param_m, 								 
					   const int param_k,	
					   double param_p,
					   double *opt_soln,
					   double *opt_obj,
					   int *lpcc_stat)
{
	//*lpcc_stat=1 means find optimal solution
	//*lpcc_stat=0 means haven't found optimal solution
	//*lpcc_stat=-1 mean lpcc is infeasible
	int status =0;
	int i,j,k;
	int lp_solnstat;
	int loop_limit;
	int cglp_solnstat;
	int cnt;
	int cut_cnt;
	int *sort_index=NULL;
	int cut_limit;
	const int initial_n_col=CPXgetnumcols(env,rx_lp);
	const int initial_n_row=CPXgetnumrows(env,rx_lp);
	double *sort_val=NULL;
	double init_lp_obj;
	const int num_var= param_n+param_m+param_m+param_k;
	double *lp_soln=NULL;
	double *alpha=NULL;
	double beta;
	double lp_obj;
	STARTINFO rx_lp_start;
	STARTINFO updated_rx_lp_start;
	STARTINFO *cg_lp_start=NULL;
	clock_t t_start,t_end;
	*lpcc_stat=LPCC_PRECONDITION_BOUNDED;
	cg_lp_start=(STARTINFO*) malloc(param_m*sizeof(STARTINFO));
	sort_val=(double*) malloc(param_m*sizeof(double));
	sort_index=(int*)malloc(param_m*sizeof(int));
	if (cg_lp_start==NULL || sort_index==NULL || sort_val==NULL)
	{
		status=-1;
		fprintf(stderr," add_disjunctivecuts(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<param_m;i++)
	{
		init_startinfo(&(cg_lp_start[i]));
	}
	init_startinfo(&rx_lp_start);
	init_startinfo(&updated_rx_lp_start);
	status=copy_startinfo(lp_start,&rx_lp_start);
	if(status) goto TERMINATE;
	lp_soln=(double*) malloc(num_var*sizeof(double));
	alpha=(double*)malloc((param_n+param_m)*sizeof(double));
	if(lp_soln==NULL || alpha==NULL)
	{
		status = NO_MEMORY;	
		fprintf (stderr, " add_disjunctivecuts(): could not allocate memory.\n");
		goto TERMINATE;
	}
	loop_limit=(int) floor(param_m*param_p);
	cut_cnt=0;
	fprintf(stdout," Add disjunctive cut at root node...(loop=%d)\n",loop_limit);
	t_start=clock();
	for (i=0;i<loop_limit;i++)
	{
		status = cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat, &updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,2);
		if (status) goto TERMINATE;
		switch(lp_solnstat)
		{
		case CPX_STAT_OPTIMAL:
			t_end=clock();
			if (i==0) init_lp_obj=lp_obj;
			else 
			{
				if (i==((int)floor(param_m/8)-1))
				{
					fprintf(stdout,"	m/8 loop:");
					fprintf(stdout," LB: %12.4f PT: %6.0fs\n",lp_obj,(double)(t_end-t_start)/CLOCKS_PER_SEC);
				}
				if (i==((int)floor(param_m/4)-1))
				{
					fprintf(stdout,"	m/4 loop:");
					fprintf(stdout," LB: %12.4f PT: %6.0fs\n",lp_obj,(double)(t_end-t_start)/CLOCKS_PER_SEC);
				}
				if (i==((int)floor(param_m/3)-1))
				{
					fprintf(stdout,"	m/3 loop:");
					fprintf(stdout," LB: %12.4f PT: %6.0fs\n",lp_obj,(double)(t_end-t_start)/CLOCKS_PER_SEC);
				}
				if (i==((int)floor(param_m/2)-1))
				{
					fprintf(stdout,"	m/2 loop:");
					fprintf(stdout," LB: %12.4f PT: %6.0fs\n",lp_obj,(double)(t_end-t_start)/CLOCKS_PER_SEC);
				}
			}
			cnt=0;
			status=copy_startinfo(&updated_rx_lp_start,&rx_lp_start);
			if(status)goto TERMINATE;
			for (j=0;j<param_m;j++)
			{
				sort_val[j]=lp_soln[param_n+j]*lp_soln[param_n+param_m+j];
				sort_index[j]=j;
				if (lp_soln[param_n+j]<=ZERO_TOLERANCE || lp_soln[param_n+param_m+j]<=ZERO_TOLERANCE)
				{
					cnt++;
				}
			}
			status=get_sort_index(param_m,sort_val,sort_index,-1);
			if(status) goto TERMINATE;
			if (cnt==param_m)
			{
				//get the feasible lpcc solution 
				*lpcc_stat=LPCC_PRECONDITION_OPT;
				for(j=0;j<num_var;j++)
				{
					opt_soln[j]=lp_soln[j];
				}
				*opt_obj=lp_obj;
			}else
			{
				cut_limit=(cnt<3)?cnt:3;
				for (j=0;j<cut_limit;j++)
				{
					//violation complementarity and generate cuts
					//fprintf(stdout,"%d: %f\n",sort_index[j],sort_val[sort_index[j]]);
					status=CGLP_cutgenerator(env,cg_lp,&(cg_lp_start[sort_index[j]]),q_coef,matrix_N,matrix_M,sort_index[j],param_k,param_m,param_n,lp_soln,alpha,&beta,&cglp_solnstat);
					if(status)goto TERMINATE;
					if (cglp_solnstat==1)
					{
						//generate disjunctive cut and add to lp
						cut_cnt++;
						status=CPXnewrows(env,rx_lp,1,&beta,NULL,NULL,NULL);
						if(status) goto TERMINATE;
						for (k=0;k<(param_n+param_m);k++)
						{
							status=CPXchgcoef(env,rx_lp,(CPXgetnumrows(env,rx_lp)-1),k,alpha[k]);
							if(status) goto TERMINATE;
						}
						status=CPXnewcols(env,rx_lp,1,NULL,NULL,NULL,NULL,NULL);
						if(status) goto TERMINATE;
						status=CPXchgcoef(env,rx_lp,(CPXgetnumrows(env,rx_lp)-1),(CPXgetnumcols(env,rx_lp)-1),-1);
						if(status) goto TERMINATE;
					}
				}
			}
			break;
		case CPX_STAT_INFEASIBLE:
			*lpcc_stat=LPCC_PRECONDITION_INFEASIBLE;
			break;
		default:
			status=-1;
			fprintf(stderr," add_disjunctivecuts(): unexpected lp_solnstat: %d\n",lp_solnstat);
			goto TERMINATE;
		}
		if(*lpcc_stat!=LPCC_PRECONDITION_BOUNDED) break;
	}
	t_end=clock();
	fprintf(stdout,"	Total time for generating cuts: %12.3fs\n",(double)(t_end-t_start)/CLOCKS_PER_SEC);
	switch(*lpcc_stat)
	{
	case LPCC_PRECONDITION_OPT:
		fprintf(stdout,"	After adding disjunctive cuts, find optimal solution %f \n",lp_obj);
		break;
	case LPCC_PRECONDITION_BOUNDED:

		status = cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat, &updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,LP_PRECONDITION_UNKNOWN);
		if (status) goto TERMINATE;
		status=clearcutpool(env,rx_lp,initial_n_row,initial_n_col,cut_cnt);
		if (status) goto TERMINATE;
		free_startinfo(&rx_lp_start);
		init_startinfo(&rx_lp_start);
		status = cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat, &updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,LP_PRECONDITION_UNKNOWN);
		if (status) goto TERMINATE;
		status=copy_startinfo(&updated_rx_lp_start,lp_start);
		if(status)goto TERMINATE;
		fprintf(stdout,"	After adding disjunctive cuts, improved lp relaxation is %f\n",lp_obj);
		break;
	case LPCC_PRECONDITION_INFEASIBLE:
		fprintf(stdout,"	After adding disjunctive cuts, find LPCC is infeasible \n");
		break;
	default:
		status=-1;
		fprintf(stderr," add_disjunctivecuts(): unexpected lpcc_stat: %d\n",*lpcc_stat);
		goto TERMINATE;
	}
TERMINATE:
	free_and_null((char**) &lp_soln);
	free_and_null((char**) &alpha);
	free_startinfo(&rx_lp_start);
	free_startinfo(&updated_rx_lp_start);
	free_and_null((char**) &sort_index);
	free_and_null((char**) &sort_val);
	if (cg_lp_start!=NULL)
	{
		for (i=0;i<param_m;i++)
		{
			free_startinfo(&(cg_lp_start[i]));
		}
	}
	free_and_null((char**) &cg_lp_start);
	return(status);
}


int add_simplecuts(CPXENVptr env, 
				  CPXLPptr rx_lp, 	
				  STARTINFO *lp_start,								 
				  const int param_n,								 
				  const int param_m, 								 
				  const int param_k,	
				  double param_p,
				  double *opt_soln,
				  double *opt_obj,
				  int *lpcc_stat)
{
	//*lpcc_stat=1 means find optimal solution
	//*lpcc_stat=0 means haven't found optimal solution
	//*lpcc_stat=-1 mean lpcc is infeasible
	int rx_lp_n_col;
	int rx_lp_n_row;
	int status =0;
	int i,j;
	int lp_solnstat;
	int loop_limit;
	int cnt;
	int cut_cnt;
	const int initial_rx_lp_col=CPXgetnumcols(env,rx_lp);
	const int initial_rx_lp_row=CPXgetnumrows(env,rx_lp);
	double init_lp_obj;
	const int num_var= param_n+param_m+param_m+param_k;
	double *lp_soln=NULL;
	double lp_obj;
	STARTINFO rx_lp_start;
	STARTINFO updated_rx_lp_start;
	CONSTRAINT_SET cut_set;
	clock_t t_start,t_end;
	*lpcc_stat=LPCC_PRECONDITION_BOUNDED;

	init_startinfo(&rx_lp_start);
	init_startinfo(&updated_rx_lp_start);
	init_constraint_set(&cut_set);
	status=copy_startinfo(lp_start,&rx_lp_start);
	if(status) goto TERMINATE;
	lp_soln=(double*) malloc(num_var*sizeof(double));
	if(lp_soln==NULL)
	{
		status = NO_MEMORY;	
		fprintf (stderr, " add_simplecuts(): could not allocate memory.\n");
		goto TERMINATE;
	}
	loop_limit=(int) floor(param_m*param_p);
	cut_cnt=0;
	fprintf(stdout," Add simple cuts at root node...(loop=%d)\n",loop_limit);
	t_start=clock();
	status = cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat, &updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,LP_PRECONDITION_UNKNOWN);
	if (status) goto TERMINATE;
	for (i=0;i<loop_limit;i++)
	{
		switch(lp_solnstat)
		{
		case CPX_STAT_OPTIMAL:
			t_end=clock();
			if (i==0) init_lp_obj=lp_obj;
			else 
			{
				if (i==((int)floor(param_m/8)-1))
				{
					fprintf(stdout,"	m/8 loop:");
					fprintf(stdout," LB: %12.4f PT: %6.0fs\n",lp_obj,(double)(t_end-t_start)/CLOCKS_PER_SEC);
				}
				if (i==((int)floor(param_m/4)-1))
				{
					fprintf(stdout,"	m/4 loop:");
					fprintf(stdout," LB: %12.4f PT: %6.0fs\n",lp_obj,(double)(t_end-t_start)/CLOCKS_PER_SEC);
				}
				if (i==((int)floor(param_m/3)-1))
				{
					fprintf(stdout,"	m/3 loop:");
					fprintf(stdout," LB: %12.4f PT: %6.0fs\n",lp_obj,(double)(t_end-t_start)/CLOCKS_PER_SEC);
				}
				if (i==((int)floor(param_m/2)-1))
				{
					fprintf(stdout,"	m/2 loop:");
					fprintf(stdout," LB: %12.4f PT: %6.0fs\n",lp_obj,(double)(t_end-t_start)/CLOCKS_PER_SEC);
				}
			}
			cnt=0;
			status=copy_startinfo(&updated_rx_lp_start,&rx_lp_start);
			if(status)goto TERMINATE;
			for (j=0;j<param_m;j++)
			{
				if (lp_soln[param_n+j]>ZERO_TOLERANCE && lp_soln[param_n+param_m+j]>ZERO_TOLERANCE)
				{
					//violation complementarity and generate cuts
					status=simplecuts_generator(env,rx_lp,param_n,param_m,param_k,j,lp_soln[param_n+j],lp_soln[param_n+param_m+j],&cut_set);
					if(status)goto TERMINATE;
					cut_cnt++;
				}else
				{
					cnt++;
				}
			}
			if (cnt==param_m)
			{
				//get the feasible lpcc solution 
				*lpcc_stat=LPCC_PRECONDITION_OPT;
				for(j=0;j<num_var;j++)
				{
					opt_soln[j]=lp_soln[j];
				}
				*opt_obj=lp_obj;
			}else
			{
				//add simple cuts to lp
				if (cut_set.cons_matrix.n_row>0)
				{	
					rx_lp_n_col=CPXgetnumcols(env,rx_lp);
					rx_lp_n_row=CPXgetnumrows(env,rx_lp);
					status = CPXaddrows(env,rx_lp,0,cut_set.cons_matrix.n_row,cut_set.cons_matrix.nnz,cut_set.rhs,NULL,cut_set.cons_matrix.matbeg,cut_set.cons_matrix.matind,cut_set.cons_matrix.matval,NULL,NULL);
					if (status) 
					{
						fprintf(stderr," add_simplecuts(): add row failed\n");
						goto TERMINATE;
					}
					status = CPXnewcols(env, rx_lp, cut_set.cons_matrix.n_row, NULL, NULL, NULL, NULL, NULL);
					if (status) goto TERMINATE;
					for (j=0; j<cut_set.cons_matrix.n_row;j++)
					{
						status = CPXchgcoef (env, rx_lp, rx_lp_n_row+j, rx_lp_n_col+j, -1.0);
						if (status) goto TERMINATE;
					}

					status = cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat, &updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,2);
					if (status) goto TERMINATE;
					status=clearcutpool(env,rx_lp,rx_lp_n_row,rx_lp_n_col,cut_set.cons_matrix.n_row);
					if (status) goto TERMINATE;
					free_startinfo(&rx_lp_start);
					init_startinfo(&rx_lp_start);
					status = cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat, &updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,2);
					if (status) goto TERMINATE;
				}
				//clear cut_set
				free_constraint_set(&cut_set);
				init_constraint_set(&cut_set);
			}
			break;
		case CPX_STAT_INFEASIBLE:
			*lpcc_stat=LPCC_PRECONDITION_INFEASIBLE;
			break;
		default:
			status=-1;
			fprintf(stderr," add_simplecuts(): unexpected lp_solnstat: %d\n",lp_solnstat);
			goto TERMINATE;
		}
		if(*lpcc_stat!=LPCC_PRECONDITION_BOUNDED) break;
	}
	free_startinfo(&rx_lp_start);
	init_startinfo(&rx_lp_start);
	status = cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat, &updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,2);
	t_end=clock();
	fprintf(stdout,"	Total time for generating cuts: %12.3f\n",(double)(t_end-t_start)/CLOCKS_PER_SEC);
	switch(*lpcc_stat)
	{
	case LPCC_PRECONDITION_OPT:
		fprintf(stdout,"	After adding simple cuts, find optimal solution %f \n",lp_obj);
		break;
	case LPCC_PRECONDITION_BOUNDED:
		status=copy_startinfo(&updated_rx_lp_start,lp_start);
		if(status)goto TERMINATE;
		fprintf(stdout,"	After adding simple cuts, improved lp relaxation is %f\n",lp_obj);
		break;
	case LPCC_PRECONDITION_INFEASIBLE:
		fprintf(stdout,"	After adding simple cuts, find LPCC is infeasible \n");
		break;
	default:
		status=-1;
		fprintf(stderr," adddisjunctivecuts(): unexpected lpcc_stat: %d\n",*lpcc_stat);
		goto TERMINATE;
	}
TERMINATE:
	free_and_null((char**) &lp_soln);
	free_startinfo(&rx_lp_start);
	free_startinfo(&updated_rx_lp_start);
	free_constraint_set(&cut_set);
	return(status);
}

//int add_flowcovercut(CPXENVptr env, 
//				   CPXLPptr rx_lp, 	
//				   STARTINFO *lp_start,								 
//				   const int param_n,								 
//				   const int param_m, 								 
//				   const int param_k,	
//				   double param_p,
//				   double *opt_soln,
//				   double *opt_obj,
//				   int *lpcc_stat)
//{
//	//*lpcc_stat=1 means find optimal solution
//	//*lpcc_stat=0 means haven't found optimal solution
//	//*lpcc_stat=-1 mean lpcc is infeasible
//	int rx_lp_n_col;
//	int rx_lp_n_row;
//	int status =0;
//	int i,j;
//	int lp_solnstat;
//	int loop_limit;
//	//int cnt;
//	int cut_cnt;
//	int violation_cnt;
//	double init_lp_obj;
//	const int num_var= param_n+param_m+param_m+param_k;
//	double *lp_soln=NULL;
//	double lp_obj;
//	STARTINFO rx_lp_start;
//	STARTINFO updated_rx_lp_start;
//	CONSTRAINT_SET cut_set;
//	clock_t t_start,t_end;
//	*lpcc_stat=0;
//
//	init_startinfo(&rx_lp_start);
//	init_startinfo(&updated_rx_lp_start);
//	init_constraint_set(&cut_set);
//	status=copy_startinfo(lp_start,&rx_lp_start);
//	if(status) goto TERMINATE;
//	lp_soln=(double*) malloc(num_var*sizeof(double));
//	if(lp_soln==NULL)
//	{
//		status = NO_MEMORY;	
//		fprintf (stderr, " add_simplecuts(): could not allocate memory.\n");
//		goto TERMINATE;
//	}
//	loop_limit=(int) floor(param_m*param_p);
//	cut_cnt=0;
//	violation_cnt=0;
//	fprintf(stdout," add flow cover cuts at root node, loop: %d\n",loop_limit);
//	t_start=clock();
//	//
//	loop_limit=1;
//	//
//	for (i=0;i<loop_limit;i++)
//	{
//		status = cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat, &updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,2);
//		if (status) goto TERMINATE;
//		switch(lp_solnstat)
//		{
//		case CPX_STAT_OPTIMAL:
//			t_end=clock();
//			if (i==((int)floor(param_m/8)-1))
//			{
//				fprintf(stdout," m/8\n");
//			}
//			if (i==((int)floor(param_m/4)-1))
//			{
//				fprintf(stdout," m/4\n");
//			}
//			if (i==((int)floor(param_m/3)-1))
//			{
//				fprintf(stdout," m/3\n");
//			}
//			if (i==((int)floor(param_m/2)-1))
//			{
//				fprintf(stdout," m/2\n");
//			}
//			if (i==0) init_lp_obj=lp_obj;
//			else fprintf(stdout," improved lb:%f, use time: %12.3f\n",lp_obj,(double)(t_end-t_start)/CLOCKS_PER_SEC);
//			status=copy_startinfo(&updated_rx_lp_start,&rx_lp_start);
//			if(status)goto TERMINATE;
//			for (j=0;j<param_m;j++)
//			{
//				if (lp_soln[param_n+j]>EP && lp_soln[param_n+param_m+j]>EP)
//				{
//					violation_cnt++;
//				}
//			}
//			//for (j=0;j<param_m;j++)
//			//{
//			//	if (lp_soln[param_n+j]>EP && lp_soln[param_n+param_m+j]>EP)
//			//	{
//			//		//violation complementarity and generate cuts
//			//		status=SimpleCuts_generator(env,rx_lp,param_n,param_m,param_k,j,lp_soln[param_n+j],lp_soln[param_n+param_m+j],&cut_set);
//			//		if(status)goto TERMINATE;
//			//		cut_cnt++;
//			//	}
//			//}
//			//status=flowcovercut_generator(env,rx_lp,rx_lp_start,param_n,param_m,param_k,lp_soln,&cut_set);
//			//status=ext_simple_cut_generator(env,rx_lp,rx_lp_start,param_n,param_m,param_k,lp_soln,2,&cut_set);
//			//status=ext_simple_cut_generator(env,rx_lp,rx_lp_start,param_n,param_m,param_k,lp_soln,3,&cut_set);
//			//status=ext_simple_cut_generator(env,rx_lp,rx_lp_start,param_n,param_m,param_k,lp_soln,4,&cut_set);
//			//status=ext_simple_cut_generator(env,rx_lp,rx_lp_start,param_n,param_m,param_k,lp_soln,5,&cut_set);
//			//status=ext_simple_cut_generator(env,rx_lp,rx_lp_start,param_n,param_m,param_k,lp_soln,6,&cut_set);
//			//status=ext_simple_cut_generator(env,rx_lp,rx_lp_start,param_n,param_m,param_k,lp_soln,7,&cut_set);
//			//status=ext_simple_cut_generator(env,rx_lp,rx_lp_start,param_n,param_m,param_k,lp_soln,8,&cut_set);
//			if (status)goto TERMINATE;
//			if (violation_cnt==0)
//			{
//				//get the feasible lpcc solution 
//				*lpcc_stat=1;
//				for(j=0;j<num_var;j++)
//				{
//					opt_soln[j]=lp_soln[j];
//				}
//				*opt_obj=lp_obj;
//			}else
//			{
//				//add simple cuts to lp
//				cut_cnt+=cut_set.cons_matrix.n_row;
//				//print_matrix(&(cut_set.cons_matrix));
//				if (cut_set.cons_matrix.n_row>0)
//				{	
//					rx_lp_n_col=CPXgetnumcols(env,rx_lp);
//					rx_lp_n_row=CPXgetnumrows(env,rx_lp);
//					status = CPXaddrows(env,rx_lp,0,cut_set.cons_matrix.n_row,cut_set.cons_matrix.nnz,cut_set.rhs,NULL,cut_set.cons_matrix.matbeg,cut_set.cons_matrix.matind,cut_set.cons_matrix.matval,NULL,NULL);
//					if (status) 
//					{
//						fprintf(stderr," add_flowcovercut(): add row failed\n");
//						goto TERMINATE;
//					}
//					status = CPXnewcols(env, rx_lp, cut_set.cons_matrix.n_row, NULL, NULL, NULL, NULL, NULL);
//					if (status) goto TERMINATE;
//					for (j=0; j<cut_set.cons_matrix.n_row;j++)
//					{
//						status = CPXchgcoef (env, rx_lp, rx_lp_n_row+j, rx_lp_n_col+j, -1.0);
//						if (status) goto TERMINATE;
//					}
//				}
//				//clear cut_set
//				free_constraint_set(&cut_set);
//				init_constraint_set(&cut_set);
//			}
//			break;
//		case CPX_STAT_INFEASIBLE:
//			*lpcc_stat=-1;
//			break;
//		default:
//			status=-1;
//			fprintf(stderr," add_flowcovercut(): unexpected lp_solnstat: %d\n",lp_solnstat);
//			goto TERMINATE;
//		}
//		if(*lpcc_stat!=0) break;
//	}
//	t_end=clock();
//	fprintf(stdout," total flow cover cuts added: %d\n",cut_cnt);
//	fprintf(stdout," total time for generating cuts: %12.3f\n",(double)(t_end-t_start)/CLOCKS_PER_SEC);
//	fprintf(stdout," before adding flow cover cuts, lp relaxation is %f\n",init_lp_obj);
//	switch(*lpcc_stat)
//	{
//	case 1:
//		fprintf(stdout," after adding flow cover cuts, find optimal solution %f \n",lp_obj);
//		break;
//	case 0:
//		status=copy_startinfo(&updated_rx_lp_start,lp_start);
//		if(status)goto TERMINATE;
//		fprintf(stdout," after adding flow cover cuts, improved lp relaxation is %f\n",lp_obj);
//		break;
//	case -1:
//		fprintf(stdout," after adding flow cover cuts, find lpcc is infeasible \n");
//		break;
//	default:
//		status=-1;
//		fprintf(stderr," add_flowcovercut(): unexpected lpcc_stat: %d\n",*lpcc_stat);
//		goto TERMINATE;
//	}
//TERMINATE:
//	free_and_null((char**) &lp_soln);
//	free_startinfo(&rx_lp_start);
//	free_startinfo(&updated_rx_lp_start);
//	free_constraint_set(&cut_set);
//	return(status);
//}

int boundcuts_generator(CPXENVptr env, 
						CPXLPptr rx_lp, 	
						STARTINFO rx_lp_start,								 
						const int param_n,								 
						const int param_m, 								 
						const int param_k,	
						const double lpcc_ub,
						const int num_var,
						int *select_index,
						int param_p,
						CONSTRAINT_SET *cut_ptr)
{
	int status =0;
	int i;
	int index_row;
	int index_col;
	double* y_ub=NULL;
	double* w_ub=NULL;
	const double one=1;
	const double zero=0;
	double rhs;
	double pre_y_ub;
	double pre_w_ub;
	y_ub=(double*) malloc(param_m*sizeof(double));
	w_ub=(double*) malloc(param_m*sizeof(double));
	if (y_ub==NULL || w_ub==NULL)
	{
		status=-1;
		fprintf(stderr," compute_yw_UB(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<param_m;i++)
	{
		if (select_index[i]==0) continue;
		//first consider computing bound on y
		status=compute_UB(env,rx_lp,rx_lp_start,param_n,param_m,param_k,lpcc_ub,num_var,'y',i,param_p,&(y_ub[i]));
		if(status) goto TERMINATE;
		//fprintf(stdout," y_%d_ub: %f\n",i,y_ub[i]);
		//then consider computing bound on w
		status=compute_UB(env,rx_lp,rx_lp_start,param_n,param_m,param_k,lpcc_ub,num_var,'w',i,param_p,&(w_ub[i]));
		if(status) goto TERMINATE;
		//fprintf(stdout," w_%d_ub: %f\n",i,w_ub[i]);
	}
	//for (i=0;i<param_m;i++)
	//{
	//	if (select_index[i]==0) continue;
	//	//fprintf(stdout,"y_ub:%f,w_ub:%f\n",y_ub[i],w_ub[i]);
	//	if (1)
	//	{
	//		if (y_ub[i]<CPX_INFBOUND && y_ub[i]>ZERO_TOLERANCE)
	//		{
	//			status=ext_bound_cuts_generator(env,rx_lp,param_n,param_m,param_k,i,y_ub[i],2,cut_ptr);
	//			if(status)goto TERMINATE;
	//		}
	//		if (w_ub[i]<CPX_INFBOUND && w_ub[i]>ZERO_TOLERANCE)
	//		{
	//			status=ext_bound_cuts_generator(env,rx_lp,param_n,param_m,param_k,i,w_ub[i],1,cut_ptr);
	//			if(status)goto TERMINATE;
	//		}
	//	}
	//}
	for (i=0;i<param_m;i++)
	{
		if (select_index[i]==0) continue;
		//modify coefficient and only modify when the bound is improved
		index_row=param_k+param_m+i;
		index_col=param_n+i;
		status=CPXgetcoef(env,rx_lp,index_row,index_col,&pre_w_ub);
		if(status) goto TERMINATE;
		index_col=param_n+param_m+i;
		status=CPXgetcoef(env,rx_lp,index_row,index_col,&pre_y_ub);
		if(status) goto TERMINATE;
		//fprintf(stdout,"pre_w_ub: %f, pre_y_ub: %f\n",pre_w_ub,pre_y_ub);
		//fprintf(stdout,"cur_w_ub: %f, cur_y_ub: %f\n",w_ub[i],y_ub[i]);
		if (pre_w_ub>1e-6 || pre_y_ub>1e-6)
		{
			if (pre_y_ub<y_ub[i])
			{
				y_ub[i]=pre_y_ub;
			}
			if (pre_w_ub<w_ub[i])
			{
				w_ub[i]=pre_w_ub;
			}
		}
		//fprintf(stdout,"cur_w_ub: %f, cur_y_ub: %f\n",w_ub[i],y_ub[i]);
		if (y_ub[i]<1e-6 && w_ub[i]<1e-6)
		{
			y_ub[i]=0;
			w_ub[i]=1;
		}
		if (y_ub[i]<CPX_INFBOUND && w_ub[i]<CPX_INFBOUND)
		{		
			index_row=param_k+param_m+i;
			rhs=y_ub[i]*w_ub[i];
			status=CPXchgrhs(env,rx_lp,1,&index_row,&rhs);
			if(status)goto TERMINATE;
			index_col=param_n+i;
			status=CPXchgcoef(env,rx_lp,index_row,index_col,w_ub[i]);
			if(status) goto TERMINATE;
			index_col=param_n+param_m+i;
			status=CPXchgcoef(env,rx_lp,index_row,index_col,y_ub[i]);
			if(status) goto TERMINATE;
			index_col=param_n+param_m+param_m+param_k+i;
			status=CPXchgcoef(env,rx_lp,index_row,index_col,one);
			if(status) goto TERMINATE;
		}
	}
	//fprintf(stdout,"bounds on y: \n");
	//for (i=0;i<param_m;i++)
	//{
	//	fprintf(stdout,"y_%d:%f\n",i,y_ub[i]);
	//}
	//fprintf(stdout,"bounds on w: \n");
	//for (i=0;i<param_m;i++)
	//{
	//	fprintf(stdout,"w_%d:%f\n",i,w_ub[i]);
	//}
TERMINATE:
	free_and_null((char**)&y_ub);
	free_and_null((char**)&w_ub);
	return (status);
}
int add_boundcuts(CPXENVptr env, 
				  CPXLPptr rx_lp, 	
				  STARTINFO *lp_start,								 
				  const int param_n,								 
				  const int param_m, 								 
				  const int param_k,	
				  const int param_p,
				  const int loop_limit,
				  const double lpcc_ub,
				  double *opt_soln,
				  double *opt_obj,
				  int *lpcc_stat)
{
	//*lpcc_stat=1 means find optimal solution
	//*lpcc_stat=0 means haven't found optimal solution
	//*lpcc_stat=-1 mean lpcc is infeasible
	int status =0;
	int i,j;
	int lp_solnstat;
	int cnt;
	int cut_limit;
	int *sort_index=NULL;
	int *select_index=NULL;
//	int rx_lp_n_col;
//	int rx_lp_n_row;
	double *sort_val=NULL;
	double init_lp_obj;
	const int num_var= param_n+param_m+param_m+param_k;
	double *lp_soln=NULL;
	double lp_obj;
	STARTINFO rx_lp_start;
	STARTINFO updated_rx_lp_start;
	clock_t t_start,t_end;
	CONSTRAINT_SET cut_set;
	init_constraint_set(&cut_set);
	sort_index=(int*)malloc(param_m*sizeof(int));
	sort_val=(double*)malloc(param_m*sizeof(double));
	select_index=(int*)calloc(param_m,sizeof(int));
	if (sort_index==NULL || sort_val==NULL || select_index==NULL)
	{
		status=-1;
		fprintf(stderr,"add_boundcuts(): unable to allocate memory\n");
		goto TERMINATE;
	}

	*lpcc_stat=LPCC_PRECONDITION_BOUNDED;
	init_startinfo(&rx_lp_start);
	init_startinfo(&updated_rx_lp_start);
	status=copy_startinfo(lp_start,&rx_lp_start);
	if(status) goto TERMINATE;
	lp_soln=(double*) malloc(num_var*sizeof(double));
	if(lp_soln==NULL)
	{
		status = NO_MEMORY;	
		fprintf (stderr, " add_boundcuts(): could not allocate memory.\n");
		goto TERMINATE;
	}
	fprintf(stdout," Add bound cut at root node...(p=%d, loop=%d)\n",param_p,loop_limit);
	t_start=clock();
	for (i=0;i<loop_limit;i++)
	{
		status = cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat, &updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,2);
		if (status) goto TERMINATE;
		//fprintf(stdout,"debug; %d\n",lp_solnstat);
		switch(lp_solnstat)
		{
		case CPX_STAT_OPTIMAL:
			t_end=clock();
			if (i==0) init_lp_obj=lp_obj;
			else fprintf(stdout,"	LB: %12.4f PT: %6.0fs\n",lp_obj,(double)(t_end-t_start)/CLOCKS_PER_SEC);
			cnt=0;
			status=copy_startinfo(&updated_rx_lp_start,&rx_lp_start);
			if(status)goto TERMINATE;
			for (j=0;j<param_m;j++)
			{
				sort_index[j]=j;
				select_index[j]=0;
				sort_val[j]=lp_soln[param_n+j]*lp_soln[param_n+param_m+j];
				if (lp_soln[param_n+j]<ZERO_TOLERANCE || lp_soln[param_n+param_m+j]<ZERO_TOLERANCE)
				{
					cnt++;
				}
			}
			status=get_sort_index(param_m,sort_val,sort_index,-1);
			if (status) goto TERMINATE;
			cut_limit=5;
			cut_limit=cut_limit<(param_m-cnt)?cut_limit:(param_m-cnt);
			for (j=0;j<cut_limit;j++)
			{
				select_index[sort_index[j]]=1;
			}
			if (cnt==param_m)
			{
				//get the feasible lpcc solution 
				*lpcc_stat=LPCC_PRECONDITION_OPT;
				for(j=0;j<num_var;j++)
				{
					opt_soln[j]=lp_soln[j];
				}
				*opt_obj=lp_obj;
			}else
			{
				//add bound cut
				status=boundcuts_generator(env,rx_lp,rx_lp_start,param_n,param_m,param_k,lpcc_ub,num_var,select_index,param_p,&cut_set);
				if(status) goto TERMINATE;
				//add ext bound cuts to lp
				//if (cut_set.cons_matrix.n_row>0 && 0)
				//{	
				//	fprintf(stdout,"add ext bound cuts: %d\n",cut_set.cons_matrix.n_row);
				//	rx_lp_n_col=CPXgetnumcols(env,rx_lp);
				//	rx_lp_n_row=CPXgetnumrows(env,rx_lp);
				//	status = CPXaddrows(env,rx_lp,0,cut_set.cons_matrix.n_row,cut_set.cons_matrix.nnz,cut_set.rhs,NULL,cut_set.cons_matrix.matbeg,cut_set.cons_matrix.matind,cut_set.cons_matrix.matval,NULL,NULL);
				//	if (status) 
				//	{
				//		fprintf(stderr," add_boundcuts(): add row failed %d\n",status);
				//		goto TERMINATE;
				//	}
				//	status = CPXnewcols(env, rx_lp, cut_set.cons_matrix.n_row, NULL, NULL, NULL, NULL, NULL);
				//	if (status) goto TERMINATE;
				//	for (j=0; j<cut_set.cons_matrix.n_row;j++)
				//	{
				//		status = CPXchgcoef (env, rx_lp, rx_lp_n_row+j, rx_lp_n_col+j, -1.0);
				//		if (status) goto TERMINATE;
				//	}

				//	status = cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat, &updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,2);
				//	if (status) goto TERMINATE;
				//	fprintf(stdout,"col: %d, row: %d, obj: %f\n",CPXgetnumcols(env,rx_lp),CPXgetnumrows(env,rx_lp),lp_obj);
				//	status=clearcutpool(env,rx_lp,rx_lp_n_row,rx_lp_n_col,cut_set.cons_matrix.n_row);
				//	if (status) goto TERMINATE;
				//	free_startinfo(&rx_lp_start);
				//	init_startinfo(&rx_lp_start);
				//	status = cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat, &updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,2);
				//	if (status) goto TERMINATE;
				//	fprintf(stdout,"col: %d, row: %d, obj: %f\n",CPXgetnumcols(env,rx_lp),CPXgetnumrows(env,rx_lp),lp_obj);
				//}
				////clear cut_set
				//free_constraint_set(&cut_set);
				//init_constraint_set(&cut_set);
			}
			break;
		case CPX_STAT_INFEASIBLE:
			*lpcc_stat=LPCC_PRECONDITION_INFEASIBLE;
			break;
		default:
			status=-1;
			fprintf(stderr," add_boundcuts(): unexpected lp_solnstat: %d\n",lp_solnstat);
			goto TERMINATE;
		}
		if(*lpcc_stat!=LPCC_PRECONDITION_BOUNDED) break;
	}
	t_end=clock();
	fprintf(stdout,"	Total time for generating cuts: %12.3f\n",(double)(t_end-t_start)/CLOCKS_PER_SEC);
	switch(*lpcc_stat)
	{
	case LPCC_PRECONDITION_OPT:
		fprintf(stdout,"	After adding bound cuts, find optimal solution %f \n",lp_obj);
		break;
	case LPCC_PRECONDITION_BOUNDED:
		status=copy_startinfo(&updated_rx_lp_start,lp_start);
		if(status)goto TERMINATE;
		fprintf(stdout,"	After adding bound cuts, improved lp relaxation is %f\n",lp_obj);
		break;
	case LPCC_PRECONDITION_INFEASIBLE:
		fprintf(stdout,"	After adding bound cuts, find lpcc is infeasible \n");
		break;
	default:
		status=-1;
		fprintf(stderr," add_boundcuts(): unexpected lpcc_stat: %d\n",*lpcc_stat);
		goto TERMINATE;
	}
TERMINATE:
	free_and_null((char**)&sort_index);
	free_and_null((char**)&sort_val);
	free_and_null((char**)&select_index);
	free_and_null((char**) &lp_soln);
	free_startinfo(&rx_lp_start);
	free_startinfo(&updated_rx_lp_start);
	free_constraint_set(&cut_set);
	return(status);
}

int get_simplex_tableau_by_violation_index(CPXENVptr env, 
										   CPXLPptr rx_lp, 
										   const int param_n,
										   const int param_m,
										   const int param_k,
										   const int violate_complementary_index,
										   double **y_tableau,
										   double **w_tableau)
{
	int status=0;
	const int num_cols=CPXgetnumcols(env,rx_lp);
	const int num_rows=CPXgetnumrows(env,rx_lp);
	int row_id;
	*y_tableau=(double*)malloc(num_cols*sizeof(double));
	*w_tableau=(double*)malloc(num_cols*sizeof(double));

	if (*y_tableau==NULL || *w_tableau==NULL)
	{
		status=NO_MEMORY;
		fprintf(stderr," get_simplex_tableau():could not allocate memory.\n");
		goto TERMINATE;
	}
	status = CPXgetijrow(env,rx_lp,-1,param_n+violate_complementary_index,&row_id);
	if (status) goto TERMINATE;
	status = CPXbinvarow(env,rx_lp,row_id,*y_tableau);
	if (status) goto TERMINATE;
	status = CPXgetijrow(env,rx_lp,-1,param_n+param_m+violate_complementary_index,&row_id);
	if (status) goto TERMINATE;
	status = CPXbinvarow(env,rx_lp,row_id,*w_tableau);
	if(status) goto TERMINATE;
TERMINATE:
	return(status);
}
int get_simplex_tableau(CPXENVptr env, 
						CPXLPptr rx_lp, 								 
						const int param_n,								 
						const int param_m, 								 
						const int param_k,	
						const int violation_cnt,
						STARTINFO lp_start,
						double *lp_soln,
						double **y_hat,
						double **w_hat,
						double **y_bar,
						double **w_bar,
						double ***y_tableau,
						double ***w_tableau,
						int *violation_index,
						const int ub_setting)
{
	int status =0;
	int i;
	int cnt;
	int condition;
	STARTINFO rx_lp_start;
	STARTINFO updated_rx_lp_start;
	init_startinfo(&rx_lp_start);
	init_startinfo(&updated_rx_lp_start);
	status=copy_startinfo(&lp_start,&rx_lp_start);
	if(status) goto TERMINATE;
	*y_hat=(double*) calloc(violation_cnt,sizeof(double));
	*w_hat=(double*) calloc(violation_cnt,sizeof(double));
	*y_bar=(double*) calloc(violation_cnt,sizeof(double));
	*w_bar=(double*) calloc(violation_cnt,sizeof(double));
	*y_tableau=(double**) calloc(violation_cnt,sizeof(double*));
	*w_tableau=(double**) calloc(violation_cnt,sizeof(double*));
	if (*y_hat==NULL || *w_hat==NULL || *y_bar==NULL || *w_bar==NULL || *y_tableau==NULL || *w_tableau==NULL)
	{
		status=-1;
		fprintf(stderr," get_simplex_tableau(): could not allocate memory\n");
		goto TERMINATE;
	}

	cnt=0;
	for (i=0;i<param_m;i++)
	{
		if (lp_soln[param_n+i]>ZERO_TOLERANCE && lp_soln[param_n+param_m+i]>ZERO_TOLERANCE)
		{
			//violation complementarity 
			//get upper bound for y and w
			if (violation_index!=NULL)
			{
				violation_index[cnt]=i;
			}
			switch(ub_setting)
			{
			case -1:
				(*y_bar)[cnt]=CPX_INFBOUND;
				(*w_bar)[cnt]=CPX_INFBOUND;
				break;
			case 1:
				(*y_bar)[cnt]=1;
				(*w_bar)[cnt]=1;
				break;
			case 0:
				status=getVARbound(env,rx_lp,&rx_lp_start,&updated_rx_lp_start,param_n+i,CPX_MAX,(*y_bar)+cnt,&condition);
				if(status) goto TERMINATE;
				status=getVARbound(env,rx_lp,&rx_lp_start,&updated_rx_lp_start,param_n+param_m+i,CPX_MAX,(*w_bar)+cnt,&condition);
				if (status) goto TERMINATE;
				break;
			default:
				status=-1;
				fprintf(stderr," get_simplex_tableau(): unexpected ub_setting: %d\n",ub_setting);
				goto TERMINATE;
			}
			(*y_hat)[cnt]=lp_soln[param_n+i];
			(*w_hat)[cnt]=lp_soln[param_n+param_m+i];
			//fprintf(stdout,"violation: y_hat: %f, w_hat: %f\n",lp_soln[param_n+i],lp_soln[param_n+param_m+i]);
			status=get_simplex_tableau_by_violation_index(env,rx_lp,param_n,param_m,param_k,i,(*y_tableau)+cnt,(*w_tableau)+cnt);
			if(status) goto TERMINATE;
			cnt++;
		}
	}
TERMINATE:
	free_startinfo(&rx_lp_start);
	free_startinfo(&updated_rx_lp_start);
	return(status);
}

//int flowcovercut_generator_1(const int num_cols,
//							 const int num_rows,
//							 const int param_n,
//							 const int param_m,
//							 const int param_k,
//							 const int violation_cnt,
//							 const int violation_index,
//							 int *cstat,
//							 double *y_hat,
//							 double *w_hat,
//							 double *y_bar,
//							 double *w_bar,
//							 double **y_tableau,
//							 double **w_tableau,
//							 CONSTRAINT_SET *cut_ptr)
//{
//	int status=0;
//	int i;
//	int nnz;
//	double sum_hat;
//	double tmp_sum_hat_1;
//	double tmp_sum_hat_2;
//	double c_1;
//	double c_2;
//	double *sum_coef=NULL;
//	double *tmp_sum_coef_1=NULL;
//	double *tmp_sum_coef_2=NULL;
//	double *cut_coef=NULL;
//	double *tmp_cut_coef_1=NULL;
//	double *tmp_cut_coef_2=NULL;
//	int *ind=NULL;
//	double *val=NULL;
//	double norm;
//	double tmp_norm_1;
//	double tmp_norm_2;
//	double max_norm;
//	double scale=1;
//	int flag;
//	sum_coef=(double*) malloc(num_cols*sizeof(double));
//	tmp_sum_coef_1=(double*)malloc(num_cols*sizeof(double));
//	tmp_sum_coef_2=(double*)malloc(num_cols*sizeof(double));
//	cut_coef=(double*) malloc(num_cols*sizeof(double));
//	tmp_cut_coef_1=(double*) malloc(num_cols*sizeof(double));
//	tmp_cut_coef_2=(double*) malloc(num_cols*sizeof(double));
//	if (sum_coef==NULL || tmp_sum_coef_1==NULL || tmp_sum_coef_2==NULL || cut_coef == NULL || tmp_cut_coef_1==NULL || tmp_cut_coef_2==NULL)
//	{
//		status=-1;
//		fprintf(stderr," flowcovercut_generator():could not allocate memory.\n");
//		goto TERMINATE;
//	}
//	//fprintf(stdout,"%f,%f\n",y_hat[violation_index],w_hat[violation_index]);
//	if (y_bar[violation_index]<1e6)
//	{
//		c_1=y_bar[violation_index]*w_hat[violation_index]/(y_bar[violation_index]-y_hat[violation_index]);
//		c_2=w_bar[violation_index];
//		sum_hat=w_hat[violation_index];
//		vector_sum(NULL,w_tableau[violation_index],num_cols,sum_coef);
//		status=compute_coef_flowcovercut_1(num_cols,cstat,w_tableau[violation_index],y_tableau[violation_index],w_hat[violation_index],y_hat[violation_index],sum_hat,y_bar[violation_index],sum_coef,cut_coef,&norm);
//		if(status) goto TERMINATE;
//		//max_norm=-CPX_INFBOUND;
//		max_norm=norm;
//		//fprintf(stdout,"norm: %f,sum_hat:%f\n",norm,sum_hat);
//		flag=0;
//		for (i=0;i<violation_cnt;i++)
//		{
//			if (i==violation_index) continue;
//			tmp_sum_hat_1=sum_hat+y_hat[i];
//			tmp_sum_hat_2=sum_hat+w_hat[i];
//			vector_sum(sum_coef,y_tableau[i],num_cols,tmp_sum_coef_1);
//			vector_sum(sum_coef,w_tableau[i],num_cols,tmp_sum_coef_2);
//			tmp_norm_1=0;
//			tmp_norm_2=0;
//			//fprintf(stdout,"tmp_sum_hat_1: %f, tmp_sum_hat_2: %f,c1: %f, c2: %f\n",tmp_sum_hat_1,tmp_sum_hat_2,c_1,c_2);
//			if (tmp_sum_hat_1<c_1-1e-6 && tmp_sum_hat_1<=c_2)
//			{
//				status=compute_coef_flowcovercut_1(num_cols,cstat,w_tableau[violation_index],y_tableau[violation_index],w_hat[violation_index],y_hat[violation_index],tmp_sum_hat_1,y_bar[violation_index],tmp_sum_coef_1,tmp_cut_coef_1,&tmp_norm_1);
//				if(status) goto TERMINATE;
//				//fprintf(stdout,"tmp_norm_1: %f\n",tmp_norm_1);
//			}
//			if (tmp_sum_hat_2<c_1-1e-6 && tmp_sum_hat_2<=c_2)
//			{
//				status=compute_coef_flowcovercut_1(num_cols,cstat,w_tableau[violation_index],y_tableau[violation_index],w_hat[violation_index],y_hat[violation_index],tmp_sum_hat_2,y_bar[violation_index],tmp_sum_coef_2,tmp_cut_coef_2,&tmp_norm_2);
//				if(status) goto TERMINATE;
//				//fprintf(stdout,"tmp_norm_2: %f\n",tmp_norm_2);
//			}
//			if (tmp_norm_1!=0 || tmp_norm_2!=0)
//			{
//				if (tmp_norm_1<=tmp_norm_2 && tmp_norm_2>max_norm*scale)
//				{
//					max_norm=tmp_norm_2;
//					sum_hat=tmp_sum_hat_2;
//					vector_sum(NULL,tmp_sum_coef_2,num_cols,sum_coef);
//					vector_sum(NULL,tmp_cut_coef_2,num_cols,cut_coef);
//					//fprintf(stdout,"improve\n");
//					flag=1;
//				}
//				if (tmp_norm_1>tmp_norm_2 && tmp_norm_1>max_norm*scale)
//				{
//					max_norm=tmp_norm_1;
//					sum_hat=tmp_sum_hat_1;
//					vector_sum(NULL,tmp_sum_coef_1,num_cols,sum_coef);
//					vector_sum(NULL,tmp_cut_coef_1,num_cols,cut_coef);
//					//fprintf(stdout,"improve\n");
//					flag=1;
//				}
//			}
//			
//		}
//		if (flag)
//		{
//			free_and_null((char**)&ind);
//			free_and_null((char**)&val);
//			status=convertArray(num_cols,cut_coef,&nnz,&ind,&val);
//			if(status)goto TERMINATE;
//			status=add_cuts(cut_ptr,num_cols,nnz,ind,val,1);
//			if(status)goto TERMINATE;
//		}
//	}
//	if (w_bar[violation_index]<1e6)
//	{
//		c_1=y_hat[violation_index]*w_bar[violation_index]/(w_bar[violation_index]-w_hat[violation_index]);
//		c_2=y_bar[violation_index];
//		sum_hat=y_hat[violation_index];
//		vector_sum(NULL,y_tableau[violation_index],num_cols,sum_coef);
//		status=compute_coef_flowcovercut_1(num_cols,cstat,y_tableau[violation_index],w_tableau[violation_index],y_hat[violation_index],w_hat[violation_index],sum_hat,w_bar[violation_index],sum_coef,cut_coef,&norm);
//		if(status) goto TERMINATE;
//		//max_norm=-CPX_INFBOUND;
//		max_norm=norm;
//		//fprintf(stdout,"norm: %f,sum_hat:%f\n",norm,sum_hat);
//		flag=0;
//		for (i=0;i<violation_cnt;i++)
//		{
//			if (i==violation_index) continue;
//			tmp_sum_hat_1=sum_hat+y_hat[i];
//			tmp_sum_hat_2=sum_hat+w_hat[i];
//			vector_sum(sum_coef,y_tableau[i],num_cols,tmp_sum_coef_1);
//			vector_sum(sum_coef,w_tableau[i],num_cols,tmp_sum_coef_2);
//			tmp_norm_1=0;
//			tmp_norm_2=0;
//			//fprintf(stdout,"tmp_sum_hat_1: %f, tmp_sum_hat_2: %f,c1: %f, c2: %f\n",tmp_sum_hat_1,tmp_sum_hat_2,c_1,c_2);
//			if (tmp_sum_hat_1<c_1-1e-6 && tmp_sum_hat_1<=c_2)
//			{
//				status=compute_coef_flowcovercut_1(num_cols,cstat,y_tableau[violation_index],w_tableau[violation_index],y_hat[violation_index],w_hat[violation_index],tmp_sum_hat_1,w_bar[violation_index],tmp_sum_coef_1,tmp_cut_coef_1,&tmp_norm_1);
//				if(status) goto TERMINATE;
//				//fprintf(stdout,"tmp_norm_1: %f\n",tmp_norm_1);
//			}
//			if (tmp_sum_hat_2<c_1-1e-6 && tmp_sum_hat_2<=c_2)
//			{
//				status=compute_coef_flowcovercut_1(num_cols,cstat,y_tableau[violation_index],w_tableau[violation_index],y_hat[violation_index],w_hat[violation_index],tmp_sum_hat_2,w_bar[violation_index],tmp_sum_coef_2,tmp_cut_coef_2,&tmp_norm_2);
//				if(status) goto TERMINATE;
//				//fprintf(stdout,"tmp_norm_2: %f\n",tmp_norm_2);
//			}
//			if (tmp_norm_1!=0 || tmp_norm_2!=0)
//			{
//				if (tmp_norm_1<=tmp_norm_2 && tmp_norm_2>max_norm*scale)
//				{
//					max_norm=tmp_norm_2;
//					sum_hat=tmp_sum_hat_2;
//					vector_sum(NULL,tmp_sum_coef_2,num_cols,sum_coef);
//					vector_sum(NULL,tmp_cut_coef_2,num_cols,cut_coef);
//					//fprintf(stdout,"improve\n");
//					flag=1;
//				}
//				if (tmp_norm_1>tmp_norm_2 && tmp_norm_1>max_norm*scale)
//				{
//					max_norm=tmp_norm_1;
//					sum_hat=tmp_sum_hat_1;
//					vector_sum(NULL,tmp_sum_coef_1,num_cols,sum_coef);
//					vector_sum(NULL,tmp_cut_coef_1,num_cols,cut_coef);
//					//fprintf(stdout,"improve\n");
//					flag=1;
//				}
//			}
//		}
//		if (flag)
//		{
//			free_and_null((char**)&ind);
//			free_and_null((char**)&val);
//			status=convertArray(num_cols,cut_coef,&nnz,&ind,&val);
//			if(status)goto TERMINATE;
//			status=add_cuts(cut_ptr,num_cols,nnz,ind,val,1);
//			if(status)goto TERMINATE;
//		}
//	}
//TERMINATE:
//	free_and_null((char**)&sum_coef);
//	free_and_null((char**)&tmp_sum_coef_1);
//	free_and_null((char**)&tmp_sum_coef_2);
//	free_and_null((char**)&cut_coef);
//	free_and_null((char**)&tmp_cut_coef_1);
//	free_and_null((char**)&tmp_cut_coef_2);
//	free_and_null((char**)&ind);
//	free_and_null((char**)&val);
//	return(status);
//}
//int flowcovercut_generator_2(const int num_cols,
//							 const int num_rows,
//							 const int param_n,
//							 const int param_m,
//							 const int param_k,
//							 const int violation_cnt,
//							 const int violation_index,
//							 int *cstat,
//							 double *y_hat,
//							 double *w_hat,
//							 double *y_bar,
//							 double *w_bar,
//							 double **y_tableau,
//							 double **w_tableau,
//							 CONSTRAINT_SET *cut_ptr)
//{
//	int status=0;
//	int i;
//	int nnz;
//	double *cut_coef=NULL;
//	double *tmp_cut_coef_1=NULL;
//	double *tmp_cut_coef_2=NULL;
//	int *ind=NULL;
//	int *V_select_1=NULL;
//	int *V_select_2=NULL;
//	double *val=NULL;
//	double norm;
//	double tmp_norm_1;
//	double tmp_norm_2;
//	double max_norm;
//	int flag;
//
//	cut_coef=(double*) malloc(num_cols*sizeof(double));
//	tmp_cut_coef_1=(double*) malloc(num_cols*sizeof(double));
//	tmp_cut_coef_2=(double*) malloc(num_cols*sizeof(double));
//	V_select_1=(int*) calloc(violation_cnt,sizeof(int));
//	V_select_2=(int*) calloc(violation_cnt,sizeof(int));
//	if (V_select_1==NULL || V_select_2==NULL || cut_coef == NULL || tmp_cut_coef_1==NULL || tmp_cut_coef_2==NULL)
//	{
//		status=-1;
//		fprintf(stderr," flowcovercut_generator():could not allocate memory.\n");
//		goto TERMINATE;
//	}
//	//fprintf(stdout,"violation: %f,%f\n",y_hat[violation_index],w_hat[violation_index]);
//	if (y_bar[violation_index]<1e6)
//	{
//		V_select_1[violation_index]=2;
//		status=compute_coef_flowcovercut_2(num_cols,cstat,violation_cnt,y_tableau,w_tableau,y_hat,w_hat,y_bar,w_bar,V_select_1,cut_coef,&norm);
//		if(status) goto TERMINATE;
//		max_norm=norm;
//		//fprintf(stdout,"norm: %f\n",norm);
//		flag=0;
//		for (i=0;i<violation_cnt;i++)
//		{
//			if (i==violation_index) continue;
//			tmp_norm_1=0;
//			tmp_norm_2=0;
//			//fprintf(stdout,"tmp_sum_hat_1: %f, tmp_sum_hat_2: %f,c1: %f, c2: %f\n",tmp_sum_hat_1,tmp_sum_hat_2,c_1,c_2);
//			V_select_1[i]=1;
//			status=compute_coef_flowcovercut_2(num_cols,cstat,violation_cnt,y_tableau,w_tableau,y_hat,w_hat,y_bar,w_bar,V_select_1,tmp_cut_coef_1,&tmp_norm_1);
//			if(status)goto TERMINATE;
//			V_select_1[i]=2;
//			status=compute_coef_flowcovercut_2(num_cols,cstat,violation_cnt,y_tableau,w_tableau,y_hat,w_hat,y_bar,w_bar,V_select_1,tmp_cut_coef_2,&tmp_norm_2);
//			if(status)goto TERMINATE;
//			V_select_1[i]=0;
//			if (tmp_norm_1!=0 || tmp_norm_2!=0)
//			{
//				if (tmp_norm_1<=tmp_norm_2 && tmp_norm_2>max_norm)
//				{
//					max_norm=tmp_norm_2;
//					V_select_1[i]=2;
//					vector_sum(NULL,tmp_cut_coef_2,num_cols,cut_coef);
//					fprintf(stdout,"improve: %f\n",max_norm);
//					flag=1;
//				}
//				if (tmp_norm_1>tmp_norm_2 && tmp_norm_1>max_norm)
//				{
//					max_norm=tmp_norm_1;
//					V_select_1[i]=1;
//					vector_sum(NULL,tmp_cut_coef_1,num_cols,cut_coef);
//					fprintf(stdout,"improve: %f\n",max_norm);
//					flag=1;
//				}
//			}
//
//		}
//		if (flag)
//		{
//			free_and_null((char**)&ind);
//			free_and_null((char**)&val);
//			status=convertArray(num_cols,cut_coef,&nnz,&ind,&val);
//			if(status)goto TERMINATE;
//			status=add_cuts(cut_ptr,num_cols,nnz,ind,val,1);
//			if(status)goto TERMINATE;
//		}
//	}
//	if (w_bar[violation_index]<1e6)
//	{
//		V_select_2[violation_index]=1;
//		status=compute_coef_flowcovercut_2(num_cols,cstat,violation_cnt,y_tableau,w_tableau,y_hat,w_hat,y_bar,w_bar,V_select_2,cut_coef,&norm);
//		if(status) goto TERMINATE;
//		max_norm=norm;
//		//fprintf(stdout,"norm: %f\n",norm);
//		flag=0;
//		for (i=0;i<violation_cnt;i++)
//		{
//			if (i==violation_index) continue;
//			tmp_norm_1=0;
//			tmp_norm_2=0;
//			//fprintf(stdout,"tmp_sum_hat_1: %f, tmp_sum_hat_2: %f,c1: %f, c2: %f\n",tmp_sum_hat_1,tmp_sum_hat_2,c_1,c_2);
//			V_select_2[i]=1;
//			status=compute_coef_flowcovercut_2(num_cols,cstat,violation_cnt,y_tableau,w_tableau,y_hat,w_hat,y_bar,w_bar,V_select_2,tmp_cut_coef_1,&tmp_norm_1);
//			if(status)goto TERMINATE;
//			V_select_2[i]=2;
//			status=compute_coef_flowcovercut_2(num_cols,cstat,violation_cnt,y_tableau,w_tableau,y_hat,w_hat,y_bar,w_bar,V_select_2,tmp_cut_coef_2,&tmp_norm_2);
//			if(status)goto TERMINATE;
//			V_select_2[i]=0;
//			if (tmp_norm_1!=0 || tmp_norm_2!=0)
//			{
//				if (tmp_norm_1<=tmp_norm_2 && tmp_norm_2>max_norm)
//				{
//					max_norm=tmp_norm_2;
//					V_select_2[i]=2;
//					vector_sum(NULL,tmp_cut_coef_2,num_cols,cut_coef);
//					fprintf(stdout,"improve: %f\n",max_norm);
//					flag=1;
//				}
//				if (tmp_norm_1>tmp_norm_2 && tmp_norm_1>max_norm)
//				{
//					max_norm=tmp_norm_1;
//					V_select_2[i]=1;
//					vector_sum(NULL,tmp_cut_coef_1,num_cols,cut_coef);
//					fprintf(stdout,"improve: %f\n",max_norm);
//					flag=1;
//				}
//			}
//		}
//		if (flag)
//		{
//			free_and_null((char**)&ind);
//			free_and_null((char**)&val);
//			status=convertArray(num_cols,cut_coef,&nnz,&ind,&val);
//			if(status)goto TERMINATE;
//			status=add_cuts(cut_ptr,num_cols,nnz,ind,val,1);
//			if(status)goto TERMINATE;
//		}
//	}
//TERMINATE:
//	free_and_null((char**)&V_select_1);
//	free_and_null((char**)&V_select_2);
//	free_and_null((char**)&cut_coef);
//	free_and_null((char**)&tmp_cut_coef_1);
//	free_and_null((char**)&tmp_cut_coef_2);
//	free_and_null((char**)&ind);
//	free_and_null((char**)&val);
//	return(status);
//}

//int compute_coef_flowcovercut_1(int num_cols,
//								int *cstat,
//								double *tableau_1,
//								double *tableau_2,
//								double hat_1,
//								double hat_2,
//								double sum_hat,
//								double ub,
//								double *sum_V_tableau,
//								double *cut_coef,
//								double *norm)
//{
//	int status=0;
//	int i;
//	double rhs;
//	double tmp;
//	double tmp_norm=0;
//	rhs=hat_1+sum_hat*hat_2/ub-sum_hat;
//	//fprintf(stdout,"%f,%f,%f,%f\n",hat_1,hat_2,sum_hat,ub);
//	if(rhs<EP)
//	{
//		status=-1;
//		fprintf(stderr," compute_coef_flowcovercut(): rhs is too small: %f\n",rhs);
//		goto TERMINATE;
//	}
//	for (i=0;i< num_cols;i++)
//	{
//		if (cstat[i]==CPX_AT_LOWER)
//		{
//			tmp=sum_V_tableau[i]<-EP?-sum_V_tableau[i]:0;
//			cut_coef[i]=(tableau_1[i]+tableau_2[i]*sum_hat/ub+tmp)/rhs;
//			//tmp_norm+=pow(tmp,2);
//			if (cut_coef[i]>EP)
//			{
//				tmp_norm+=pow(cut_coef[i],2);
//				
//			}
//		}else if(cstat[i]==CPX_BASIC)
//		{
//			cut_coef[i]=0;
//		}else
//		{
//			status=-1;
//			fprintf(stderr," compute_coef_flowcovercut(): unexpected base status: %d\n",cstat[i]);
//			goto TERMINATE;
//		}
//	}
//	tmp_norm=pow(tmp_norm,0.5);
//	tmp_norm=1/tmp_norm;
//	*norm=tmp_norm;
//TERMINATE:
//	return(status);
//}
//
//int compute_coef_flowcovercut_2(int num_cols,
//								int *cstat,
//								const int violation_cnt,
//								double **y_tableau,
//								double **w_tableau,
//								double *y_hat,
//								double *w_hat,
//								double *y_bar,
//								double *w_bar,
//								int *V_select,
//								double *cut_coef,
//								double *norm)
//{
//	int status=0;
//	int i,j;
//	double rhs;
//	double tmp;
//	double tmp_sum_1;
//	double tmp_sum_2;
//	double tmp_norm=0;
//	double lambda;
//	lambda=0;
//	for (i=0;i<violation_cnt;i++)
//	{
//		switch(V_select[i])
//		{
//		case 0:
//			break;
//		case 1:
//			lambda+=y_bar[i]-y_hat[i];
//			//fprintf(stdout,"ybar: %f, yhat: %f,lembda: %f\n",y_bar[i],y_hat[i],lambda);
//			break;
//		case 2:
//			lambda+=w_bar[i]-w_hat[i];
//			//fprintf(stdout,"wbar: %f, what: %f,lembda: %f\n",w_bar[i],w_hat[i],lambda);
//			break;
//		default:
//			status=-1;
//			fprintf(stderr," compute_coef_flowcovercut_2(): unexpected V_select:%d\n",V_select[i]);
//			goto TERMINATE;
//		}
//	}
//	rhs=0;
//	for (i=0;i<violation_cnt;i++)
//	{
//		switch(V_select[i])
//		{
//		case 0:
//			break;
//		case 1:
//			tmp=(y_bar[i]-lambda)>0?(y_bar[i]-lambda):0;
//			//fprintf(stdout,"%f,yhat:%f\n",y_bar[i]-lambda,y_hat[i]);
//			rhs+=tmp*w_hat[i]/w_bar[i];
//			break;
//		case 2:
//			tmp=(w_bar[i]-lambda)>0?(w_bar[i]-lambda):0;
//			//fprintf(stdout,"%f,what:%f\n",w_bar[i]-lambda,w_hat[i]);
//			rhs+=tmp*y_hat[i]/y_bar[i];
//			break;
//		default:
//			status=-1;
//			fprintf(stderr," compute_coef_flowcovercut_2(): unexpected V_select:%d\n",V_select[i]);
//			goto TERMINATE;
//		}
//	}
//	//fprintf(stdout,"lembda: %f\n",lambda);
//	if(rhs<1e-6)
//	{
//		*norm=0;
//		goto TERMINATE;
//	}
//
//	for (i=0;i< num_cols;i++)
//	{
//		if (cstat[i]==CPX_AT_LOWER)
//		{
//			tmp_sum_1=0;
//			tmp_sum_2=0;
//			for (j=0;j<violation_cnt;j++)
//			{
//				switch(V_select[j])
//				{
//				case 0:
//					break;
//				case 1:
//					tmp=(y_bar[j]-lambda)>0?(y_bar[j]-lambda):0;
//					tmp_sum_1+=y_tableau[j][i];
//					tmp_sum_2+=tmp*w_tableau[j][i]/w_bar[j];
//					break;
//				case 2:
//					tmp=(w_bar[j]-lambda)>0?(w_bar[j]-lambda):0;
//					tmp_sum_1+=w_tableau[j][i];
//					tmp_sum_2+=tmp*y_tableau[j][i]/y_bar[j];
//					break;
//				default:
//					status=-1;
//					fprintf(stderr," compute_coef_flowcovercut_2(): unexpected V_select:%d\n",V_select[i]);
//					goto TERMINATE;
//				}
//			}
//			tmp_sum_1=tmp_sum_1>0?tmp_sum_1:0;
//			cut_coef[i]=(tmp_sum_1+tmp_sum_2)/rhs;
//			//tmp_norm+=pow(tmp_sum_1,2);
//			if (cut_coef[i]>EP)
//			{
//				tmp_norm+=pow(cut_coef[i],2);
//
//			}
//		}else if(cstat[i]==CPX_BASIC)
//		{
//			cut_coef[i]=0;
//		}else
//		{
//			status=-1;
//			fprintf(stderr," compute_coef_flowcovercut(): unexpected base status: %d\n",cstat[i]);
//			goto TERMINATE;
//		}
//	}
//	tmp_norm=pow(tmp_norm,0.5);
//	tmp_norm=tmp_norm>1e-8?tmp_norm:1e-8;
//	tmp_norm=1/tmp_norm;
//	*norm=tmp_norm;
//TERMINATE:
//	return(status);
//}

int compute_coef_ext_simplecut(int num_cols,
								int *cstat,
								const int violation_cnt,
								double **y_tableau,
								double **w_tableau,
								double *y_hat,
								double *w_hat,
								double *y_bar,
								double *w_bar,
								int **V_select,
								const int candidate_size,
								const int V_size,
								double *cut_coef,
								double *norm)
{
	int status=0;
	int i,j,k;
	int cnt;
	double *coef=NULL;
	double tmp_norm;
	coef=(double*)calloc(num_cols,sizeof(double));
	if (coef==NULL)
	{
		status=-1;
		fprintf(stderr," compute_coef_ext_simplecut(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<num_cols;i++)
	{
		switch(cstat[i])
		{
		case CPX_BASIC:
			cut_coef[i]=0;
			break;
		case CPX_AT_LOWER:
			cut_coef[i]=-CPX_INFBOUND;
			break;
		default:
			status=-1;
			fprintf(stderr," compute_coef_ext_simplecut(): unexpected cstat: %d\n",cstat[i]);
			goto TERMINATE;
		}
	}
	fprintf(stdout,"candidate: %d\n",candidate_size);
	for (i=0;i<candidate_size;i++)
	{
		tmp_norm=0;
		for (j=0;j<num_cols;j++)
		{
			switch(cstat[j])
			{
			case CPX_BASIC:
				break;
			case CPX_AT_LOWER:
				coef[j]=0;
				cnt=0;
				for (k=0;k<violation_cnt;k++)
				{
					switch(V_select[i][k])
					{
					case 0:
						break;
					case 1:
						coef[j]+=y_tableau[k][j]/y_hat[k];
						cnt++;
						break;
					case 2:
						coef[j]+=w_tableau[k][j]/w_hat[k];
						cnt++;
						break;
					default:
						status=-1;
						fprintf(stderr," compute_coef_ext_simplecut(): unexpected V_select: %d\n",V_select[i][k]);
						goto TERMINATE;
					}
				}
				if (cnt!=V_size)
				{
					status=-1;
					fprintf(stderr," compute_coef_ext_simplecut(): unexpected condition:%d,%d\n",cnt,V_size);
					goto TERMINATE;
				}
				coef[j]=coef[j]/cnt;
				if (cut_coef[j]<coef[j])
				{
					cut_coef[j]=coef[j];
				}
				break;
			default:
				status=-1;
				fprintf(stderr," compute_coef_ext_simplecut(): unexpected cstat: %d\n",cstat[j]);
				goto TERMINATE;
			}
			tmp_norm+=pow(max(cut_coef[j],0),2);
			//fprintf(stdout,"%f,",coef[j]);
		}
		//fprintf(stdout,"\n");
	}

	tmp_norm=pow(tmp_norm,0.5);
	fprintf(stdout,"norm:%f\n",tmp_norm);
	*norm=tmp_norm;
TERMINATE:
	free_and_null((char**)&coef);
	return(status);
}

int generate_candidate_V_set(int num_cols,
							 const int violation_cnt,
							 int *cstat,
							 double **y_tableau,
							 double **w_tableau,
							 double *y_hat,
							 double *w_hat,
							 double *y_bar,
							 double *w_bar,
							 const int V_size,
							 int ***V_set,
							 int *cand_cnt)
{
	int status=0;
	int i,j;
	int tmp_n;
	const int candidate_cnt=(int) pow(2,V_size);
	int *sort_index=NULL;
	double *prod_violation=NULL;
	double norm_1,norm_2;
	if (violation_cnt<V_size)
	{
		status=-1;
		fprintf(stderr," generate_candidate_V_set(): violation_cnt is smaller than V_size: violation_cnt: %d, V_size: %d\n",violation_cnt,V_size);
		goto TERMINATE;
	}
	*V_set=(int**)malloc(candidate_cnt*sizeof(int*));
	sort_index=(int*)malloc(violation_cnt*sizeof(int));
	prod_violation=(double*)malloc(violation_cnt*sizeof(double));
	if ((*V_set)==NULL || sort_index==NULL || prod_violation==NULL)
	{
		status=-1;
		fprintf(stderr," generate_candidate_V_set(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<candidate_cnt;i++)
	{
		(*V_set)[i]=(int*)calloc(violation_cnt,sizeof(int));
		if ((*V_set)[i]==NULL)
		{
			status=-1;
			fprintf(stderr," generate_candidate_V_set(): unable to allocate memory\n");
			goto TERMINATE;
		}
	}
	for (i=0;i<violation_cnt;i++)
	{
		sort_index[i]=i;
		norm_1=0;
		norm_2=0;
		for (j=0;j<num_cols;j++)
		{
			if (cstat[j]==CPX_BASIC || cstat[j]==CPX_FREE_SUPER || cstat[j]==CPX_AT_UPPER) continue;
			if(y_tableau[i][j]>ZERO_TOLERANCE)norm_1+=pow(y_tableau[i][j],2);
			if(w_tableau[i][j]>ZERO_TOLERANCE)norm_2+=pow(w_tableau[i][j],2);
		}
		//prod_violation[i]=pow(y_hat[i]*w_hat[i],1)/(pow(norm_1,0.5)*pow(norm_2,0.5));
		prod_violation[i]=pow(y_hat[i]*w_hat[i],1);
	}
	status=get_sort_index(violation_cnt,prod_violation,sort_index,-1);
	if(status) goto TERMINATE;
	for (i=0;i<candidate_cnt;i++)
	{
		tmp_n=i;
		//fprintf(stdout,"\ni:%d\n",i);
		for (j=0;j<V_size;j++)
		{

			(*V_set)[i][sort_index[j]]=1+(tmp_n%2);
			//fprintf(stdout,"%d",(*V_set)[i][sort_index[j]]-1);
			tmp_n/=2;
		}
	}
	*cand_cnt=candidate_cnt;
TERMINATE:
	free_and_null((char**)&sort_index);
	free_and_null((char**)&prod_violation);
	return(status);
}

//int flowcovercut_generator(CPXENVptr env, 
//						   CPXLPptr rx_lp, 
//						   STARTINFO lp_start,
//						   const int param_n,
//						   const int param_m,
//						   const int param_k,
//						   double *lp_soln,
//						   CONSTRAINT_SET *cut_ptr)
//{
//	int status=0;
//	int i;
//	int violation_cnt;
//	int *cstat=NULL;
//	int *rstat=NULL;
//	const int num_cols=CPXgetnumcols(env,rx_lp);
//	const int num_rows=CPXgetnumrows(env,rx_lp);
//	double *y_hat=NULL;
//	double *w_hat=NULL;
//	double *y_bar=NULL;
//	double *w_bar=NULL;
//	double **y_tableau=NULL;
//	double **w_tableau=NULL;
//	int **V_set=NULL;
//	cstat=(int*) malloc(num_cols*sizeof(int));
//	rstat=(int*) malloc(num_rows*sizeof(int));
//	if (cstat==NULL || rstat==NULL)
//	{
//		status=-1;
//		fprintf(stderr," flowcovercut_generator(): can not allocate memory\n");
//		goto TERMINATE;
//	}
//	//count the violation complemnetarity
//	violation_cnt=0;
//	for (i=0;i<param_m;i++)
//	{
//		if (lp_soln[param_n+i]>EP && lp_soln[param_n+param_m+i]>EP)
//		{
//			violation_cnt++;
//		}
//	}
//	if (violation_cnt==0) goto TERMINATE;
//	//fprintf(stdout,"call flow cover generator, violation cnt: %d\n",violation_cnt);
//	//get simplex information
//	//fprintf(stdout,"get simplex tableau information\n");
//	status=get_simplex_tableau(env,rx_lp,param_n,param_m,param_k,violation_cnt,lp_start,lp_soln,&y_hat,&w_hat,&y_bar,&w_bar,&y_tableau,&w_tableau,NULL,0);
//	if (status) goto TERMINATE;
//	//get base information
//	status=CPXgetbase(env,rx_lp,cstat,rstat);
//	if(status) goto TERMINATE;
//	//fprintf(stdout,"generate flow cut\n");
//	for (i=0;i<violation_cnt;i++)
//	{
//		//fprintf(stdout,"%d\n",i);
//		//fprintf(stdout,"check %f,%f\n",y_hat[i],w_hat[i]);
//		//status=flowcovercut_generator_1(num_cols,num_rows,param_n,param_m,param_k,violation_cnt,i,cstat,y_hat,w_hat,y_bar,w_bar,y_tableau,w_tableau,cut_ptr);
//		//if(status) goto TERMINATE;
//		status=flowcovercut_generator_2(num_cols,num_rows,param_n,param_m,param_k,violation_cnt,i,cstat,y_hat,w_hat,y_bar,w_bar,y_tableau,w_tableau,cut_ptr);
//		if(status) goto TERMINATE;
//	}
//TERMINATE:
//	free_and_null((char**)&cstat);
//	free_and_null((char**)&rstat);
//	free_and_null((char**)&y_hat);
//	free_and_null((char**)&w_hat);
//	free_and_null((char**)&y_bar);
//	free_and_null((char**)&w_bar);
//	for (i=0;i<violation_cnt;i++)
//	{
//		free_and_null((char**)&(y_tableau[i]));
//		free_and_null((char**)&(w_tableau[i]));
//	}
//	free_and_null((char**)&y_tableau);
//	free_and_null((char**)&w_tableau);
//	return(status);
//}
int ext_simple_cut_generator(CPXENVptr env, 
						   CPXLPptr rx_lp, 
						   STARTINFO lp_start,
						   const int param_n,
						   const int param_m,
						   const int param_k,
						   double *lp_soln,
						   const int init_v_size,
						   CONSTRAINT_SET *cut_ptr)
{
	int status=0;
	int i;
	int violation_cnt;
	int V_size;
	int *cstat=NULL;
	int *rstat=NULL;
	int candidate_cnt;
	int *ind=NULL;
	int nnz;
	const int num_cols=CPXgetnumcols(env,rx_lp);
	const int num_rows=CPXgetnumrows(env,rx_lp);
	double *y_hat=NULL;
	double *w_hat=NULL;
	double *y_bar=NULL;
	double *w_bar=NULL;
	double **y_tableau=NULL;
	double **w_tableau=NULL;
	double *val=NULL;
	double *cut_coef=NULL;
	double norm;
	int **V_set=NULL;
	cstat=(int*) malloc(num_cols*sizeof(int));
	rstat=(int*) malloc(num_rows*sizeof(int));
	cut_coef=(double*) malloc(num_cols*sizeof(double));
	if (cstat==NULL || rstat==NULL || cut_coef==NULL)
	{
		status=-1;
		fprintf(stderr," flowcovercut_generator(): can not allocate memory\n");
		goto TERMINATE;
	}
	//count the violation complemnetarity
	violation_cnt=0;
	for (i=0;i<param_m;i++)
	{
		if (lp_soln[param_n+i]>ZERO_TOLERANCE && lp_soln[param_n+param_m+i]>ZERO_TOLERANCE)
		{
			violation_cnt++;
		}
	}
	if (violation_cnt==0) goto TERMINATE;
	//fprintf(stdout,"call flow cover generator, violation cnt: %d\n",violation_cnt);
	//get simplex information
	//fprintf(stdout,"get simplex tableau information\n");
	status=get_simplex_tableau(env,rx_lp,param_n,param_m,param_k,violation_cnt,lp_start,lp_soln,&y_hat,&w_hat,&y_bar,&w_bar,&y_tableau,&w_tableau,NULL,-1);
	if (status) goto TERMINATE;
	//get base information
	status=CPXgetbase(env,rx_lp,cstat,rstat);
	if(status) goto TERMINATE;
	V_size=min(violation_cnt,init_v_size);
	status=generate_candidate_V_set(num_cols,violation_cnt,cstat,y_tableau,w_tableau,y_hat,w_hat,y_bar,w_bar,V_size,&V_set,&candidate_cnt);
	//fprintf(stdout,"generate flow cut\n");
	status=compute_coef_ext_simplecut(num_cols,cstat,violation_cnt,y_tableau,w_tableau,y_hat,w_hat,y_bar,w_bar,V_set,candidate_cnt,V_size,cut_coef,&norm);
	if(status) goto TERMINATE;
	free_and_null((char**)&ind);
	free_and_null((char**)&val);
	status=convertArray(num_cols,cut_coef,&nnz,&ind,&val);
	if(status)goto TERMINATE;
	status=add_cuts(cut_ptr,num_cols,nnz,ind,val,1);
	if(status)goto TERMINATE;
TERMINATE:
	free_and_null((char**)&ind);
	free_and_null((char**)&val);
	free_and_null((char**)&cut_coef);
	for (i=0;i<candidate_cnt;i++)
	{
		free_and_null((char**)&(V_set[i]));
	}
	free_and_null((char**)&V_set);
	free_and_null((char**)&cstat);
	free_and_null((char**)&rstat);
	free_and_null((char**)&y_hat);
	free_and_null((char**)&w_hat);
	free_and_null((char**)&y_bar);
	free_and_null((char**)&w_bar);
	for (i=0;i<violation_cnt;i++)
	{
		free_and_null((char**)&(y_tableau[i]));
		free_and_null((char**)&(w_tableau[i]));
	}
	free_and_null((char**)&y_tableau);
	free_and_null((char**)&w_tableau);
	return(status);
}

int add_mccormickcut(CPXENVptr env,
					CPXLPptr lp, 
					const int param_n,
					const int param_m, 
					const int param_k, 
					MATRIX matrix_N, 
					MATRIX matrix_M, 
					double *q_coef,
					const int refinement_cnt,
					int *processing_stat,
					int *initial_subproblem_cnt,
					int **initial_subproblem_stat,
					int *x_start_index,
					int *ybar_start_index,
					int *mcc_start_row,
					int *mcc_start_col,
					double ***p_x_lb,
					double ***p_x_ub,
					double ***p_y_lb,
					double ***p_y_ub,
					const int	   x_p_cnt,
					const int	   y_p_cnt)
{
	//matrix need to be grouped by row
	//min c*x+d*y
	//st.  A*x+B*y		-s=b
	//     N*x+M*y-w		=-q
	//	   y>=0;w>=0;s>=0
	//	   (x>=0)
	//the order of var in LP is x,y,w,s
	//*processing_stat=-1: problem infeasible
	int status=0;
	int i,j,k;
	int lp_solnstat;
	int lp_num_row;
	int lp_num_col;
	int linnzcnt;
	int quadnzcnt;
	int nzcnt;
	int mccboundconstraint_position_row;
	int mccboundconstraint_position_col;
	int row_indice[4];
	int col_indice[4];
	//int *qmatbeg=NULL;
	//int *qmatcnt=NULL;
	//int *qmatind=NULL;
	int *linind=NULL;
	int *quadrow=NULL;
	int *quadcol=NULL;
	int matrix_M_status;
	int condition;
	int partition_cnt=0;
	int *partition_x_cnt=NULL;
	int *partition_y_cnt=NULL;
	const int n_col=CPXgetnumcols(env,lp);
	const int n_row=CPXgetnumrows(env,lp);
	//double *qmatval=NULL;
	double *linval=NULL;
	double *quadval=NULL;
	double *ybar_lb=NULL;
	double *ybar_ub=NULL;
	double *sigma_lb=NULL;
	double *sigma_ub=NULL;
	double *x_lb=NULL;
	double *x_ub=NULL;
	double **partition_x_lb=NULL;
	double **partition_x_ub=NULL;
	double **partition_ybar_lb=NULL;
	double **partition_ybar_ub=NULL;
	double **y_star=NULL;
	double *My_star=NULL;
	double y_starMy_star;
	double lbmcc;
	double rhs[4];
	double tmplb;
	double tmpub;
	double initial_lbmcc;
	double violation_level;
	double unbounded_flag;
	double tmp_val;
	const char l_sense='L';
	const char u_sense='U';
	MATRIX matrix_N_t;
	MATRIX matrix_M_t;
	MATRIX matrix_M_bar;// symmetrized M
	CPXLPptr copy_lp = NULL;
	STARTINFO lp_start;
	STARTINFO update_lp_start;
	clock_t t_start,t_end;
	*processing_stat=1;
	*x_start_index=-1;
	*ybar_start_index=-1;
	*mcc_start_col=-1;
	*mcc_start_row=-1;
	fprintf(stdout," Mccormick refinements (%d refinements)\n",refinement_cnt);
	//fprintf(stdout,"param_n:%d,param_m:%d,param_k:%d",param_n,param_m,param_k);
	init_startinfo(&lp_start);
	init_startinfo(&update_lp_start);
	//initialize
	copy_lp = CPXcloneprob(env,lp,&status);
	if (status) goto TERMINATE;
	init_matrix(&matrix_N_t);
	init_matrix(&matrix_M_t);
	init_matrix(&matrix_M_bar);
	status=copy_matrix(&matrix_N,&matrix_N_t);
	if (status) goto TERMINATE;
	status=transpose_matrix(&matrix_N_t);
	if (status) goto TERMINATE;
	status=copy_matrix(&matrix_M,&matrix_M_t);
	if (status) goto TERMINATE;
	status=transpose_matrix(&matrix_M_t);
	if (status) goto TERMINATE;
	//add var y_bar into lp
	ybar_lb=(double*) malloc(param_n*sizeof(double));
	ybar_ub=(double*) malloc(param_n*sizeof(double));
	if ( ybar_lb==NULL || ybar_ub==NULL)	
	{
		status=-1;
		fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<param_n;i++)
	{
		ybar_lb[i]=-CPX_INFBOUND;
		ybar_ub[i]=CPX_INFBOUND;
	}
	status=CPXnewcols(env,copy_lp,param_n,NULL,ybar_lb,ybar_ub,NULL,NULL);
	if (status) goto TERMINATE; 
	status=CPXnewrows(env,copy_lp,param_n,NULL,NULL,NULL,NULL);
	if (status) goto TERMINATE;
	for (i=0;i<param_n;i++)
	{
		for (j=0;j<matrix_N_t.matcnt[i];j++)
		{
			status=CPXchgcoef(env,copy_lp,n_row+i,param_n+matrix_N_t.matind[matrix_N_t.matbeg[i]+j],matrix_N_t.matval[matrix_N_t.matbeg[i]+j]);
			if (status) goto TERMINATE;
		}
		status=CPXchgcoef(env,copy_lp,n_row+i,n_col+i,-1.0);
		if (status) goto TERMINATE;
	}
	//add var sigma into lp (sigma_i=ybar_i*x_i)
	sigma_lb=(double*) malloc(param_n*sizeof(double));
	sigma_ub=(double*) malloc(param_n*sizeof(double));
	if(sigma_lb==NULL || sigma_ub==NULL)
	{
		status=-1;
		fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<param_n;i++)
	{
		sigma_lb[i]=-CPX_INFBOUND;
		sigma_ub[i]=CPX_INFBOUND;
	}
	status=CPXnewcols(env,copy_lp,param_n,NULL,sigma_lb,sigma_ub,NULL,NULL);
	if(status) goto TERMINATE;

	t_start=clock();
	// get x,ybar's lowerbound and upperbound
	x_ub=(double*) malloc(param_n*sizeof(double));
	x_lb=(double*) malloc(param_n*sizeof(double));
	if (x_ub==NULL || x_lb==NULL )
	{
		status=-1;
		fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<param_n;i++)
	{
		x_lb[i]=-CPX_INFBOUND;
		x_ub[i]=CPX_INFBOUND;
		ybar_lb[i]=-CPX_INFBOUND;
		ybar_ub[i]=CPX_INFBOUND;
		status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,i,CPX_MIN,&tmplb,&condition);
		if(status)goto TERMINATE;
		if (condition==1)
		{
			x_lb[i]=tmplb-10*ZERO_TOLERANCE;
		}
		status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,i,CPX_MAX,&tmpub,&condition);
		if(status)goto TERMINATE;
		if (condition==1)
		{
			x_ub[i]=tmpub+10*ZERO_TOLERANCE;
		}	
		status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,n_col+i,CPX_MIN,&tmplb,&condition);
		if(status)goto TERMINATE;
		if (condition==1)
		{
			ybar_lb[i]=tmplb-10*ZERO_TOLERANCE;
		}
		status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,n_col+i,CPX_MAX,&tmpub,&condition);
		if(status)goto TERMINATE;
		if (condition==1)
		{
			ybar_ub[i]=tmpub+10*ZERO_TOLERANCE;
		}		
	}
	//print the upper and lower bound for x_i, ybar_i and y_i
	//fprintf(stdout,"\nprint the initial upper and lower bound of x\n");
	//for (i=0;i<param_n;i++)
	//{
	//	fprintf(stdout,"x[%d]: ub:%f, lb:%f\n",i,x_ub[i],x_lb[i]);
	//}
	//fprintf(stdout,"\nprint the initial upper and lower bound of ybar\n");
	//for (i=0;i<param_n;i++)
	//{
	//	fprintf(stdout,"ybar[%d]: ub:%f, lb:%f\n",i,ybar_ub[i],ybar_lb[i]);
	//}
	//add the McCormick cuts from upper and lower bound of x[i] and ybar[i]
	violation_level=0;
	unbounded_flag=1;
	for (j=0;j<param_n;j++)
	{
		if (x_lb[j]>-CPX_INFBOUND && x_ub[j]<CPX_INFBOUND && ybar_lb[j] >-CPX_INFBOUND && ybar_ub[j]<CPX_INFBOUND)
		{
			violation_level+=0.5*fabs(x_ub[j]-x_lb[j])*fabs(ybar_ub[j]-ybar_lb[j]);
			//fprintf(stdout,"%d: %f\n",j,0.5*fabs(x_ub[j]-x_lb[j])*fabs(ybar_ub[j]-ybar_lb[j]));
		}else unbounded_flag=0;
	}
	//if (unbounded_flag)
	//{
	//	fprintf(stdout,"total violation level:%f\n",violation_level);
	//}
	mccboundconstraint_position_row=CPXgetnumrows(env,copy_lp);
	mccboundconstraint_position_col=CPXgetnumcols(env,copy_lp);
	for (i=0;i<param_n;i++)
	{
		lp_num_row=CPXgetnumrows(env,copy_lp);
		lp_num_col=CPXgetnumcols(env,copy_lp);
		if(x_lb[i]>-CPX_INFBOUND && ybar_lb[i]>-CPX_INFBOUND)rhs[0]=x_lb[i]*ybar_lb[i];
		else rhs[0]=0;
		if(x_ub[i]<CPX_INFBOUND && ybar_ub[i]<CPX_INFBOUND)rhs[1]=x_ub[i]*ybar_ub[i];
		else rhs[1]=0;
		if(x_lb[i]>-CPX_INFBOUND && ybar_ub[i]<CPX_INFBOUND)rhs[2]=x_lb[i]*ybar_ub[i];
		else rhs[2]=0;
		if(x_ub[i]<CPX_INFBOUND && ybar_lb[i]>-CPX_INFBOUND)rhs[3]=x_ub[i]*ybar_lb[i];
		else rhs[3]=0;

		status=CPXnewrows(env,copy_lp,4,rhs,NULL,NULL,NULL);
		if(status) goto TERMINATE;
		status=CPXnewcols(env,copy_lp,4,NULL,NULL,NULL,NULL,NULL);
		if(status) goto TERMINATE;
		if (x_lb[i]>-CPX_INFBOUND && ybar_lb[i]>-CPX_INFBOUND)
		{
			status=CPXchgcoef(env,copy_lp,lp_num_row,lp_num_col,1);
			if (status) goto TERMINATE;
			status=CPXchgcoef(env,copy_lp,lp_num_row,i,ybar_lb[i]);
			if (status) goto TERMINATE;
			status=CPXchgcoef(env,copy_lp,lp_num_row,n_col+i,x_lb[i]);
			if (status) goto TERMINATE;
			status=CPXchgcoef(env,copy_lp,lp_num_row,n_col+param_n+i,-1);
			if (status) goto TERMINATE;
		}
		if (x_ub[i]<CPX_INFBOUND && ybar_ub[i]<CPX_INFBOUND)
		{
			status=CPXchgcoef(env,copy_lp,lp_num_row+1,lp_num_col+1,1);
			if (status) goto TERMINATE;
			status=CPXchgcoef(env,copy_lp,lp_num_row+1,i,ybar_ub[i]);
			if (status) goto TERMINATE;
			status=CPXchgcoef(env,copy_lp,lp_num_row+1,n_col+i,x_ub[i]);
			if (status) goto TERMINATE;
			status=CPXchgcoef(env,copy_lp,lp_num_row+1,n_col+param_n+i,-1);
			if (status) goto TERMINATE;
		}
		if (x_lb[i]>-CPX_INFBOUND && ybar_ub[i]<CPX_INFBOUND)
		{
			status=CPXchgcoef(env,copy_lp,lp_num_row+2,lp_num_col+2,-1);
			if (status) goto TERMINATE;
			status=CPXchgcoef(env,copy_lp,lp_num_row+2,i,ybar_ub[i]);
			if (status) goto TERMINATE;
			status=CPXchgcoef(env,copy_lp,lp_num_row+2,n_col+i,x_lb[i]);
			if (status) goto TERMINATE;
			status=CPXchgcoef(env,copy_lp,lp_num_row+2,n_col+param_n+i,-1);
			if (status) goto TERMINATE;
		}
		if (x_ub[i]<CPX_INFBOUND && ybar_lb[i]>-CPX_INFBOUND)
		{
			status=CPXchgcoef(env,copy_lp,lp_num_row+3,lp_num_col+3,-1);
			if (status) goto TERMINATE;
			status=CPXchgcoef(env,copy_lp,lp_num_row+3,i,ybar_lb[i]);
			if (status) goto TERMINATE;
			status=CPXchgcoef(env,copy_lp,lp_num_row+3,n_col+i,x_ub[i]);
			if (status) goto TERMINATE;
			status=CPXchgcoef(env,copy_lp,lp_num_row+3,n_col+param_n+i,-1);
			if (status) goto TERMINATE;
		}
	}
	//get symmetrized M 
	status=symmetrizing_matrix(&matrix_M_t,&matrix_M_bar,&matrix_M_status);
    if(status) printf(" status after symmetrizing_matrix %d\n",status);
	if(status) goto TERMINATE;
	if (matrix_M_status==0)
	{
		////if M is skewed symmetric matrix, then we don't need to setup QCP to do refinements
		//fprintf(stdout," matrix M is skew symmetric matrix\n");
		//add sum{j in 1..m} q[j]*y[j] + sum{i in 1..n} sigma[i] = 0;
		status=CPXnewrows(env,copy_lp,1,NULL,NULL,NULL,NULL);
		if(status) goto TERMINATE;
		lp_num_row=CPXgetnumrows(env,copy_lp);
		for (i=0;i<param_m;i++)
		{
			status=CPXchgcoef(env,copy_lp,lp_num_row-1,param_n+i,q_coef[i]);
			if(status) goto TERMINATE;
		}
		for (i=0;i<param_n;i++)
		{
			status=CPXchgcoef(env,copy_lp,lp_num_row-1,n_col+param_n+i,1);
			if(status) goto TERMINATE;
		}
		status = CPXlpopt(env, copy_lp);
		if ( status ) goto TERMINATE;
		status = CPXsolution (env, copy_lp, &lp_solnstat, &lbmcc, NULL, NULL, NULL, NULL);
		if ( status ) goto TERMINATE;
		t_end=clock();
		//fprintf(stdout,"solnstat:%d\n",lp_solnstat);
		fprintf(stdout,"	LB: %12.4f PT: %12.3f\n",lbmcc,(double)(t_end-t_start)/CLOCKS_PER_SEC);
		//refinement
		t_start=clock();
		for (i=0;i<refinement_cnt;i++)
		{
			//fprintf(stdout,"%d refinement:\n",i);
			//refine bound
			for (j=0;j<param_n;j++)
			{
				status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,j,CPX_MIN,&tmplb,&condition);
				if(status)goto TERMINATE;
				if(x_lb[j]<tmplb && condition==1)x_lb[j]=tmplb-10*ZERO_TOLERANCE;
				if(condition==-1)
				{
					*processing_stat=-1;
					goto TERMINATE;
				}
				status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,j,CPX_MAX,&tmpub,&condition);
				if(status)goto TERMINATE;
				if(x_ub[j]>tmpub && condition==1)x_ub[j]=tmpub+10*ZERO_TOLERANCE;
				if(condition==-1)
				{
					*processing_stat=-1;
					goto TERMINATE;
				}
				status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,n_col+j,CPX_MIN,&tmplb,&condition);
				if(status)goto TERMINATE;
				if(ybar_lb[j]<tmplb && condition==1)ybar_lb[j]=tmplb-10*ZERO_TOLERANCE;
				if(condition==-1)
				{
					*processing_stat=-1;
					goto TERMINATE;
				}
				status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,n_col+j,CPX_MAX,&tmpub,&condition);
				if(status)goto TERMINATE;
				if(ybar_ub[j]>tmpub && condition==1)ybar_ub[j]=tmpub+10*ZERO_TOLERANCE;
				if(condition==-1)
				{
					*processing_stat=-1;
					goto TERMINATE;
				}
				if (x_lb[j]>x_ub[j])
				{
					x_lb[j]=x_ub[j];
				}
				if (ybar_lb[j]>ybar_ub[j])
				{
					ybar_lb[j]=ybar_ub[j];
				}
			}
			//modify the coefficent in copy_lp
			for (j=0;j<param_n;j++)
			{
				row_indice[0]=mccboundconstraint_position_row+4*j;
				row_indice[1]=mccboundconstraint_position_row+4*j+1;
				row_indice[2]=mccboundconstraint_position_row+4*j+2;
				row_indice[3]=mccboundconstraint_position_row+4*j+3;
				col_indice[0]=mccboundconstraint_position_col+4*j;
				col_indice[1]=mccboundconstraint_position_col+4*j+1;
				col_indice[2]=mccboundconstraint_position_col+4*j+2;
				col_indice[3]=mccboundconstraint_position_col+4*j+3;
				if (x_lb[j]>-CPX_INFBOUND && ybar_lb[j]>-CPX_INFBOUND)
				{
					rhs[0]=x_lb[j]*ybar_lb[j];
					status=CPXchgcoef(env,copy_lp,row_indice[0],j,ybar_lb[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[0],n_col+j,x_lb[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[0],col_indice[0],1);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[0],n_col+param_n+j,-1);
					if (status) goto TERMINATE;
				}else
				{
					rhs[0]=0;
				}

				if (x_ub[j]<CPX_INFBOUND && ybar_ub[j]<CPX_INFBOUND)
				{
					rhs[1]=x_ub[j]*ybar_ub[j];
					status=CPXchgcoef(env,copy_lp,row_indice[1],j,ybar_ub[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[1],n_col+j,x_ub[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[1],col_indice[1],1);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[1],n_col+param_n+j,-1);
					if (status) goto TERMINATE;
				}else
				{
					rhs[1]=0;
				}

				if (x_lb[j]>-CPX_INFBOUND && ybar_ub[j]<CPX_INFBOUND)
				{
					rhs[2]=x_lb[j]*ybar_ub[j];
					status=CPXchgcoef(env,copy_lp,row_indice[2],j,ybar_ub[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[2],n_col+j,x_lb[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[2],col_indice[2],-1);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[2],n_col+param_n+j,-1);
					if (status) goto TERMINATE;
				}else
				{
					rhs[2]=0;
				}

				if (x_ub[j]<CPX_INFBOUND && ybar_lb[j]>-CPX_INFBOUND)
				{
					rhs[3]=x_ub[j]*ybar_lb[j];
					status=CPXchgcoef(env,copy_lp,row_indice[3],j,ybar_lb[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[3],n_col+j,x_ub[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[3],col_indice[3],-1);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[3],n_col+param_n+j,-1);
					if (status) goto TERMINATE;
				}else
				{
					rhs[3]=0;
				}
				status=CPXchgrhs(env,copy_lp,4,row_indice,rhs);
				if(status) goto TERMINATE;
			}
			unbounded_flag=1;
			violation_level=0;
			for (j=0;j<param_n;j++)
			{
				if (x_lb[j]>-CPX_INFBOUND && x_ub[j]<CPX_INFBOUND && ybar_lb[j] >-CPX_INFBOUND && ybar_ub[j]<CPX_INFBOUND)
				{
					violation_level+=0.5*fabs(x_ub[j]-x_lb[j])*fabs(ybar_ub[j]-ybar_lb[j]);
					//fprintf(stdout,"%d: %f\n",j,0.5*fabs(x_ub[j]-x_lb[j])*fabs(ybar_ub[j]-ybar_lb[j]));
				}else unbounded_flag=0;
			}
			if (unbounded_flag)
			{
				fprintf(stdout,"total violation level:%f\n",violation_level);
			}
		}		
		unbounded_flag=1;
		for (i=0;i<param_n;i++)
		{
			if (x_lb[i]==-CPX_INFBOUND || x_ub[i]==CPX_INFBOUND && ybar_lb[i]==-CPX_INFBOUND && ybar_ub[i]==CPX_INFBOUND) unbounded_flag=0;
			//x_lb[i]=x_lb[i]-10*EP;
			//x_ub[i]=x_ub[i]+10*EP;
			//ybar_lb[i]=ybar_lb[i]-10*EP;
			//ybar_ub[i]=ybar_ub[i]+10*EP;
			//fprintf(stdout,"x[%d]_lb:%f,x[%d]_ub:%f,y[%d]_bar_lb:%f,y[%d]_bar_ub:%f\n",i,x_lb[i],i,x_ub[i],i,ybar_lb[i],i,ybar_ub[i]);
		}
		status = CPXlpopt (env, copy_lp);
		if ( status ) goto TERMINATE;
		status = CPXsolution (env, copy_lp, &lp_solnstat, &lbmcc, NULL, NULL, NULL, NULL);
		if ( status ) goto TERMINATE;
		t_end=clock();
		//fprintf(stdout,"solnstat:%d\n",lp_solnstat);
		fprintf(stdout,"	LB: %12.4f PT: %12.3f\n",lbmcc,(double)(t_end-t_start)/CLOCKS_PER_SEC);
		if (unbounded_flag && (param_n < 3))
		{
			//partition refinement:
			if ((partition_x_cnt=(int*)malloc(param_n*sizeof(int)))==NULL ||
				(partition_y_cnt=(int*)malloc(param_n*sizeof(int)))==NULL)
			{
				status=-1;
				goto TERMINATE;
			}
			for (i=0;i<param_n;i++)
			{
				partition_x_cnt[i]=1;
				partition_y_cnt[i]=1;
			}
			partition_x_cnt[0]=x_p_cnt;
			partition_x_cnt[1]=x_p_cnt;
			partition_y_cnt[0]=y_p_cnt;
			partition_y_cnt[1]=y_p_cnt;
			status=partition_range(x_lb,x_ub,ybar_lb,ybar_ub,&partition_x_lb,&partition_x_ub,&partition_ybar_lb,&partition_ybar_ub,param_n,partition_x_cnt,partition_y_cnt,&partition_cnt);
			if(status) goto TERMINATE;
			//fprintf(stdout,"%d\n",partition_cnt);
			*initial_subproblem_cnt=partition_cnt;
			*initial_subproblem_stat=(int*)malloc(partition_cnt*sizeof(int));
			y_star=(double**)malloc(partition_cnt*sizeof(double*));
			if (y_star==NULL || *initial_subproblem_stat==NULL)
			{
				status=-1;
				fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
				goto TERMINATE;
			}
			for (i=0;i<partition_cnt;i++)
			{
				fprintf(stdout,"	Further Mccormick refinement on partition %d\n",i);
				status=partition_McCormickRefinement(env,copy_lp,param_n,param_m,param_k,&((*initial_subproblem_stat)[i]),partition_x_lb[i],partition_x_ub[i],partition_ybar_lb[i],partition_ybar_ub[i],&(y_star[i]),3,n_col,mccboundconstraint_position_row,mccboundconstraint_position_col);
				if(status) goto TERMINATE;	
				if ((*initial_subproblem_stat)[i]==-1)
				{
					fprintf(stdout,"		this partition is pruned after refinement\n");
				}
			}	
			//add mccormick cut to the orignal LP
			for (i=0;i<param_n;i++)
			{
				status=CPXchgbds(env,lp,1,&i,&l_sense,&(x_lb[i]));
				if(status) goto TERMINATE;
				status=CPXchgbds(env,lp,1,&i,&u_sense,&(x_ub[i]));
				if(status) goto TERMINATE;
			}
			*x_start_index=0;
			//add y_bar
			*ybar_start_index=CPXgetnumcols(env,lp);
			status=CPXnewcols(env,lp,param_n,NULL,ybar_lb,ybar_ub,NULL,NULL);
			if (status) goto TERMINATE; 
			status=CPXnewrows(env,lp,param_n,NULL,NULL,NULL,NULL);
			if (status) goto TERMINATE;
			for (i=0;i<param_n;i++)
			{
				for (j=0;j<matrix_N_t.matcnt[i];j++)
				{
					status=CPXchgcoef(env,lp,n_row+i,param_n+matrix_N_t.matind[matrix_N_t.matbeg[i]+j],matrix_N_t.matval[matrix_N_t.matbeg[i]+j]);
					if (status) goto TERMINATE;
				}
				status=CPXchgcoef(env,lp,n_row+i,n_col+i,-1.0);
				if (status) goto TERMINATE;
			}
			//add var sigma into lp (sigma_i=ybar_i*x_i)
			status=CPXnewcols(env,lp,param_n,NULL,sigma_lb,sigma_ub,NULL,NULL);
			if(status) goto TERMINATE;
			//add the McCormick cuts from upper and lower bound of x[i] and ybar[i]
			*mcc_start_row=CPXgetnumrows(env,lp);
			*mcc_start_col=CPXgetnumcols(env,lp);
			for (i=0;i<param_n;i++)
			{
				lp_num_row=CPXgetnumrows(env,lp);
				lp_num_col=CPXgetnumcols(env,lp);
				rhs[0]=x_lb[i]*ybar_lb[i];
				rhs[1]=x_ub[i]*ybar_ub[i];
				rhs[2]=x_lb[i]*ybar_ub[i];
				rhs[3]=x_ub[i]*ybar_lb[i];

				status=CPXnewrows(env,lp,4,rhs,NULL,NULL,NULL);
				if(status) goto TERMINATE;
				status=CPXnewcols(env,lp,4,NULL,NULL,NULL,NULL,NULL);
				if(status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,lp_num_col,1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,lp_num_col+1,1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,lp_num_col+2,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,lp_num_col+3,-1);
				if (status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,i,ybar_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,i,ybar_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,i,ybar_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,i,ybar_lb[i]);
				if (status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,n_col+i,x_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,n_col+i,x_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,n_col+i,x_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,n_col+i,x_ub[i]);
				if (status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
			}
			//add the constraint related with mccormic cut
			//sum{j in 1..m} q[j]*y[j] +  sum{i in 1..n} sigma[i] = 0;
			status=CPXnewrows(env,lp,1,NULL,NULL,NULL,NULL);
			if(status) goto TERMINATE;
			lp_num_row=CPXgetnumrows(env,lp);
			for (i=0;i<param_m;i++)
			{
				status=CPXchgcoef(env,lp,lp_num_row-1,param_n+i,q_coef[i]);
				if(status) goto TERMINATE;
			}
			for (i=0;i<param_n;i++)
			{
				status=CPXchgcoef(env,lp,lp_num_row-1,n_col+param_n+i,1);
				if(status) goto TERMINATE;
			}
		}else
		{
			partition_cnt=1;
			*initial_subproblem_cnt=partition_cnt;
			*initial_subproblem_stat=(int*)malloc(partition_cnt*sizeof(int));
			y_star=(double**)malloc(partition_cnt*sizeof(double*));
			y_star[0]=(double*)malloc(param_m*sizeof(double));
			if (y_star==NULL || 
				*initial_subproblem_stat==NULL ||
				y_star[0]==NULL)
			{
				status=-1;
				fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
				goto TERMINATE;
			}
			//add mccormick cut to the orignal LP
			for (i=0;i<param_n;i++)
			{
				status=CPXchgbds(env,lp,1,&i,&l_sense,&(x_lb[i]));
				if(status) goto TERMINATE;
				status=CPXchgbds(env,lp,1,&i,&u_sense,&(x_ub[i]));
				if(status) goto TERMINATE;
			}
			//add y_bar
			status=CPXnewcols(env,lp,param_n,NULL,ybar_lb,ybar_ub,NULL,NULL);
			if (status) goto TERMINATE; 
			status=CPXnewrows(env,lp,param_n,NULL,NULL,NULL,NULL);
			if (status) goto TERMINATE;
			for (i=0;i<param_n;i++)
			{
				for (j=0;j<matrix_N_t.matcnt[i];j++)
				{
					status=CPXchgcoef(env,lp,n_row+i,param_n+matrix_N_t.matind[matrix_N_t.matbeg[i]+j],matrix_N_t.matval[matrix_N_t.matbeg[i]+j]);
					if (status) goto TERMINATE;
				}
				status=CPXchgcoef(env,lp,n_row+i,n_col+i,-1.0);
				if (status) goto TERMINATE;
			}
			//add var sigma into lp (sigma_i=ybar_i*x_i)
			status=CPXnewcols(env,lp,param_n,NULL,sigma_lb,sigma_ub,NULL,NULL);
			if(status) goto TERMINATE;
			//add the McCormick cuts from upper and lower bound of x[i] and ybar[i]
			for (i=0;i<param_n;i++)
			{
				lp_num_row=CPXgetnumrows(env,lp);
				lp_num_col=CPXgetnumcols(env,lp);
				rhs[0]=x_lb[i]*ybar_lb[i];
				rhs[1]=x_ub[i]*ybar_ub[i];
				rhs[2]=x_lb[i]*ybar_ub[i];
				rhs[3]=x_ub[i]*ybar_lb[i];

				status=CPXnewrows(env,lp,4,rhs,NULL,NULL,NULL);
				if(status) goto TERMINATE;
				status=CPXnewcols(env,lp,4,NULL,NULL,NULL,NULL,NULL);
				if(status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,lp_num_col,1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,lp_num_col+1,1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,lp_num_col+2,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,lp_num_col+3,-1);
				if (status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,i,ybar_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,i,ybar_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,i,ybar_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,i,ybar_lb[i]);
				if (status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,n_col+i,x_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,n_col+i,x_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,n_col+i,x_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,n_col+i,x_ub[i]);
				if (status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
			}
			//add the constraint related with mccormic cut
			//sum{j in 1..m} q[j]*y[j] +  sum{i in 1..n} sigma[i] = 0;
			status=CPXnewrows(env,lp,1,NULL,NULL,NULL,NULL);
			if(status) goto TERMINATE;
			lp_num_row=CPXgetnumrows(env,lp);
			for (i=0;i<param_m;i++)
			{
				status=CPXchgcoef(env,lp,lp_num_row-1,param_n+i,q_coef[i]);
				if(status) goto TERMINATE;
			}
			for (i=0;i<param_n;i++)
			{
				status=CPXchgcoef(env,lp,lp_num_row-1,n_col+param_n+i,1);
				if(status) goto TERMINATE;
			}
		}
	}else if(matrix_M_status==1)
	{
		//add the quadratic constraint related with mccormic cut
		//sum{j in 1..m} q[j]*y[j] + sum{i in 1..m, j in 1..m} Mt[i,j] * y[i] * y[j]+ sum{i in 1..n} sigma[i] <= 0;
/*		linnzcnt=param_m+param_n;
		linind=(int*) malloc(linnzcnt*sizeof(int));
		linval=(double*) malloc(linnzcnt*sizeof(double));
		if (linind==NULL || linval==NULL)
		{
			status=-1;
			fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
			goto TERMINATE;
		}
		for (i=0;i<param_m;i++)
		{
			linind[i]=param_n+i;
			linval[i]=q_coef[i];	
		}
		for (i=0;i<param_n;i++)
		{
			linind[param_m+i]=n_col+param_n+i;
			linval[param_m+i]=1;
		}
		quadnzcnt=matrix_M_bar.nnz;
		quadcol=(int*)malloc(quadnzcnt*sizeof(int));
		quadrow=(int*)malloc(quadnzcnt*sizeof(int));
		quadval=(double*)malloc(quadnzcnt*sizeof(double));
		if (quadcol==NULL || quadrow==NULL || quadval==NULL)
		{
			status=-1;
			fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
			goto TERMINATE;
		}
		if (matrix_M_bar.param==0)
		{
			//group by column
			nzcnt=0;
			for (i=0;i<param_m;i++)
			{
				for (j=0;j<matrix_M_bar.matcnt[i];j++)
				{
					quadcol[nzcnt]=param_n+i;
					quadrow[nzcnt]=param_n+matrix_M_bar.matind[matrix_M_bar.matbeg[i]+j];
					quadval[nzcnt]=matrix_M_bar.matval[matrix_M_bar.matbeg[i]+j];
					nzcnt++;
				}
			}
		}else if (matrix_M_bar.param==1)
		{
			//group by row
			nzcnt=0;
			for (i=0;i<param_m;i++)
			{
				for (j=0;j<matrix_M_bar.matcnt[i];j++)
				{
					quadrow[nzcnt]=param_n+i;
					quadcol[nzcnt]=param_n+matrix_M_bar.matind[matrix_M_bar.matbeg[i]+j];
					quadval[nzcnt]=matrix_M_bar.matval[matrix_M_bar.matbeg[i]+j];
					nzcnt++;
				}
			}
		}else
		{
			status=-1;
			fprintf(stderr," add_mccormickcut(): unexpected matrix param=%d",matrix_M_bar.param);
			goto TERMINATE;
		}
		//status = CPXsetintparam (env, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
		//if ( status )goto TERMINATE;
		status=CPXaddqconstr(env,copy_lp,linnzcnt,quadnzcnt,0.0,'L',linind,linval,quadrow,quadcol,quadval,NULL);
		if (status) goto TERMINATE;
		status = CPXbaropt (env, copy_lp);	
		if ( status ) 
		{
			//fprintf(stdout,"status:%d",status);
			goto TERMINATE;
		}
		
		status = CPXsolution (env, copy_lp, &lp_solnstat, &lbmcc, NULL, NULL, NULL, NULL);
		if ( status ) goto TERMINATE;
		t_end=clock();
		//fprintf(stdout,"%d\n",lp_solnstat);
		fprintf(stdout,"	LB: %12.4f PT: %12.3f\n",lbmcc,(double)(t_end-t_start)/CLOCKS_PER_SEC);
		initial_lbmcc=lbmcc;
		//refinement
		t_start=clock();
		for (i=0;i<refinement_cnt;i++)
		{
			//refine bound
			//fprintf(stdout,"refinement %d:\n",i);
			for (j=0;j<param_n;j++)
			{
				status=getVARbound_baropt(env,copy_lp,j,CPX_MIN,&tmplb,&condition);
				if(status)goto TERMINATE;
				if(x_lb[j]<tmplb &&condition==1)x_lb[j]=tmplb-10*ZERO_TOLERANCE;
				if(condition==-1)
				{
					*processing_stat=-1;
					goto TERMINATE;
				}
				status=getVARbound_baropt(env,copy_lp,j,CPX_MAX,&tmpub,&condition);
				if(status)goto TERMINATE;
				if(x_ub[j]>tmpub && condition==1)x_ub[j]=tmpub+10*ZERO_TOLERANCE;
				if(condition==-1)
				{
					*processing_stat=-1;
					goto TERMINATE;
				}
				status=getVARbound_baropt(env,copy_lp,n_col+j,CPX_MIN,&tmplb,&condition);
				if(status)goto TERMINATE;
				if(ybar_lb[j]<tmplb && condition==1)ybar_lb[j]=tmplb-10*ZERO_TOLERANCE;
				if(condition==-1)
				{
					*processing_stat=-1;
					goto TERMINATE;
				}
				status=getVARbound_baropt(env,copy_lp,n_col+j,CPX_MAX,&tmpub,&condition);
				if(status)goto TERMINATE;
				if(ybar_ub[j]>tmpub && condition==1)ybar_ub[j]=tmpub+10*ZERO_TOLERANCE;
				if(condition==-1)
				{
					*processing_stat=-1;
					goto TERMINATE;
				}
				if (x_lb[j]>x_ub[j])
				{
					x_lb[j]=x_ub[j];
				}
				if (ybar_lb[j]>ybar_ub[j])
				{
					ybar_lb[j]=ybar_ub[j];
				}
			}		
			//fprintf(stdout,"\nprint the upper and lower bound of x\n");
			//for (j=0;j<param_n;j++)
			//{
			//	fprintf(stdout,"x[%d]: ub:%f, lb:%f\n",j,x_ub[j],x_lb[j]);
			//}
			////fprintf(stdout,"\nprint the upper and lower bound of ybar\n");
			//for (j=0;j<param_n;j++)
			//{
			//	fprintf(stdout,"ybar[%d]: ub:%f, lb:%f\n",j,ybar_ub[j],ybar_lb[j]);
			//}
			//modify the coefficent in copy_lp
			for (j=0;j<param_n;j++)
			{				
				row_indice[0]=mccboundconstraint_position_row+4*j;
				row_indice[1]=mccboundconstraint_position_row+4*j+1;
				row_indice[2]=mccboundconstraint_position_row+4*j+2;
				row_indice[3]=mccboundconstraint_position_row+4*j+3;
				col_indice[0]=mccboundconstraint_position_col+4*j;
				col_indice[1]=mccboundconstraint_position_col+4*j+1;
				col_indice[2]=mccboundconstraint_position_col+4*j+2;
				col_indice[3]=mccboundconstraint_position_col+4*j+3;
				if (x_lb[j]>-CPX_INFBOUND && ybar_lb[j]>-CPX_INFBOUND)
				{
					rhs[0]=x_lb[j]*ybar_lb[j];
					status=CPXchgcoef(env,copy_lp,row_indice[0],j,ybar_lb[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[0],n_col+j,x_lb[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[0],col_indice[0],1);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[0],n_col+param_n+j,-1);
					if (status) goto TERMINATE;
				}else
				{
					rhs[0]=0;
				}

				if (x_ub[j]<CPX_INFBOUND && ybar_ub[j]<CPX_INFBOUND)
				{
					rhs[1]=x_ub[j]*ybar_ub[j];
					status=CPXchgcoef(env,copy_lp,row_indice[1],j,ybar_ub[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[1],n_col+j,x_ub[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[1],col_indice[1],1);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[1],n_col+param_n+j,-1);
					if (status) goto TERMINATE;
				}else
				{
					rhs[1]=0;
				}

				if (x_lb[j]>-CPX_INFBOUND && ybar_ub[j]<CPX_INFBOUND)
				{
					rhs[2]=x_lb[j]*ybar_ub[j];
					status=CPXchgcoef(env,copy_lp,row_indice[2],j,ybar_ub[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[2],n_col+j,x_lb[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[2],col_indice[2],-1);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[2],n_col+param_n+j,-1);
					if (status) goto TERMINATE;
				}else
				{
					rhs[2]=0;
				}

				if (x_ub[j]<CPX_INFBOUND && ybar_lb[j]>-CPX_INFBOUND)
				{
					rhs[3]=x_ub[j]*ybar_lb[j];
					status=CPXchgcoef(env,copy_lp,row_indice[3],j,ybar_lb[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[3],n_col+j,x_ub[j]);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[3],col_indice[3],-1);
					if (status) goto TERMINATE;
					status=CPXchgcoef(env,copy_lp,row_indice[3],n_col+param_n+j,-1);
					if (status) goto TERMINATE;
				}else
				{
					rhs[3]=0;
				}
				status=CPXchgrhs(env,copy_lp,4,row_indice,rhs);
				if(status) goto TERMINATE;
			}
			unbounded_flag=1;
			violation_level=0;
			for (j=0;j<param_n;j++)
			{
				if (x_lb[j]>-CPX_INFBOUND && x_ub[j]<CPX_INFBOUND && ybar_lb[j] >-CPX_INFBOUND && ybar_ub[j]<CPX_INFBOUND)
				{
					violation_level+=0.5*fabs(x_ub[j]-x_lb[j])*fabs(ybar_ub[j]-ybar_lb[j]);
					//fprintf(stdout,"%d: %f\n",j,0.5*fabs(x_ub[j]-x_lb[j])*fabs(ybar_ub[j]-ybar_lb[j]));
				}else unbounded_flag=0;
			}
			//if (unbounded_flag)
			//{
			//	fprintf(stdout,"total violation level:%f\n",violation_level);
			//}
		}
		unbounded_flag=1;
		for (i=0;i<param_n;i++)
		{
			if (x_lb[i]==-CPX_INFBOUND || x_ub[i]==CPX_INFBOUND && ybar_lb[i]==-CPX_INFBOUND && ybar_ub[i]==CPX_INFBOUND) unbounded_flag=0;
			if (x_lb[i]>x_ub[i])
			{
				x_lb[i]=x_ub[i];
			}
			if (ybar_lb[i]>ybar_ub[i])
			{
				ybar_lb[i]=ybar_ub[i];
			}
			x_lb[i]=x_lb[i]-10*ZERO_TOLERANCE;
			x_ub[i]=x_ub[i]+10*ZERO_TOLERANCE;
			ybar_lb[i]=ybar_lb[i]-10*ZERO_TOLERANCE;
			ybar_ub[i]=ybar_ub[i]+10*ZERO_TOLERANCE;
			//fprintf(stdout,"x[%d]_lb:%f,x[%d]_ub:%f,y[%d]_bar_lb:%f,y[%d]_bar_ub:%f\n",i,x_lb[i],i,x_ub[i],i,ybar_lb[i],i,ybar_ub[i]);
		}
		status = CPXbaropt (env, copy_lp);
		if ( status ) goto TERMINATE;
		status = CPXsolution (env, copy_lp, &lp_solnstat, &lbmcc, NULL, NULL, NULL, NULL);
		//fprintf(stdout,"solution status: %d\n",lp_solnstat);
		if ( status ) goto TERMINATE;
		t_end=clock();
		fprintf(stdout,"	LB: %12.4f PT: %12.3f\n",lbmcc,(double)(t_end-t_start)/CLOCKS_PER_SEC);
		if (unbounded_flag)
		{
			//partition refinement:
			if ((partition_x_cnt=(int*)malloc(param_n*sizeof(int)))==NULL ||
				(partition_y_cnt=(int*)malloc(param_n*sizeof(int)))==NULL)
			{
				status=-1;
				goto TERMINATE;
			}
			for (i=0;i<param_n;i++)
			{
				partition_x_cnt[i]=1;
				partition_y_cnt[i]=1;
			}
			partition_x_cnt[0]=x_p_cnt;
			partition_x_cnt[1]=x_p_cnt;
			partition_y_cnt[0]=y_p_cnt;
			partition_y_cnt[1]=y_p_cnt;
			status=partition_range(x_lb,x_ub,ybar_lb,ybar_ub,&partition_x_lb,&partition_x_ub,&partition_ybar_lb,&partition_ybar_ub,param_n,partition_x_cnt,partition_y_cnt,&partition_cnt);
			if(status) goto TERMINATE;
			//fprintf(stdout,"%d\n",partition_cnt);
			*initial_subproblem_cnt=partition_cnt;
			*initial_subproblem_stat=(int*)malloc(partition_cnt*sizeof(int));
			y_star=(double**)malloc(partition_cnt*sizeof(double*));
			if (y_star==NULL || *initial_subproblem_stat==NULL)
			{
				status=-1;
				fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
				goto TERMINATE;
			}
			for (i=0;i<partition_cnt;i++)
			{
				fprintf(stdout,"	Further Mccormick refinement on partition %d\n",i);
				status=partition_McCormickRefinement(env,copy_lp,param_n,param_m,param_k,&((*initial_subproblem_stat)[i]),partition_x_lb[i],partition_x_ub[i],partition_ybar_lb[i],partition_ybar_ub[i],&(y_star[i]),3,n_col,mccboundconstraint_position_row,mccboundconstraint_position_col);
				if(status) goto TERMINATE;	
				if ((*initial_subproblem_stat)[i]==-1)
				{
					fprintf(stdout,"			this partition is pruned after refinement\n");
				}
			}	
	
			//add mccormick cut to the orignal LP
			for (i=0;i<param_n;i++)
			{
				status=CPXchgbds(env,lp,1,&i,&l_sense,&(x_lb[i]));
				if(status) goto TERMINATE;
				status=CPXchgbds(env,lp,1,&i,&u_sense,&(x_ub[i]));
				if(status) goto TERMINATE;
			}
			*x_start_index=0;
			//add y_bar
			*ybar_start_index=CPXgetnumcols(env,lp);
			status=CPXnewcols(env,lp,param_n,NULL,ybar_lb,ybar_ub,NULL,NULL);
			if (status) goto TERMINATE; 
			status=CPXnewrows(env,lp,param_n,NULL,NULL,NULL,NULL);
			if (status) goto TERMINATE;
			for (i=0;i<param_n;i++)
			{
				for (j=0;j<matrix_N_t.matcnt[i];j++)
				{
					status=CPXchgcoef(env,lp,n_row+i,param_n+matrix_N_t.matind[matrix_N_t.matbeg[i]+j],matrix_N_t.matval[matrix_N_t.matbeg[i]+j]);
					if (status) goto TERMINATE;
				}
				status=CPXchgcoef(env,lp,n_row+i,n_col+i,-1.0);
				if (status) goto TERMINATE;
			}
			////add var sigma into lp (sigma_i=ybar_i*x_i)
			status=CPXnewcols(env,lp,param_n,NULL,sigma_lb,sigma_ub,NULL,NULL);
			if(status) goto TERMINATE;
			//add the McCormick cuts from upper and lower bound of x[i] and ybar[i]
			*mcc_start_row=CPXgetnumrows(env,lp);
			*mcc_start_col=CPXgetnumcols(env,lp);
			for (i=0;i<param_n;i++)
			{
				lp_num_row=CPXgetnumrows(env,lp);
				lp_num_col=CPXgetnumcols(env,lp);
				rhs[0]=x_lb[i]*ybar_lb[i];
				rhs[1]=x_ub[i]*ybar_ub[i];
				rhs[2]=x_lb[i]*ybar_ub[i];
				rhs[3]=x_ub[i]*ybar_lb[i];
				status=CPXnewrows(env,lp,4,rhs,NULL,NULL,NULL);
				if(status) goto TERMINATE;
				status=CPXnewcols(env,lp,4,NULL,NULL,NULL,NULL,NULL);
				if(status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,lp_num_col,1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,lp_num_col+1,1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,lp_num_col+2,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,lp_num_col+3,-1);
				if (status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,i,ybar_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,i,ybar_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,i,ybar_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,i,ybar_lb[i]);
				if (status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,n_col+i,x_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,n_col+i,x_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,n_col+i,x_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,n_col+i,x_ub[i]);
				if (status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
			}
			//add subgradient approximation
			My_star=(double*)malloc(param_m*sizeof(double));
			if (My_star==NULL)
			{
				status=-1;
				fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
				goto TERMINATE;
			}
			for (i=0;i<partition_cnt;i++)
			{
				if ((*initial_subproblem_stat)[i]!=-1)
				{
					y_starMy_star=0;
					for (j=0;j<param_m;j++)
					{
						My_star[j]=0;
						for (k=0;k<matrix_M_bar.matcnt[j];k++)
						{
							My_star[j]+=matrix_M_bar.matval[matrix_M_bar.matbeg[j]+k]*y_star[i][matrix_M_bar.matind[matrix_M_bar.matbeg[j]+k]];
							y_starMy_star+=matrix_M_bar.matval[matrix_M_bar.matbeg[j]+k]*y_star[i][matrix_M_bar.matind[matrix_M_bar.matbeg[j]+k]]*y_star[i][j];
						}
					}
					status=CPXnewrows(env,lp,1,&y_starMy_star,&l_sense,NULL,NULL);
					if(status) goto TERMINATE;
					lp_num_row=CPXgetnumrows(env,lp);
					for (j=0;j<param_m;j++)
					{
						status=CPXchgcoef(env,lp,lp_num_row-1,param_n+j,q_coef[j]+2*My_star[j]);
						if(status) goto TERMINATE;
					}
					for (j=0;j<param_n;j++)
					{
						status=CPXchgcoef(env,lp,lp_num_row-1,n_col+param_n+j,1);
						if(status) goto TERMINATE;
					}
				}
			}
		}else
		{
			partition_cnt=1;
			*initial_subproblem_cnt=partition_cnt;
			*initial_subproblem_stat=(int*)malloc(partition_cnt*sizeof(int));
			y_star=(double**)malloc(partition_cnt*sizeof(double*));
			y_star[0]=(double*)malloc(param_m*sizeof(double));
			if (y_star==NULL || 
				*initial_subproblem_stat==NULL ||
				y_star[0]==NULL)
			{
				status=-1;
				fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
				goto TERMINATE;
			}
			
			status=CPXgetx(env,copy_lp,y_star[0],param_n,param_n+param_m-1);
			//add mccormick cut to the orignal LP
			for (i=0;i<param_n;i++)
			{
				tmp_val=x_lb[i];
				status=CPXchgbds(env,lp,1,&i,&l_sense,&tmp_val);
				if(status) goto TERMINATE;
				tmp_val=x_ub[i];
				status=CPXchgbds(env,lp,1,&i,&u_sense,&tmp_val);
				if(status) goto TERMINATE;
			}
			//add y_bar
			status=CPXnewcols(env,lp,param_n,NULL,ybar_lb,ybar_ub,NULL,NULL);
			if (status) goto TERMINATE; 
			status=CPXnewrows(env,lp,param_n,NULL,NULL,NULL,NULL);
			if (status) goto TERMINATE;
			for (i=0;i<param_n;i++)
			{
				for (j=0;j<matrix_N_t.matcnt[i];j++)
				{
					status=CPXchgcoef(env,lp,n_row+i,param_n+matrix_N_t.matind[matrix_N_t.matbeg[i]+j],matrix_N_t.matval[matrix_N_t.matbeg[i]+j]);
					if (status) goto TERMINATE;
				}
				status=CPXchgcoef(env,lp,n_row+i,n_col+i,-1.0);
				if (status) goto TERMINATE;
			}
			////add var sigma into lp (sigma_i=ybar_i*x_i)
			status=CPXnewcols(env,lp,param_n,NULL,sigma_lb,sigma_ub,NULL,NULL);
			if(status) goto TERMINATE;
			//add the McCormick cuts from upper and lower bound of x[i] and ybar[i]
			for (i=0;i<param_n;i++)
			{
				lp_num_row=CPXgetnumrows(env,lp);
				lp_num_col=CPXgetnumcols(env,lp);
				rhs[0]=x_lb[i]*ybar_lb[i];
				rhs[1]=x_ub[i]*ybar_ub[i];
				rhs[2]=x_lb[i]*ybar_ub[i];
				rhs[3]=x_ub[i]*ybar_lb[i];
				status=CPXnewrows(env,lp,4,rhs,NULL,NULL,NULL);
				if(status) goto TERMINATE;
				status=CPXnewcols(env,lp,4,NULL,NULL,NULL,NULL,NULL);
				if(status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,lp_num_col,1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,lp_num_col+1,1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,lp_num_col+2,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,lp_num_col+3,-1);
				if (status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,i,ybar_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,i,ybar_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,i,ybar_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,i,ybar_lb[i]);
				if (status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,n_col+i,x_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,n_col+i,x_ub[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,n_col+i,x_lb[i]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,n_col+i,x_ub[i]);
				if (status) goto TERMINATE;

				status=CPXchgcoef(env,lp,lp_num_row,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+1,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+2,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,lp,lp_num_row+3,n_col+param_n+i,-1);
				if (status) goto TERMINATE;
			}
			//add subgradient approximation
			My_star=(double*)malloc(param_m*sizeof(double));
			if (My_star==NULL)
			{
				status=-1;
				fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
				goto TERMINATE;
			}
			for (i=0;i<partition_cnt;i++)
			{
				if ((*initial_subproblem_stat)[i]!=-1)
				{
					y_starMy_star=0;
					for (j=0;j<param_m;j++)
					{
						My_star[j]=0;
						for (k=0;k<matrix_M_bar.matcnt[j];k++)
						{
							My_star[j]+=matrix_M_bar.matval[matrix_M_bar.matbeg[j]+k]*y_star[i][matrix_M_bar.matind[matrix_M_bar.matbeg[j]+k]];
							y_starMy_star+=matrix_M_bar.matval[matrix_M_bar.matbeg[j]+k]*y_star[i][matrix_M_bar.matind[matrix_M_bar.matbeg[j]+k]]*y_star[i][j];
						}
					}
					status=CPXnewrows(env,lp,1,&y_starMy_star,&l_sense,NULL,NULL);
					if(status) goto TERMINATE;
					lp_num_row=CPXgetnumrows(env,lp);
					for (j=0;j<param_m;j++)
					{
						status=CPXchgcoef(env,lp,lp_num_row-1,param_n+j,q_coef[j]+2*My_star[j]);
						if(status) goto TERMINATE;
					}
					for (j=0;j<param_n;j++)
					{
						status=CPXchgcoef(env,lp,lp_num_row-1,n_col+param_n+j,1);
						if(status) goto TERMINATE;
					}
				}
			}
		}     */
	}else
	{
		status=-1;
		fprintf(stderr," add_mccormickcut(): unexpected matrix_M_status: %d\n",matrix_M_status);
		goto TERMINATE;
	}
	t_end=clock();
	fprintf(stdout,"	Total time for Mccormick refinement: %f\n",(double)(t_end-t_start)/CLOCKS_PER_SEC);
	*p_x_lb=partition_x_lb;
	*p_x_ub=partition_x_ub;
	*p_y_lb=partition_ybar_lb;
	*p_y_ub=partition_ybar_ub;

TERMINATE:
	if ( copy_lp != NULL ) CPXfreeprob (env, &copy_lp);
	free_startinfo(&lp_start);
	free_startinfo(&update_lp_start);
	free_matrix(&matrix_N_t);
	free_matrix(&matrix_M_t);
	free_matrix(&matrix_M_bar);
	free_and_null((char**)&linind);
	free_and_null((char**)&quadrow);
	free_and_null((char**)&quadcol);
	free_and_null((char**)&linval);
	free_and_null((char**)&quadval);
	free_and_null((char**)&sigma_lb);
	free_and_null((char**)&sigma_ub);
	free_and_null((char**)&x_lb);
	free_and_null((char**)&x_ub);
	free_and_null((char**) &ybar_lb);
	free_and_null((char**) &ybar_ub);
	free_and_null((char**) &My_star);
	if (y_star!=NULL)
	{
		for (i=0;i<partition_cnt;i++)
		{
			free_and_null((char**) &y_star[i]);
		}
		free_and_null((char**) &y_star);
	}
	free_and_null((char**)&partition_x_cnt);
	free_and_null((char**)&partition_y_cnt);
	return(status);
}


//int add_mccormickcut(CPXENVptr env,
//					 CPXLPptr lp, 
//					 const int param_n,
//					 const int param_m, 
//					 const int param_k, 
//					 MATRIX matrix_N, 
//					 MATRIX matrix_M, 
//					 double *q_coef,
//					 const int refinement_cnt)
//{
//	//matrix need to be grouped by row
//	//min c*x+d*y
//	//st.  A*x+B*y		-s=b
//	//     N*x+M*y-w		=-q
//	//	   y>=0;w>=0;s>=0
//	//	   (x>=0)
//	//the order of var in LP is x,y,w,s
//	int status=0;
//	int i,j,k;
//	int lp_solnstat;
//	int lp_num_row;
//	int lp_num_col;
//	int linnzcnt;
//	int quadnzcnt;
//	int nzcnt;
//	int mccboundconstraint_position_row;
//	int mccboundconstraint_position_col;
//	int row_indice[4];
//	int col_indice[4];
//	int *qmatbeg=NULL;
//	int *qmatcnt=NULL;
//	int *qmatind=NULL;
//	int *linind=NULL;
//	int *quadrow=NULL;
//	int *quadcol=NULL;
//	int matrix_M_status;
//	int flag_1;
//	int flag_2;
//	int flag_3;
//	int flag_4;
//	int n;
//	int condition;
//	int partition_cnt=0;
//	int *partition_x_cnt=NULL;
//	int *partition_y_cnt=NULL;
//	const int n_col=CPXgetnumcols(env,lp);
//	const int n_row=CPXgetnumrows(env,lp);
//	double *qmatval=NULL;
//	double *linval=NULL;
//	double *quadval=NULL;
//	double *ybar_lb=NULL;
//	double *ybar_ub=NULL;
//	double *sigma_lb=NULL;
//	double *sigma_ub=NULL;
//	double *x_lb=NULL;
//	double *x_ub=NULL;
//
//	double **partition_x_lb=NULL;
//	double **partition_x_ub=NULL;
//	double **partition_ybar_lb=NULL;
//	double **partition_ybar_ub=NULL;
//
//	double **y_star=NULL;
//	//double *lp_soln=NULL;
//	double *My_star=NULL;
//	double y_starMy_star;
//	//	double lb_yMy;
//	double lbmcc;
//	double rhs[4];
//	double tmplb;
//	double tmpub;
//	double initial_lbmcc;
//	double violation_level;
//	double unbounded_flag;
//	double sum_y_ub;
//	double sum_y_lb;
//	double sum_w_ub;
//	double sum_w_lb;
//	const char l_sense='L';
//	const char u_sense='U';
//	MATRIX matrix_N_t;
//	MATRIX matrix_M_t;
//	MATRIX matrix_M_bar;// symmetrized M
//	CPXLPptr copy_lp = NULL;
//	STARTINFO lp_start;
//	STARTINFO update_lp_start;
//	clock_t t_start,t_end;
//	fprintf(stdout," add Mccormick cuts (%d refinements)\n",refinement_cnt);
//	//double branch_range;
//	init_startinfo(&lp_start);
//	init_startinfo(&update_lp_start);
//	//initialize
//	copy_lp = CPXcloneprob(env,lp,&status);
//	if (status) goto TERMINATE;
//	init_matrix(&matrix_N_t);
//	init_matrix(&matrix_M_t);
//	init_matrix(&matrix_M_bar);
//	status=copy_matrix(&matrix_N,&matrix_N_t);
//	if (status) goto TERMINATE;
//	status=transpose_matrix(&matrix_N_t);
//	if (status) goto TERMINATE;
//	status=copy_matrix(&matrix_M,&matrix_M_t);
//	if (status) goto TERMINATE;
//	//print_full_matrix(&matrix_M_t);
//	status=transpose_matrix(&matrix_M_t);
//	if (status) goto TERMINATE;
//	//print_full_matrix(&matrix_M_t);
//	//add var y_bar into lp
//	ybar_lb=(double*) malloc(param_n*sizeof(double));
//	ybar_ub=(double*) malloc(param_n*sizeof(double));
//	if ( ybar_lb==NULL || ybar_ub==NULL)	
//	{
//		status=-1;
//		fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
//		goto TERMINATE;
//	}
//	for (i=0;i<param_n;i++)
//	{
//		ybar_lb[i]=-CPX_INFBOUND;
//		ybar_ub[i]=CPX_INFBOUND;
//	}
//	status=CPXnewcols(env,copy_lp,param_n,NULL,ybar_lb,ybar_ub,NULL,NULL);
//	if (status) goto TERMINATE; 
//	status=CPXnewrows(env,copy_lp,param_n,NULL,NULL,NULL,NULL);
//	if (status) goto TERMINATE;
//	for (i=0;i<param_n;i++)
//	{
//		for (j=0;j<matrix_N_t.matcnt[i];j++)
//		{
//			status=CPXchgcoef(env,copy_lp,n_row+i,param_n+matrix_N_t.matind[matrix_N_t.matbeg[i]+j],matrix_N_t.matval[matrix_N_t.matbeg[i]+j]);
//			if (status) goto TERMINATE;
//		}
//		status=CPXchgcoef(env,copy_lp,n_row+i,n_col+i,-1.0);
//		if (status) goto TERMINATE;
//	}
//
//	//add var sigma into lp (sigma_i=ybar_i*x_i)
//	sigma_lb=(double*) malloc(param_n*sizeof(double));
//	sigma_ub=(double*) malloc(param_n*sizeof(double));
//	if(sigma_lb==NULL || sigma_ub==NULL)
//	{
//		status=-1;
//		fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
//		goto TERMINATE;
//	}
//	for (i=0;i<param_n;i++)
//	{
//		sigma_lb[i]=-CPX_INFBOUND;
//		sigma_ub[i]=CPX_INFBOUND;
//	}
//	status=CPXnewcols(env,copy_lp,param_n,NULL,sigma_lb,sigma_ub,NULL,NULL);
//	if(status) goto TERMINATE;
//
//	t_start=clock();
//	// get x,y,ybar's lowerbound and upperbound
//	x_ub=(double*) malloc(param_n*sizeof(double));
//	x_lb=(double*) malloc(param_n*sizeof(double));
//	//y_ub=(double*) malloc(param_m*sizeof(double));
//	//y_lb=(double*) malloc(param_m*sizeof(double));
//	if (x_ub==NULL || x_lb==NULL )//|| y_ub== NULL || y_lb==NULL)
//	{
//		status=-1;
//		fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
//		goto TERMINATE;
//	}
//	for (i=0;i<param_n;i++)
//	{
//		status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,i,CPX_MIN,&tmplb,&condition);
//		if(status)goto TERMINATE;
//		status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,i,CPX_MAX,&tmpub,&condition);
//		if(status)goto TERMINATE;
//		x_lb[i]=tmplb;
//		x_ub[i]=tmpub;
//		status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,n_col+i,CPX_MIN,&tmplb,&condition);
//		if(status)goto TERMINATE;
//		status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,n_col+i,CPX_MAX,&tmpub,&condition);
//		if(status)goto TERMINATE;
//		ybar_lb[i]=tmplb;
//		ybar_ub[i]=tmpub;
//
//		////branching test
//		//if(i==0)
//		//{
//		//	//branch_range=fabs(x_lb[i]-x_ub[i])/10;
//		//	//x_lb[i]=x_ub[i]-branch_range*6;
//		//	//x_ub[i]=x_ub[i]-branch_range*5;
//		//	branch_range=fabs(ybar_ub[i]-ybar_lb[i])/10;
//		//	ybar_lb[i]=ybar_ub[i]-branch_range*5;
//		//	ybar_ub[i]=ybar_ub[i]-branch_range*4;
//		//}
//		//if(i==1)
//		//{
//		//	//branch_range=fabs(x_lb[i]-x_ub[i])/10;
//		//	//x_lb[i]=x_ub[i]-branch_range*4;
//		//	//x_ub[i]=x_ub[i]-branch_range*3;
//		//	branch_range=fabs(ybar_ub[i]-ybar_lb[i])/10;
//		//	ybar_lb[i]=ybar_ub[i]-branch_range*10;
//		//	ybar_ub[i]=ybar_ub[i]-branch_range*5;
//		//}
//	}
//	//print the upper and lower bound for x_i, ybar_i and y_i
//	fprintf(stdout,"\nprint the upper and lower bound of x\n");
//	for (i=0;i<param_n;i++)
//	{
//		fprintf(stdout,"x[%d]: ub:%f, lb:%f\n",i,x_ub[i],x_lb[i]);
//	}
//	fprintf(stdout,"\nprint the upper and lower bound of ybar\n");
//	for (i=0;i<param_n;i++)
//	{
//		fprintf(stdout,"ybar[%d]: ub:%f, lb:%f\n",i,ybar_ub[i],ybar_lb[i]);
//	}
//	//add the McCormick cuts from upper and lower bound of x[i] and ybar[i]
//	mccboundconstraint_position_row=CPXgetnumrows(env,copy_lp);
//	mccboundconstraint_position_col=CPXgetnumcols(env,copy_lp);
//	for (i=0;i<param_n;i++)
//	{
//		lp_num_row=CPXgetnumrows(env,copy_lp);
//		lp_num_col=CPXgetnumcols(env,copy_lp);
//		flag_1=1;
//		flag_2=1;
//		flag_3=1;
//		flag_4=1;
//		if(x_lb[i]==-CPX_INFBOUND)
//		{
//			flag_1=0;
//			flag_3=0;
//		}
//		if (ybar_lb[i]==-CPX_INFBOUND)
//		{
//			flag_1=0;
//			flag_4=0;
//		}
//		if (x_ub[i]==CPX_INFBOUND)
//		{
//			flag_2=0;
//			flag_4=0;
//		}
//		if (ybar_ub[i]==CPX_INFBOUND)
//		{
//			flag_2=0;
//			flag_3=0;
//		}
//		if(flag_1)rhs[0]=x_lb[i]*ybar_lb[i];
//		else rhs[0]=0;
//		if(flag_2)rhs[1]=x_ub[i]*ybar_ub[i];
//		else rhs[1]=0;
//		if(flag_3)rhs[2]=x_lb[i]*ybar_ub[i];
//		else rhs[2]=0;
//		if(flag_4)rhs[3]=x_ub[i]*ybar_lb[i];
//		else rhs[3]=0;
//
//		status=CPXnewrows(env,copy_lp,4,rhs,NULL,NULL,NULL);
//		if(status) goto TERMINATE;
//		status=CPXnewcols(env,copy_lp,4,NULL,NULL,NULL,NULL,NULL);
//		if(status) goto TERMINATE;
//		//fprintf(stdout,"%d,%d,%d,%d\n",flag_1,flag_2,flag_3,flag_4);
//		if (flag_1)
//		{
//			status=CPXchgcoef(env,copy_lp,lp_num_row,lp_num_col,1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,copy_lp,lp_num_row,i,ybar_lb[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,copy_lp,lp_num_row,n_col+i,x_lb[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,copy_lp,lp_num_row,n_col+param_n+i,-1);
//			if (status) goto TERMINATE;
//		}
//		if (flag_2)
//		{
//			status=CPXchgcoef(env,copy_lp,lp_num_row+1,lp_num_col+1,1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,copy_lp,lp_num_row+1,i,ybar_ub[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,copy_lp,lp_num_row+1,n_col+i,x_ub[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,copy_lp,lp_num_row+1,n_col+param_n+i,-1);
//			if (status) goto TERMINATE;
//		}
//		if (flag_3)
//		{
//			status=CPXchgcoef(env,copy_lp,lp_num_row+2,lp_num_col+2,-1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,copy_lp,lp_num_row+2,i,ybar_ub[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,copy_lp,lp_num_row+2,n_col+i,x_lb[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,copy_lp,lp_num_row+2,n_col+param_n+i,-1);
//			if (status) goto TERMINATE;
//		}
//		if (flag_4)
//		{
//			status=CPXchgcoef(env,copy_lp,lp_num_row+3,lp_num_col+3,-1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,copy_lp,lp_num_row+3,i,ybar_lb[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,copy_lp,lp_num_row+3,n_col+i,x_ub[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,copy_lp,lp_num_row+3,n_col+param_n+i,-1);
//			if (status) goto TERMINATE;
//		}
//	}
//	//get symmetrized M 
//	status=symmetrizing_matrix(&matrix_M_t,&matrix_M_bar,&matrix_M_status);
//	if(status) goto TERMINATE;
//	//print_full_matrix(&matrix_M_bar);
//	if (matrix_M_status==0)
//	{
//		//if M is skewed symmetric matrix, then we don't need to setup QCP to do refinements
//		fprintf(stdout," matrix M is skew symmetric matrix\n");
//		//add sum{j in 1..m} q[j]*y[j] + sum{i in 1..n} sigma[i] = 0;
//		status=CPXnewrows(env,copy_lp,1,NULL,NULL,NULL,NULL);
//		if(status) goto TERMINATE;
//		lp_num_row=CPXgetnumrows(env,copy_lp);
//		for (i=0;i<param_m;i++)
//		{
//			status=CPXchgcoef(env,copy_lp,lp_num_row-1,param_n+i,q_coef[i]);
//			if(status) goto TERMINATE;
//		}
//		for (i=0;i<param_n;i++)
//		{
//			status=CPXchgcoef(env,copy_lp,lp_num_row-1,n_col+param_n+i,1);
//			if(status) goto TERMINATE;
//		}
//
//		status = CPXlpopt(env, copy_lp);
//		if ( status ) goto TERMINATE;
//
//		status = CPXsolution (env, copy_lp, &lp_solnstat, &lbmcc, NULL, NULL, NULL, NULL);
//		if ( status ) goto TERMINATE;
//		t_end=clock();
//		fprintf(stdout,"solnstat:%d\n",lp_solnstat);
//		fprintf(stdout," %f, %12.3f\n",lbmcc,(double)(t_end-t_start)/CLOCKS_PER_SEC);
//		//refinement
//		t_start=clock();
//		y_star=(double**) calloc(refinement_cnt,sizeof(double*));
//		if (y_star==NULL)
//		{
//			status=-1;
//			fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
//			goto TERMINATE;
//		}
//		for (i=0;i<refinement_cnt;i++)
//		{
//			fprintf(stdout,"%d refinement:\n",i);
//			//refine bound
//			for (j=0;j<param_n;j++)
//			{
//
//				//status=get_refined_VAR_bound(env,copy_lp,param_n,param_m,param_k,j,CPX_MIN,2,&tmplb);
//				//if(status)goto TERMINATE;
//				//if(x_lb[j]<tmplb)x_lb[j]=tmplb;
//				//status=get_refined_VAR_bound(env,copy_lp,param_n,param_m,param_k,j,CPX_MAX,2,&tmpub);
//				//if(status)goto TERMINATE;
//				//if(x_ub[j]>tmpub)x_ub[j]=tmpub;
//				//status=get_refined_VAR_bound(env,copy_lp,param_n,param_m,param_k,n_col+j,CPX_MIN,2,&tmplb);
//				//if(status)goto TERMINATE;
//				//if(ybar_lb[j]<tmplb)ybar_lb[j]=tmplb;
//				//status=get_refined_VAR_bound(env,copy_lp,param_n,param_m,param_k,n_col+j,CPX_MAX,2,&tmpub);
//				//if(status)goto TERMINATE;
//				//if(ybar_ub[j]>tmpub)ybar_ub[j]=tmpub;
//
//
//				status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,j,CPX_MIN,&tmplb,&condition);
//				if(status)goto TERMINATE;
//				if(x_lb[j]<tmplb)x_lb[j]=tmplb;
//				status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,j,CPX_MAX,&tmpub,&condition);
//				if(status)goto TERMINATE;
//				if(x_ub[j]>tmpub)x_ub[j]=tmpub;
//				status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,n_col+j,CPX_MIN,&tmplb,&condition);
//				if(status)goto TERMINATE;
//				if(ybar_lb[j]<tmplb)ybar_lb[j]=tmplb;
//				status=getVARbound(env,copy_lp,&lp_start,&update_lp_start,n_col+j,CPX_MAX,&tmpub,&condition);
//				if(status)goto TERMINATE;
//				if(ybar_ub[j]>tmpub)ybar_ub[j]=tmpub;
//			}
//			//fprintf(stdout,"\nprint the upper and lower bound of x\n");
//			for (j=0;j<param_n;j++)
//			{
//				//fprintf(stdout,"x[%d]: ub:%f, lb:%f\n",j,x_ub[j],x_lb[j]);
//				//x_ub[j]=x_ub[j]+fabs(x_ub[j])*1e-5;
//				//x_lb[j]=x_lb[j]-fabs(x_lb[j])*1e-5;
//
//			}
//			//fprintf(stdout,"\nprint the upper and lower bound of ybar\n");
//			for (j=0;j<param_n;j++)
//			{
//				//fprintf(stdout,"ybar[%d]: ub:%f, lb:%f\n",j,ybar_ub[j],ybar_lb[j]);
//				//ybar_ub[j]=ybar_ub[j]+fabs(ybar_ub[j])*1e-5;
//				//ybar_lb[j]=ybar_lb[j]-fabs(ybar_lb[j])*1e-5;
//			}
//			//modify the coefficent in copy_lp
//			for (j=0;j<param_n;j++)
//			{
//				rhs[0]=x_lb[j]*ybar_lb[j];
//				rhs[1]=x_ub[j]*ybar_ub[j];
//				rhs[2]=x_lb[j]*ybar_ub[j];
//				rhs[3]=x_ub[j]*ybar_lb[j];
//				row_indice[0]=mccboundconstraint_position_row+4*j;
//				row_indice[1]=mccboundconstraint_position_row+4*j+1;
//				row_indice[2]=mccboundconstraint_position_row+4*j+2;
//				row_indice[3]=mccboundconstraint_position_row+4*j+3;
//				status=CPXchgrhs(env,copy_lp,4,row_indice,rhs);
//				if(status) goto TERMINATE;
//
//				status=CPXchgcoef(env,copy_lp,row_indice[0],j,ybar_lb[j]);
//				if (status) goto TERMINATE;
//				status=CPXchgcoef(env,copy_lp,row_indice[1],j,ybar_ub[j]);
//				if (status) goto TERMINATE;
//				status=CPXchgcoef(env,copy_lp,row_indice[2],j,ybar_ub[j]);
//				if (status) goto TERMINATE;
//				status=CPXchgcoef(env,copy_lp,row_indice[3],j,ybar_lb[j]);
//				if (status) goto TERMINATE;
//
//				status=CPXchgcoef(env,copy_lp,row_indice[0],n_col+j,x_lb[j]);
//				if (status) goto TERMINATE;
//				status=CPXchgcoef(env,copy_lp,row_indice[1],n_col+j,x_ub[j]);
//				if (status) goto TERMINATE;
//				status=CPXchgcoef(env,copy_lp,row_indice[2],n_col+j,x_lb[j]);
//				if (status) goto TERMINATE;
//				status=CPXchgcoef(env,copy_lp,row_indice[3],n_col+j,x_ub[j]);
//				if (status) goto TERMINATE;
//			}
//			//solve 
//			status = CPXlpopt (env, copy_lp);
//			if ( status ) goto TERMINATE;
//			//get y_star
//			y_star[i]=(double*)malloc(param_m*sizeof(double));
//			if (y_star[i]==NULL)
//			{
//				status=-1;
//				fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
//				goto TERMINATE;
//			}
//			status = CPXsolution (env, copy_lp, &lp_solnstat, &lbmcc, NULL, NULL, NULL, NULL);
//			if ( status ) goto TERMINATE;
//			status = CPXgetx(env,copy_lp,y_star[i],param_n,param_n+param_m-1);
//			if(status) goto TERMINATE;
//			t_end=clock();
//			fprintf(stdout,"solnstat:%d\n",lp_solnstat);
//			fprintf(stdout," %f,%12.3f\n",lbmcc,(double)(t_end-t_start)/CLOCKS_PER_SEC);
//		}	
//		//add mccormick cut to the orignal LP
//		for (i=0;i<param_n;i++)
//		{
//			fprintf(stdout,"x_lb:%f,x_ub:%f,y_bar_lb:%f,y_bar_ub:%f\n",x_lb[i],x_ub[i],ybar_lb[i],ybar_ub[i]);
//		}
//		for (i=0;i<param_n;i++)
//		{
//			status=CPXchgbds(env,lp,1,&i,&l_sense,&(x_lb[i]));
//			if(status) goto TERMINATE;
//			status=CPXchgbds(env,lp,1,&i,&u_sense,&(x_ub[i]));
//			if(status) goto TERMINATE;
//		}
//		//add y_bar
//		status=CPXnewcols(env,lp,param_n,NULL,ybar_lb,ybar_ub,NULL,NULL);
//		if (status) goto TERMINATE; 
//		status=CPXnewrows(env,lp,param_n,NULL,NULL,NULL,NULL);
//		if (status) goto TERMINATE;
//		for (i=0;i<param_n;i++)
//		{
//			for (j=0;j<matrix_N_t.matcnt[i];j++)
//			{
//				status=CPXchgcoef(env,lp,n_row+i,param_n+matrix_N_t.matind[matrix_N_t.matbeg[i]+j],matrix_N_t.matval[matrix_N_t.matbeg[i]+j]);
//				if (status) goto TERMINATE;
//			}
//			status=CPXchgcoef(env,lp,n_row+i,n_col+i,-1.0);
//			if (status) goto TERMINATE;
//		}
//		//add var sigma into lp (sigma_i=ybar_i*x_i)
//		status=CPXnewcols(env,lp,param_n,NULL,sigma_lb,sigma_ub,NULL,NULL);
//		if(status) goto TERMINATE;
//
//		//add the McCormick cuts from upper and lower bound of x[i] and ybar[i]
//		for (i=0;i<param_n;i++)
//		{
//			lp_num_row=CPXgetnumrows(env,lp);
//			lp_num_col=CPXgetnumcols(env,lp);
//			rhs[0]=x_lb[i]*ybar_lb[i];
//			rhs[1]=x_ub[i]*ybar_ub[i];
//			rhs[2]=x_lb[i]*ybar_ub[i];
//			rhs[3]=x_ub[i]*ybar_lb[i];
//
//			status=CPXnewrows(env,lp,4,rhs,NULL,NULL,NULL);
//			if(status) goto TERMINATE;
//			status=CPXnewcols(env,lp,4,NULL,NULL,NULL,NULL,NULL);
//			if(status) goto TERMINATE;
//
//			status=CPXchgcoef(env,lp,lp_num_row,lp_num_col,1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+1,lp_num_col+1,1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+2,lp_num_col+2,-1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+3,lp_num_col+3,-1);
//			if (status) goto TERMINATE;
//
//			status=CPXchgcoef(env,lp,lp_num_row,i,ybar_lb[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+1,i,ybar_ub[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+2,i,ybar_ub[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+3,i,ybar_lb[i]);
//			if (status) goto TERMINATE;
//
//			status=CPXchgcoef(env,lp,lp_num_row,n_col+i,x_lb[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+1,n_col+i,x_ub[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+2,n_col+i,x_lb[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+3,n_col+i,x_ub[i]);
//			if (status) goto TERMINATE;
//
//			status=CPXchgcoef(env,lp,lp_num_row,n_col+param_n+i,-1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+1,n_col+param_n+i,-1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+2,n_col+param_n+i,-1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+3,n_col+param_n+i,-1);
//			if (status) goto TERMINATE;
//		}
//		//add the constraint related with mccormic cut
//		//sum{j in 1..m} q[j]*y[j] +  sum{i in 1..n} sigma[i] = 0;
//		status=CPXnewrows(env,lp,1,NULL,NULL,NULL,NULL);
//		if(status) goto TERMINATE;
//		lp_num_row=CPXgetnumrows(env,lp);
//		for (i=0;i<param_m;i++)
//		{
//			status=CPXchgcoef(env,lp,lp_num_row-1,param_n+i,q_coef[i]);
//			if(status) goto TERMINATE;
//		}
//		for (i=0;i<param_n;i++)
//		{
//			status=CPXchgcoef(env,lp,lp_num_row-1,n_col+param_n+i,1);
//			if(status) goto TERMINATE;
//		}
//	}else if(matrix_M_status==1)
//	{
//		//add the quadratic constraint related with mccormic cut
//		//sum{j in 1..m} q[j]*y[j] + sum{i in 1..m, j in 1..m} Mt[i,j] * y[i] * y[j]+ sum{i in 1..n} sigma[i] <= 0;
//		linnzcnt=param_m+param_n;
//		linind=(int*) malloc(linnzcnt*sizeof(int));
//		linval=(double*) malloc(linnzcnt*sizeof(double));
//		if (linind==NULL || linval==NULL)
//		{
//			status=-1;
//			fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
//			goto TERMINATE;
//		}
//		for (i=0;i<param_m;i++)
//		{
//			linind[i]=param_n+i;
//			linval[i]=q_coef[i];	
//		}
//		for (i=0;i<param_n;i++)
//		{
//			linind[param_m+i]=n_col+param_n+i;
//			linval[param_m+i]=1;
//		}
//		quadnzcnt=matrix_M_bar.nnz;
//		quadcol=(int*)malloc(quadnzcnt*sizeof(int));
//		quadrow=(int*)malloc(quadnzcnt*sizeof(int));
//		quadval=(double*)malloc(quadnzcnt*sizeof(double));
//		if (quadcol==NULL || quadrow==NULL || quadval==NULL)
//		{
//			status=-1;
//			fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
//			goto TERMINATE;
//		}
//		if (matrix_M_bar.param==0)
//		{
//			//group by column
//			nzcnt=0;
//			for (i=0;i<param_m;i++)
//			{
//				for (j=0;j<matrix_M_bar.matcnt[i];j++)
//				{
//					quadcol[nzcnt]=param_n+i;
//					quadrow[nzcnt]=param_n+matrix_M_bar.matind[matrix_M_bar.matbeg[i]+j];
//					quadval[nzcnt]=matrix_M_bar.matval[matrix_M_bar.matbeg[i]+j];
//					nzcnt++;
//				}
//			}
//		}else if (matrix_M_bar.param==1)
//		{
//			//group by row
//			nzcnt=0;
//			for (i=0;i<param_m;i++)
//			{
//				for (j=0;j<matrix_M_bar.matcnt[i];j++)
//				{
//					quadrow[nzcnt]=param_n+i;
//					quadcol[nzcnt]=param_n+matrix_M_bar.matind[matrix_M_bar.matbeg[i]+j];
//					quadval[nzcnt]=matrix_M_bar.matval[matrix_M_bar.matbeg[i]+j];
//					nzcnt++;
//				}
//			}
//		}else
//		{
//			status=-1;
//			fprintf(stderr," add_mccormickcut(): unexpected matrix param=%d",matrix_M_bar.param);
//			goto TERMINATE;
//		}
//		status=CPXaddqconstr(env,copy_lp,linnzcnt,quadnzcnt,0.0,'L',linind,linval,quadrow,quadcol,quadval,NULL);
//		if (status) goto TERMINATE;
//		//status = CPXsetintparam (env, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
//		//if ( status )goto TERMINATE;
//		status = CPXbaropt (env, copy_lp);
//		if ( status ) 
//		{
//			fprintf(stdout,"status: %d\n",status);
//			goto TERMINATE;
//		}
//		status = CPXsolution (env, copy_lp, &lp_solnstat, &lbmcc, NULL, NULL, NULL, NULL);
//		if ( status ) goto TERMINATE;
//		t_end=clock();
//		fprintf(stdout,"%d\n",lp_solnstat);
//		fprintf(stdout,"initial %f, %12.3f\n",lbmcc,(double)(t_end-t_start)/CLOCKS_PER_SEC);
//		initial_lbmcc=lbmcc;
//		//refinement
//		t_start=clock();
//		y_star=(double**) malloc(refinement_cnt*sizeof(double*));
//		if (y_star==NULL)
//		{
//			status=-1;
//			fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
//			goto TERMINATE;
//		}
//
//		for (i=0;i<refinement_cnt;i++)
//		{
//			//refine bound
//			fprintf(stdout,"refinement %d:\n",i);
//			for (j=0;j<param_n;j++)
//			{
//				status=getVARbound_baropt(env,copy_lp,j,CPX_MIN,&tmplb,&condition);
//				if(status)goto TERMINATE;
//				if(x_lb[j]<tmplb)x_lb[j]=tmplb;
//				status=getVARbound_baropt(env,copy_lp,j,CPX_MAX,&tmpub,&condition);
//				if(status)goto TERMINATE;
//				if(x_ub[j]>tmpub)x_ub[j]=tmpub;
//				status=getVARbound_baropt(env,copy_lp,n_col+j,CPX_MIN,&tmplb,&condition);
//				if(status)goto TERMINATE;
//				if(ybar_lb[j]<tmplb)ybar_lb[j]=tmplb;
//				status=getVARbound_baropt(env,copy_lp,n_col+j,CPX_MAX,&tmpub,&condition);
//				if(status)goto TERMINATE;
//				if(ybar_ub[j]>tmpub)ybar_ub[j]=tmpub;
//			}
//
//			//fprintf(stdout,"\nprint the upper and lower bound of x\n");
//			for (j=0;j<param_n;j++)
//			{
//				fprintf(stdout,"x[%d]: ub:%f, lb:%f\n",j,x_ub[j],x_lb[j]);
//			}
//			//fprintf(stdout,"\nprint the upper and lower bound of ybar\n");
//			for (j=0;j<param_n;j++)
//			{
//				fprintf(stdout,"ybar[%d]: ub:%f, lb:%f\n",j,ybar_ub[j],ybar_lb[j]);
//			}
//			//modify the coefficent in copy_lp
//			for (j=0;j<param_n;j++)
//			{				
//				flag_1=1;
//				flag_2=1;
//				flag_3=1;
//				flag_4=1;
//				row_indice[0]=mccboundconstraint_position_row+4*j;
//				row_indice[1]=mccboundconstraint_position_row+4*j+1;
//				row_indice[2]=mccboundconstraint_position_row+4*j+2;
//				row_indice[3]=mccboundconstraint_position_row+4*j+3;
//				col_indice[0]=mccboundconstraint_position_col+4*j;
//				col_indice[1]=mccboundconstraint_position_col+4*j+1;
//				col_indice[2]=mccboundconstraint_position_col+4*j+2;
//				col_indice[3]=mccboundconstraint_position_col+4*j+3;
//				if (x_lb[j]>-CPX_INFBOUND && ybar_lb[j]>-CPX_INFBOUND)
//				{
//					rhs[0]=x_lb[j]*ybar_lb[j];
//					status=CPXchgcoef(env,copy_lp,row_indice[0],j,ybar_lb[j]);
//					if (status) goto TERMINATE;
//					status=CPXchgcoef(env,copy_lp,row_indice[0],n_col+j,x_lb[j]);
//					if (status) goto TERMINATE;
//					status=CPXchgcoef(env,copy_lp,row_indice[0],col_indice[0],1);
//					if (status) goto TERMINATE;
//					status=CPXchgcoef(env,copy_lp,row_indice[0],n_col+param_n+j,-1);
//					if (status) goto TERMINATE;
//				}else
//				{
//					rhs[0]=0;
//					//flag_1=0;
//				}
//
//				if (x_ub[j]<CPX_INFBOUND && ybar_ub[j]<CPX_INFBOUND)
//				{
//					rhs[1]=x_ub[j]*ybar_ub[j];
//					status=CPXchgcoef(env,copy_lp,row_indice[1],j,ybar_ub[j]);
//					if (status) goto TERMINATE;
//					status=CPXchgcoef(env,copy_lp,row_indice[1],n_col+j,x_ub[j]);
//					if (status) goto TERMINATE;
//					status=CPXchgcoef(env,copy_lp,row_indice[1],col_indice[1],1);
//					if (status) goto TERMINATE;
//					status=CPXchgcoef(env,copy_lp,row_indice[1],n_col+param_n+j,-1);
//					if (status) goto TERMINATE;
//				}else
//				{
//					rhs[1]=0;
//					//flag_2=0;
//				}
//
//				if (x_lb[j]>-CPX_INFBOUND && ybar_ub[j]<CPX_INFBOUND)
//				{
//					rhs[2]=x_lb[j]*ybar_ub[j];
//					status=CPXchgcoef(env,copy_lp,row_indice[2],j,ybar_ub[j]);
//					if (status) goto TERMINATE;
//					status=CPXchgcoef(env,copy_lp,row_indice[2],n_col+j,x_lb[j]);
//					if (status) goto TERMINATE;
//					status=CPXchgcoef(env,copy_lp,row_indice[2],col_indice[2],-1);
//					if (status) goto TERMINATE;
//					status=CPXchgcoef(env,copy_lp,row_indice[2],n_col+param_n+j,-1);
//					if (status) goto TERMINATE;
//				}else
//				{
//					rhs[2]=0;
//					//flag_3=0;
//				}
//
//				if (x_ub[j]<CPX_INFBOUND && ybar_lb[j]>-CPX_INFBOUND)
//				{
//					rhs[3]=x_ub[j]*ybar_lb[j];
//					status=CPXchgcoef(env,copy_lp,row_indice[3],j,ybar_lb[j]);
//					if (status) goto TERMINATE;
//					status=CPXchgcoef(env,copy_lp,row_indice[3],n_col+j,x_ub[j]);
//					if (status) goto TERMINATE;
//					status=CPXchgcoef(env,copy_lp,row_indice[3],col_indice[3],-1);
//					if (status) goto TERMINATE;
//					status=CPXchgcoef(env,copy_lp,row_indice[3],n_col+param_n+j,-1);
//					if (status) goto TERMINATE;
//				}else
//				{
//					rhs[3]=0;
//					//flag_4=0;
//				}
//				status=CPXchgrhs(env,copy_lp,4,row_indice,rhs);
//				if(status) goto TERMINATE;
//				//fprintf(stdout,"%d,%d,%d,%d\n",flag_1,flag_2,flag_3,flag_4);
//			}
//			unbounded_flag=1;
//			violation_level=0;
//			for (j=0;j<param_n;j++)
//			{
//				if (x_lb[j]>-CPX_INFBOUND && x_ub[j]<CPX_INFBOUND && ybar_lb[j] >-CPX_INFBOUND && ybar_ub[j]<CPX_INFBOUND)
//				{
//					violation_level+=0.5*fabs(x_ub[j]-x_lb[j])*fabs(ybar_ub[j]-ybar_lb[j]);
//					fprintf(stdout,"%d: %f\n",j,0.5*fabs(x_ub[j]-x_lb[j])*fabs(ybar_ub[j]-ybar_lb[j]));
//				}else unbounded_flag=0;
//			}
//			if (unbounded_flag)
//			{
//				fprintf(stdout,"total violation level:%f\n",violation_level);
//			}
//			y_star[i]=NULL;
//		}
//
//
//		status = CPXbaropt (env, copy_lp);
//		if ( status ) goto TERMINATE;
//		//get y_star
//		y_star[0]=(double*)malloc(param_m*sizeof(double));
//		if (y_star[0]==NULL)
//		{
//			status=-1;
//			fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
//			goto TERMINATE;
//		}
//		status = CPXsolution (env, copy_lp, &lp_solnstat, &lbmcc, NULL, NULL, NULL, NULL);
//		fprintf(stdout,"%d\n",lp_solnstat);
//		if ( status ) goto TERMINATE;
//		status = CPXgetx(env,copy_lp,y_star[0],param_n,param_n+param_m-1);
//		if(status) goto TERMINATE;
//		t_end=clock();
//		fprintf(stdout," %f,%12.3f\n",lbmcc,(double)(t_end-t_start)/CLOCKS_PER_SEC);
//
//		for (i=0;i<param_n;i++)
//		{
//			if (x_lb[i]>x_ub[i])
//			{
//				x_lb[i]=x_ub[i];
//			}
//			if (ybar_lb[i]>ybar_ub[i])
//			{
//				ybar_lb[i]=ybar_ub[i];
//			}
//			fprintf(stdout,"x[%d]_lb:%f,x[%d]_ub:%f,y[%d]_bar_lb:%f,y[%d]_bar_ub:%f\n",i,x_lb[i],i,x_ub[i],i,ybar_lb[i],i,ybar_ub[i]);
//		}
//		//partition refinement:
//		//partition_cnt=4;
//		//partition_x_ub=(double**) malloc(partition_cnt*sizeof(double*));
//		//partition_x_lb=(double**) malloc(partition_cnt*sizeof(double*));
//		//partition_ybar_ub=(double**) malloc(partition_cnt*sizeof(double*));
//		//partition_ybar_lb=(double**) malloc(partition_cnt*sizeof(double*));
//		//for (i=0;i<partition_cnt;i++)
//		//{
//		//	partition_x_lb[i]=(double*)malloc(param_n*sizeof(double));
//		//	partition_x_ub[i]=(double*)malloc(param_n*sizeof(double));
//		//	partition_ybar_lb[i]=(double*)malloc(param_n*sizeof(double));
//		//	partition_ybar_ub[i]=(double*)malloc(param_n*sizeof(double));
//		//	n=i;
//		//	for (j=0;j<param_n;j++)
//		//	{
//		//		//partition_x_lb[i][j]=x_lb[j];
//		//		//partition_x_ub[i][j]=x_ub[j];
//		//		//if (fmod(n,2)==1)
//		//		//{
//		//		//	partition_ybar_lb[i][j]=(ybar_ub[j]+ybar_lb[j])/2;
//		//		//	partition_ybar_ub[i][j]=ybar_ub[j];
//		//		//}else
//		//		//{
//		//		//	partition_ybar_lb[i][j]=ybar_lb[j];
//		//		//	partition_ybar_ub[i][j]=(ybar_ub[j]+ybar_lb[j])/2;
//		//		//}
//		//		partition_ybar_lb[i][j]=ybar_lb[j];
//		//		partition_ybar_ub[i][j]=ybar_ub[j];
//		//		if (j<2)
//		//		{
//		//			if (fmod(n,2)==1)
//		//			{
//		//				partition_x_lb[i][j]=(x_lb[j]+x_ub[j])/2;
//		//				partition_x_ub[i][j]=x_ub[j];
//		//			}else
//		//			{
//		//				partition_x_lb[i][j]=x_lb[j];
//		//				partition_x_ub[i][j]=(x_lb[j]+x_ub[j])/2;
//		//			}
//		//			n=n/2;
//		//		}else
//		//		{
//		//			partition_x_lb[i][j]=x_lb[j];
//		//			partition_x_ub[i][j]=x_ub[j];
//		//		}
//
//		//	}
//		//}
//		if ((partition_x_cnt=(int*)malloc(param_n*sizeof(int)))==NULL ||
//			(partition_y_cnt=(int*)malloc(param_n*sizeof(int)))==NULL)
//		{
//			status=-1;
//			goto TERMINATE;
//		}
//		for (i=0;i<param_n;i++)
//		{
//			partition_x_cnt[i]=1;
//			partition_y_cnt[i]=1;
//		}
//		partition_x_cnt[0]=2;
//		partition_x_cnt[1]=2;
//		partition_y_cnt[0]=2;
//		partition_y_cnt[1]=2;
//		status=partition_range(x_lb,x_ub,ybar_lb,ybar_ub,&partition_x_lb,&partition_x_ub,&partition_ybar_lb,&partition_ybar_ub,param_n,partition_x_cnt,partition_y_cnt,&partition_cnt);
//		if(status) goto TERMINATE;
//		fprintf(stdout,"%d\n",partition_cnt);
//		for (i=0;i<partition_cnt;i++)
//		{
//			if(i!=13)continue;
//			fprintf(stdout,"partition %d\n",i);
//			status=partition_McCormickRefinement(env,copy_lp,param_n,param_m,param_k,partition_x_lb[i],partition_x_ub[i],partition_ybar_lb[i],partition_ybar_ub[i],10,n_col,mccboundconstraint_position_row,mccboundconstraint_position_col);
//			fprintf(stdout,"%d\n\n\n",status);
//			//if(status) goto TERMINATE;
//			fprintf(stdout,"%f,%f;    %f,%f\n",partition_x_lb[i],partition_x_ub[i],partition_ybar_lb[i],partition_ybar_ub[i]);		
//			for (j=0;j<param_n;j++)
//			{
//				x_lb[j]=partition_x_lb[i][j];
//				x_ub[j]=partition_x_ub[i][j];
//				ybar_lb[j]=partition_ybar_lb[i][j];
//				ybar_ub[j]=partition_ybar_ub[i][j];
//			}
//		}
//
//		//get up and lower bound on sum of y and w
//		//y_ub
//		//status=zeroobjectivefunc(env,copy_lp);
//		//if(status)goto TERMINATE;
//		//for (i=0;i<param_m;i++)
//		//{
//		//	status=CPXchgcoef(env,copy_lp,-1,param_n+i,-1);
//		//	if (status) goto TERMINATE;
//		//}
//		//status =CPXbaropt(env,copy_lp);
//		//if(status) goto TERMINATE;
//		//status = CPXsolution (env, copy_lp, &lp_solnstat, &sum_y_ub, NULL, NULL, NULL, NULL);
//		//if ( status ) goto TERMINATE;
//		////fprintf(stdout,"%d\n",lp_solnstat);
//		//if(lp_solnstat!=1) sum_y_ub=CPX_INFBOUND;
//		//else sum_y_ub=-sum_y_ub;
//		////y_lb
//		//status=zeroobjectivefunc(env,copy_lp);
//		//if(status)goto TERMINATE;
//		//for (i=0;i<param_m;i++)
//		//{
//		//	status=CPXchgcoef(env,copy_lp,-1,param_n+i,1);
//		//	if (status) goto TERMINATE;
//		//}
//		//status =CPXbaropt(env,copy_lp);
//		//if(status) goto TERMINATE;
//		//status = CPXsolution (env, copy_lp, &lp_solnstat, &sum_y_lb, NULL, NULL, NULL, NULL);
//		//if ( status ) goto TERMINATE;
//		//fprintf(stdout,"%d\n",lp_solnstat);
//		//if(lp_solnstat!=1) sum_y_lb=-CPX_INFBOUND;
//		////w_ub
//		//status=zeroobjectivefunc(env,copy_lp);
//		//if(status)goto TERMINATE;
//		//for (i=0;i<param_m;i++)
//		//{
//		//	status=CPXchgcoef(env,copy_lp,-1,param_n+param_m+i,-1);
//		//	if (status) goto TERMINATE;
//		//}
//		//status =CPXbaropt(env,copy_lp);
//		//if(status) goto TERMINATE;
//		//status = CPXsolution (env, copy_lp, &lp_solnstat, &sum_w_ub, NULL, NULL, NULL, NULL);
//		//if ( status ) goto TERMINATE;
//		//fprintf(stdout,"%d\n",lp_solnstat);
//		//if(lp_solnstat!=1) sum_w_ub=CPX_INFBOUND;
//		//else sum_w_ub=-sum_w_ub;
//		////w_lb
//		//status=zeroobjectivefunc(env,copy_lp);
//		//if(status)goto TERMINATE;
//		//for (i=0;i<param_m;i++)
//		//{
//		//	status=CPXchgcoef(env,copy_lp,-1,param_n+param_m+i,1);
//		//	if (status) goto TERMINATE;
//		//}
//		//status =CPXbaropt(env,copy_lp);
//		//if(status) goto TERMINATE;
//		//status = CPXsolution (env, copy_lp, &lp_solnstat, &sum_w_lb, NULL, NULL, NULL, NULL);
//		//if ( status ) goto TERMINATE;
//		//fprintf(stdout,"%d\n",lp_solnstat);
//		//if(lp_solnstat!=1) sum_w_lb=-CPX_INFBOUND;
//		////end
//		//fprintf(stdout,"sum y ub: %f, lb: %f\n",sum_y_ub,sum_y_lb);
//		//fprintf(stdout,"sum w ub: %f, lb: %f\n",sum_w_ub,sum_w_lb);
//
//
//
//		//zero copy_lp objective
//		//status=zeroobjectivefunc(env,copy_lp);
//		//if(status) goto TERMINATE;
//		////add quadratic objective function
//		////set up Q matrix
//		//lp_num_col=CPXgetnumcols(env,copy_lp);
//		//qmatbeg=(int*) calloc(lp_num_col,sizeof(int));
//		//qmatcnt=(int*) calloc(lp_num_col,sizeof(int));
//		//qmatind=(int*) calloc(matrix_M_bar.nnz,sizeof(int));
//		//qmatval=(double*) calloc(matrix_M_bar.nnz,sizeof(double));
//		//if (qmatbeg==NULL || qmatcnt==NULL || qmatind==NULL || qmatval==NULL)
//		//{
//		//	status=-1;
//		//	fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
//		//	goto TERMINATE;
//		//}
//		//for (i=0;i<lp_num_col;i++)
//		//{
//		//	if (i<param_n)qmatcnt[i]=0;
//		//	else if (i>=param_n && i<param_n+param_m)
//		//	{
//		//		qmatcnt[i]=matrix_M_bar.matcnt[i-param_n];
//		//	}else qmatcnt[i]=0;
//		//	if (i==0) qmatbeg[i]=0;
//		//	else qmatbeg[i]=qmatbeg[i-1]+qmatcnt[i-1];
//		//}
//		//for (i=0;i<matrix_M_bar.nnz;i++)
//		//{
//		//	qmatval[i]=matrix_M_bar.matval[i];
//		//	qmatind[i]=matrix_M_bar.matind[i]+param_n;
//		//}
//		//status=CPXcopyquad(env,copy_lp,qmatbeg,qmatcnt,qmatind,qmatval);
//		//if(status) goto TERMINATE;
//		//status = CPXbaropt (env, copy_lp);
//		//if (status==CPXERR_Q_NOT_POS_DEF)
//		//{
//		//	fprintf(stdout," matrix M is not psd\n");
//		//	status=0;
//		//	goto TERMINATE;
//		//}else if(status) goto TERMINATE;
//
//		//status = CPXsolution (env, copy_lp, &lp_solnstat, &lb_yMy, NULL, NULL, NULL, NULL);
//		//if ( status ) goto TERMINATE;
//		//lb_yMy*=-2;
//		//fprintf(stdout," ytMty's lb is %f\n",-lb_yMy);
//		//add mccormick cut to the orignal LP
//
//
//
//		for (i=0;i<param_n;i++)
//		{
//			status=CPXchgbds(env,lp,1,&i,&l_sense,&(x_lb[i]));
//			if(status) goto TERMINATE;
//			status=CPXchgbds(env,lp,1,&i,&u_sense,&(x_ub[i]));
//			if(status) goto TERMINATE;
//		}
//		//add y_bar
//		status=CPXnewcols(env,lp,param_n,NULL,ybar_lb,ybar_ub,NULL,NULL);
//		if (status) goto TERMINATE; 
//		status=CPXnewrows(env,lp,param_n,NULL,NULL,NULL,NULL);
//		if (status) goto TERMINATE;
//		for (i=0;i<param_n;i++)
//		{
//			for (j=0;j<matrix_N_t.matcnt[i];j++)
//			{
//				status=CPXchgcoef(env,lp,n_row+i,param_n+matrix_N_t.matind[matrix_N_t.matbeg[i]+j],matrix_N_t.matval[matrix_N_t.matbeg[i]+j]);
//				if (status) goto TERMINATE;
//			}
//			status=CPXchgcoef(env,lp,n_row+i,n_col+i,-1.0);
//			if (status) goto TERMINATE;
//		}
//		////add var sigma into lp (sigma_i=ybar_i*x_i)
//		status=CPXnewcols(env,lp,param_n,NULL,sigma_lb,sigma_ub,NULL,NULL);
//		if(status) goto TERMINATE;
//
//		//add the McCormick cuts from upper and lower bound of x[i] and ybar[i]
//
//		for (i=0;i<param_n;i++)
//		{
//			lp_num_row=CPXgetnumrows(env,lp);
//			lp_num_col=CPXgetnumcols(env,lp);
//			rhs[0]=x_lb[i]*ybar_lb[i];
//			rhs[1]=x_ub[i]*ybar_ub[i];
//			rhs[2]=x_lb[i]*ybar_ub[i];
//			rhs[3]=x_ub[i]*ybar_lb[i];
//			status=CPXnewrows(env,lp,4,rhs,NULL,NULL,NULL);
//			if(status) goto TERMINATE;
//			status=CPXnewcols(env,lp,4,NULL,NULL,NULL,NULL,NULL);
//			if(status) goto TERMINATE;
//
//			status=CPXchgcoef(env,lp,lp_num_row,lp_num_col,1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+1,lp_num_col+1,1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+2,lp_num_col+2,-1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+3,lp_num_col+3,-1);
//			if (status) goto TERMINATE;
//
//			status=CPXchgcoef(env,lp,lp_num_row,i,ybar_lb[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+1,i,ybar_ub[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+2,i,ybar_ub[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+3,i,ybar_lb[i]);
//			if (status) goto TERMINATE;
//
//			status=CPXchgcoef(env,lp,lp_num_row,n_col+i,x_lb[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+1,n_col+i,x_ub[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+2,n_col+i,x_lb[i]);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+3,n_col+i,x_ub[i]);
//			if (status) goto TERMINATE;
//
//			status=CPXchgcoef(env,lp,lp_num_row,n_col+param_n+i,-1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+1,n_col+param_n+i,-1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+2,n_col+param_n+i,-1);
//			if (status) goto TERMINATE;
//			status=CPXchgcoef(env,lp,lp_num_row+3,n_col+param_n+i,-1);
//			if (status) goto TERMINATE;
//		}
//		//add the constraint related with mccormic cut
//		//sum{j in 1..m} q[j]*y[j] +  sum{i in 1..n} sigma[i]  0;
//		//status=CPXnewrows(env,lp,1,&lb_yMy,&l_sense,NULL,NULL);
//		//if(status) goto TERMINATE;
//		//lp_num_row=CPXgetnumrows(env,lp);
//		//for (i=0;i<param_m;i++)
//		//{
//		//	status=CPXchgcoef(env,lp,lp_num_row-1,param_n+i,q_coef[i]);
//		//	if(status) goto TERMINATE;
//		//}
//		//for (i=0;i<param_n;i++)
//		//{
//		//	status=CPXchgcoef(env,lp,lp_num_row-1,n_col+param_n+i,1);
//		//	if(status) goto TERMINATE;
//		//}
//		//add subgradient approximation
//		My_star=(double*)malloc(param_m*sizeof(double));
//		if (My_star==NULL)
//		{
//			status=-1;
//			fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
//			goto TERMINATE;
//		}
//		//print_full_matrix(&matrix_M_bar);
//		for (i=0;i<1;i++)
//		{
//			y_starMy_star=0;
//			for (j=0;j<param_m;j++)
//			{
//				My_star[j]=0;
//				for (k=0;k<matrix_M_bar.matcnt[j];k++)
//				{
//					My_star[j]+=matrix_M_bar.matval[matrix_M_bar.matbeg[j]+k]*y_star[i][matrix_M_bar.matind[matrix_M_bar.matbeg[j]+k]];
//					y_starMy_star+=matrix_M_bar.matval[matrix_M_bar.matbeg[j]+k]*y_star[i][matrix_M_bar.matind[matrix_M_bar.matbeg[j]+k]]*y_star[i][j];
//				}
//			}
//			//fprintf(stdout,"y_star %d:\n",i);
//			//for (j=0;j<param_m;j++)
//			//{
//			//	fprintf(stdout,"%f\n",y_star[i][j]);
//			//}
//			//fprintf(stdout,"My_star %d:\n",i);
//			//for (j=0;j<param_m;j++)
//			//{
//			//	fprintf(stdout,"%f\n",My_star[j]);
//			//}
//			//fprintf(stdout,"%f\n",y_starMy_star);
//			status=CPXnewrows(env,lp,1,&y_starMy_star,&l_sense,NULL,NULL);
//			if(status) goto TERMINATE;
//			lp_num_row=CPXgetnumrows(env,lp);
//			for (j=0;j<param_m;j++)
//			{
//				status=CPXchgcoef(env,lp,lp_num_row-1,param_n+j,q_coef[j]+2*My_star[j]);
//				if(status) goto TERMINATE;
//			}
//			for (j=0;j<param_n;j++)
//			{
//				status=CPXchgcoef(env,lp,lp_num_row-1,n_col+param_n+j,1);
//				if(status) goto TERMINATE;
//			}
//		}
//		////add sum y and w ub and lb
//		//if (sum_y_lb>-CPX_INFBOUND)
//		//{
//		//	status=CPXnewrows(env,lp,1,&sum_y_lb,NULL,NULL,NULL);
//		//	if(status) goto TERMINATE;
//		//	for (i=0;i<param_m;i++)
//		//	{
//		//		status=CPXchgcoef(env,lp,CPXgetnumrows(env,lp)-1,param_n+i,1);
//		//		if(status) goto TERMINATE;
//		//	}
//		//	status=CPXnewcols(env,lp,1,NULL,NULL,NULL,NULL,NULL);
//		//	if(status) goto TERMINATE;
//		//	status=CPXchgcoef(env,lp,CPXgetnumrows(env,lp)-1,CPXgetnumcols(env,lp)-1,-1);
//		//	if(status) goto TERMINATE;
//		//} 
//		//if (sum_y_ub<CPX_INFBOUND)
//		//{
//		//	status=CPXnewrows(env,lp,1,&sum_y_ub,NULL,NULL,NULL);
//		//	if(status) goto TERMINATE;
//		//	for (i=0;i<param_m;i++)
//		//	{
//		//		status=CPXchgcoef(env,lp,CPXgetnumrows(env,lp)-1,param_n+i,1);
//		//		if(status) goto TERMINATE;
//		//	}
//		//	status=CPXnewcols(env,lp,1,NULL,NULL,NULL,NULL,NULL);
//		//	if(status) goto TERMINATE;
//		//	status=CPXchgcoef(env,lp,CPXgetnumrows(env,lp)-1,CPXgetnumcols(env,lp)-1,1);
//		//	if(status) goto TERMINATE;
//		//}
//		//if (sum_w_lb>-CPX_INFBOUND)
//		//{
//		//	status=CPXnewrows(env,lp,1,&sum_w_lb,NULL,NULL,NULL);
//		//	if(status) goto TERMINATE;
//		//	for (i=0;i<param_m;i++)
//		//	{
//		//		status=CPXchgcoef(env,lp,CPXgetnumrows(env,lp)-1,param_n+param_m+i,1);
//		//		if(status) goto TERMINATE;
//		//	}
//		//	status=CPXnewcols(env,lp,1,NULL,NULL,NULL,NULL,NULL);
//		//	if(status) goto TERMINATE;
//		//	status=CPXchgcoef(env,lp,CPXgetnumrows(env,lp)-1,CPXgetnumcols(env,lp)-1,-1);
//		//	if(status) goto TERMINATE;
//		//}
//		//if (sum_w_ub<CPX_INFBOUND)
//		//{
//		//	status=CPXnewrows(env,lp,1,&sum_w_ub,NULL,NULL,NULL);
//		//	if(status) goto TERMINATE;
//		//	for (i=0;i<param_m;i++)
//		//	{
//		//		status=CPXchgcoef(env,lp,CPXgetnumrows(env,lp)-1,param_n+param_m+i,1);
//		//		if(status) goto TERMINATE;
//		//	}
//		//	status=CPXnewcols(env,lp,1,NULL,NULL,NULL,NULL,NULL);
//		//	if(status) goto TERMINATE;
//		//	status=CPXchgcoef(env,lp,CPXgetnumrows(env,lp)-1,CPXgetnumcols(env,lp)-1,1);
//		//	if(status) goto TERMINATE;
//		//}
//
//		//status = CPXsetintparam (env, CPX_PARAM_NUMERICALEMPHASIS, CPX_OFF);
//		//if ( status )goto TERMINATE;
//	}else
//	{
//		status=-1;
//		fprintf(stderr," add_mccormickcut(): unexpected matrix_M_status: %d\n",matrix_M_status);
//		goto TERMINATE;
//	}
//
//
//TERMINATE:
//	if ( copy_lp != NULL ) CPXfreeprob (env, &copy_lp);
//	free_startinfo(&lp_start);
//	free_startinfo(&update_lp_start);
//	free_matrix(&matrix_N_t);
//	free_matrix(&matrix_M_t);
//	free_matrix(&matrix_M_bar);
//	free_and_null((char**)&qmatbeg);
//	free_and_null((char**)&qmatcnt);
//	free_and_null((char**)&qmatind);
//	free_and_null((char**)&qmatval);
//	free_and_null((char**)&linind);
//	free_and_null((char**)&quadrow);
//	free_and_null((char**)&quadcol);
//	free_and_null((char**)&linval);
//	free_and_null((char**)&quadval);
//	free_and_null((char**)&sigma_lb);
//	free_and_null((char**)&sigma_ub);
//	free_and_null((char**)&x_lb);
//	free_and_null((char**)&x_ub);
//	//free_and_null((char**)&y_lb);
//	//free_and_null((char**)&y_ub);
//	free_and_null((char**) &ybar_lb);
//	free_and_null((char**) &ybar_ub);
//	free_and_null((char**) &My_star);
//	if (partition_x_lb!=NULL)
//	{
//		for (i=0;i<partition_cnt;i++)
//		{
//			free_and_null((char**)&(partition_x_lb[i]));
//		}
//		free_and_null((char**) &partition_x_lb);
//	}
//	if (partition_x_ub!=NULL)
//	{
//		for (i=0;i<partition_cnt;i++)
//		{
//			free_and_null((char**) &(partition_x_ub[i]));
//		}
//		free_and_null((char**) &partition_x_ub);
//	}
//	if (partition_ybar_lb!=NULL)
//	{
//		for (i=0;i<partition_cnt;i++)
//		{
//			free_and_null((char**) &(partition_ybar_lb[i]));
//		}
//		free_and_null((char**) &partition_ybar_lb);
//	}
//	if (partition_ybar_ub!=NULL)
//	{
//		for (i=0;i<partition_cnt;i++)
//		{
//			free_and_null((char**) &partition_ybar_ub[i]);
//		}
//		free_and_null((char**) &partition_ybar_ub);
//	}
//	if(y_star!=NULL)
//	{
//		for (i=0;i<refinement_cnt;i++)
//		{
//			free_and_null((char**) &(y_star[i]));
//		}
//		free_and_null((char**) &y_star);
//	}
//	return(status);
//}
int get_refined_VAR_bound(CPXENVptr env, 
						  CPXLPptr lp, 									 
						  const int param_n,								 
						  const int param_m, 								 
						  const int param_k,	
						  const int var_index,
						  const int maxormin,
						  const int param_p,
						  double *var_bound)
{
	int status =0;
	int i,j;
	int lp_solnstat;
	int tmp_index;
	int violation_index;
	int divided_value;
	int* select_indice=NULL;
	char lu;
	double bd;
	double* lp_soln=NULL;
	double* select_violation=NULL;
	double* initial_bd=NULL;
	double tmp_violation;
	double violation;
	double lp_obj;
	double max_bound;
	const double one=1;
	const double zero=0;
	const double negative_one=-1;
	CPXLPptr copy_lp = NULL;
	STARTINFO lp_start_1;
	STARTINFO lp_start_2;
	const int num_var=CPXgetnumcols(env,lp);
	lp_soln=(double*) malloc(num_var*sizeof(double));
	if (lp_soln==NULL)
	{
		status=-1;
		fprintf(stderr," get_refined_VAR_bound(): unable to allocate memory\n");
		goto TERMINATE;
	}
	if (param_p>=1)
	{
		select_indice=(int*) calloc(param_p,sizeof(int));
		select_violation=(double*) malloc(param_p*sizeof(double));
		initial_bd=(double*) malloc(param_p*sizeof(double));
		if (select_indice==NULL || select_violation==NULL || initial_bd==NULL)
		{
			status=-1;
			fprintf(stderr," get_refined_VAR_bound(): unable to allocate memory\n");
			goto TERMINATE;
		}
		for (i=0;i<param_p;i++)
		{
			select_violation[i]=-CPX_INFBOUND;
		}
	}

	init_startinfo(&lp_start_1);
	init_startinfo(&lp_start_2);
	copy_lp = CPXcloneprob(env,lp,&status);
	if (status) 
	{
		status=-1;
		fprintf(stderr," get_refined_VAR_bound(): unable to copy lp \n");
		goto TERMINATE;
	}
	if (maxormin==CPX_MAX || maxormin==CPX_MIN) CPXchgobjsen(env,copy_lp,maxormin);
	else
	{
		status=-1;
		fprintf(stderr," get_refined_VAR_bound(): unexpected maxormin:%d\n",maxormin);
		goto TERMINATE;
	}
	status=zeroobjectivefunc(env,copy_lp);
	if(status) goto TERMINATE;
	status=CPXchgobj(env,copy_lp,1,&var_index,&one);
	if(status) goto TERMINATE;
	status=cplexSolveLP(env,copy_lp,lp_start_1,num_var,lp_soln,&lp_obj,&lp_solnstat,&lp_start_2,CPX_ALG_AUTOMATIC,1,2);
	if(status) goto TERMINATE;
	if (lp_solnstat==CPX_STAT_OPTIMAL)
	{
		*var_bound=lp_obj;

		if (param_p>=1)
		{
			//select the most param_p violation
			for (i=0;i<param_m;i++)
			{
				violation=lp_soln[param_n+i]*lp_soln[param_n+param_m+i];
				violation_index=i;
				//fprintf(stdout,"%d: %f\n",i,violation);
				for (j=0;j<param_p;j++)
				{
					if(violation>select_violation[j])
					{
						tmp_violation=select_violation[j];
						select_violation[j]=violation;
						violation=tmp_violation;
						tmp_index=select_indice[j];
						select_indice[j]=violation_index;
						violation_index=tmp_index;
					}
				}
			}
			//fprintf(stdout," biggest %d violation:\n",param_p);
			//for (i=0;i<param_p;i++)
			//{
			//	fprintf(stdout," %d: %f\n",select_indice[i],select_violation[i]);
			//}
			if (maxormin==CPX_MAX)
			{
				max_bound=-CPX_INFBOUND;
			}else 
			{
				max_bound=CPX_INFBOUND;
			}
			for (i=0;i<pow(2,param_p);i++)
			{
				divided_value=i;
				//fprintf(stdout,"%d: ",i);
				for (j=0;j<param_p;j++)
				{
					if ((divided_value%2)==0)
					{
						lu='U';
						bd=0.0;
						tmp_index=param_n+select_indice[j];
						status=CPXgetub(env,copy_lp,&(initial_bd[j]),tmp_index,tmp_index);
						if(status) goto TERMINATE;
						status=CPXtightenbds(env,copy_lp,1,&tmp_index,&lu,&bd);
						if(status) goto TERMINATE;
					}else
					{
						lu='U';
						bd=0.0;
						tmp_index=param_n+param_m+select_indice[j];
						status=CPXgetub(env,copy_lp,&(initial_bd[j]),tmp_index,tmp_index);
						if(status) goto TERMINATE;
						status=CPXtightenbds(env,copy_lp,1,&tmp_index,&lu,&bd);
						if(status) goto TERMINATE;
					}
					//fprintf(stdout," %d",divided_value%2);
					divided_value/=2;
				}
				//fprintf(stdout,"\n");
				status=cplexSolveLP(env,copy_lp,lp_start_1,num_var,lp_soln,&lp_obj,&lp_solnstat,&lp_start_2,CPX_ALG_AUTOMATIC,1,2);
				if(status) goto TERMINATE;
				if (lp_solnstat==CPX_STAT_OPTIMAL)
				{
					status=copy_startinfo(&lp_start_2,&lp_start_1);
					if(status) goto TERMINATE;
					if (maxormin==CPX_MAX)
					{
						if (max_bound<lp_obj)
						{
							max_bound=lp_obj;
						}
					}else 
					{
						if (max_bound>lp_obj)
						{
							max_bound=lp_obj;
						}
					}

				}else if (lp_solnstat==CPX_STAT_INFEASIBLE)
				{

				}else
				{
					status=-1;
					fprintf(stderr," compute_yw_UB(): unexpected solnstat: %d",lp_solnstat);
					goto TERMINATE;
				}
				//recover lp
				divided_value=i;
				//fprintf(stdout,"%d: ",i);
				for (j=0;j<param_p;j++)
				{
					if ((divided_value%2)==0)
					{
						lu='U';
						bd=initial_bd[j];
						tmp_index=param_n+select_indice[j];
						status=CPXtightenbds(env,copy_lp,1,&tmp_index,&lu,&bd);
					}else
					{
						lu='U';
						bd=initial_bd[j];
						tmp_index=param_n+param_m+select_indice[j];
						status=CPXtightenbds(env,copy_lp,1,&tmp_index,&lu,&bd);
					}
					//fprintf(stdout," %d",divided_value%2);
					divided_value/=2;
				}
			}
			if (maxormin==CPX_MAX)
			{
				if (max_bound>-CPX_INFBOUND)
				{
					//fprintf(stdout," %f, imroved bd: %f\n",*UB,max_UB);
					*var_bound=max_bound;
				}else
				{
					*var_bound=-CPX_INFBOUND;
					//fprintf(stdout," %f\n",*UB);
				}
			}else 
			{
				if (max_bound<CPX_INFBOUND)
				{
					//fprintf(stdout," %f, imroved bd: %f\n",*UB,max_UB);
					*var_bound=max_bound;
				}else
				{
					*var_bound=CPX_INFBOUND;
					//fprintf(stdout," %f\n",*UB);
				}
			}
		}
	}else if (lp_solnstat==CPX_STAT_UNBOUNDED)
	{
		if (maxormin==CPX_MAX)
		{
			*var_bound=CPX_INFBOUND;
		}else 
		{
			*var_bound=-CPX_INFBOUND;
		}
		if (param_p>=1)
		{
			//select the most param_p violation
			for (i=0;i<param_m;i++)
			{
				violation=lp_soln[param_n+i]*lp_soln[param_n+param_m+i];
				violation_index=i;
				//fprintf(stdout,"%d: %f\n",i,violation);
				for (j=0;j<param_p;j++)
				{
					if(violation>select_violation[j])
					{
						tmp_violation=select_violation[j];
						select_violation[j]=violation;
						violation=tmp_violation;
						tmp_index=select_indice[j];
						select_indice[j]=violation_index;
						violation_index=tmp_index;
					}
				}
			}
			//fprintf(stdout," biggest %d violation:\n",param_p);
			//for (i=0;i<param_p;i++)
			//{
			//	fprintf(stdout," %d: %f\n",select_indice[i],select_violation[i]);
			//}
			if (maxormin==CPX_MAX)
			{
				max_bound=-CPX_INFBOUND;
			}else 
			{
				max_bound=CPX_INFBOUND;
			}
			for (i=0;i<pow(2,param_p);i++)
			{
				divided_value=i;
				//fprintf(stdout,"%d: ",i);
				for (j=0;j<param_p;j++)
				{
					if ((divided_value%2)==0)
					{
						lu='U';
						bd=0.0;
						tmp_index=param_n+select_indice[j];
						status=CPXgetub(env,copy_lp,&(initial_bd[j]),tmp_index,tmp_index);
						if(status) goto TERMINATE;
						status=CPXtightenbds(env,copy_lp,1,&tmp_index,&lu,&bd);
						if(status) goto TERMINATE;
					}else
					{
						lu='U';
						bd=0.0;
						tmp_index=param_n+param_m+select_indice[j];
						status=CPXgetub(env,copy_lp,&(initial_bd[j]),tmp_index,tmp_index);
						if(status) goto TERMINATE;
						status=CPXtightenbds(env,copy_lp,1,&tmp_index,&lu,&bd);
						if(status) goto TERMINATE;
					}
					//fprintf(stdout," %d",divided_value%2);
					divided_value/=2;
				}
				//fprintf(stdout,"\n");
				status=cplexSolveLP(env,copy_lp,lp_start_1,num_var,lp_soln,&lp_obj,&lp_solnstat,&lp_start_2,CPX_ALG_AUTOMATIC,1,2);
				if(status) goto TERMINATE;
				if (lp_solnstat==CPX_STAT_OPTIMAL)
				{
					status=copy_startinfo(&lp_start_2,&lp_start_1);
					if(status) goto TERMINATE;
					if (maxormin==CPX_MAX)
					{
						if (max_bound<lp_obj)
						{
							max_bound=lp_obj;
						}
					}else 
					{
						if (max_bound>lp_obj)
						{
							max_bound=lp_obj;
						}
					}

				}else if (lp_solnstat==CPX_STAT_INFEASIBLE)
				{

				}else if (lp_solnstat==CPX_STAT_UNBOUNDED)
				{
					if (maxormin==CPX_MAX)
					{
						max_bound=CPX_INFBOUND;
					}else 
					{
						max_bound=-CPX_INFBOUND;
					}
				}
				else
				{
					status=-1;
					fprintf(stderr," compute_yw_UB(): unexpected solnstat: %d",lp_solnstat);
					goto TERMINATE;
				}
				//recover lp
				divided_value=i;
				//fprintf(stdout,"%d: ",i);
				for (j=0;j<param_p;j++)
				{
					if ((divided_value%2)==0)
					{
						lu='U';
						bd=initial_bd[j];
						tmp_index=param_n+select_indice[j];
						status=CPXtightenbds(env,copy_lp,1,&tmp_index,&lu,&bd);
					}else
					{
						lu='U';
						bd=initial_bd[j];
						tmp_index=param_n+param_m+select_indice[j];
						status=CPXtightenbds(env,copy_lp,1,&tmp_index,&lu,&bd);
					}
					//fprintf(stdout," %d",divided_value%2);
					divided_value/=2;
				}
			}
			if (maxormin==CPX_MAX)
			{
				if (max_bound>-CPX_INFBOUND)
				{
					//fprintf(stdout," %f, imroved bd: %f\n",*UB,max_UB);
					*var_bound=max_bound;
				}else
				{
					*var_bound=-CPX_INFBOUND;
					//fprintf(stdout," %f\n",*UB);
				}
			}else 
			{
				if (max_bound<CPX_INFBOUND)
				{
					//fprintf(stdout," %f, imroved bd: %f\n",*UB,max_UB);
					*var_bound=max_bound;
				}else
				{
					*var_bound=CPX_INFBOUND;
					//fprintf(stdout," %f\n",*UB);
				}
			}
		}
	}else
	{
		status=-1;
		fprintf(stderr," compute_yw_UB(): unexpected solnstat: %d",lp_solnstat);
		goto TERMINATE;
	}
TERMINATE:
	free_and_null((char**)&select_indice);
	free_and_null((char**)&lp_soln);
	free_and_null((char**)&select_violation);
	free_and_null((char**)&initial_bd);
	free_startinfo(&lp_start_1);
	free_startinfo(&lp_start_2);
	if ( copy_lp != NULL ) CPXfreeprob (env, &copy_lp);
	return (status);
}

int clearcutpool(CPXENVptr env,
				 CPXLPptr lp,
				 int start_row,
				 int start_col,
				 int cut_cnt)
{
	int status=0;
	int i;
	const int lp_n_col=CPXgetnumcols(env,lp);
	const int lp_n_row=CPXgetnumrows(env,lp);
	double *x_soln=NULL;
	int *delcolset=NULL;
	int *delrowset=NULL;
	int rm_cut_cnt=0;
	x_soln=(double*)malloc(lp_n_col*sizeof(double));
	delcolset =(int*) calloc(lp_n_col,sizeof(int));
	delrowset =(int*) calloc(lp_n_row,sizeof(int));
	if (x_soln==NULL || delcolset==NULL || delrowset==NULL )
	{
		status=-1;
		fprintf(stderr,"clearcutpool(): unable to allocate memory\n");
		goto TERMINATE;
	}
	status=CPXgetx(env,lp,x_soln,0,lp_n_col-1);
	if(status) goto TERMINATE;
	for (i=0;i<cut_cnt;i++)
	{
		if (x_soln[start_col+i]>1e-6)
		{
			delcolset[start_col+i]=1;
			delrowset[start_row+i]=1;
			rm_cut_cnt++;
		}
	}
	status=CPXdelsetcols(env,lp,delcolset);
	if(status)goto TERMINATE;
	status=CPXdelsetrows(env,lp,delrowset);
	if(status)goto TERMINATE;
	//for (i=0;i<cut_cnt;i++)
	//{
	//	fprintf(stdout,"%d\n",delrowset[start_row+i]);
	//}
	//fprintf(stdout,"total remove cuts: %d\n",rm_cut_cnt);
	//fprintf(stdout,"initial lp col: %d, row: %d\n",lp_n_col,lp_n_row);
	//fprintf(stdout,"final lp col: %d, row: %d\n",CPXgetnumcols(env,lp),CPXgetnumrows(env,lp));
TERMINATE:
	free_and_null((char**)&x_soln);
	free_and_null((char**)&delcolset);
	free_and_null((char**)&delrowset);
	return(status);
}

int ext_bound_cuts_generator(CPXENVptr env, 
							 CPXLPptr rx_lp, 
							 const int param_n,
							 const int param_m,
							 const int param_k,
							 const int violate_complementary_index,
							 const double ub,
							 const int direction,
							 CONSTRAINT_SET *cut_ptr)
{
	int status=0;
	const int num_cols=CPXgetnumcols(env,rx_lp);
	const int num_rows=CPXgetnumrows(env,rx_lp);
	int i;
	int row_id;
	int *ind = NULL;
	int intzero=0;
	int nnz=0;
	double one=1.0;
	double *val = NULL;
	double *tableau=NULL;
	double violation_val;
	//double tmp_lb;
	int index_1;
	int index_2;
	switch(direction)
	{
	case 1:
		index_1=param_n+violate_complementary_index;
		index_2=param_n+param_m+violate_complementary_index;
		break;
	case 2:
		index_1=param_n+param_m+violate_complementary_index;
		index_2=param_n+violate_complementary_index;
		break;
	default:
		status=-1;
		goto TERMINATE;
	}
	status=CPXgetx(env,rx_lp,&violation_val,index_1,index_1);
	if(status) goto TERMINATE;
	//fprintf(stdout,"v: %f\n",violation_val);
	tableau = (double*) malloc(num_cols*sizeof(double));
	ind = (int*) malloc(num_cols*sizeof(int));
	val = (double*) malloc(num_cols*sizeof(double));
	if (tableau==NULL || ind ==NULL || val== NULL)
	{
		status=NO_MEMORY;
		fprintf(stderr," ext_bound_cuts_generator():could not allocate memory.\n");
		goto TERMINATE;
	}
	status = CPXgetijrow(env,rx_lp,-1,index_1,&row_id);
	if (status) goto TERMINATE;
	status = CPXbinvarow(env,rx_lp,row_id,tableau);
	if (status) goto TERMINATE;
	for (i=0;i<num_cols;i++)
	{
		if (i==index_1) 
		{
			//fprintf(stdout,"%f\n",tableau[i]);
			continue;
		}
		if (i==index_2)
		{
			ind[nnz]=i;
			val[nnz]=-(violation_val/ub);
			//fprintf(stdout,"%f\n",tableau[i]);
			nnz++;
		}else if (tableau[i]>0)
		{
			ind[nnz]=i;
			val[nnz]=tableau[i];
			nnz++;	
		}
	}
	//fprintf(stdout,"cut coef:\n");
	//for (i=0;i<nnz;i++)
	//{
	//	fprintf(stdout," %f\n",val[i]);
	//}
	status=add_cuts(cut_ptr,num_cols,nnz,ind,val,-0.001);
	if(status) goto TERMINATE;
TERMINATE:
	free_and_null((char**)&tableau);
	free_and_null((char**)&ind);
	free_and_null((char**)&val);
	return(status);
}

int partition_McCormickRefinement(CPXENVptr env,
								 CPXLPptr lp,
								 int param_n,
								 int param_m,
								 int param_k,
								 int *subproblem_stat,
								 double *x_lb,
								 double *x_ub,
								 double *ybar_lb,
								 double *ybar_ub,
								 double **ystar,
								 int refinement_cnt,
								 int y_bar_start_index,
								 int mccboundconstraint_position_row,
								 int mccboundconstraint_position_col)
{
	int status=0;
	CPXLPptr copy_lp = NULL;
	int i,j;
	int tmp_index;
	int unbounded_flag;
	int condition;
	int solnstat;
	int err_cnt;
	int cnt;
	double pre_violation_level;
	double violation_level;
	double tmplb;
	double tmpub;
	double obj;
	double rhs[4];
    char lu;
	int row_indice[4],col_indice[4];
	*subproblem_stat=1;
	*ystar=NULL;
	copy_lp = CPXcloneprob(env,lp,&status);
	if (status) goto TERMINATE;
	//fprintf(stdout,"partition refinement:\n");
	//fprintf(stdout,"initial bound:\n");
	//for (i=0;i<param_n;i++)
	//{
	//	fprintf(stdout,"x[%d]: lb: %f, ub: %f\n",i,x_lb[i],x_ub[i]);
	//	fprintf(stdout,"y_bar[%d]: lb: %f, ub: %f\n",i,ybar_lb[i],ybar_ub[i]);
	//}
	for (i=0;i<param_n;i++)
	{
		lu='L';
		status=CPXchgbds(env,copy_lp,1,&i,&lu,&(x_lb[i]));
		if(status) goto TERMINATE;
		lu='U';
		status=CPXchgbds(env,copy_lp,1,&i,&lu,&(x_ub[i]));
		if(status) goto TERMINATE;
		tmp_index=y_bar_start_index+i;
		lu='L';
		status=CPXchgbds(env,copy_lp,1,&tmp_index,&lu,&(ybar_lb[i]));
		if(status) goto TERMINATE;
		lu='U';
		status=CPXchgbds(env,copy_lp,1,&tmp_index,&lu,&(ybar_ub[i]));
		if(status) goto TERMINATE;
	}
	pre_violation_level=CPX_INFBOUND;
	for (i=0;i<refinement_cnt;i++)
	{
		err_cnt=0;
		cnt=0;
		for (j=0;j<param_n;j++)
		{
			if (fabs(x_ub[j]-x_lb[j])>ZERO_TOLERANCE)
			{
				status=getVARbound_baropt(env,copy_lp,j,CPX_MIN,&tmplb,&condition);
				if(status)goto TERMINATE;
				if(x_lb[j]<tmplb && condition==1)x_lb[j]=tmplb;
				if (condition==-1)
				{
					*subproblem_stat=-1;
					goto TERMINATE;
				}
				if (condition==0)
				{
					err_cnt++;
				}
				cnt++;
				status=getVARbound_baropt(env,copy_lp,j,CPX_MAX,&tmpub,&condition);
				if(status)goto TERMINATE;
				if(x_ub[j]>tmpub)x_ub[j]=tmpub;
				if (condition==-1)
				{
					*subproblem_stat=-1;
					goto TERMINATE;
				}
				if (condition==0)
				{
					err_cnt++;
				}
				cnt++;
			}
			if (fabs(ybar_ub[j]-ybar_lb[j])>ZERO_TOLERANCE)
			{
				status=getVARbound_baropt(env,copy_lp,y_bar_start_index+j,CPX_MIN,&tmplb,&condition);
				if(status)goto TERMINATE;
				if(ybar_lb[j]<tmplb)ybar_lb[j]=tmplb;
				if (condition==-1)
				{
					*subproblem_stat=-1;
					goto TERMINATE;
				}
				if (condition==0)
				{
					err_cnt++;
				}
				cnt++;
				status=getVARbound_baropt(env,copy_lp,y_bar_start_index+j,CPX_MAX,&tmpub,&condition);
				if(status)goto TERMINATE;
				if(ybar_ub[j]>tmpub)ybar_ub[j]=tmpub;
				if (condition==-1)
				{
					*subproblem_stat=-1;
					goto TERMINATE;
				}
				if (condition==0)
				{
					err_cnt++;
				}
				cnt++;
			}
		}


		//modify the coefficent in copy_lp
		for (j=0;j<param_n;j++)
		{				
			row_indice[0]=mccboundconstraint_position_row+4*j;
			row_indice[1]=mccboundconstraint_position_row+4*j+1;
			row_indice[2]=mccboundconstraint_position_row+4*j+2;
			row_indice[3]=mccboundconstraint_position_row+4*j+3;
			col_indice[0]=mccboundconstraint_position_col+4*j;
			col_indice[1]=mccboundconstraint_position_col+4*j+1;
			col_indice[2]=mccboundconstraint_position_col+4*j+2;
			col_indice[3]=mccboundconstraint_position_col+4*j+3;
			if (x_lb[j]>-CPX_INFBOUND && ybar_lb[j]>-CPX_INFBOUND)
			{
				rhs[0]=x_lb[j]*ybar_lb[j];
				status=CPXchgcoef(env,copy_lp,row_indice[0],j,ybar_lb[j]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,copy_lp,row_indice[0],y_bar_start_index+j,x_lb[j]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,copy_lp,row_indice[0],col_indice[0],1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,copy_lp,row_indice[0],y_bar_start_index+param_n+j,-1);
				if (status) goto TERMINATE;
			}else
			{
				rhs[0]=0;
			}
			if (x_ub[j]<CPX_INFBOUND && ybar_ub[j]<CPX_INFBOUND)
			{
				rhs[1]=x_ub[j]*ybar_ub[j];
				status=CPXchgcoef(env,copy_lp,row_indice[1],j,ybar_ub[j]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,copy_lp,row_indice[1],y_bar_start_index+j,x_ub[j]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,copy_lp,row_indice[1],col_indice[1],1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,copy_lp,row_indice[1],y_bar_start_index+param_n+j,-1);
				if (status) goto TERMINATE;
			}else
			{
				rhs[1]=0;
			}
			if (x_lb[j]>-CPX_INFBOUND && ybar_ub[j]<CPX_INFBOUND)
			{
				rhs[2]=x_lb[j]*ybar_ub[j];
				status=CPXchgcoef(env,copy_lp,row_indice[2],j,ybar_ub[j]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,copy_lp,row_indice[2],y_bar_start_index+j,x_lb[j]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,copy_lp,row_indice[2],col_indice[2],-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,copy_lp,row_indice[2],y_bar_start_index+param_n+j,-1);
				if (status) goto TERMINATE;
			}else
			{
				rhs[2]=0;
			}

			if (x_ub[j]<CPX_INFBOUND && ybar_lb[j]>-CPX_INFBOUND)
			{
				rhs[3]=x_ub[j]*ybar_lb[j];
				status=CPXchgcoef(env,copy_lp,row_indice[3],j,ybar_lb[j]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,copy_lp,row_indice[3],y_bar_start_index+j,x_ub[j]);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,copy_lp,row_indice[3],col_indice[3],-1);
				if (status) goto TERMINATE;
				status=CPXchgcoef(env,copy_lp,row_indice[3],y_bar_start_index+param_n+j,-1);
				if (status) goto TERMINATE;
			}else
			{
				rhs[3]=0;
			}
			status=CPXchgrhs(env,copy_lp,4,row_indice,rhs);
			if(status) goto TERMINATE;
		}
		unbounded_flag=1;
		violation_level=0;
		for (j=0;j<param_n;j++)
		{
			//fprintf(stdout,"x[%d]: ub:%f, lb:%f\n",j,x_ub[j],x_lb[j]);
			//fprintf(stdout,"ybar[%d]: ub:%f, lb:%f\n",j,ybar_ub[j],ybar_lb[j]);
			if (x_lb[j]>-CPX_INFBOUND && x_ub[j]<CPX_INFBOUND && ybar_lb[j] >-CPX_INFBOUND && ybar_ub[j]<CPX_INFBOUND)
			{
				violation_level+=0.5*fabs(x_ub[j]-x_lb[j])*fabs(ybar_ub[j]-ybar_lb[j]);
			}else unbounded_flag=0;
		}
		if (unbounded_flag)
		{
			//fprintf(stdout,"violation level:%f\n",violation_level);
			//if (fabs(pre_violation_level-violation_level)<0.5) break;
			pre_violation_level=violation_level;
		}
		if (cnt==err_cnt) break;
	}
	//fprintf(stdout,"final bounds:\n");
	//for (j=0;j<param_n;j++)
	//{
	//	fprintf(stdout,"x[%d]: ub:%f, lb:%f\n",j,x_ub[j],x_lb[j]);
	//	fprintf(stdout,"ybar[%d]: ub:%f, lb:%f\n",j,ybar_ub[j],ybar_lb[j]);
	//}
	status = CPXbaropt (env, copy_lp);
	if ( status ) goto TERMINATE;
	//get y_star
	(*ystar)=(double*)malloc(param_m*sizeof(double));
	if ((*ystar)==NULL)
	{
		status=-1;
		fprintf(stderr," unable to allocate memory\n");
		goto TERMINATE;
	}
	status = CPXsolution (env, copy_lp, &solnstat, &obj, NULL, NULL, NULL, NULL);	
	if ( status ) goto TERMINATE;
	//fprintf(stdout,"%d; %f\n",solnstat,obj);
	status = CPXgetx(env,copy_lp,(*ystar),param_n,param_n+param_m-1);
	if(status) goto TERMINATE;
	if (solnstat==CPX_STAT_INFEASIBLE)
	{
		*subproblem_stat=-1;
	}
TERMINATE:
	if(copy_lp!=NULL) CPXfreeprob(env,&copy_lp);
	return(status);
}

