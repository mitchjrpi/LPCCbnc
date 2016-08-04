#define _CRT_SECURE_NO_DEPRECATE
#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "fcgimisc.h"
#include "LPCC.h"
#include "LPCCsolver.h"
#include "utilities.h"
#include "instancegenerator.h"
#include "CPLEXcontrol.h"
#include "cutsgenerator.h"
#include <time.h>

extern int nodecnt;
extern int leftnodes;

int solveLPCC (PARAM parameter)
{
	int status = 0;
	int param_n;
	int param_m;
	int param_k;
	int    i;
	int LPCC_status;
	double LPCC_lb=-INF_BOUND;
	double LPCC_ub=INF_BOUND;
	double *LPCC_soln =NULL;
	double *c_coef	=NULL;
	double *d_coef=NULL;
	double *b_coef=NULL;
	double *q_coef=NULL;
	double **matrix_A= NULL;
	double **matrix_B=NULL;
	double **matrix_N=NULL;
	double **matrix_M=NULL;
	MATRIX _matrix_A;
	MATRIX _matrix_B;
	MATRIX _matrix_N;
	MATRIX _matrix_M;
	init_matrix(&_matrix_A);
	init_matrix(&_matrix_B);
	init_matrix(&_matrix_N);
	init_matrix(&_matrix_M);
	fprintf(stdout,"Reading inputdata...\n");
	switch (parameter.format_mode)
	{
	case 1:
		//fprintf(stdout,"full matrix (separate files) format\n");
		status = readdata_format_1(parameter.inputfile,&param_n,&param_m,&param_k,&c_coef,&d_coef,&b_coef,&q_coef,&matrix_A,&matrix_B,&matrix_N,&matrix_M);
		if ( status ) goto TERMINATE;
		status=get_matrix(matrix_A,param_k,param_n,1,&_matrix_A);
		if (status) goto TERMINATE;
		status=get_matrix(matrix_B,param_k,param_m,1,&_matrix_B);
		if (status) goto TERMINATE;
		status=get_matrix(matrix_N,param_m,param_n,1,&_matrix_N);
		if (status) goto TERMINATE;
		status=get_matrix(matrix_M,param_m,param_m,1,&_matrix_M);
		if (status) goto TERMINATE;
		break;
	case 2:
		//fprintf(stdout,"full matrix (single file) format\n");
		status = readdata_format_2(parameter.inputfile,&param_n,&param_m,&param_k,&c_coef,&d_coef,&b_coef,&q_coef,&matrix_A,&matrix_B,&matrix_N,&matrix_M);
		if ( status ) goto TERMINATE;
		status=get_matrix(matrix_A,param_k,param_n,1,&_matrix_A);
		if (status) goto TERMINATE;
		status=get_matrix(matrix_B,param_k,param_m,1,&_matrix_B);
		if (status) goto TERMINATE;
		status=get_matrix(matrix_N,param_m,param_n,1,&_matrix_N);
		if (status) goto TERMINATE;
		status=get_matrix(matrix_M,param_m,param_m,1,&_matrix_M);
		if (status) goto TERMINATE;
		break;
	case 3:
		//fprintf(stdout,"compact matrix format\n");
		status = readdata_format_3(parameter.inputfile,&param_n,&param_m,&param_k,&c_coef,&d_coef,&b_coef,&q_coef,&_matrix_A,&_matrix_B,&_matrix_N,&_matrix_M);
		if ( status ) goto TERMINATE;
		break;
	case 4:
		//fprintf(stdout,"BIP .mps format\n");
		status = convertMPS(parameter.inputfile,
			&param_n,
			&param_m,
			&param_k,
			&c_coef,
			&d_coef,
			&b_coef,
			&q_coef,
			&matrix_A,
			&matrix_B,
			&matrix_N,
			&matrix_M);
		if ( status ) goto TERMINATE;
		status=get_matrix(matrix_A,param_k,param_n,1,&_matrix_A);
		if (status) goto TERMINATE;
		status=get_matrix(matrix_B,param_k,param_m,1,&_matrix_B);
		if (status) goto TERMINATE;
		status=get_matrix(matrix_N,param_m,param_n,1,&_matrix_N);
		if (status) goto TERMINATE;
		status=get_matrix(matrix_M,param_m,param_m,1,&_matrix_M);
		if (status) goto TERMINATE;
		break;
	default:
		status=-1;
		fprintf(stderr,"Unknown format, please check the input file again\n");
		goto TERMINATE;
	}

	status = LPCCSolver(param_n,
						param_m,
						param_k,
						c_coef,
						d_coef,
						b_coef,
						q_coef,
						_matrix_A,
						_matrix_B,
						_matrix_M,
						_matrix_N,
						parameter,
						&LPCC_soln,
						&LPCC_status,
						&LPCC_lb,
						&LPCC_ub);
	if(status) goto TERMINATE;
TERMINATE:
	/* Free pointer */
	free_matrix(&_matrix_A);
	free_matrix(&_matrix_B);
	free_matrix(&_matrix_N);
	free_matrix(&_matrix_M);
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
	free_and_null ((char **) &LPCC_soln);
	return (status);
}  
int LPCCSolver(const int param_n, 
			   const int param_m,
			   const int param_k,
			   double *c_coef,
			   double *d_coef,
			   double *b_coef,
			   double *q_coef,
			   MATRIX _matrix_A,
			   MATRIX _matrix_B,
			   MATRIX _matrix_M,
			   MATRIX _matrix_N,
			   PARAM parameter,
			   double **LPCC_soln,
			   int	*LPCC_stat,
			   double *LPCC_lb,
			   double *LPCC_ub)
{
	int status = 0;
	double used_time=0;
	STARTINFO init_rx_lp_start;
	STARTINFO updated_rx_lp_start;
	CPXENVptr env = NULL;
	CPXLPptr rx_lp = NULL;
	CPXLPptr cg_lp = NULL;
	NODE *headnodeptr=NULL;

	int solnstat;
	int LPCC_precondition;
	double init_lb=-INF_BOUND;
	double best_lb=-INF_BOUND;
	double improved_lb=-INF_BOUND;	
	clock_t solve_start,solve_end,preprocess_start,preprocess_end;
	time_t solve_startt,solve_endt,preprocess_startt,preprocess_endt;
        double proc_time_wall,solve_time_wall;
	*LPCC_lb=-INF_BOUND;
	*LPCC_ub=INF_BOUND;


	//initialize the global variable
	nodecnt = 0;
	leftnodes = 0;
	init_startinfo(&init_rx_lp_start);
	init_startinfo(&updated_rx_lp_start);
	*LPCC_soln = (double*) malloc((param_n+param_m+param_m+param_k)*sizeof(double));
	if (*LPCC_soln == NULL) 
	{
		status = NO_MEMORY;	
		fprintf (stderr, "Could not allocate memory.\n");
		goto TERMINATE;
	}
	///* Initialize the CPLEX environment */
	env = CPXopenCPLEX (&status);
	if ( env == NULL ) {
		char  errmsg[1024];
		fprintf (stderr, "Could not open CPLEX environment.\n");
		CPXgeterrorstring (env, status, errmsg);
		fprintf (stderr, "%s", errmsg);
		goto TERMINATE;
	}
	rx_lp = CPXcreateprob (env, &status, "rx_lp");
	if ( rx_lp == NULL ) {
		fprintf (stderr, "Failed to create LP.\n");
		goto TERMINATE;
	}
	cg_lp = CPXcreateprob (env, &status, "cg_lp");
	if ( cg_lp == NULL ) {
		fprintf (stderr, "Failed to create LP.\n");
		goto TERMINATE;
	}
	status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
	if ( status ) {
		fprintf (stderr, 
			"Failure to turn on screen indicator, error %d.\n", status);
		goto TERMINATE;
	}
	status = CPXsetintparam (env, CPX_PARAM_DATACHECK, CPX_ON);
	if ( status ) {
		fprintf (stderr, 
			"Failure to turn on data checking, error %d.\n", status);
		goto TERMINATE;
	}
	//setup problem
	status=setupGeneralLPCC_lp(env, rx_lp,c_coef,d_coef,b_coef,q_coef,_matrix_A,_matrix_B,_matrix_N,_matrix_M,parameter.format_mode);
	if(status) goto TERMINATE;
	// fprintf(stdout," setup initial CGLP...\n");
	status=setup_initial_CGLP(env,cg_lp,param_k,param_m,param_n,b_coef,q_coef,_matrix_A,_matrix_B,_matrix_N,_matrix_M);
	if(status) goto TERMINATE;
	//solve problem
	fprintf(stdout,"Start solving...\n");
	solve_start=clock();
	solve_startt=time(NULL);
	preprocess_start=clock();
	preprocess_startt=time(NULL);
	// printf("%d preprocess_startt\n",(int)preprocess_startt);
	status=preprocess(env,rx_lp,cg_lp,_matrix_N,_matrix_M,q_coef,param_n,param_m,param_k,&headnodeptr,&init_lb,&improved_lb,LPCC_ub,*LPCC_soln,&LPCC_precondition,parameter);
        // printf(" status after preprocess %d\n",status);
	if (status) goto TERMINATE;
	preprocess_end = clock();
	preprocess_endt = time(NULL);
	fprintf(stdout," Total preprocessing time: %f\n",(double)(preprocess_end-preprocess_start)/CLOCKS_PER_SEC);
	proc_time_wall = difftime(preprocess_endt, preprocess_startt);
	fprintf(stdout," Total preprocessing wall clock time: %.0f\n",proc_time_wall);
	switch (LPCC_precondition)
	{
	case LPCC_PRECONDITION_BOUNDED:
		switch(parameter.solve_mode)
		{
		case 'c':
			status=CPXsetdblparam(env,CPX_PARAM_EPGAP,parameter.opt_tolerance);
			if(status) goto TERMINATE;
			status= CPXsetdblparam(env,CPX_PARAM_TILIM,parameter.time_limit);
			if(status) goto TERMINATE;
			status=CPX_LPCCSolver(env,rx_lp,param_n,param_m,param_k,*LPCC_soln);
			if (status) goto TERMINATE;
			solve_end = clock();
			solve_endt = time(NULL);
                        // printf("%s", asctime(localtime(&preprocess_endt)));
                        // printf("%s", asctime(localtime(&solve_endt)));
	                solve_time_wall = difftime(solve_endt, preprocess_endt);
			nodecnt = CPXgetnodecnt(env,rx_lp);
			solnstat=CPXgetstat (env, rx_lp);
			switch(solnstat)
			{
			case CPXMIP_INFEASIBLE:
				*LPCC_stat=LPCC_INFEASIBLE;
				*LPCC_lb=INF_BOUND;
				*LPCC_ub=INF_BOUND;
				break;
			case CPXMIP_INForUNBD:
				*LPCC_stat=LPCC_INForUNBD;
				break;
			case CPXMIP_OPTIMAL_TOL:
			case CPXMIP_OPTIMAL:
				status=CPXgetobjval(env,rx_lp,LPCC_ub);
				if(status) goto TERMINATE;
				status=CPXgetx(env,rx_lp,*LPCC_soln,0,(param_n+param_m+param_m+param_k-1));
				if(status) goto TERMINATE;
				*LPCC_stat=LPCC_OPTIMAL;
				*LPCC_lb=*LPCC_ub;
				break;
			case CPXMIP_TIME_LIM_FEAS:
				status=CPXgetobjval(env,rx_lp,LPCC_ub);
				if(status) goto TERMINATE;
				status=CPXgetx(env,rx_lp,*LPCC_soln,0,(param_n+param_m+param_m+param_k-1));
				if(status) goto TERMINATE;
				status=CPXgetbestobjval(env,rx_lp,&best_lb);
				if(status) goto TERMINATE;
				*LPCC_stat=LPCC_TERMINATE_WITH_SOLN;
				*LPCC_lb=best_lb;
				break;
			case CPXMIP_TIME_LIM_INFEAS:
				status=CPXgetbestobjval(env,rx_lp,&best_lb);
				if(status) goto TERMINATE;
				*LPCC_stat=LPCC_TERMINATE_WITHOUT_SOLN;
				*LPCC_lb=best_lb;
				*LPCC_ub=INF_BOUND;
				break;
            case CPXMIP_UNBOUNDED:
                *LPCC_stat=LPCC_UNBOUNDED;
                *LPCC_lb=-INF_BOUND;
                *LPCC_ub=INF_BOUND;
                break;
            /* case CPX_STAT_OPTIMAL:
                status=CPXgetobjval(env,rx_lp,LPCC_ub);
                if(status) goto TERMINATE;
                status=CPXgetx(env,rx_lp,*LPCC_soln,0,(param_n+param_m+param_m+param_k-1));
                if(status) goto TERMINATE;
                *LPCC_stat=LPCC_OPTIMAL;
                *LPCC_lb=*LPCC_ub;
                break; */
			default:
				fprintf(stderr," unexpected MIP CPLEX solution stauts 1: %d\n",solnstat);
				status=-1;
				goto TERMINATE;
			}
			break;
		case 'b':
			used_time=0;
			status=LPCCSolver_bounded(env,rx_lp,param_n,param_m,param_k,&headnodeptr,LPCC_lb,LPCC_ub,*LPCC_soln,LPCC_stat,parameter,&used_time);
			if (status) goto TERMINATE;
			solve_end = clock();
			solve_endt = time(NULL);
                        // printf("%s", asctime(localtime(&preprocess_endt)));
                        // printf("%s", asctime(localtime(&solve_endt)));
	                solve_time_wall = difftime(solve_endt, preprocess_endt);
			break;
		default:
			status=-1;
			fprintf(stderr,"parameter should be c or b\n");
			goto TERMINATE;
		}
		break;
	case LPCC_PRECONDITION_INFEASIBLE:
		*LPCC_stat=LPCC_INFEASIBLE;
		*LPCC_lb=INF_BOUND;
		*LPCC_ub=INF_BOUND;
		solve_end = clock();
		solve_endt = time(NULL);
                        // printf("%s", asctime(localtime(&preprocess_endt)));
                        // printf("%s", asctime(localtime(&solve_endt)));
	                solve_time_wall = difftime(solve_endt, preprocess_endt);
		break;
	case LPCC_PRECONDITION_OPT:
		*LPCC_lb=*LPCC_ub;
		solve_end = clock();
		solve_endt = time(NULL);
                        // printf("%s", asctime(localtime(&preprocess_endt)));
                        // printf("%s", asctime(localtime(&solve_endt)));
	                solve_time_wall = difftime(solve_endt, preprocess_endt);
		break;
	case LPCC_PRECONDITION_UNBOUNDED:
		switch(parameter.solve_mode)
		{
		case 'c':
			status=CPXsetdblparam(env,CPX_PARAM_EPGAP,parameter.opt_tolerance);
			if(status) goto TERMINATE;
			status= CPXsetdblparam(env,CPX_PARAM_TILIM,parameter.time_limit);
			if(status) goto TERMINATE;
			status=CPX_LPCCSolver(env,rx_lp,param_n,param_m,param_k,*LPCC_soln);
			if (status) goto TERMINATE;
			solve_end = clock();
		        solve_endt = time(NULL);
                        // printf("%s", asctime(localtime(&preprocess_endt)));
                        // printf("%s", asctime(localtime(&solve_endt)));
	                solve_time_wall = difftime(solve_endt, preprocess_endt);
			nodecnt = CPXgetnodecnt(env,rx_lp);
			solnstat=CPXgetstat (env, rx_lp);
			switch(solnstat)
			{
			case CPXMIP_INFEASIBLE:
				*LPCC_stat=LPCC_INFEASIBLE;
				break;
			case CPXMIP_INForUNBD:
				*LPCC_stat=LPCC_INForUNBD;
				break;
			case CPXMIP_OPTIMAL_TOL:
			case CPXMIP_OPTIMAL:
				status=CPXgetobjval(env,rx_lp,LPCC_ub);
				if(status) goto TERMINATE;
				status=CPXgetx(env,rx_lp,*LPCC_soln,0,(param_n+param_m+param_m+param_k-1));
				if(status) goto TERMINATE;
				*LPCC_stat=LPCC_OPTIMAL;
				*LPCC_lb=*LPCC_ub;
				break;
			case CPXMIP_TIME_LIM_FEAS:
				status=CPXgetobjval(env,rx_lp,LPCC_ub);
				if(status) goto TERMINATE;
				status=CPXgetx(env,rx_lp,*LPCC_soln,0,(param_n+param_m+param_m+param_k-1));
				if(status) goto TERMINATE;
				status=CPXgetbestobjval(env,rx_lp,&best_lb);
				if(status) goto TERMINATE;
				*LPCC_stat=LPCC_TERMINATE_WITH_SOLN;
				*LPCC_lb=best_lb;
				break;
			case CPXMIP_TIME_LIM_INFEAS:
				status=CPXgetbestobjval(env,rx_lp,&best_lb);
				if(status) goto TERMINATE;
				*LPCC_stat=LPCC_TERMINATE_WITHOUT_SOLN;
				*LPCC_lb=best_lb;
				*LPCC_ub=INF_BOUND;
				break;
			case CPXMIP_UNBOUNDED:
				*LPCC_stat=LPCC_UNBOUNDED;
				*LPCC_lb=-INF_BOUND;
				*LPCC_ub=INF_BOUND;
				break;
            /* case CPX_STAT_OPTIMAL:
                status=CPXgetobjval(env,rx_lp,LPCC_ub);
                if(status) goto TERMINATE;
                status=CPXgetx(env,rx_lp,*LPCC_soln,0,(param_n+param_m+param_m+param_k-1));
                if(status) goto TERMINATE;
                *LPCC_stat=LPCC_OPTIMAL;
                *LPCC_lb=*LPCC_ub;
                break; */
            default:
				fprintf(stderr," unexpected MIP CPLEX solution stauts 2: %d\n",solnstat);
				status=-1;
				goto TERMINATE;
			}
			break;
		case 'b':
			used_time=0;
			status=LPCCSolver_unbounded(env,rx_lp,param_n,param_m,param_k,&headnodeptr,LPCC_lb,LPCC_ub,*LPCC_soln,LPCC_stat,parameter,&used_time);
			if (status) goto TERMINATE;
			solve_end = clock();
		        solve_endt = time(NULL);
                        // printf("%s", asctime(localtime(&preprocess_endt)));
                        // printf("%s", asctime(localtime(&solve_endt)));
	                solve_time_wall = difftime(solve_endt, preprocess_endt);
			break;
		default:
			status=-1;
			fprintf(stderr,"parameter should be c or b\n");
			goto TERMINATE;
		}
		break;
	}
	printResult(*LPCC_stat,param_n,param_m,param_k,nodecnt,parameter,(double)(preprocess_end-preprocess_start)/CLOCKS_PER_SEC,(double)(solve_end-preprocess_end)/CLOCKS_PER_SEC,proc_time_wall,solve_time_wall,init_lb,improved_lb,*LPCC_lb,*LPCC_ub);
	writeSolution(*LPCC_stat,param_n,param_m,param_k,nodecnt,parameter,(double)(preprocess_end-preprocess_start)/CLOCKS_PER_SEC,(double)(solve_end-preprocess_end)/CLOCKS_PER_SEC,proc_time_wall,solve_time_wall,init_lb,improved_lb,*LPCC_lb,*LPCC_ub,*LPCC_soln);
	
TERMINATE:
        if ( status != 0 ) fprintf(stderr," status %i\n",status);
	if ( rx_lp != NULL ) CPXfreeprob (env, &rx_lp);
	if ( cg_lp != NULL ) CPXfreeprob (env, &cg_lp);
	if ( env != NULL ) CPXcloseCPLEX (&env);
	free_startinfo(&init_rx_lp_start);
	free_startinfo(&updated_rx_lp_start);
	return (status);
}

int preprocess(CPXENVptr env, 
			   CPXLPptr rx_lp, 
			   CPXLPptr cg_lp,
			   MATRIX matrix_N,
			   MATRIX matrix_M,
			   double *q_coef,
			   const int param_n,
			   const int param_m, 
			   const int param_k, 
			   NODE **headnodeptr, 
			   double *init_lb, 
			   double *improve_lb,
			   double *ub, 
			   double *feasible_x, 
			   int *LPCC_precondition,
			   PARAM preprocess_param)
{
	int status =0;
	int i,j;
	int lp_solnstat;
	int recovery_stat;
	int recovery_refinement_cnt;
	const int num_var= param_n+param_m+param_m+param_k;
	int mcc_processing_stat=0;
	int initial_subproblem_cnt=1;
	int *initial_subproblem_stat=NULL;
	int x_start_index=-1;
	int ybar_start_index=-1;
	int mcc_start_row=-1;
	int mcc_start_col=-1;
	int initial_problem_cnt=0;
	double **p_x_lb=NULL;
	double **p_x_ub=NULL;
	double **p_y_lb=NULL;
	double **p_y_ub=NULL;
	double initial_range;
	double tmp_range;
	double refinement_range_lb;
	double refinement_range_ub;
	double *lp_soln=NULL;
	double lp_obj;
	STARTINFO rx_lp_start;
	STARTINFO updated_rx_lp_start;
	clock_t t_start,t_end;
	NODE *tmp_node=NULL;
	init_startinfo(&rx_lp_start);
	init_startinfo(&updated_rx_lp_start);
	lp_soln=(double*) malloc(num_var*sizeof(double));
	if(lp_soln==NULL)
	{
		status = NO_MEMORY;	
		fprintf (stderr, " preprocess(): could not allocate memory.\n");
		goto TERMINATE;
	}
	*headnodeptr=(NODE*) malloc(sizeof(NODE));
	if (*headnodeptr==NULL)
	{
		status = NO_MEMORY;
		fprintf(stderr," preprocess(): could not allocate memory.\n");
		goto TERMINATE;
	}
	init_node(*headnodeptr);
	(*headnodeptr)->link = NULL;
	(*headnodeptr)->n=param_n;
	(*headnodeptr)->m=param_m;
	(*headnodeptr)->k=param_k;
	(*headnodeptr)->x_start_index=-1;
	(*headnodeptr)->y_bar_start_index=-1;
	(*headnodeptr)->mcc_index_col=-1;
	(*headnodeptr)->mcc_index_row=-1;
	(*headnodeptr)->x_lb= (double*) malloc(param_n*sizeof(double));
	(*headnodeptr)->x_ub= (double*) malloc(param_n*sizeof(double));
	(*headnodeptr)->y_bar_lb=(double*) malloc(param_n*sizeof(double));
	(*headnodeptr)->y_bar_ub=(double*) malloc(param_n*sizeof(double));
	(*headnodeptr)->nodebranchinfo = (int*) malloc(param_m*sizeof(int));
	(*headnodeptr)->nodex = (double*) malloc(num_var*sizeof(double));
	if ((*headnodeptr)->nodebranchinfo==NULL ||
		(*headnodeptr)->nodex==NULL ||
		(*headnodeptr)->x_lb==NULL ||
		(*headnodeptr)->x_ub==NULL ||
		(*headnodeptr)->y_bar_lb==NULL ||
		(*headnodeptr)->y_bar_ub==NULL)
	{
		status = NO_MEMORY;
		fprintf(stderr," preprocess(): could not allocate memory.\n");
		goto TERMINATE;
	}
	for (j=0;j<param_m;j++)
	{
		(*headnodeptr)->nodebranchinfo[j]=0;
	}
	fprintf(stdout,"Preprocessing...\n");
	//solve LP
	status = cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat, &updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,LP_PRECONDITION_UNKNOWN);//m
	if (status) goto TERMINATE;
	switch(lp_solnstat)
	{
	case CPX_STAT_OPTIMAL:
		*LPCC_precondition=LPCC_PRECONDITION_BOUNDED;
		fprintf(stdout," LPCC relaxation: %f\n",lp_obj);
		*init_lb=lp_obj;
		*improve_lb=lp_obj;	
		if (preprocess_param.feasibility_recovery_breath>0)
		{
			fprintf(stdout," Feasibility recovery...\n");
			t_start=clock();
			status=LPCC_feasible_recovery(env,rx_lp,updated_rx_lp_start,param_n,param_m,param_k,(**headnodeptr),-CPX_INFBOUND,CPX_INFBOUND,num_var,lp_soln,ub,feasible_x,&recovery_stat,(int)ceil(param_m*preprocess_param.feasibility_recovery_breath),preprocess_param.feasibility_recovery_depth);
			if(status) goto TERMINATE;
			if (recovery_stat==1)
			{
				//refine 
				if (preprocess_param.feasibility_recovery_refinement==1)
				{
					fprintf(stdout,"	Feasibility recovery refinement...\n");
					recovery_refinement_cnt=0;
					refinement_range_lb=lp_obj;
					refinement_range_ub=lp_obj+(*ub-lp_obj)/2;
					initial_range=*ub-lp_obj;
					if(status) goto TERMINATE;
					for (;;)
					{
						if ((fabs(refinement_range_lb-refinement_range_ub)/initial_range)<0.01 || (fabs(*ub-lp_obj)/max(1,lp_obj))<0.05 || refinement_range_lb>=refinement_range_ub )break;//
						fprintf(stdout,"		Search range  %f - %f\n",refinement_range_lb,refinement_range_ub);
						status=LPCC_feasible_recovery(env,rx_lp,updated_rx_lp_start,param_n,param_m,param_k,(**headnodeptr),refinement_range_lb,refinement_range_ub,num_var,feasible_x,ub,feasible_x,&recovery_stat,1,preprocess_param.feasibility_recovery_depth);
						if(status) goto TERMINATE;
						if (recovery_stat==1)
						{
							refinement_range_ub=refinement_range_lb+(*ub-refinement_range_lb)/2;
						}else
						{
							tmp_range=refinement_range_ub-refinement_range_lb;
							refinement_range_lb=refinement_range_ub;
							refinement_range_ub=refinement_range_lb+tmp_range/2;
						}
					}
				}
			}else
			{
				for (j=0;j<10;j++)
				{
					status=LPCC_feasible_recover_M_psd_fix_x(env,rx_lp,param_n,param_m,param_k,(param_n+param_m+param_m),matrix_N,matrix_M,q_coef,(**headnodeptr),lp_soln,ub,lp_soln,&recovery_stat);
					if(status) goto TERMINATE;	
					if (recovery_stat==-1)break;
					status=LPCC_feasible_recover_M_psd_fix_y_bar(env,rx_lp,param_n,param_m,param_k,(param_n+param_m+param_m),matrix_N,matrix_M,q_coef,(**headnodeptr),lp_soln,ub,lp_soln,&recovery_stat);
					if(status) goto TERMINATE;
					status=LPCC_feasible_recover_by_rounding(env,rx_lp,rx_lp_start,param_n,param_m,param_k,num_var,(**headnodeptr),lp_soln,ub,feasible_x,&recovery_stat);
					if(status) goto TERMINATE;
					if (recovery_stat==1) break;
				}
			}
			if (recovery_stat==1)
			{
				fprintf(stdout,"	Final find the best ub: %f\n",*ub);
				status=lp_ub_control(env,rx_lp,*ub,0);
				if(status)goto TERMINATE;
			}else
			{
				fprintf(stdout,"	Unable to recover a feasible solution\n");
			}
			t_end=clock();
			fprintf(stdout,"	Total feasibility recovery time: %12.3f\n",(double)(t_end-t_start)/CLOCKS_PER_SEC);
		}
		if (preprocess_param.disjunctivecut_loop_param>0)
		{
			status=add_disjunctivecuts(env,rx_lp,cg_lp,q_coef,matrix_N,matrix_M,&updated_rx_lp_start,param_n,param_m,param_k,preprocess_param.disjunctivecut_loop_param,feasible_x,ub,LPCC_precondition);
			if(status) goto TERMINATE;
		}
		if((*LPCC_precondition)==LPCC_PRECONDITION_OPT || (*LPCC_precondition)==LPCC_PRECONDITION_INFEASIBLE) goto TERMINATE; 
		if (preprocess_param.simplecut_loop_param>0) 
		{
			status=add_simplecuts(env,rx_lp,&updated_rx_lp_start,param_n,param_m,param_k,preprocess_param.simplecut_loop_param,feasible_x,ub,LPCC_precondition);
			if(status)goto TERMINATE;
		}
		if((*LPCC_precondition)==LPCC_PRECONDITION_OPT || (*LPCC_precondition)==LPCC_PRECONDITION_INFEASIBLE) goto TERMINATE; 

		if (preprocess_param.boundcut_p>=0) 
		{
			status=add_boundcuts(env,rx_lp,&updated_rx_lp_start,param_n,param_m,param_k,preprocess_param.boundcut_p,preprocess_param.boundcut_loop,*ub,feasible_x,ub,LPCC_precondition);
			if(status) goto TERMINATE;
		}
		if((*LPCC_precondition)==LPCC_PRECONDITION_OPT || (*LPCC_precondition)==LPCC_PRECONDITION_INFEASIBLE) goto TERMINATE; 

		/*  if (preprocess_param.MccormickRefineCnt>=0)
		{
			status=add_mccormickcut(env,
									rx_lp,
									param_n,
									param_m,
									param_k,
									matrix_N,
									matrix_M,
									q_coef,
									preprocess_param.MccormickRefineCnt,
									&mcc_processing_stat,
									&initial_subproblem_cnt,
									&initial_subproblem_stat,
									&x_start_index,
									&ybar_start_index,
									&mcc_start_row,
									&mcc_start_col,
									&p_x_lb,
									&p_x_ub,
									&p_y_lb,
									&p_y_ub,
									preprocess_param.partition_x,
									preprocess_param.partition_y);
                        printf(" status after add_mccormickcut %d\n",status);
			if(status)goto TERMINATE;
			if (mcc_processing_stat==-1)
			{
				*LPCC_precondition=LPCC_PRECONDITION_INFEASIBLE;
				goto TERMINATE;
			}
		}   */
		if (initial_subproblem_cnt==1)
		{
			status=copy_startinfo(&updated_rx_lp_start,&rx_lp_start);
			if(status) goto TERMINATE;
			status=cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat,&updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,2);
			if(status) goto TERMINATE;
			if (lp_solnstat!=CPX_STAT_OPTIMAL)
			{
				status=-1;
				fprintf(stderr," preprocess(): unexpected solnstat: %d\n",lp_solnstat);
				goto TERMINATE;
			}
			(*headnodeptr)->nodevalue=lp_obj;
			(*headnodeptr)->nodelb = lp_obj;
			status = copy_startinfo(&updated_rx_lp_start,&(*headnodeptr)->nodestartinfo);
			if(status) goto TERMINATE;
			for (j=0; j<num_var;j++)
			{
				(*headnodeptr)->nodex[j]=lp_soln[j];
			}
			*improve_lb = lp_obj;
			initial_problem_cnt=1;
		}else
		{
			popnode(headnodeptr,&tmp_node);
			for (i=0;i<initial_subproblem_cnt;i++)
			{
				if (initial_subproblem_stat[i]!=-1)
				{
					tmp_node->x_start_index=x_start_index;
					tmp_node->y_bar_start_index=ybar_start_index;
					tmp_node->mcc_index_col=mcc_start_col;
					tmp_node->mcc_index_row=mcc_start_row;
					for (j=0;j<param_n;j++)
					{
						tmp_node->x_lb[j]=p_x_lb[i][j];
						tmp_node->x_ub[j]=p_x_ub[i][j];
						tmp_node->y_bar_lb[j]=p_y_lb[i][j];
						tmp_node->y_bar_ub[j]=p_y_ub[i][j];
					}
					status=setupnode_rx_lp(env,
										   rx_lp,
										   param_n,
										   param_m,
										   param_k,
										   *tmp_node);
					if(status) goto TERMINATE;
					status=copy_startinfo(&updated_rx_lp_start,&rx_lp_start);
					if(status) goto TERMINATE;
					status=cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat,&updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,2);
					if(status) goto TERMINATE;
					if (lp_solnstat!=CPX_STAT_OPTIMAL && lp_solnstat!=CPX_STAT_INFEASIBLE)
					{
						status=-1;
						fprintf(stderr," preprocess(): unexpected solnstat: %d\n",lp_solnstat);
						goto TERMINATE;
					}
					tmp_node->nodevalue=lp_obj;
					tmp_node->nodelb = lp_obj;
					status = copy_startinfo(&updated_rx_lp_start,&tmp_node->nodestartinfo);
					if(status) goto TERMINATE;
					for (j=0; j<num_var;j++)
					{
						tmp_node->nodex[j]=lp_soln[j];
					}
					if (lp_solnstat==CPX_STAT_OPTIMAL)
					{
						if (*improve_lb<lp_obj)*improve_lb=lp_obj;
						status=pushnode_ordered(headnodeptr,*tmp_node);
						if(status) goto TERMINATE;
						initial_problem_cnt++;
					}
					
				}
			}		
		}
		fprintf(stdout," Initial subproblem generated: %d\n",initial_problem_cnt);
		break;
	case CPX_STAT_INFEASIBLE:
		*LPCC_precondition=LPCC_PRECONDITION_INFEASIBLE;
		break;
	case CPX_STAT_UNBOUNDED:
		*LPCC_precondition=LPCC_PRECONDITION_UNBOUNDED;
		(*headnodeptr)->nodevalue=-CPX_INFBOUND;
		(*headnodeptr)->nodelb = -CPX_INFBOUND;
		status = copy_startinfo(&updated_rx_lp_start,&(*headnodeptr)->nodestartinfo);
		if(status) goto TERMINATE;
		for (j=0; j<num_var;j++)
		{
			(*headnodeptr)->nodex[j]=lp_soln[j];
		}
		break;
	default:
		status=-1;
		fprintf(stderr," preprocess(): unexpected lp_solnstat: %d\n",lp_solnstat);
		goto TERMINATE;
	}
TERMINATE:
        printf(" terminating status in preprocess %d\n",status);
	free_and_null((char**) &lp_soln);
	free_node(tmp_node);
	free_startinfo(&rx_lp_start);
	free_startinfo(&updated_rx_lp_start);
	free_and_null((char**)&initial_subproblem_stat);
	return(status);
}


int conditionaltest(CPXENVptr env, 
					CPXLPptr rx_lp, 
					const int param_n,
					const int param_m, 
					const int param_k,
					const int violate_complementary_index, 
					int *fix_signal)
{
	int status=0;
	const int num_cols=CPXgetnumcols(env,rx_lp);
	int i;
	int row_id;
	int fix_1=1;
	int fix_2=1;
	double *tableau_1=NULL;
	double *tableau_2=NULL;
	tableau_1 = (double*) malloc(num_cols*sizeof(double));
	tableau_2 = (double*) malloc(num_cols*sizeof(double));
	if (tableau_1==NULL || tableau_2==NULL)
	{
		status=NO_MEMORY;
		fprintf(stderr," conditionaltest():could not allocate memory.\n");
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
		
		if (tableau_1[i]>10e-10)
		{
			fix_2=0;
		}
		if (tableau_2[i]>10e-10)
		{
			fix_1=0;
		}
	}

	if (fix_1==1 && fix_2==1)
	{
		*fix_signal=3;
	}else if (fix_1==1)
	{
		*fix_signal=1;
	}else if (fix_2==1)
	{
		*fix_signal=2;
	}else
	{
		*fix_signal=0;
	}
TERMINATE:
	free_and_null((char**)&tableau_1);
	free_and_null((char**)&tableau_2);
	return(status);
}

int Branch_Select_bounded_1(CPXENVptr env,
						  CPXLPptr rx_lp, 
						  const int param_n,
						  const int param_m, 
						  const int param_k, 
						  NODE activenode, 
						  double *ub, 
						  NODE ***branchnode_ptr, 
						  int **branch_indicator_p,
						  double cutofflevel)
{
	//branch_indicator 1: infeasible or overub node 2: bounded feasible LPCC node 3: bounded infeasible LPCC node
	int status = 0;
	int i,j,k;
	int tmp_infeasible_node_cnt;
	int tmp_overub_node_cnt;
	int tmp_feasible_soln_node_cnt;
	int tmp_active_node_cnt;
	int max_infeasible_overub_node_cnt=-1;
	//all tmp var are used to store information for strong branching
	int tmp_fixed_index;
	int tmp_branchIndex;
	int branchIndex;
	int *tmp_solnstat=NULL;
	int *tmp_branch_indicator_p=NULL; // 1: infeasible or overub 2: feasible 3: lb
	int *tmp_branchnode_direction=NULL;
	int *branchnode_direction=NULL;
	double tmp_active_node_lb;
	double max_active_node_lb=-CPX_INFBOUND;
	double *tmp_branch_objval=NULL;
	double *branch_objval=NULL;
	double **tmp_branch_x = NULL;
	double **branch_x=NULL;
	STARTINFO *tmp_branchnode_start=NULL;
	STARTINFO *branchnode_start=NULL;
	//initialize
	tmp_solnstat = (int*) malloc(2*sizeof(int));
	tmp_branch_indicator_p = (int*) malloc(2*sizeof(int));
	tmp_branchnode_direction = (int*) malloc(2*sizeof(int*));
	branchnode_direction=(int*) malloc(2*sizeof(int*));
	tmp_branch_objval=(double*) malloc(2*sizeof(double));
	branch_objval=(double*) malloc(2*sizeof(double));
	tmp_branch_x = (double**) malloc(2*sizeof(double*));
	branch_x=(double**) malloc(2*sizeof(double*));
	tmp_branchnode_start = (STARTINFO*) malloc(2*sizeof(STARTINFO));
	branchnode_start= (STARTINFO*) malloc(2*sizeof(STARTINFO));
	*branchnode_ptr = (NODE**) calloc(2,sizeof(NODE*));
	*branch_indicator_p = (int*) malloc(2*sizeof(int));
	if (tmp_solnstat == NULL 
		|| tmp_branch_indicator_p == NULL
		|| tmp_branchnode_direction == NULL || branchnode_direction==NULL
		|| tmp_branch_objval == NULL || branch_objval == NULL
		|| tmp_branch_x == NULL || branch_x == NULL
		|| tmp_branchnode_start == NULL || branchnode_start == NULL
		|| *branchnode_ptr == NULL || *branch_indicator_p == NULL)
	{
		status = NO_MEMORY;
		fprintf(stderr,"could not allocate memory.\n");
		goto TERMINATE;
	}
	for (i=0;i<2;i++)
	{
		init_startinfo(&tmp_branchnode_start[i]);
		init_startinfo(&branchnode_start[i]);
		tmp_branch_x[i]=(double*) malloc((param_n+param_m+param_m+param_k)*sizeof(double));
		branch_x[i]=(double*) malloc((param_n+param_m+param_m+param_k)*sizeof(double));
		if (tmp_branch_x[i]==NULL || branch_x[i]==NULL)
		{
			status=NO_MEMORY;
			fprintf(stderr,"could not allocate memory.\n");
			goto TERMINATE;
		}
		status=init_branch_node(activenode,&((*branchnode_ptr)[i]));
		if(status) goto TERMINATE;
		(*branch_indicator_p)[i]=-1;
	}
	tmp_branchIndex=-1;
	branchIndex=-1;
	status=setupnode_rx_lp(env,rx_lp,param_n,param_m,param_k,activenode);
	if(status) goto TERMINATE;
	nodecnt=nodecnt+2;
	//end initialize
	if(checkallComplementary(param_n,param_m,param_k,activenode.nodex)==0)
	{
		status=-1;
		fprintf(stderr,"Branch_Select_bounded_1(): unexpected activenode condition\n");
		goto TERMINATE;
	}
	//strong branching
	for ( i = 0; i < param_m; i++) 
	{
		if (activenode.nodebranchinfo[i]!=0 
			||(checkComplementary(param_n+i,activenode.nodex) != 1 || checkComplementary(param_n+param_m+i,activenode.nodex) !=1)) continue;

		tmp_infeasible_node_cnt=0;
		tmp_overub_node_cnt=0;
		tmp_feasible_soln_node_cnt=0;
		tmp_active_node_cnt=0;
		tmp_active_node_lb=-CPX_INFBOUND;
		tmp_branchIndex=i;
		for(j=0;j<2;j++)
		{
			tmp_fixed_index=param_n+i+param_m*j;
			tmp_branchnode_direction[j]=j+1;
			status=relaxLPSolver(env,rx_lp,&(activenode.nodecuts),(param_n+param_m+param_m+param_k),1,&tmp_fixed_index,&(activenode.nodestartinfo),tmp_branch_x[j],&tmp_branch_objval[j],&tmp_solnstat[j],&tmp_branchnode_start[j],0);
			if(status) goto TERMINATE;
			switch(tmp_solnstat[j])
			{
			case CPX_STAT_OPTIMAL : 
				if (tmp_branch_objval[j]>strengthen_value(*ub))
				{
					tmp_branch_indicator_p[j]= 1;
					tmp_overub_node_cnt++;
				}else
				{
					if (checkallComplementary(param_n,param_m,param_k,tmp_branch_x[j])==0)
					{
						tmp_branch_indicator_p[j]=2;
						tmp_feasible_soln_node_cnt++;
					}else{
						tmp_branch_indicator_p[j]=3;
						tmp_active_node_cnt++;
						if(tmp_active_node_lb<tmp_branch_objval[j]) tmp_active_node_lb=tmp_branch_objval[j];
					}
				}
				break;
			case CPX_STAT_INFEASIBLE: 
				tmp_branch_indicator_p[j]=1;
				tmp_infeasible_node_cnt++;
				break;
			default:
				status = -1;
				fprintf(stderr," Branch_Select_bounded(): unexpected solnstat: %d\n",tmp_solnstat[j]);
				goto TERMINATE;
			}
		}	
		if (tmp_feasible_soln_node_cnt>0 || max_active_node_lb<tmp_active_node_lb || tmp_active_node_cnt<2)
		{
			max_active_node_lb = tmp_active_node_lb;
			branchIndex=tmp_branchIndex;
			for (j=0;j<2;j++)
			{
				switch(tmp_branch_indicator_p[j])
				{
				case 1 : //infesible or overub
					(*branch_indicator_p)[j]=1;
					break;
				case 2 : //feasible
					(*branch_indicator_p)[j]=2;
					for (k=0;k<(param_n+param_m+param_m+param_k);k++)
					{
						branch_x[j][k]=tmp_branch_x[j][k];
					}
					branch_objval[j]=tmp_branch_objval[j];
					break;
				case 3 : //normal lb
					(*branch_indicator_p)[j]=3;
					branchnode_direction[j]=tmp_branchnode_direction[j];
					for (k=0;k<(param_n+param_m+param_m+param_k);k++)
					{
						branch_x[j][k]=tmp_branch_x[j][k];
					}
					branch_objval[j]=tmp_branch_objval[j];
					status=copy_startinfo(&tmp_branchnode_start[j],&branchnode_start[j]);
					if(status) goto TERMINATE;
					break;
				default:
					status = -1;
					fprintf(stderr," Branch_Select_bounded_1(): unexpected tmp_branch_indicator_p: %d\n",tmp_branch_indicator_p[j]);
					goto TERMINATE;
				}	
			}
			if (tmp_feasible_soln_node_cnt>0 || max_active_node_lb>=cutofflevel || tmp_active_node_cnt<2) break;
		}
	}
	//generating branching nodes
	for (j=0;j<2;j++)
	{
		switch((*branch_indicator_p)[j])
		{
		case 1 : //infesible or overub
			break;
		case 2 : //feasible
			for (k=0;k<(param_n+param_m+param_m+param_k);k++)
			{
				(*branchnode_ptr)[j]->nodex[k]=branch_x[j][k];
			}
			(*branchnode_ptr)[j]->nodelb=branch_objval[j];
			(*branchnode_ptr)[j]->nodevalue=branch_objval[j];
			break;
		case 3 : //normal lb
			(*branchnode_ptr)[j]->nodebranchinfo[branchIndex]=branchnode_direction[j];
			for (k=0;k<(param_n+param_m+param_m+param_k);k++)
			{
				(*branchnode_ptr)[j]->nodex[k]=branch_x[j][k];
			}
			(*branchnode_ptr)[j]->nodelb=branch_objval[j];
			(*branchnode_ptr)[j]->nodevalue=branch_objval[j];
			status=copy_startinfo(&branchnode_start[j],&((*branchnode_ptr)[j]->nodestartinfo));
			if(status) goto TERMINATE;
			break;
		default:
			status = -1;
			fprintf(stderr," Branch_Select_bounded(): unexpected branch_indicator_p: %d\n",(*branch_indicator_p)[j]);
			goto TERMINATE;
		}	
	}
TERMINATE:
	free(tmp_solnstat);
	free(tmp_branch_indicator_p);
	free(tmp_branch_objval);
	free(branch_objval);
	for (i=0;i<2;i++)
	{
		free(tmp_branch_x[i]);
		free(branch_x[i]);
		free_startinfo(&(tmp_branchnode_start[i]));
		free_startinfo(&(branchnode_start[i]));
	}
	free(tmp_branchnode_direction);
	free(branchnode_direction);
	free(tmp_branch_x);
	free(branch_x);
	free(tmp_branchnode_start);
	free(branchnode_start);
	return (status);
}

int Branch_Select_bounded_2(CPXENVptr env,
							CPXLPptr rx_lp, 
							const int param_n,
							const int param_m, 
							const int param_k, 
							NODE activenode, 
							double *ub, 
							NODE ***branchnode_ptr, 
							int **branch_indicator_p,
							BRANCHHIST *branch_hist)
{
	//branch_indicator 1: infeasible or overub node 2: bounded feasible LPCC node 3: bounded infeasible LPCC node
	int status = 0;
	int i,j,k;
	int branch_index_by_score;
	int branchIndex;
	int tmp_fixed_index;
	int *tmp_solnstat=NULL;
	int *tmp_branch_indicator_p=NULL; // 1: infeasible or overub 2: feasible 3: lb
	int *tmp_branchnode_direction=NULL;
	int *branchnode_direction=NULL;
	double select_violation_1,select_violation_2;
	double *tmp_branch_objval=NULL;
	double *branch_objval=NULL;
	double **tmp_branch_x = NULL;
	double **branch_x=NULL;
	STARTINFO *tmp_branchnode_start=NULL;
	STARTINFO *branchnode_start=NULL;
	//initialize
	tmp_solnstat = (int*) malloc(2*sizeof(int));
	tmp_branch_indicator_p = (int*) malloc(2*sizeof(int));
	tmp_branchnode_direction = (int*) malloc(2*sizeof(int));
	branchnode_direction=(int*) malloc(2*sizeof(int));
	tmp_branch_objval=(double*) malloc(2*sizeof(double));
	branch_objval=(double*) malloc(2*sizeof(double));
	tmp_branch_x = (double**) malloc(2*sizeof(double*));
	branch_x=(double**) malloc(2*sizeof(double*));
	tmp_branchnode_start = (STARTINFO*) malloc(2*sizeof(STARTINFO));
	branchnode_start= (STARTINFO*) malloc(2*sizeof(STARTINFO));
	*branchnode_ptr = (NODE**) calloc(2,sizeof(NODE*));
	*branch_indicator_p = (int*) malloc(2*sizeof(int));
	if (tmp_solnstat == NULL 
		|| tmp_branch_indicator_p == NULL
		|| tmp_branchnode_direction == NULL || branchnode_direction==NULL
		|| tmp_branch_objval == NULL || branch_objval == NULL
		|| tmp_branch_x == NULL || branch_x == NULL
		|| tmp_branchnode_start == NULL || branchnode_start == NULL
		|| *branchnode_ptr == NULL || *branch_indicator_p == NULL)
	{
		status = NO_MEMORY;
		fprintf(stderr,"could not allocate memory.\n");
		goto TERMINATE;
	}
	for (i=0;i<2;i++)
	{
		init_startinfo(&tmp_branchnode_start[i]);
		init_startinfo(&branchnode_start[i]);
		tmp_branch_x[i]=(double*) malloc((param_n+param_m+param_m+param_k)*sizeof(double));
		branch_x[i]=(double*) malloc((param_n+param_m+param_m+param_k)*sizeof(double));
		if (tmp_branch_x[i]==NULL || branch_x[i]==NULL)
		{
			status=NO_MEMORY;
			fprintf(stderr," Branch_Select_bounded_2(): could not allocate memory.\n");
			goto TERMINATE;
		}
		status=init_branch_node(activenode,&((*branchnode_ptr)[i]));
		if(status) goto TERMINATE;
		(*branch_indicator_p)[i]=-1;
	}
	status=setupnode_rx_lp(env,rx_lp,param_n,param_m,param_k,activenode);
	if(status) goto TERMINATE;
	nodecnt=nodecnt+2;
	//end initialize
	if(checkallComplementary(param_n,param_m,param_k,activenode.nodex)==0)
	{
		status=-1;
		fprintf(stderr,"Branch_Select_bounded(): unexpected activenode condition\n");
		goto TERMINATE;
	}
	status=select_branch_by_score(env,rx_lp,&(activenode.nodestartinfo),param_n,param_m,param_k,&branch_index_by_score,branch_hist,&activenode,&select_violation_1,&select_violation_2);
	if (status) goto TERMINATE;
	branchIndex=branch_index_by_score;
	if (activenode.nodebranchinfo[branch_index_by_score]!=0 
		|| (checkComplementary(param_n+branch_index_by_score,activenode.nodex) != 1 || checkComplementary(param_n+param_m+branch_index_by_score,activenode.nodex) !=1))
	{
		status=-1;
		goto TERMINATE;
	}
	for(j=0;j<2;j++)
	{
		tmp_fixed_index=param_n+branch_index_by_score+j*param_m;
		tmp_branchnode_direction[j]=j+1;
		status=relaxLPSolver(env,rx_lp,&(activenode.nodecuts),(param_n+param_m+param_m+param_k),1,&tmp_fixed_index,&(activenode.nodestartinfo),tmp_branch_x[j],&tmp_branch_objval[j],&tmp_solnstat[j],&tmp_branchnode_start[j],0);
		if(status) goto TERMINATE;
		switch(tmp_solnstat[j])
		{
		case CPX_STAT_OPTIMAL : 
			if (tmp_branch_objval[j]>strengthen_value(*ub))
			{
				tmp_branch_indicator_p[j]= 1;
			}else
			{
				if (checkallComplementary(param_n,param_m,param_k,tmp_branch_x[j])==0)
				{
					tmp_branch_indicator_p[j]=2;
				}else{
					tmp_branch_indicator_p[j]=3;
				}
			}
			break;
		case CPX_STAT_INFEASIBLE: 
			tmp_branch_indicator_p[j]=1;
			break;
		default:
			status = -1;
			fprintf(stderr," Branch_Select_bounded(): unexpected solnstat: %d\n",tmp_solnstat[j]);
			goto TERMINATE;
		}
	}	
	for (j=0;j<2;j++)
	{
		switch(tmp_branch_indicator_p[j])
		{
		case 1 : //infesible or overub
			(*branch_indicator_p)[j]=1;
			break;
		case 2 : //feasible
			(*branch_indicator_p)[j]=2;
			for (k=0;k<(param_n+param_m+param_m+param_k);k++)
			{
				branch_x[j][k]=tmp_branch_x[j][k];
			}
			branch_objval[j]=tmp_branch_objval[j];
			break;
		case 3 : //normal lb
			(*branch_indicator_p)[j]=3;
			branchnode_direction[j]=tmp_branchnode_direction[j];
			for (k=0;k<(param_n+param_m+param_m+param_k);k++)
			{
				branch_x[j][k]=tmp_branch_x[j][k];
			}
			branch_objval[j]=tmp_branch_objval[j];
			status=copy_startinfo(&tmp_branchnode_start[j],&branchnode_start[j]);
			if(status) goto TERMINATE;
			break;
		default:
			status = -1;
			fprintf(stderr," Branch_Select_bounded(): unexpected tmp_branch_indicator_p: %d\n",tmp_branch_indicator_p[j]);
			goto TERMINATE;
		}	
	}
	//generating branching nodes
	for (j=0;j<2;j++)
	{
		switch((*branch_indicator_p)[j])
		{
		case 1 : //infesible or overub
			break;
		case 2 : //feasible
			for (k=0;k<(param_n+param_m+param_m+param_k);k++)
			{
				(*branchnode_ptr)[j]->nodex[k]=branch_x[j][k];
			}
			(*branchnode_ptr)[j]->nodelb=branch_objval[j];
			(*branchnode_ptr)[j]->nodevalue=branch_objval[j];
			break;
		case 3 : //normal lb
			(*branchnode_ptr)[j]->nodebranchinfo[branchIndex]=branchnode_direction[j];
			(*branchnode_ptr)[j]->branchIndexofParentNode=branchIndex;
			(*branchnode_ptr)[j]->branchDirectionofParentNode=branchnode_direction[j];
			if(branchnode_direction[j]==1)(*branchnode_ptr)[j]->parentNodeBranchViolation=select_violation_1;
			else (*branchnode_ptr)[j]->parentNodeBranchViolation=select_violation_2;
			for (k=0;k<(param_n+param_m+param_m+param_k);k++)
			{
				(*branchnode_ptr)[j]->nodex[k]=branch_x[j][k];
			}
			(*branchnode_ptr)[j]->nodelb=branch_objval[j];
			(*branchnode_ptr)[j]->nodevalue=branch_objval[j];
			status=copy_startinfo(&branchnode_start[j],&((*branchnode_ptr)[j]->nodestartinfo));
			if(status) goto TERMINATE;
			break;
		default:
			status = -1;
			fprintf(stderr," Branch_Select_bounded(): unexpected branch_indicator_p: %d\n",(*branch_indicator_p)[j]);
			goto TERMINATE;
		}	
	}
TERMINATE:
	free(tmp_solnstat);
	free(tmp_branch_indicator_p);
	free(tmp_branch_objval);
	free(branch_objval);
	for (i=0;i<2;i++)
	{
		free(tmp_branch_x[i]);
		free(branch_x[i]);
		free_startinfo(&(tmp_branchnode_start[i]));
		free_startinfo(&(branchnode_start[i]));
	}
	free(branchnode_direction);
	free(tmp_branchnode_direction);
	free(tmp_branch_x);
	free(branch_x);
	free(tmp_branchnode_start);
	free(branchnode_start);
	return (status);
}

int select_branch_by_score(CPXENVptr env, 
				   CPXLPptr rx_lp, 	
				   STARTINFO *lp_start,								 
				   const int param_n,								 
				   const int param_m, 								 
				   const int param_k,	
				   int *branch_index,
				   BRANCHHIST *branch_hist,
				   NODE *activenode,
				   double *violation_1,
				   double *violation_2)
{
	int status =0;
	int j;
	int lp_solnstat;
	int select_index;
	
	int flag;
	int cnt=0;
	double highest_score;
	
	double tmp_violation_1;
	double tmp_violation_2;
	double tmp_coef_norm_1;
	double tmp_coef_norm_2;
	double tmp_reduce_norm_1;
	double tmp_reduce_norm_2;
	double tmp_parallelism_1;
	double tmp_parallelism_2;

	double tmp_violation_level_1;
	double tmp_violation_level_2;
	double tmp_psuedo_cost_1;
	double tmp_psuedo_cost_2;
	double tmp_score;

	double *coef_norm_product=NULL;
	double *reduce_norm_product=NULL;
	double *violation_product=NULL;
	double *violation_level_product=NULL;
	double *psuedo_cost_product=NULL;
	double *parallelism_product=NULL;
	
	double coef_norm_product_norm=0;
	double reduce_norm_product_norm=0;
	double violation_product_norm=0;
	double violation_level_product_norm=0;
	double psuedo_cost_product_norm=0;
	double parallelism_product_norm=0;

	double w1,w2,w3,w4,w5,w6;

	const int num_var= param_n+param_m+param_m+param_k;
	double *lp_soln=NULL;
	double lp_obj;

	double violation_cnt;
	double unfix_cnt;
	double violationlevel;


	STARTINFO rx_lp_start;
	STARTINFO updated_rx_lp_start;


	init_startinfo(&rx_lp_start);
	init_startinfo(&updated_rx_lp_start);
	status=copy_startinfo(lp_start,&rx_lp_start);
	if(status) goto TERMINATE;
	lp_soln=(double*) malloc(num_var*sizeof(double));
	coef_norm_product=(double*) calloc(param_m,sizeof(double));
	reduce_norm_product=(double*) calloc(param_m,sizeof(double));
	violation_product=(double*) calloc(param_m,sizeof(double));
	violation_level_product=(double*) calloc(param_m,sizeof(double));
	psuedo_cost_product=(double*) calloc(param_m,sizeof(double));
	parallelism_product=(double*) calloc(param_m,sizeof(double));
	if(lp_soln==NULL || coef_norm_product==NULL || reduce_norm_product==NULL || violation_product==NULL || violation_level_product==NULL || psuedo_cost_product==NULL || parallelism_product==NULL)
	{
		status = NO_MEMORY;	
		fprintf (stderr, " select_branch_by_score(): could not allocate memory.\n");
		goto TERMINATE;
	}

	status = cplexSolveLP(env,rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat, &updated_rx_lp_start,CPX_ALG_AUTOMATIC,1,2);
	if (status) goto TERMINATE;
	select_index=-1;
	highest_score=-CPX_INFBOUND;
	if (lp_solnstat==CPX_STAT_OPTIMAL)
	{
		violation_cnt=0;
		unfix_cnt=0;
		for (j=0;j<param_m;j++)
		{
			if (lp_soln[param_n+j]>ZERO_TOLERANCE && lp_soln[param_n+param_m+j]>ZERO_TOLERANCE)
			{
				violation_cnt++;
				//average_violation+=lp_soln[param_n+j]*lp_soln[param_n+param_m+j];
			}
			if (activenode->nodebranchinfo[j]==0)
			{
				unfix_cnt++;
			}
		}
		if (violation_cnt==0)
		{
			status=-1;
			goto TERMINATE;
		}
		if (activenode->branchIndexofParentNode==-1 || activenode->branchDirectionofParentNode==-1)
		{
			//fprintf(stdout,"headnode\n");
		}else
		{
			//fprintf(stdout,"parent branch index:%d: parent branch direction:%d\n",activenode->branchIndexofParentNode,activenode->branchDirectionofParentNode);
			//fprintf(stdout,"parent obj: %f; cur obj: %f\n",activenode->parentNodeValue,lp_obj);
			//fprintf(stdout,"violation level: %f; unfixcnt: %f; param_m: %d\n",violation_cnt/param_m,unfix_cnt,param_m);
			violationlevel=1-violation_cnt/param_m;
			status=update_branchhist_info(branch_hist,activenode->branchIndexofParentNode,activenode->branchDirectionofParentNode,violationlevel,fabs(lp_obj-activenode->parentNodeValue)/(activenode->parentNodeBranchViolation));
			if(status) goto TERMINATE;
		}
		flag=1;
		for (j=0;j<param_m;j++)
		{
			if (lp_soln[param_n+j]>ZERO_TOLERANCE && lp_soln[param_n+param_m+j]>ZERO_TOLERANCE)
			{
				//violation complementarity
				tmp_violation_1=lp_soln[param_n+j];
				tmp_violation_2=lp_soln[param_n+param_m+j];
				status=branch_score(env,rx_lp,param_n,param_m,param_k,j,lp_soln[param_n+j],lp_soln[param_n+param_m+j],&tmp_coef_norm_1,&tmp_coef_norm_2,&tmp_reduce_norm_1,&tmp_reduce_norm_2,&tmp_parallelism_1,&tmp_parallelism_2);
				if(status)goto TERMINATE;
				//fprintf(stdout,"%d; score: %f\n",j,tmp_score);
				if (tmp_coef_norm_1==-1 || tmp_coef_norm_2==-1)
				{
					select_index=j;
					break;
				}
				status=get_branchhist_info(branch_hist,j,1,&tmp_violation_level_1,&tmp_psuedo_cost_1);
				if(status) goto TERMINATE;
				status=get_branchhist_info(branch_hist,j,2,&tmp_violation_level_2,&tmp_psuedo_cost_2);
				if(status) goto TERMINATE;

				coef_norm_product[j]=pow(tmp_violation_1*tmp_violation_2/(tmp_coef_norm_1*tmp_coef_norm_2),0.5);
				reduce_norm_product[j]=pow(tmp_reduce_norm_1*tmp_reduce_norm_2,0.5);
				violation_product[j]=pow(tmp_violation_1*tmp_violation_2,0.5);
				violation_level_product[j]=pow(tmp_violation_level_1*tmp_violation_level_2,0.5);
				psuedo_cost_product[j]=pow(tmp_violation_1*tmp_violation_2*tmp_psuedo_cost_1*tmp_psuedo_cost_2,0.5);
				parallelism_product[j]=pow(tmp_parallelism_1*tmp_parallelism_2,0.5);

				coef_norm_product_norm+=pow(coef_norm_product[j],2);
				reduce_norm_product_norm+=pow(reduce_norm_product[j],2);
				violation_product_norm+=pow(violation_product[j],2);
				violation_level_product_norm+=pow(violation_level_product[j],2);
				psuedo_cost_product_norm+=pow(psuedo_cost_product[j],2);
				parallelism_product_norm+=pow(parallelism_product[j],2);
				//coef_norm_product_norm+=pow(coef_norm_product[j],1);
				//reduce_norm_product_norm+=pow(reduce_norm_product[j],1);
				//violation_product_norm+=pow(violation_product[j],1);
				//violation_level_product_norm+=pow(violation_level_product[j],1);
				//psuedo_cost_product_norm+=pow(psuedo_cost_product[j],1);
				//parallelism_product_norm+=pow(parallelism_product[j],1);

				cnt++;
				if (min(branch_hist->branch_1_cnt[j],branch_hist->branch_2_cnt[j])<=8)
				{
					flag=0;
				}

			}
		}

	} else
	{
		status=-1;
		fprintf(stderr," select_branch_by_score(): unexpected lp_solnstat: %d\n",lp_solnstat);
		goto TERMINATE;
	}
	if (select_index==-1)
	{
		coef_norm_product_norm=pow(coef_norm_product_norm,0.5);
		reduce_norm_product_norm=pow(reduce_norm_product_norm,0.5);
		violation_product_norm=pow(violation_product_norm,0.5);
		violation_level_product_norm=pow(violation_level_product_norm,0.5);
		psuedo_cost_product_norm=pow(psuedo_cost_product_norm,0.5);
		parallelism_product_norm=pow(parallelism_product_norm,0.5);
		//coef_norm_product_norm=pow(coef_norm_product_norm,1)/cnt;
		//reduce_norm_product_norm=pow(reduce_norm_product_norm,1)/cnt;
		//violation_product_norm=pow(violation_product_norm,1)/cnt;
		//violation_level_product_norm=pow(violation_level_product_norm,1)/cnt;
		//psuedo_cost_product_norm=pow(psuedo_cost_product_norm,1)/cnt;
		//parallelism_product_norm=pow(parallelism_product_norm,1)/cnt;
		w1=0;
		w2=2;
		w3=0;
		w4=0;
		w5=0;
		w6=0;
		//w1=1;
		//w2=1;
		//w3=0;
		//w4=2;
		//w5=2;
		//w6=0;
		for (j=0;j<param_m;j++)
		{
			if (lp_soln[param_n+j]>ZERO_TOLERANCE && lp_soln[param_n+param_m+j]>ZERO_TOLERANCE)
			{

				tmp_score=w1*coef_norm_product[j]/coef_norm_product_norm;
				tmp_score+=w2*violation_product[j]/violation_product_norm;	
				tmp_score+=w3*reduce_norm_product[j]/reduce_norm_product_norm;
				tmp_score+=w4*violation_level_product[j]/violation_level_product_norm;
				tmp_score+=w5*psuedo_cost_product[j]/psuedo_cost_product_norm;
				tmp_score+=w6*parallelism_product[j]/parallelism_product_norm;	
				//tmp_score=w1*coef_norm_product[j]/coef_norm_product_norm/(coef_norm_product[j]/coef_norm_product_norm+1);
				//tmp_score+=w2*violation_product[j]/violation_product_norm/(violation_product[j]/violation_product_norm+1);		
				//tmp_score+=w3*reduce_norm_product[j]/reduce_norm_product_norm/(reduce_norm_product[j]/reduce_norm_product_norm+1);
				//tmp_score+=w4*violation_level_product[j]/violation_level_product_norm/(violation_level_product[j]/violation_level_product_norm+1);
				//tmp_score+=w5*psuedo_cost_product[j]/psuedo_cost_product_norm/(psuedo_cost_product[j]/psuedo_cost_product_norm+1);
				//tmp_score+=w6*parallelism_product[j]/parallelism_product_norm/(parallelism_product[j]/parallelism_product_norm+1);	
				//fprintf(stdout,"%f; %f; %f; %f; %f; %f;\n",coef_norm_product[j]/coef_norm_product_norm,violation_product[j]/violation_product_norm,reduce_norm_product[j]/reduce_norm_product_norm,violation_level_product[j]/violation_level_product_norm,psuedo_cost_product[j]/psuedo_cost_product_norm,parallelism_product[j]/parallelism_product_norm);
				//fprintf(stdout,"%f,%f, %d, %d,\n",tmp_score,balance_score,branch_hist->branch_problem_cnt_1[j],branch_hist->branch_problem_cnt_2[j]);
				if (highest_score<tmp_score)
				{
					highest_score=tmp_score;
					select_index=j;
				}
				//fprintf(stdout,"tmp_score:%f\n",tmp_score);
			}
		}
		//fprintf(stdout," %d,%f,%f\n",select_index,highest_score,psuedo_cost_product_norm);
	}else
	{
		//fprintf(stdout,"fix complementarity\n");
	}

	*violation_1=lp_soln[param_n+select_index];
	*violation_2=lp_soln[param_n+param_m+select_index];
	//fprintf(stdout,"select idnex: %d, %f,%f\n",select_index,pow(violation_product[select_index]/coef_norm_product[select_index],2),pow(violation_product[select_index],2));
	*branch_index=select_index;
	//fprintf(stdout,"selected index:%d\n",select_index);
TERMINATE:
	free_and_null((char**) &lp_soln);
	free_and_null((char**) &coef_norm_product);
	free_and_null((char**) &reduce_norm_product);
	free_and_null((char**) &violation_product);
	free_and_null((char**) &violation_level_product);
	free_and_null((char**) &psuedo_cost_product);
	free_and_null((char**) &parallelism_product);
	free_startinfo(&rx_lp_start);
	free_startinfo(&updated_rx_lp_start);
	return(status);
}

int branch_score(CPXENVptr env, 
				 CPXLPptr rx_lp, 
				 const int param_n,
				 const int param_m,
				 const int param_k,
				 const int violate_complementary_index,
				 const double violation_complementary_1, 
				 const double violation_complementary_2, 
				 double *coef_norm_1,
				 double *coef_norm_2,
				 double *reduce_norm_1,
				 double *reduce_norm_2,
				 double *parallelism_1,
				 double *parallelism_2)
{
	int status=0;
	const int num_cols=CPXgetnumcols(env,rx_lp);
	const int num_rows=CPXgetnumrows(env,rx_lp);
	int i;
	int row_id;
	int fix_check_1;
	int fix_check_2;
	double coef_1_norm=0;
	double coef_2_norm=0;
	double reduce_cost_1_norm=0;
	double reduce_cost_2_norm=0;
	double parallelism_cost_1=0;
	double parallelism_cost_2=0;
	double *val = NULL;
	double *tableau_1=NULL;
	double *tableau_2=NULL;
	int *cstat=NULL;
	int *rstat=NULL;
	double *reduce_cost=NULL;
	double coef_1;
	double coef_2;
	double ub;
	tableau_1 = (double*) malloc(num_cols*sizeof(double));
	tableau_2 = (double*) malloc(num_cols*sizeof(double));
	cstat=(int*) malloc(num_cols*sizeof(int));
	rstat=(int*) malloc(num_rows*sizeof(int));
	reduce_cost=(double*) malloc(num_cols*sizeof(double));
	val = (double*) malloc(num_cols*sizeof(double));
	if (tableau_1==NULL || tableau_2==NULL || cstat==NULL || rstat==NULL || reduce_cost==NULL || val== NULL)
	{
		status=NO_MEMORY;
		fprintf(stderr," branch_score():could not allocate memory.\n");
		goto TERMINATE;
	}
	status=CPXgetbase(env,rx_lp,cstat,rstat);
	if (status) goto TERMINATE;
	status = CPXgetijrow(env,rx_lp,-1,param_n+violate_complementary_index,&row_id);
	if (status) goto TERMINATE;
	status = CPXbinvarow(env,rx_lp,row_id,tableau_1);
	if (status) goto TERMINATE;
	status = CPXgetijrow(env,rx_lp,-1,param_n+param_m+violate_complementary_index,&row_id);
	if (status) goto TERMINATE;
	status = CPXbinvarow(env,rx_lp,row_id,tableau_2);
	if(status) goto TERMINATE;
	status=CPXgetdj(env,rx_lp,reduce_cost,0,num_cols-1);
	if(status) goto TERMINATE;
	fix_check_1=1;
	fix_check_2=1;
	for (i=0;i<num_cols;i++)
	{
		//if (i>=param_n && i<(param_n+param_m))
		//{
		//	if (cstat[i]!=CPX_BASIC && cstat[i+param_m]!=CPX_BASIC)
		//	{
		//		
		//		if ((tableau_1[i]>EP && tableau_1[i+param_m]>EP) || (tableau_2[i]>EP && tableau_2[i+param_m]>EP))
		//		{
		//			fprintf(stdout,"one complementarity is all zero!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		//			fprintf(stdout,"%f,%f,%f,%f\n",tableau_1[i],tableau_2[i],tableau_1[i+param_m],tableau_2[i+param_m]);
		//		}
		//		
		//	}
		//}
		if (cstat[i]==CPX_BASIC) 
		{
			val[i]=0;
		}else
		{
			coef_1= tableau_1[i];//violation_complementary_1;
			coef_2= tableau_2[i];//violation_complementary_2;
			val[i]=coef_1>coef_2 ? coef_1:coef_2;
			if (cstat[i]==CPX_AT_LOWER || cstat[i]==CPX_FREE_SUPER || cstat[i]==CPX_AT_UPPER)
			{
				if (i>=param_n && i<(param_n+2*param_m))
				{
					CPXgetub(env,rx_lp,&ub,i,i);
					if (ub<1e-6) continue;
				}

				if (cstat[i]==CPX_AT_LOWER)
				{
					if (coef_1>ZERO_TOLERANCE) //
					{
						fix_check_1=0;	
						coef_1_norm+=pow(coef_1,2);
						reduce_cost_1_norm+=pow(reduce_cost[i],2);
						parallelism_cost_1+=reduce_cost[i]*coef_1;
					}
					if (coef_2>ZERO_TOLERANCE) // 
					{
						fix_check_2=0;	
						coef_2_norm+=pow(coef_2,2);
						reduce_cost_2_norm+=pow(reduce_cost[i],2);
						parallelism_cost_2+=reduce_cost[i]*coef_2;
					}
				}
				//if (cstat[i]==CPX_AT_UPPER)
				//{
				//	if (coef_1<-EP)
				//	{
				//		fix_check_1=0;	
				//		coef_1_norm+=pow(coef_1,2);
				//		reduce_cost_1_norm+=pow(reduce_cost[i],2);
				//		parallelism_cost_1+=reduce_cost[i]*coef_1;
				//	}
				//	if (coef_2<-EP)
				//	{
				//		fix_check_2=0;	
				//		coef_2_norm+=pow(coef_2,2);
				//		reduce_cost_2_norm+=pow(reduce_cost[i],2);
				//		parallelism_cost_2+=reduce_cost[i]*coef_2;
				//	}
				//}
				if (cstat[i]==CPX_FREE_SUPER)
				{
					if (fabs(coef_1)>ZERO_TOLERANCE)
					{
						fix_check_1=0;	
						coef_1_norm+=pow(coef_1,2);
						reduce_cost_1_norm+=pow(reduce_cost[i],2);
						parallelism_cost_1+=reduce_cost[i]*coef_1;
					}
					if (fabs(coef_2)>ZERO_TOLERANCE)
					{
						fix_check_2=0;
						coef_2_norm+=pow(coef_2,2);
						reduce_cost_2_norm+=pow(reduce_cost[i],2);
						parallelism_cost_2+=reduce_cost[i]*coef_2;
					}
				}
			}
		}
	}
	coef_1_norm=max(pow(coef_1_norm,0.5),ZERO_TOLERANCE);
	coef_2_norm=max(pow(coef_2_norm,0.5),ZERO_TOLERANCE);
	reduce_cost_1_norm=max(pow(reduce_cost_1_norm,0.5),ZERO_TOLERANCE);
	reduce_cost_2_norm=max(pow(reduce_cost_2_norm,0.5),ZERO_TOLERANCE);
	parallelism_cost_1=max(fabs(parallelism_cost_1)/(reduce_cost_1_norm*coef_1_norm),ZERO_TOLERANCE);
	parallelism_cost_2=max(fabs(parallelism_cost_2)/(reduce_cost_2_norm*coef_2_norm),ZERO_TOLERANCE);
	
	*coef_norm_1=coef_1_norm;
	*coef_norm_2=coef_2_norm;
	*reduce_norm_1=reduce_cost_1_norm;
	*reduce_norm_2=reduce_cost_2_norm;
	*parallelism_1=parallelism_cost_1;
	*parallelism_2=parallelism_cost_2;

	if (fix_check_1)
	{
		*coef_norm_1=-1;
		*reduce_norm_1=-1;
		*parallelism_1=-1;
	}
	if (fix_check_2)
	{
		*coef_norm_2=-1;
		*reduce_norm_2=-1;
		*parallelism_2=-1;
	}
TERMINATE:
	free_and_null((char**)&tableau_1);
	free_and_null((char**)&tableau_2);
	free_and_null((char**)&cstat);
	free_and_null((char**)&rstat);
	free_and_null((char**)&reduce_cost);
	free_and_null((char**)&val);
	return(status);
}

int LPCCSolver_bounded(CPXENVptr env,
					   CPXLPptr rx_lp,
					   const int param_n,
					   const int param_m, 
					   const int param_k,
					   NODE **headnodeptr, 
					   double *lb,
					   double *ub,
					   double *feasible_x, 
					   int *LPCC_stat,
					   PARAM parameter,
					   double *used_time)
{
	int status = 0;
	int i,j;
	int *branch_indicator=NULL;
	int fix_count=1;
	int node_count=2;
	int update_ub_flag=0;
	NODE **branchnode_ptr = NULL;
	NODE *activenode_ptr =NULL;
	double cutofflevel;
	double total_bc_time;
	BRANCHHIST branch_hist;
	clock_t t_start,t_end,t_current;
	t_start=clock();
	init_branchhist(&branch_hist);
	init_branchhist_info(&branch_hist,param_m);
	
	while(*headnodeptr!=NULL)
	{
        t_current=clock();
        //fprintf(stdout,"LB: %12.4f UB:%12.4f EN: %8d LN:%8d PT: %6.0fs\n",(*headnodeptr)->nodelb,*ub,nodecnt,leftnodes,(double)(t_current-t_start+(*used_time))/CLOCKS_PER_SEC);
        update_ub_flag=0;
		if(fabs(*ub-(*headnodeptr)->nodelb)/max(1e-10,fabs(*ub))<parameter.opt_tolerance || *ub<weaken_value((*headnodeptr)->nodelb)) 
		{
			popnode(headnodeptr,&activenode_ptr);
			free_node(activenode_ptr);
			free(activenode_ptr);
			activenode_ptr=NULL;
			continue;
		}
		if (fmod(nodecnt,500)==0)
		{
			t_current=clock();
			if(*ub<CPX_INFBOUND)fprintf(stdout,"mod500: LB: %12.4f UB:%12.4f EN: %8d LN:%8d PT: %6.0fs\n",(*headnodeptr)->nodelb,*ub,nodecnt,leftnodes,(double)(t_current-t_start+(*used_time))/CLOCKS_PER_SEC);
			else fprintf(stdout,"LB: %12.4f UB: -- EN: %8d, LN:%8d, PT: %6.0fs\n",(*headnodeptr)->nodelb,nodecnt,leftnodes,(double)(t_current-t_start+(*used_time))/CLOCKS_PER_SEC);
		}
		if (*ub<CPX_INFBOUND)
		{
			cutofflevel=(*ub);
		}else
		{
			cutofflevel=CPX_INFBOUND;
		}
		popnode(headnodeptr,&activenode_ptr);
		if (fixedcnt_node(activenode_ptr)<=7)
		{
			status=Branch_Select_bounded_1(env,rx_lp,param_n,param_m,param_k,*activenode_ptr,ub,&branchnode_ptr,&branch_indicator,cutofflevel);
			if (status) goto TERMINATE;
		}else
		{
			status=Branch_Select_bounded_2(env,rx_lp,param_n,param_m,param_k,*activenode_ptr,ub,&branchnode_ptr,&branch_indicator,&branch_hist);
			if (status) goto TERMINATE;
		}
		

		// get sort index 
		
		//
		for (i=0;i<node_count;i++)
		{	
			switch(branch_indicator[i])
			{
			case 1 : //infesible or overub
				break;
			case 2 : //feasible
				if ((*(branchnode_ptr[i])).nodelb<weaken_value(*ub))
				{
					*ub=(*(branchnode_ptr[i])).nodelb;
					update_ub_flag=1;
					for( j= 0; j< (param_n+param_m+param_m+param_k); j++)
					{
						feasible_x[j] = (*(branchnode_ptr[i])).nodex[j];
					}
				}
				break;
			case 3 : //normal lb
				pushnode_ordered(headnodeptr,*(branchnode_ptr[i]));
				break;
			default:
				status = -1;
				fprintf(stderr," LPCC_solver(): unexpected branch_indicator: %d\n",branch_indicator[i]);
				goto TERMINATE;
			}	
		}
		if (update_ub_flag)
		{
			fprintf(stdout,"	Find better upperbound: %f\n", *ub);
			status=lp_ub_control(env,rx_lp,*ub,0);
			if (status) goto TERMINATE;
			filternode(headnodeptr,*ub);
		}
		free_node(activenode_ptr);	
		free(activenode_ptr);
		activenode_ptr=NULL;
		
		free(branch_indicator);
		branch_indicator=NULL;
		
		if (branchnode_ptr!=NULL)
		{
			for (i=0;i<node_count;i++)
			{
				free_node(branchnode_ptr[i]);
				free(branchnode_ptr[i]);
				branchnode_ptr[i]=NULL;
			}
			free(branchnode_ptr);
			branchnode_ptr=NULL;
		}
		t_end=clock();
		total_bc_time=(double)(t_end-t_start+(*used_time))/CLOCKS_PER_SEC;
		if (total_bc_time>parameter.time_limit)
		{
			fprintf(stdout,"Time limit hit\n");
			break;
		}
	}
	if (*ub<CPX_INFBOUND && *headnodeptr==NULL)
	{
		*LPCC_stat=LPCC_OPTIMAL;
		*lb=*ub;
	}
	if (*ub<CPX_INFBOUND && *headnodeptr!=NULL)
	{
		*LPCC_stat=LPCC_TERMINATE_WITH_SOLN;
		*lb=(*headnodeptr)->nodelb;
	}
	if (*ub==CPX_INFBOUND && *headnodeptr!=NULL)
	{
		*LPCC_stat=LPCC_TERMINATE_WITHOUT_SOLN;
		*lb=(*headnodeptr)->nodelb;
	}
	if (*ub==CPX_INFBOUND && *headnodeptr==NULL)
	{
		*LPCC_stat=LPCC_INFEASIBLE;
	}
	
TERMINATE:
	//fprintf(stdout,"UB:%12.4f EN: %8d LN:%8d PT: %6.0fs\n",*ub,nodecnt,leftnodes,(double)(t_current-t_start+(*used_time))/CLOCKS_PER_SEC);
    free_branchhist(&branch_hist);
    free_node(activenode_ptr);
	free(activenode_ptr);
	free(branch_indicator);
	if (branchnode_ptr !=NULL)
	{
		for (i=0;i<node_count;i++)
		{
			free_node(branchnode_ptr[i]);
		}
		free(branchnode_ptr);
		branchnode_ptr=NULL;
	}
	return(status);
}

int Branch_Select_unbounded(CPXENVptr env, 
						   CPXLPptr rx_lp,
						   const int param_n,
						   const int param_m, 
						   const int param_k,
						   NODE activenode, 
						   double *ub,
						   const int fix_count, 
						   const int nodes_count,
						   int *unbounded_stat, 
						   double *feasible_ray, 
						   NODE ***branchnode_ptr, 
						   int **branch_indicator_p)
{
	//branch_indicator 1: infeasible or overub node 2: bounded feasible LPCC node 3: bounded infeasible LPCC node 4: unbounded node(unbounded LPCC ray) 5: unbounded node(no unbounded LPCC ray)
	int status = 0;
	int i,j,k,l,m,n;
	int tmp_infeasible_node_cnt;
	int tmp_overub_node_cnt;
	int tmp_bounded_feasible_soln_node_cnt;
	int tmp_bounded_active_node_cnt;
	int tmp_unbounded_ray_node_cnt;
	int tmp_unbounded_active_node_cnt;
	const int num_var=param_n+param_m+param_m+param_k;
	//int max_infeasible_overub_node_cnt=-1;

	//all tmp var are used to store information for strong branching
	int *tmp_fixed_index=NULL;
	int *tmp_solnstat=NULL;

	int *tmp_branchIndex=NULL;
	int *branchIndex=NULL;

	int *tmp_branch_indicator_p=NULL; // 1: infeasible or overub 2: feasible 3: lb

	int **tmp_branchnode_direction=NULL;
	int **branchnode_direction=NULL;

	int checkflag=0;
	int recover_solnstat;
	double recover_objval;
	double *recover_x=NULL;

	//double tmp_average_active_node_lb;
	//double max_average_active_node_lb=-CPX_INFBOUND;

	double *tmp_branch_objval=NULL;
	double *branch_objval=NULL;
	
	double *tmp_voilation_cnt=NULL;
	double violation_cnt;
	double max_violation_cnt=CPX_INFBOUND;

	double **tmp_branch_x = NULL;
	double **branch_x=NULL;

	STARTINFO *tmp_branchnode_start=NULL;
	STARTINFO *branchnode_start=NULL;
	//initialize
	tmp_fixed_index = (int*) malloc(fix_count*sizeof(int));
	tmp_solnstat = (int*) malloc(nodes_count*sizeof(int));

	tmp_branchIndex = (int*) malloc(fix_count*sizeof(int));
	branchIndex=(int*) malloc(fix_count*sizeof(int));

	tmp_branch_indicator_p = (int*) malloc(nodes_count*sizeof(int));

	tmp_branchnode_direction = (int**) malloc(nodes_count*sizeof(int*));
	branchnode_direction=(int**) malloc(nodes_count*sizeof(int*));
	
	recover_x=(double*)malloc(num_var*sizeof(double));

	tmp_branch_objval=(double*) malloc(nodes_count*sizeof(double));
	branch_objval=(double*) malloc(nodes_count*sizeof(double));

	tmp_voilation_cnt=(double*) malloc(nodes_count*sizeof(double));

	tmp_branch_x = (double**) malloc(nodes_count*sizeof(double*));
	branch_x=(double**) malloc(nodes_count*sizeof(double*));

	tmp_branchnode_start = (STARTINFO*) malloc(nodes_count*sizeof(STARTINFO));
	branchnode_start= (STARTINFO*) malloc(nodes_count*sizeof(STARTINFO));

	*branchnode_ptr = (NODE**) calloc(nodes_count,sizeof(NODE*));
	*branch_indicator_p = (int*) malloc(nodes_count*sizeof(int));

	if (tmp_fixed_index == NULL ||  tmp_solnstat == NULL 
		|| tmp_branchIndex == NULL || branchIndex==NULL 
		|| tmp_branch_indicator_p == NULL
		|| *tmp_branchnode_direction == NULL || *branchnode_direction==NULL
		|| recover_x==NULL
		|| tmp_branch_objval == NULL || branch_objval == NULL
		|| tmp_voilation_cnt==NULL
		|| *tmp_branch_x == NULL || branch_x == NULL
		|| tmp_branchnode_start == NULL || branchnode_start == NULL
		|| *branchnode_ptr == NULL || *branch_indicator_p == NULL)
	{
		status = NO_MEMORY;
		fprintf(stderr," Branch_Select_unbounded(): could not allocate memory.\n");
		goto TERMINATE;
	}
	for (i=0;i<nodes_count;i++)
	{
		init_startinfo(&tmp_branchnode_start[i]);
		init_startinfo(&branchnode_start[i]);
		tmp_branchnode_direction[i]= (int*) malloc(fix_count*sizeof(int));
		branchnode_direction[i]=(int*) malloc(fix_count*sizeof(int));
		tmp_branch_x[i]=(double*) malloc(num_var*sizeof(double));
		branch_x[i]=(double*) malloc(num_var*sizeof(double));
		if (tmp_branch_x[i]==NULL || branch_x[i]==NULL || tmp_branchnode_direction[i]==NULL || branchnode_direction[i]==NULL)
		{
			status=NO_MEMORY;
			fprintf(stderr," Branch_Select_unbounded(): could not allocate memory.\n");
			goto TERMINATE;
		}
		status=init_branch_node(activenode,&((*branchnode_ptr)[i]));
		if(status) goto TERMINATE;
		(*branch_indicator_p)[i]=-1;
	}
	for (i=0;i<fix_count;i++)
	{
		tmp_branchIndex[i]=-1;
		branchIndex[i]=-1;
	}
	status=setupnode_rx_lp(env,rx_lp,param_n,param_m,param_k,activenode);
	if(status) goto TERMINATE;
	nodecnt=nodecnt+2;
	// we need to check whether n_disj can be divided by fix_count
	if (fmod(param_m,fix_count)>0.001)
	{
		status=-1;
		fprintf(stderr,"Branch_Select_unbounded(): n_disj is not a multipler of fix_count\n");
		goto TERMINATE;
	}
	//end of initialize
	if(checkallComplementary(param_n,param_m,param_k,activenode.nodex)==0)
	{
		fprintf(stdout," find potential unbounded ray,");
		status=LPCC_feasible_recover_by_rounding(env,rx_lp,activenode.nodestartinfo,param_n,param_m,param_k,num_var,activenode,activenode.nodex,&recover_objval,recover_x,&recover_solnstat);
		if (status) goto TERMINATE;
		if(recover_solnstat==CPX_STAT_UNBOUNDED)
		{
			fprintf(stdout," recover unbounded ray\n");
			for( j= 0; j< num_var; j++)
			{
				feasible_ray[j] = recover_x[j];
			}
			*unbounded_stat=-1;
			goto TERMINATE;
		}else if(recover_solnstat==CPX_STAT_INFEASIBLE)
		{
			fprintf(stdout," unable to recover unbounded ray\n");
			checkflag=1;
		}else
		{
			fprintf(stderr,"Branch_Select_unbounded(): unexpected recover_solnstat: %d\n",recover_solnstat );
			status=-1;
			goto TERMINATE;
		}
	}
	//branch select
	for ( i = 0; i < param_m/fix_count; i++) 
	{
		k=0;
		for(j=0;j<fix_count;j++)
		{
			if(activenode.nodebranchinfo[i*fix_count+j]==0) k++;
		}
		if(k>0 && k<fix_count)
		{
			status = -1;
			fprintf(stderr,"Branch_Select_unbounded(): unexpected nodebranchinfo\n");
			goto TERMINATE;
		}
		if(k==fix_count)
		{
			l=0;
			for(j=0;j<fix_count;j++)
			{
				if ((checkComplementary(param_n+i*fix_count+j,activenode.nodex) == 1 && checkComplementary(param_n+param_m+i*fix_count+j,activenode.nodex) ==1) || checkflag)
				{
					l++;
					break;
				}
			}
			
			if(l>0)
			{
				tmp_infeasible_node_cnt=0;
				tmp_overub_node_cnt=0;
				tmp_bounded_feasible_soln_node_cnt=0;
				tmp_bounded_active_node_cnt=0;
				tmp_unbounded_ray_node_cnt=0;
				tmp_unbounded_active_node_cnt=0;
				violation_cnt=CPX_INFBOUND;
				//tmp_average_active_node_lb=-CPX_INFBOUND;
				for(j=0;j<fix_count;j++)
				{
					tmp_branchIndex[j]=i*fix_count+j;
				}
				for(j=0;j<nodes_count;j++)
				{
					//get the fix index for particular node
					n=j;
					for (m=0;m<fix_count;m++)
					{
						if(fmod(n,2)>0.001)
						{
							tmp_fixed_index[m]=param_n+fix_count*i+m;
							tmp_branchnode_direction[j][m]=1;
						}else
						{
							tmp_fixed_index[m]=param_n+param_m+fix_count*i+m;
							tmp_branchnode_direction[j][m]=2;
						}
						n=n/2;
					}
					status=relaxLPSolver(env,rx_lp,&(activenode.nodecuts),num_var,fix_count,tmp_fixed_index,&(activenode.nodestartinfo),tmp_branch_x[j],&tmp_branch_objval[j],&tmp_solnstat[j],&tmp_branchnode_start[j],1);
					if(status) goto TERMINATE;
					switch(tmp_solnstat[j])
					{
					case CPX_STAT_OPTIMAL : 
						if (tmp_branch_objval[j]>strengthen_value(*ub))
						{
							tmp_branch_indicator_p[j]= 1;
							tmp_overub_node_cnt++;
						}else
						{
							if (checkallComplementary(param_n,param_m,param_k,tmp_branch_x[j])==0)
							{
								tmp_branch_indicator_p[j]=2;
								tmp_bounded_feasible_soln_node_cnt++;
							}else
							{
								tmp_branch_indicator_p[j]=3;
								tmp_bounded_active_node_cnt++;
								//tmp_average_active_node_lb+=tmp_branch_objval[j];
								//if(tmp_average_active_node_lb<tmp_branch_objval[j]) tmp_average_active_node_lb=tmp_branch_objval[j];
							}
						}
						break;
					case CPX_STAT_INFEASIBLE: 
						tmp_branch_indicator_p[j]=1;
						tmp_infeasible_node_cnt++;
						break;
					case CPX_STAT_UNBOUNDED:
						if (checkallComplementary(param_n,param_m,param_k,tmp_branch_x[j])==0)
						{
							tmp_branch_indicator_p[j]=4;
							tmp_unbounded_ray_node_cnt++;
							//tmp_voilation_cnt[j]=countallHomoComplementary(&D1_p->cons_matrix,&D2_p->cons_matrix,tmp_branch_x[j]);
							//fprintf(stdout,"violation_cnt:%f\n",tmp_voilation_cnt[j]);
						}else
						{
							tmp_branch_indicator_p[j]=5;
							tmp_unbounded_active_node_cnt++;
							tmp_voilation_cnt[j]=countallComplementary(param_n,param_m,param_k,tmp_branch_x[j]);
							if (violation_cnt>tmp_voilation_cnt[j])
							{
								violation_cnt=tmp_voilation_cnt[j];
							}			
						}
						break;
					default:
						status = -1;
						fprintf(stderr," Branch_Select_unbounded(): unexpected solnstat: %d\n",tmp_solnstat[j]);
						goto TERMINATE;
					}
				}	
				//fprintf(stdout,"%d\n",tmp_bounded_active_node_cnt+tmp_overub_node_cnt+tmp_infeasible_node_cnt);
				if (tmp_unbounded_ray_node_cnt>0 || tmp_bounded_feasible_soln_node_cnt>0|| tmp_infeasible_node_cnt>0 || tmp_overub_node_cnt>0 || tmp_bounded_active_node_cnt>0 || max_violation_cnt>violation_cnt)
				{
					max_violation_cnt=violation_cnt;
					//fprintf(stdout,"%f\n",max_violation_cnt);
					for (j=0;j<fix_count;j++)
					{
						branchIndex[j]=tmp_branchIndex[j];
					}
					for (j=0;j<nodes_count;j++)
					{
						switch(tmp_branch_indicator_p[j])
						{
						case 1 : 
							(*branch_indicator_p)[j]=1;
							break;
						case 2 : 
							(*branch_indicator_p)[j]=2;
							for (m=0;m<num_var;m++)
							{
								branch_x[j][m]=tmp_branch_x[j][m];
							}
							branch_objval[j]=tmp_branch_objval[j];
							break;
						case 3 : 
							(*branch_indicator_p)[j]=3;
							for (m=0;m<fix_count;m++)
							{
								branchnode_direction[j][m]=tmp_branchnode_direction[j][m];
							}
							for (m=0;m<num_var;m++)
							{
								branch_x[j][m]=tmp_branch_x[j][m];
							}
							branch_objval[j]=tmp_branch_objval[j];
							status=copy_startinfo(&tmp_branchnode_start[j],&branchnode_start[j]);
							if(status) goto TERMINATE;
							break;
						case 4:
							(*branch_indicator_p)[j]=4;
							for (m=0;m<fix_count;m++)
							{
								branchnode_direction[j][m]=tmp_branchnode_direction[j][m];
							}
							for (m=0;m<num_var;m++)
							{
								branch_x[j][m]=tmp_branch_x[j][m];
							}
							branch_objval[j]=-CPX_INFBOUND;
							status=copy_startinfo(&tmp_branchnode_start[j],&branchnode_start[j]);
							if(status) goto TERMINATE;
							break;
						case 5:
							(*branch_indicator_p)[j]=5;
							for (m=0;m<fix_count;m++)
							{
								branchnode_direction[j][m]=tmp_branchnode_direction[j][m];
							}
							for (m=0;m<num_var;m++)
							{
								branch_x[j][m]=tmp_branch_x[j][m];
							}
							branch_objval[j]=-CPX_INFBOUND;
							status=copy_startinfo(&tmp_branchnode_start[j],&branchnode_start[j]);
							if(status) goto TERMINATE;
							break;
						default:
							status = -1;
							fprintf(stderr," Branch_Select_unbounded(): unexpected tmp_branch_indicator_p: %d\n",tmp_branch_indicator_p[j]);
							goto TERMINATE;
						}	
					}
					if (tmp_unbounded_ray_node_cnt>0 || tmp_bounded_feasible_soln_node_cnt>0|| tmp_infeasible_node_cnt>0 ||  tmp_overub_node_cnt>0 || tmp_bounded_active_node_cnt>0 ) break;
				}
			}
		}
	}
	//generating branching nodes
	for (j=0;j<nodes_count;j++)
	{
		switch((*branch_indicator_p)[j])
		{
		case 1: 
			break;
		case 2: 
			for (m=0;m<num_var;m++)
			{
				(*branchnode_ptr)[j]->nodex[m]=branch_x[j][m];
			}
			(*branchnode_ptr)[j]->nodelb=branch_objval[j];
			(*branchnode_ptr)[j]->nodevalue=branch_objval[j];
			break;
		case 3: 
		case 4:
		case 5:
			for (m=0;m<fix_count;m++)
			{
				(*branchnode_ptr)[j]->nodebranchinfo[branchIndex[m]]=branchnode_direction[j][m];
			}
			for (m=0;m<num_var;m++)
			{
				(*branchnode_ptr)[j]->nodex[m]=branch_x[j][m];
			}
			(*branchnode_ptr)[j]->nodelb=branch_objval[j];
			(*branchnode_ptr)[j]->nodevalue=branch_objval[j];
			status=copy_startinfo(&branchnode_start[j],&((*branchnode_ptr)[j]->nodestartinfo));
			if(status) goto TERMINATE;
			break;
		default:
			status = -1;
			fprintf(stderr," Branch_Select_unbounded(): unexpected branch_indicator_p: %d\n",(*branch_indicator_p)[j]);
			goto TERMINATE;
		}	
	}
	
TERMINATE:
	free_and_null((char**) &tmp_fixed_index);
	free_and_null((char**) &tmp_solnstat);
	free_and_null((char**) &tmp_branchIndex);
	free_and_null((char**) &branchIndex);
	free_and_null((char**) &tmp_branch_indicator_p);
	free_and_null((char**) &tmp_branch_objval);
	free_and_null((char**) &branch_objval);
	free_and_null((char**) &tmp_voilation_cnt);
	free_and_null((char**) &recover_x);
	for (i=0;i<nodes_count;i++)
	{
		free_and_null((char**) &tmp_branchnode_direction[i]);
		free_and_null((char**) &branchnode_direction[i]);
		free_and_null((char**) &tmp_branch_x[i]);
		free_and_null((char**) &branch_x[i]);
		free_startinfo(&tmp_branchnode_start[i]);
		free_startinfo(&branchnode_start[i]);
	}

	free_and_null((char**) &tmp_branchnode_direction);
	free_and_null((char**) &branchnode_direction);
	free_and_null((char**) &tmp_branch_x);
	free_and_null((char**) &branch_x);
	free_and_null((char**) &tmp_branchnode_start);
	free_and_null((char**) &branchnode_start);

	return(status);
}
int LPCCSolver_unbounded(CPXENVptr env, 
						 CPXLPptr rx_lp, 
						 const int param_n,
						 const int param_m, 
						 const int param_k,
						 NODE **headnodeptr, 
						 double *lb,
						 double *ub,
						 double *feasible_x, 
						 int *LPCC_stat,
						 PARAM parameter,
						 double *used_time)
{
	int status = 0;
	int i,j;
	int *branch_indicator=NULL;
	int unbounded_node_cnt=1;
	int solnstat=0;
	int fix_count=1;
	int node_count=2;
	int update_ub_flag=0;
	const int num_var=param_n+param_m+param_m+param_k;
	NODE **branchnode_ptr = NULL;
	NODE *activenode_ptr =NULL;
	NODE *bounded_headnode_ptr=NULL;
	BRANCHHIST branch_hist;
	BRANCHHIST bounded_branch_hist;
	clock_t t_start,t_end,t_current;
	double total_bc_time;
	*used_time=0;
	t_start=clock();
	init_branchhist(&branch_hist);
	init_branchhist_info(&branch_hist,param_m);
	init_branchhist(&bounded_branch_hist);
	init_branchhist_info(&bounded_branch_hist,param_m);
	while(*headnodeptr!=NULL)
	{
		//fprintf(stdout," unbounded nodes left: %d\n",unbounded_node_cnt);
		if (fmod(nodecnt,500)==0)
		{
			t_current=clock();
			if(*ub<CPX_INFBOUND)fprintf(stdout,"LB: -- UB:%12.4f EN: %8d LN:%8d PT: %6.0fs\n",(*headnodeptr)->nodelb,*ub,nodecnt,leftnodes,(double)(t_current-t_start)/CLOCKS_PER_SEC);
			else fprintf(stdout,"LB: -- UB: -- EN: %8d, LN:%8d, PT: %6.0fs\n",(*headnodeptr)->nodelb,nodecnt,leftnodes,(double)(t_current-t_start)/CLOCKS_PER_SEC);
		}
		popnode(headnodeptr,&activenode_ptr);
		unbounded_node_cnt--;
		update_ub_flag=0;
		status=Branch_Select_unbounded(env,rx_lp,param_n,param_m,param_k,*activenode_ptr,ub,fix_count,node_count,&solnstat,feasible_x,&branchnode_ptr,&branch_indicator);
		if (status) goto TERMINATE;
		if (solnstat==-1)
		{
			*LPCC_stat=LPCC_UNBOUNDED;
			deletenode(headnodeptr);
			deletenode(&bounded_headnode_ptr);
			goto TERMINATE;
		}
		for (i=0;i<node_count;i++)
		{	
			switch(branch_indicator[i])
			{
			case 1 : 
				break;
			case 2 : 
				if ((*(branchnode_ptr[i])).nodelb<=strengthen_value(*ub))
				{
					*ub=(*(branchnode_ptr[i])).nodelb;
					update_ub_flag=1;
					for( j= 0; j< num_var; j++)
					{
						feasible_x[j] = (*(branchnode_ptr[i])).nodex[j];
					}
				}
				break;
			case 3 : 
				pushnode_ordered(&bounded_headnode_ptr,*(branchnode_ptr[i]));
				break;
			case 4:
				pushnode_unordered(headnodeptr,*(branchnode_ptr[i]));
				unbounded_node_cnt++;
				break;
			case 5:
				pushnode_unordered(headnodeptr,*(branchnode_ptr[i]));
				unbounded_node_cnt++;
				break;
			default:
				status = -1;
				fprintf(stderr,"LPCCSolver_unbounded(): unexpected branch_indicator: %d\n",branch_indicator[i]);
				goto TERMINATE;
			}	
		}
		if (update_ub_flag)
		{
			fprintf(stdout,"	Find better upperbound: %f\n", *ub);
			filternode(&bounded_headnode_ptr,*ub);
		}
		free_node(activenode_ptr);
		free(activenode_ptr);
		activenode_ptr=NULL;
		free_and_null((char**)&branch_indicator);
		if (branchnode_ptr!=NULL)
		{
			for (i=0;i<node_count;i++)
			{
				free_node(branchnode_ptr[i]);
			}
			free(branchnode_ptr);
			branchnode_ptr=NULL;
		}
		t_end=clock();
		total_bc_time=(double)(t_end-t_start)/CLOCKS_PER_SEC;
		if (total_bc_time>parameter.time_limit)
		{
			fprintf(stdout,"Time limit hit\n");
			break;
		}
	}
	if (total_bc_time>parameter.time_limit)
	{
		if (*ub<CPX_INFBOUND && *headnodeptr!=NULL)
		{
			*lb=-CPX_INFBOUND;
			*LPCC_stat=LPCC_TERMINATE_WITH_SOLN;
		}
		if (*ub==CPX_INFBOUND && *headnodeptr!=NULL)
		{
			*lb=-CPX_INFBOUND;
			*LPCC_stat=LPCC_TERMINATE_WITHOUT_SOLN;
		}
	}else
	{
		if(bounded_headnode_ptr==NULL)
		{
			if (*headnodeptr==NULL && *ub<CPX_INFBOUND)
			{
				*LPCC_stat=LPCC_OPTIMAL;
				*lb=*ub;
			}
			if (*ub==CPX_INFBOUND && *headnodeptr==NULL)
			{
				*LPCC_stat=LPCC_INFEASIBLE;
				*lb=CPX_INFBOUND;

			}
		}else
		{
			*used_time=total_bc_time;
			//fprintf(stdout,"LPCC is bounded, continue to solve bounded problem\n" );
			status=LPCCSolver_bounded(env,rx_lp,param_n,param_m,param_k,&bounded_headnode_ptr,lb,ub,feasible_x,LPCC_stat,parameter,used_time);
			if (status) goto TERMINATE;
		}
	}
TERMINATE:
	free_node(activenode_ptr);
	free_branchhist(&branch_hist);
	free_branchhist(&bounded_branch_hist);
	free(activenode_ptr);
	free_and_null((char**)&branch_indicator);
	if (*branchnode_ptr !=NULL)
	{
		for (i=0;i<node_count;i++)
		{
			free_node(branchnode_ptr[i]);
		}
		free(branchnode_ptr);
	}
	return(status);

}
int LPCC_feasible_recover_by_rounding(CPXENVptr env, 
									  CPXLPptr rx_lp, 
									  STARTINFO rx_lp_start,
									  const int param_n,
									  const int param_m, 
									  const int param_k,
									  const int num_var,
									  NODE active_node, 
									  double *rounding_point,
									  double *recover_objval, 
									  double *recover_x, 
									  int *recover_solnstat)
{
	int status =0;
	int indices;
	int i;
	int *activeset=active_node.nodebranchinfo;
	char lu;
	double bd;
	double gap1;
	double gap2;
	CPXLPptr copy_rx_lp = NULL;
	STARTINFO tmp;
	init_startinfo(&tmp);
	copy_rx_lp = CPXcloneprob(env,rx_lp,&status);
	if (status) {
		fprintf(stderr," LPCC_feasible_recover_by_rounding(): unable to copy rx_lp \n");
		goto TERMINATE;
	}
	status=setupnode_rx_lp(env,copy_rx_lp,param_n,param_m,param_k,active_node);
	for (i=0; i< param_m; i++) {	
		if ( activeset[i] ==0) {	
			gap1 = rounding_point[param_n+i];
			gap2 =rounding_point[param_n+param_m+i] ;
			if (gap1<=gap2){
				indices = param_n+i;	
			} else {		
				indices = param_n+param_m + i;
			}
			lu='U';
			bd=0.0;
			status=CPXtightenbds(env,copy_rx_lp,1,&indices,&lu,&bd);
			if (status) {
				fprintf(stderr," LPCC_feasible_recover_by_rounding(): unable to tighten bound\n");
				goto TERMINATE;
			}
		}
	}
	status = cplexSolveLP(env,copy_rx_lp,rx_lp_start,param_n+param_m+param_m+param_k,recover_x,recover_objval,recover_solnstat,&tmp,CPX_ALG_AUTOMATIC,1,2);
	if (status) goto TERMINATE;
TERMINATE:
	free_startinfo(&tmp);
	if ( copy_rx_lp != NULL ) {
		status = CPXfreeprob (env, &copy_rx_lp);
		if ( status ) {  
			fprintf (stderr, " LPCC_feasible_recover_by_rounding(): CPXfreeprob failed\n");
		}
	}
	return (status);
}

int LPCC_feasible_recover_M_psd_fix_x(CPXENVptr env, 
									  CPXLPptr rx_lp, 
									  const int param_n,
									  const int param_m, 
									  const int param_k,
									  const int num_var,
									  MATRIX matrix_N, 
									  MATRIX matrix_M, 
									  double *q_coef,
									  NODE active_node,
									  double *x_soln,
									  double *recover_objval, 
									  double *recover_x, 
									  int *recover_solnstat)
{
	int status =0;
	int i,j;
	int matrix_M_status;


	int lp_num_col;
	int lp_num_row;

	int lp_solnstat;
	int indices;
	int *qmatcnt=NULL;
	int *qmatbeg=NULL;
	int *qmatind=NULL;
	char lu;

	double *qmatval=NULL;
	double tmp_val;
	double *y_tilde_lb=NULL;
	CPXLPptr copy_rx_lp = NULL;
	MATRIX matrix_M_bar;
	MATRIX matrix_N_t;
	*recover_solnstat=0;
	init_matrix(&matrix_M_bar);
	init_matrix(&matrix_N_t);
	status=copy_matrix(&matrix_N,&matrix_N_t);
	if(status) goto TERMINATE;
	status=transpose_matrix(&matrix_N_t);
	if(status) goto TERMINATE;
	status=symmetrizing_matrix(&matrix_M,&matrix_M_bar,&matrix_M_status);
	if(status) goto TERMINATE;
	copy_rx_lp = CPXcloneprob(env,rx_lp,&status);
	if (status) {
		fprintf(stderr," LPCC_feasible_recover_by_rounding(): unable to copy rx_lp \n");
		goto TERMINATE;
	}
	status=setupnode_rx_lp(env,copy_rx_lp,param_n,param_m,param_k,active_node);

	status=zeroobjectivefunc(env,copy_rx_lp);
	if(status) goto TERMINATE;
	//add quadratic term
	lp_num_col=CPXgetnumcols(env,copy_rx_lp);
	qmatbeg=(int*) calloc(lp_num_col,sizeof(int));
	qmatcnt=(int*) calloc(lp_num_col,sizeof(int));
	qmatind=(int*) calloc(matrix_M_bar.nnz,sizeof(int));
	qmatval=(double*) calloc(matrix_M_bar.nnz,sizeof(double));
	if (qmatbeg==NULL || qmatcnt==NULL || qmatind==NULL || qmatval==NULL)
	{
		status=-1;
		fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<lp_num_col;i++)
	{
		if (i<param_n)qmatcnt[i]=0;
		else if (i>=param_n && i<param_n+param_m)
		{
			qmatcnt[i]=matrix_M_bar.matcnt[i-param_n];
		}else qmatcnt[i]=0;
		if (i==0) qmatbeg[i]=0;
		else qmatbeg[i]=qmatbeg[i-1]+qmatcnt[i-1];
	}
	for (i=0;i<matrix_M_bar.nnz;i++)
	{
		qmatval[i]=matrix_M_bar.matval[i];
		qmatind[i]=matrix_M_bar.matind[i]+param_n;
	}
	status=CPXcopyquad(env,copy_rx_lp,qmatbeg,qmatcnt,qmatind,qmatval);
	if(status) goto TERMINATE;	
	//add constraint y_tilde=Nt*y
	y_tilde_lb=(double*) malloc(param_n*sizeof(double));
	if(y_tilde_lb==NULL)
	{
		status=-1;
		fprintf(stderr," unable to allocate memory\n");
		goto TERMINATE;
	}
	
	for(i=0;i<param_n;i++)
	{
		y_tilde_lb[i]=-CPX_INFBOUND;
	}
	lp_num_col=CPXgetnumcols(env,copy_rx_lp);
	lp_num_row=CPXgetnumrows(env,copy_rx_lp);
	status=CPXnewcols(env,copy_rx_lp,param_n,NULL,y_tilde_lb,NULL,NULL,NULL);
	if(status) goto TERMINATE;
	status=CPXnewrows(env,copy_rx_lp,param_n,NULL,NULL,NULL,NULL);
	if (status) goto TERMINATE;
	for (i=0;i<param_n;i++)
	{
		for (j=0;j<matrix_N_t.matcnt[i];j++)
		{
			status=CPXchgcoef(env,copy_rx_lp,lp_num_row+i,param_n+matrix_N_t.matind[matrix_N_t.matbeg[i]+j],matrix_N_t.matval[matrix_N_t.matbeg[i]+j]);
			if (status) goto TERMINATE;
		}
		status=CPXchgcoef(env,copy_rx_lp,lp_num_row+i,lp_num_col+i,-1.0);
		if (status) goto TERMINATE;
	}
	
	for (i=0;i<param_n;i++)
	{
		lu='U';
		tmp_val=x_soln[i];
		//fprintf(stdout,"tmp_val:%f\n",tmp_val);
		indices=i;
		status=CPXtightenbds(env,copy_rx_lp,1,&indices,&lu,&tmp_val);
		if(status) goto TERMINATE;
		lu='L';
		status=CPXtightenbds(env,copy_rx_lp,1,&indices,&lu,&tmp_val);
		if(status) goto TERMINATE;
		indices=lp_num_col+i;
		status=CPXchgcoef(env,copy_rx_lp,-1,indices,tmp_val/2);
		if(status) goto TERMINATE;
	}

	for (i=0;i<param_m;i++)
	{
		status=CPXchgcoef(env,copy_rx_lp,-1,param_n+i,q_coef[i]/2);
		if(status) goto TERMINATE;
	}

	status = CPXlpopt (env, copy_rx_lp);
	if (status==CPXERR_Q_NOT_POS_DEF)
	{
		fprintf(stdout," Matrix M is not psd\n");
		*recover_solnstat=-1;
		status=0;
		goto TERMINATE;
	}else if(status) goto TERMINATE;
	status = CPXsolution (env, copy_rx_lp, &lp_solnstat, &tmp_val, NULL, NULL, NULL, NULL);
	if ( status ) goto TERMINATE;
	status=CPXgetx(env,copy_rx_lp,recover_x,0,num_var-1);
TERMINATE:

	if ( copy_rx_lp != NULL ) {
		status = CPXfreeprob (env, &copy_rx_lp);
		if ( status ) {  
			fprintf (stderr, " LPCC_feasible_recover_by_rounding(): CPXfreeprob failed\n");
		}
	}

	free_and_null((char**) &qmatbeg);
	free_and_null((char**) &qmatcnt);
	free_and_null((char**) &qmatind);
	free_and_null((char**) &qmatval);
	free_and_null((char**) &y_tilde_lb);
	free_matrix(&matrix_N_t);
	free_matrix(&matrix_M_bar);
	return (status);
}

int LPCC_feasible_recover_M_psd_fix_y_bar(CPXENVptr env, 
										  CPXLPptr rx_lp, 
										  const int param_n,
										  const int param_m, 
										  const int param_k,
										  const int num_var,
										  MATRIX matrix_N, 
										  MATRIX matrix_M, 
										  double *q_coef,
										  NODE active_node, 
										  double *x_soln,
										  double *recover_objval, 
										  double *recover_x, 
										  int *recover_solnstat)
{
	int status =0;
	int i,j;
	int matrix_M_status;


	int lp_num_col;
	int lp_num_row;

	int lp_solnstat;
	int indices;
	int *qmatcnt=NULL;
	int *qmatbeg=NULL;
	int *qmatind=NULL;
	char lu;

	double *qmatval=NULL;
	double tmp_val;
	double *y_tilde_lb=NULL;
	CPXLPptr copy_rx_lp = NULL;
	MATRIX matrix_M_bar;
	MATRIX matrix_N_t;
	*recover_solnstat=0;
	init_matrix(&matrix_M_bar);
	init_matrix(&matrix_N_t);
	status=copy_matrix(&matrix_N,&matrix_N_t);
	if(status) goto TERMINATE;
	status=transpose_matrix(&matrix_N_t);
	if(status) goto TERMINATE;
	status=symmetrizing_matrix(&matrix_M,&matrix_M_bar,&matrix_M_status);
	if(status) goto TERMINATE;
	copy_rx_lp = CPXcloneprob(env,rx_lp,&status);
	if (status) {
		fprintf(stderr," LPCC_feasible_recover_by_rounding(): unable to copy rx_lp \n");
		goto TERMINATE;
	}
	status=setupnode_rx_lp(env,copy_rx_lp,param_n,param_m,param_k,active_node);

	status=zeroobjectivefunc(env,copy_rx_lp);
	if(status) goto TERMINATE;

	lp_num_col=CPXgetnumcols(env,copy_rx_lp);
	qmatbeg=(int*) calloc(lp_num_col,sizeof(int));
	qmatcnt=(int*) calloc(lp_num_col,sizeof(int));
	qmatind=(int*) calloc(matrix_M_bar.nnz,sizeof(int));
	qmatval=(double*) calloc(matrix_M_bar.nnz,sizeof(double));
	if (qmatbeg==NULL || qmatcnt==NULL || qmatind==NULL || qmatval==NULL)
	{
		status=-1;
		fprintf(stderr," add_mccormickcut(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<lp_num_col;i++)
	{
		if (i<param_n)qmatcnt[i]=0;
		else if (i>=param_n && i<param_n+param_m)
		{
			qmatcnt[i]=matrix_M_bar.matcnt[i-param_n];
		}else qmatcnt[i]=0;
		if (i==0) qmatbeg[i]=0;
		else qmatbeg[i]=qmatbeg[i-1]+qmatcnt[i-1];
	}
	for (i=0;i<matrix_M_bar.nnz;i++)
	{
		qmatval[i]=matrix_M_bar.matval[i];
		qmatind[i]=matrix_M_bar.matind[i]+param_n;
	}
	status=CPXcopyquad(env,copy_rx_lp,qmatbeg,qmatcnt,qmatind,qmatval);
	if(status) goto TERMINATE;
	//add constraint y_tilde=Nt*y
	y_tilde_lb=(double*) malloc(param_n*sizeof(double));
	if(y_tilde_lb==NULL)
	{
		status=-1;
		fprintf(stderr," unable to allocate memory\n");
		goto TERMINATE;
	}

	for(i=0;i<param_n;i++)
	{
		y_tilde_lb[i]=-CPX_INFBOUND;
	}
	lp_num_col=CPXgetnumcols(env,copy_rx_lp);
	lp_num_row=CPXgetnumrows(env,copy_rx_lp);
	status=CPXnewcols(env,copy_rx_lp,param_n,NULL,y_tilde_lb,NULL,NULL,NULL);
	if(status) goto TERMINATE;
	status=CPXnewrows(env,copy_rx_lp,param_n,NULL,NULL,NULL,NULL);
	if (status) goto TERMINATE;
	for (i=0;i<param_n;i++)
	{
		for (j=0;j<matrix_N_t.matcnt[i];j++)
		{
			status=CPXchgcoef(env,copy_rx_lp,lp_num_row+i,param_n+matrix_N_t.matind[matrix_N_t.matbeg[i]+j],matrix_N_t.matval[matrix_N_t.matbeg[i]+j]);
			if (status) goto TERMINATE;
		}
		status=CPXchgcoef(env,copy_rx_lp,lp_num_row+i,lp_num_col+i,-1.0);
		if (status) goto TERMINATE;
	}

	for (i=0;i<param_n;i++)
	{
		tmp_val=0;
		for (j=0;j<matrix_N_t.matcnt[i];j++)
		{
			tmp_val+=matrix_N_t.matval[matrix_N_t.matbeg[i]+j]*x_soln[param_n+matrix_N_t.matind[matrix_N_t.matbeg[i]+j]];
		}
		lu='U';
		indices=lp_num_col+i;
		status=CPXtightenbds(env,copy_rx_lp,1,&indices,&lu,&tmp_val);
		if(status) goto TERMINATE;
		lu='L';
		status=CPXtightenbds(env,copy_rx_lp,1,&indices,&lu,&tmp_val);
		if(status) goto TERMINATE;
		indices=i;
		status=CPXchgcoef(env,copy_rx_lp,-1,indices,tmp_val/2);
		if(status) goto TERMINATE;
	}

	for (i=0;i<param_m;i++)
	{
		status=CPXchgcoef(env,copy_rx_lp,-1,param_n+i,q_coef[i]/2);
		if(status) goto TERMINATE;
	}

	status = CPXlpopt (env, copy_rx_lp);
	if (status==CPXERR_Q_NOT_POS_DEF)
	{
		fprintf(stdout," Matrix M is not psd\n");
		*recover_solnstat=-1;
		status=0;
		goto TERMINATE;
	}else if(status) goto TERMINATE;
	status = CPXsolution (env, copy_rx_lp, &lp_solnstat, &tmp_val, NULL, NULL, NULL, NULL);

	if ( status ) goto TERMINATE;
	//fprintf(stdout," solnstat; %d\n",lp_solnstat);
	//fprintf(stdout," gap: %f\n",tmp_val);
	status=CPXgetx(env,copy_rx_lp,recover_x,0,num_var-1);
TERMINATE:

	if ( copy_rx_lp != NULL ) {
		status = CPXfreeprob (env, &copy_rx_lp);
		if ( status ) {  
			fprintf (stderr, " LPCC_feasible_recover_by_rounding(): CPXfreeprob failed\n");
		}
	}

	free_and_null((char**) &qmatbeg);
	free_and_null((char**) &qmatcnt);
	free_and_null((char**) &qmatind);
	free_and_null((char**) &qmatval);
	free_and_null((char**) &y_tilde_lb);
	free_matrix(&matrix_N_t);
	free_matrix(&matrix_M_bar);
	return (status);
}

int LPCC_feasible_recover_by_fixing_variable(CPXENVptr env, 
											 CPXLPptr rx_lp, 
											 STARTINFO rx_lp_start,
											 const int param_n,
											 const int param_m, 
											 const int param_k,
											 const int num_var,
											 NODE active_node, 
											 double *pump_objective_coef,							
											 double *recover_objval, 								 
											 double *recover_x, 								 
											 int *recover_solnstat)
{
	int status =0;
	int indices;
	int i;
	int *activeset=active_node.nodebranchinfo;
	char lu;
	double bd;
	//double gap1;
	//double gap2;
	CPXLPptr copy_rx_lp = NULL;
	STARTINFO tmp;
	init_startinfo(&tmp);
	copy_rx_lp = CPXcloneprob(env,rx_lp,&status);
	if (status) {
		fprintf(stderr," LPCC_feasible_recover_by_fixing_variable(): unable to copy rx_lp \n");
		goto TERMINATE;
	}
	status=setupnode_rx_lp(env,copy_rx_lp,param_n,param_m,param_k,active_node);
	for (i=0; i< param_m; i++) {	
		if ( activeset[i] ==0) {	
			if (fabs(pump_objective_coef[i])<1e-6 && fabs(pump_objective_coef[param_m+i]-1)<1e-6)
			{
				indices=param_n+param_m+i;
			}else if (fabs(pump_objective_coef[i]-1)<1e-6 && fabs(pump_objective_coef[param_m+i])<1e-6)
			{
				indices=param_n+i;
			}else
			{
				status=-1;
				fprintf(stderr," LPCC_feasible_recover_by_fixing_variable(): unexpected pump_objective_coef: %f, %f\n",pump_objective_coef[i],pump_objective_coef[param_m+i]);
				goto TERMINATE;
			}
			lu='U';
			bd=0.0;
			status=CPXtightenbds(env,copy_rx_lp,1,&indices,&lu,&bd);
			if (status) {
				fprintf(stderr," LPCC_feasible_recover_by_fixing_variable(): unable to tighten bound\n");
				goto TERMINATE;
			}
		}
	}
	status = cplexSolveLP(env,copy_rx_lp,rx_lp_start,param_n+param_m+param_m+param_k,recover_x,recover_objval,recover_solnstat,&tmp,CPX_ALG_AUTOMATIC,1,2);
	if (status) goto TERMINATE;
TERMINATE:
	free_startinfo(&tmp);
	if ( copy_rx_lp != NULL ) {
		status = CPXfreeprob (env, &copy_rx_lp);
		if ( status ) {  
			fprintf (stderr, " LPCC_feasible_recover_by_fixing_variable(): CPXfreeprob failed\n");
		}
	}
	return (status);
}
int feasiblity_local_search(CPXENVptr env, 
						    CPXLPptr rx_lp, 	
						    STARTINFO rx_lp_start,								 
						    const int param_n,								 
						    const int param_m, 								 
						    const int param_k,								 
						    NODE active_node,		
							const int initial_index,
							double *initial_pump_objective_coef,
						    const double initial_gap,
						    const int num_var,								 
							double *feasible_pump_objective_coef,
						    int *local_search_stat,
							int control_param)
{
	int status =0;
	int i,j,k;
	int selectedindex;
	int lp_solnstat;
	int repeat=0;
	int* indices=NULL;
	int limit;
	int *activeset=active_node.nodebranchinfo;
	double* pump_objective_coef=NULL;
	double* tmp_pump_objective_coef=NULL;
	double* lp_soln=NULL;
	double lp_obj;
	double sum_gap;
	CPXLPptr copy_rx_lp = NULL;
	STARTINFO lp_start_1,lp_start_2;
	indices=(int*)malloc(2*param_m*sizeof(int));
	pump_objective_coef=(double*)calloc(2*param_m,sizeof(double));
	tmp_pump_objective_coef=(double*)calloc(2*param_m,sizeof(double));
	lp_soln=(double*)malloc(num_var*sizeof(double));
	if (indices==NULL || pump_objective_coef==NULL || tmp_pump_objective_coef==NULL || lp_soln==NULL)
	{
		status=-1;
		fprintf(stderr," feasiblity_local_search(): unable to allocate memory\n");
		goto TERMINATE;
	}
	init_startinfo(&lp_start_1);
	init_startinfo(&lp_start_2);
	copy_rx_lp = CPXcloneprob(env,rx_lp,&status);
	if (status) {
		fprintf(stderr," feasiblity_local_search(): unable to copy rx_lp \n");
		goto TERMINATE;
	}
	status=copy_startinfo(&rx_lp_start,&lp_start_1);
	if(status) goto TERMINATE;
	for (i=0; i< param_m; i++) {	
		indices[i]=param_n+i;
		indices[param_m+i]=param_n+param_m+i;
		pump_objective_coef[i]=initial_pump_objective_coef[i];
		pump_objective_coef[param_m+i]=initial_pump_objective_coef[param_m+i];
	}
	status=zeroobjectivefunc(env,copy_rx_lp);
	if(status) goto TERMINATE;
	sum_gap=initial_gap;
	if (fabs(sum_gap)<ZERO_TOLERANCE)
	{
		status=-1;
		fprintf(stderr," feasibility_local_search(): initial_gap %f is too small",initial_gap);
		goto TERMINATE;
	}else
	{
		//start heuristic method to find feasible solution
		repeat=1;
		selectedindex=initial_index;
		limit=0;
		while (repeat)
		{
			if (control_param!=-1 && limit>control_param)break;
			limit++;
			repeat=0;
			for (i=0;i<param_m;i++)
			{
				if (activeset[i]==0 && i!=selectedindex)
				{
					for (j=0;j<param_m;j++)
					{
						if (j==i)
						{
							tmp_pump_objective_coef[j]=1-pump_objective_coef[j];
							tmp_pump_objective_coef[param_m+j]=1-pump_objective_coef[param_m+j];
						} 
						else
						{
							tmp_pump_objective_coef[j]=pump_objective_coef[j];
							tmp_pump_objective_coef[param_m+j]=pump_objective_coef[param_m+j];
						}

					}
					status=CPXchgobj(env,copy_rx_lp,2*param_m,indices,tmp_pump_objective_coef);
					if(status) goto TERMINATE;
					status=cplexSolveLP(env,copy_rx_lp,lp_start_1,num_var,lp_soln,&lp_obj,&lp_solnstat,&lp_start_2,CPX_ALG_AUTOMATIC,1,2);
					if(status) goto TERMINATE;
					if (lp_solnstat!=CPX_STAT_OPTIMAL)
					{
						status=-1;
						fprintf(stderr," feasiblity_local_search(): unexpected lp solnstat: %d",lp_solnstat);
					}
					status=copy_startinfo(&lp_start_2,&lp_start_1);
					if(status) goto TERMINATE;
					if (fabs(lp_obj)<ZERO_TOLERANCE)
					{
						for (k=0;k<2*param_m;k++)
						{
							feasible_pump_objective_coef[k]=tmp_pump_objective_coef[k];
						}
						*local_search_stat=1;
						goto TERMINATE;
					}else if (sum_gap>lp_obj)
					{
						//get a better sum_gap
						sum_gap=lp_obj;
						selectedindex=i;
						pump_objective_coef[selectedindex]=1-pump_objective_coef[selectedindex];
						pump_objective_coef[param_m+selectedindex]=1-pump_objective_coef[param_m+selectedindex];
						//fprintf(stdout,"feasible gap: %f\n",sum_gap);
						repeat=1;
					}
				}
			}
		}
		*local_search_stat=0;
		//fprintf(stdout," recovery process failed\n");
	}
TERMINATE:
	free_and_null((char**)&indices);
	free_and_null((char**)&pump_objective_coef);
	free_and_null((char**)&tmp_pump_objective_coef);
	free_and_null((char**)&lp_soln);
	free_startinfo(&lp_start_1);
	free_startinfo(&lp_start_2);
	if ( copy_rx_lp != NULL ) {
		status = CPXfreeprob (env, &copy_rx_lp);
		if ( status ) {  
			fprintf (stderr, " LPCC_feasible_recover_by_ray(): CPXfreeprob failed\n");
		}
	}

	return (status);
}

int LPCC_feasible_recovery(CPXENVptr env, 
						   CPXLPptr rx_lp, 	
						   STARTINFO rx_lp_start,								 
						   const int param_n,								 
						   const int param_m, 								 
						   const int param_k,								 
						   NODE active_node,	
						   const double recovery_range_lb,
						   const double recovery_range_ub,
						   const int num_var,
						   double *initial_soln,								 
						   double *recover_objval, 
						   double *recover_soln,
						   int *recover_stat,
						   int control_param_1,
						   int control_param_2)
{
	//*recover_stat=0: recovery failed
	//*recover_stat=1: recovery succeed
	int status =0;
	int i,j;
	int lp_solnstat;
	int local_search_stat;
	int roundingstat;
	int cnt;
	int tmp_index;
	int switch_index;
	int *activeset=active_node.nodebranchinfo;
	int* indices=NULL;
	int* sorted_gap_array_index=NULL;
	double* pump_objective_coef=NULL;
	double* tmp_pump_objective_coef=NULL;
	double* lp_soln=NULL;
	double* gap_array=NULL;
	double lp_obj;
	double sum_gap;
	double y_gap;
	double w_gap;
	CPXLPptr copy_rx_lp = NULL;
	STARTINFO lp_start_1,lp_start_2;
	if (control_param_1<1)
	{
		control_param_1=1;
	}
	indices=(int*)malloc(2*param_m*sizeof(int));
	pump_objective_coef=(double*)calloc(2*param_m,sizeof(double));
	tmp_pump_objective_coef=(double*)calloc(2*param_m,sizeof(double));
	lp_soln=(double*)malloc(num_var*sizeof(double));
	gap_array=(double*)malloc(param_m*sizeof(double));
	sorted_gap_array_index=(int*)malloc(param_m*sizeof(int));
	if (indices==NULL || pump_objective_coef==NULL || tmp_pump_objective_coef==NULL || lp_soln==NULL || gap_array==NULL || sorted_gap_array_index==NULL)
	{
		status=-1;
		fprintf(stderr," LPCC_feasible_recovery(): unable to allocate memory\n");
		goto TERMINATE;
	}
	//fprintf(stdout," Feasible recovery...\n");
	init_startinfo(&lp_start_1);
	init_startinfo(&lp_start_2);
	copy_rx_lp = CPXcloneprob(env,rx_lp,&status);
	if (status) {
		fprintf(stderr," LPCC_feasible_recovery(): unable to copy rx_lp \n");
		goto TERMINATE;
	}
	status=setupnode_rx_lp(env,copy_rx_lp,param_n,param_m,param_k,active_node);
	if(status) goto TERMINATE;
	for (i=0; i< param_m; i++) {	
		indices[i]=param_n+i;
		indices[param_m+i]=param_n+param_m+i;
		gap_array[i]=-1;
		sorted_gap_array_index[i]=-1;
		switch(activeset[i])
		{
		case 0:
			y_gap = initial_soln[param_n+i];
			w_gap =initial_soln[param_n+param_m+i] ;
			if (y_gap<=w_gap){
				pump_objective_coef[i]=1;
				pump_objective_coef[param_m+i]=0;
			} else {		
				pump_objective_coef[i]=0;
				pump_objective_coef[param_m+i]=1;
			}
			break;
		case 1:
			pump_objective_coef[i]=1;
			pump_objective_coef[param_m+i]=0;
			break;
		case 2:
			pump_objective_coef[i]=0;
			pump_objective_coef[param_m+i]=1;
			break;
		case 3:
			//either y or w is zero, and here we choose y
			pump_objective_coef[i]=1;
			pump_objective_coef[param_m+i]=0;
			break;
		default:
			status=-1;
			fprintf(stderr," LPCC_feasible_recovery(): unexpected activeset value:%d\n",activeset[i]);
			goto TERMINATE;
		}
	}
	if (recovery_range_lb>-CPX_INFBOUND)
	{
		status=add_lb_constraint(env,copy_rx_lp,recovery_range_lb);
		if(status) goto TERMINATE;
	}
	if (recovery_range_ub<CPX_INFBOUND)
	{
		status=add_ub_constraint(env,copy_rx_lp,recovery_range_ub);
		if(status) goto TERMINATE;
	}
	//fprintf(stdout,"add lb:%f, add ub:%f\n",recovery_range_lb,recovery_range_ub);
	status=zeroobjectivefunc(env,copy_rx_lp);
	if(status) goto TERMINATE;

	
	status=CPXchgobj(env,copy_rx_lp,2*param_m,indices,pump_objective_coef);
	if(status) goto TERMINATE;
	status=cplexSolveLP(env,copy_rx_lp,rx_lp_start,num_var,lp_soln,&lp_obj,&lp_solnstat,&lp_start_1,CPX_ALG_AUTOMATIC,1,2);
	if (status) goto TERMINATE;
	if (lp_solnstat!=CPX_STAT_OPTIMAL) 
	{
		status=-1;
		fprintf(stderr," LPCC_feasible_recovery(): unexpected lp_solnstat: %d",lp_solnstat);
		goto TERMINATE;
	}
	if (fabs(lp_obj)<ZERO_TOLERANCE)
	{
		//recover a feasible solution
		status=LPCC_feasible_recover_by_fixing_variable(env,rx_lp,rx_lp_start,param_n,param_m,param_k,num_var,active_node,pump_objective_coef,recover_objval,recover_soln,&roundingstat);
		if(status) goto TERMINATE;
		if (roundingstat!=CPX_STAT_OPTIMAL)
		{
			status=-1;
			fprintf(stderr," LPCC_feasible_recovery(): unexpected rounding status:%d\n",roundingstat);
			goto TERMINATE;
		}
		*recover_stat=1;
		fprintf(stdout,"	Find a feasible solution with objective:%f by directly rounding\n",*recover_objval);
		goto TERMINATE;
	}else
	{
		sum_gap=lp_obj;
		//fprintf(stdout," initial feasible gap: %f\n",sum_gap);
		//start heuristic method to find feasible solution
		cnt=0;
		for (i=0;i<param_m;i++)
		{
			if (activeset[i]==0)
			{
				for (j=0;j<param_m;j++)
				{
					if (j==i)
					{
						tmp_pump_objective_coef[j]=1-pump_objective_coef[j];
						tmp_pump_objective_coef[param_m+j]=1-pump_objective_coef[param_m+j];
					} 
					else
					{
						tmp_pump_objective_coef[j]=pump_objective_coef[j];
						tmp_pump_objective_coef[param_m+j]=pump_objective_coef[param_m+j];
					}

				}
				status=CPXchgobj(env,copy_rx_lp,2*param_m,indices,tmp_pump_objective_coef);
				if(status) goto TERMINATE;
				status=cplexSolveLP(env,copy_rx_lp,lp_start_1,num_var,lp_soln,&lp_obj,&lp_solnstat,&lp_start_2,CPX_ALG_AUTOMATIC,1,2);
				if(status) goto TERMINATE;
				status=copy_startinfo(&lp_start_2,&lp_start_1);
				if(status) goto TERMINATE;
				if (fabs(lp_obj)<ZERO_TOLERANCE)
				{
					//recover a feasible solution
					status=LPCC_feasible_recover_by_fixing_variable(env,rx_lp,rx_lp_start,param_n,param_m,param_k,num_var,active_node,tmp_pump_objective_coef,recover_objval,recover_soln,&roundingstat);
					if(status) goto TERMINATE;
					if (roundingstat!=CPX_STAT_OPTIMAL)
					{
						status=-1;
						fprintf(stderr," LPCC_feasible_recovery(): unexpected rounding status:%d\n",roundingstat);
						goto TERMINATE;
					}
					*recover_stat=1;
					fprintf(stdout,"	Find a feasible solution with objective:%f\n",*recover_objval);
					goto TERMINATE;
				}else if (sum_gap>lp_obj)
				{
					//get a better sum_gap
					cnt++;
					gap_array[i]=lp_obj;
					tmp_index=i;
					//fprintf(stdout,"cnt: %d, gap: %f\n",cnt,lp_obj);
					for (j=0;j<cnt;j++)
					{
						if (sorted_gap_array_index[j]<0)
						{
							sorted_gap_array_index[j]=tmp_index;
						}else
						{
							if (gap_array[sorted_gap_array_index[j]]>gap_array[tmp_index])
							{
								switch_index=sorted_gap_array_index[j];
								sorted_gap_array_index[j]=tmp_index;
								tmp_index=switch_index;
							}
						}
					}
				}
			}
		}
		//fprintf(stdout," improved feasible gap: (sorted)\n");
		//for (i=0;i<cnt;i++)
		//{
		//	fprintf(stdout,"cnt: %d, feasible gap: %f\n",i,gap_array[sorted_gap_array_index[i]]);
		//}
		//fprintf(stdout," further heuristic search by exploiting each improved feasible gap:\n");
		for (i=0;i<cnt;i++)
		{
			if (i>control_param_1)break;
			//fprintf(stdout,"cnt: %d, exploiting feasible gap: %f\n",i,gap_array[sorted_gap_array_index[i]]);
			pump_objective_coef[sorted_gap_array_index[i]]=1-pump_objective_coef[sorted_gap_array_index[i]];
			pump_objective_coef[param_m+sorted_gap_array_index[i]]=1-pump_objective_coef[param_m+sorted_gap_array_index[i]];
			feasiblity_local_search(env,copy_rx_lp,lp_start_1,param_n,param_m,param_k,active_node,sorted_gap_array_index[i],pump_objective_coef,gap_array[sorted_gap_array_index[i]],num_var,tmp_pump_objective_coef,&local_search_stat,control_param_2);
			if (local_search_stat==1)
			{
				//recover a feasible solution
				status=LPCC_feasible_recover_by_fixing_variable(env,rx_lp,rx_lp_start,param_n,param_m,param_k,num_var,active_node,tmp_pump_objective_coef,recover_objval,recover_soln,&roundingstat);
				if(status) goto TERMINATE;
				if (roundingstat!=CPX_STAT_OPTIMAL)
				{
					status=-1;
					fprintf(stderr," LPCC_feasible_recovery(): unexpected rounding status:%d\n",roundingstat);
					goto TERMINATE;
				}
				*recover_stat=1;
				fprintf(stdout,"	Find a feasible solution with objective:%f\n",*recover_objval);
				goto TERMINATE;
			}
			pump_objective_coef[sorted_gap_array_index[i]]=1-pump_objective_coef[sorted_gap_array_index[i]];
			pump_objective_coef[param_m+sorted_gap_array_index[i]]=1-pump_objective_coef[param_m+sorted_gap_array_index[i]];
		}
	}
	*recover_stat=0;
	//fprintf(stdout," recovery process failed\n");
TERMINATE:
	free_and_null((char**)&indices);
	free_and_null((char**)&pump_objective_coef);
	free_and_null((char**)&tmp_pump_objective_coef);
	free_and_null((char**)&lp_soln);
	free_startinfo(&lp_start_1);
	free_startinfo(&lp_start_2);
	if ( copy_rx_lp != NULL ) {
		status = CPXfreeprob (env, &copy_rx_lp);
		if ( status ) {  
			fprintf (stderr, " LPCC_feasible_recover_by_ray(): CPXfreeprob failed\n");
		}
	}

	return (status);
}

int fixing_direction(CPXENVptr env, 
				 CPXLPptr rx_lp, 
				 const int param_n,
				 const int param_m,
				 const int param_k,
				 const int violate_complementary_index,
				 const double violation_complementary_1, 
				 const double violation_complementary_2, 
				 double *direction)
{
	int status=0;
	const int num_cols=CPXgetnumcols(env,rx_lp);
	const int num_rows=CPXgetnumrows(env,rx_lp);
	int i;
	int row_id;
	int fix_check_1;
	int fix_check_2;
	int coef_1_pos_cnt;
	int coef_1_neg_cnt;
	int coef_1_zero_cnt;
	int coef_2_pos_cnt;
	int coef_2_neg_cnt;
	int coef_2_zero_cnt;
	double reduce_cost_sum_1;
	double reduce_cost_sum_2;
	double coef_1_max_value;
	double coef_2_max_value;
	double coef_1_geom;
	double coef_2_geom;
	double score_1;
	double score_2;
	double coef_1_norm_score;
	double coef_1_pos_norm_score;
	double coef_1_neg_norm_score;
	double coef_2_norm_score;
	double coef_2_pos_norm_score;
	double coef_2_neg_norm_score;
	double *val = NULL;
	double *tableau_1=NULL;
	double *tableau_2=NULL;
	int *cstat=NULL;
	int *rstat=NULL;
	double *reduce_cost=NULL;
	double coef_1;
	double coef_2;
	*direction=0;
	tableau_1 = (double*) malloc(num_cols*sizeof(double));
	tableau_2 = (double*) malloc(num_cols*sizeof(double));
	cstat=(int*) malloc(num_cols*sizeof(int));
	rstat=(int*) malloc(num_cols*sizeof(int));
	reduce_cost=(double*) malloc(num_cols*sizeof(double));
	val = (double*) malloc(num_cols*sizeof(double));
	if (tableau_1==NULL || tableau_2==NULL || cstat==NULL || rstat==NULL || reduce_cost==NULL || val== NULL)
	{
		status=NO_MEMORY;
		fprintf(stderr," branch_score():could not allocate memory.\n");
		goto TERMINATE;
	}
	status=CPXgetbase(env,rx_lp,cstat,rstat);
	if (status) goto TERMINATE;
	status = CPXgetijrow(env,rx_lp,-1,param_n+violate_complementary_index,&row_id);
	if (status) goto TERMINATE;
	status = CPXbinvarow(env,rx_lp,row_id,tableau_1);
	if (status) goto TERMINATE;
	status = CPXgetijrow(env,rx_lp,-1,param_n+param_m+violate_complementary_index,&row_id);
	if (status) goto TERMINATE;
	status = CPXbinvarow(env,rx_lp,row_id,tableau_2);
	if(status) goto TERMINATE;
	status=CPXgetdj(env,rx_lp,reduce_cost,0,num_cols-1);
	if(status) goto TERMINATE;
	coef_1_pos_cnt=0;
	coef_1_neg_cnt=0;
	coef_1_zero_cnt=0;
	coef_2_pos_cnt=0;
	coef_2_neg_cnt=0;
	coef_2_zero_cnt=0;
	coef_1_norm_score=0;
	coef_1_pos_norm_score=0;
	coef_1_neg_norm_score=0;
	coef_2_norm_score=0;
	coef_2_pos_norm_score=0;
	coef_2_neg_norm_score=0;
	coef_1_geom=1;
	coef_2_geom=1;
	reduce_cost_sum_1=0;
	reduce_cost_sum_2=0;
	coef_1_max_value=-CPX_INFBOUND;
	coef_2_max_value=-CPX_INFBOUND;
	fix_check_1=1;
	fix_check_2=1;
	for (i=0;i<num_cols;i++)
	{

		if (cstat[i]==CPX_BASIC) 
		{
			val[i]=0;
		}else
		{
			coef_1= tableau_1[i]/violation_complementary_1;
			coef_2= tableau_2[i]/violation_complementary_2;
			val[i]=coef_1>coef_2 ? coef_1:coef_2;
			//fprintf(stdout,"val: %f\n",val[i]);
			if (cstat[i]==CPX_AT_LOWER || cstat[i]==CPX_FREE_SUPER)
			{
				if (cstat[i]==CPX_AT_LOWER)
				{
					//fprintf(stdout,"coef_1: %f, coef_2: %f\n",coef_1,coef_2);
					coef_1_norm_score+=pow(coef_1,2);
					coef_2_norm_score+=pow(coef_2,2);
					if (coef_1>ZERO_TOLERANCE)
					{
						fix_check_1=0;
						coef_1_pos_norm_score+=pow(coef_1,2);
						coef_1_geom*=1/coef_1;
						coef_1_pos_cnt++;
						reduce_cost_sum_1+=reduce_cost[i];
						if (coef_1>coef_1_max_value)
						{
							coef_1_max_value=coef_1;
						}
					}else if (coef_1<-ZERO_TOLERANCE)
					{
						coef_1_neg_norm_score+=pow(coef_1,2);
						coef_1_neg_cnt++;
					}else
					{
						coef_1_zero_cnt++;
					}
					if (coef_2>ZERO_TOLERANCE)
					{
						fix_check_2=0;
						coef_2_pos_norm_score+=pow(coef_2,2);
						coef_2_geom*=1/coef_2;
						coef_2_pos_cnt++;
						reduce_cost_sum_2+=reduce_cost[i];
						if (coef_2>coef_2_max_value)
						{
							coef_2_max_value=coef_2;
						}
					}else if (coef_2<-ZERO_TOLERANCE)
					{
						coef_2_neg_norm_score+=pow(coef_2,2);
						coef_2_neg_cnt++;
					}else
					{
						coef_2_zero_cnt++;
					}
				}
				if (cstat[i]==CPX_FREE_SUPER)
				{
					coef_1_norm_score+=pow(coef_1,2);
					coef_2_norm_score+=pow(coef_2,2);
					if (fabs(coef_1)>ZERO_TOLERANCE)
					{
						fix_check_1=0;
						coef_1_pos_norm_score+=pow(coef_1,2);
						coef_1_neg_norm_score+=pow(coef_1,2);
						coef_1_geom*=1/fabs(coef_1);
						coef_1_pos_cnt++;
						coef_1_neg_cnt++;
						reduce_cost_sum_1+=reduce_cost[i];
						if (fabs(coef_1)>coef_1_max_value)
						{
							coef_1_max_value=fabs(coef_1);
						}

					}else
					{
						coef_1_zero_cnt++;
					}
					if (fabs(coef_2)>ZERO_TOLERANCE)
					{
						fix_check_2=0;
						coef_2_pos_norm_score+=pow(coef_2,2);
						coef_2_neg_norm_score+=pow(coef_2,2);
						coef_2_geom*=1/fabs(coef_2);
						coef_2_pos_cnt++;
						coef_2_neg_cnt++;
						reduce_cost_sum_2+=reduce_cost[i];
						if (fabs(coef_2)>coef_2_max_value)
						{
							coef_2_max_value=fabs(coef_2);
						}
					}else
					{
						coef_2_zero_cnt++;
					}
				}
			}
		}
	}
	if (fix_check_1 || fix_check_2)
	{
		if (fix_check_1)
		{
			*direction=2;
		}
		if (fix_check_2)
		{
			*direction=1;
		}
	}else
	{
		coef_1_pos_norm_score=coef_1_pos_norm_score>1e-6?coef_1_pos_norm_score:1e-6;
		coef_1_neg_norm_score=coef_1_neg_norm_score>1e-6?coef_1_neg_norm_score:1e-6;
		coef_2_pos_norm_score=coef_2_pos_norm_score>1e-6?coef_2_pos_norm_score:1e-6;
		coef_2_neg_norm_score=coef_2_neg_norm_score>1e-6?coef_2_neg_norm_score:1e-6;
		coef_1_max_value=coef_1_max_value>1e-6?coef_1_max_value:1e-6;
		coef_2_max_value=coef_2_max_value>1e-6?coef_2_max_value:1e-6;
		coef_1_norm_score=pow(coef_1_norm_score,0.5);
		coef_2_norm_score=pow(coef_2_norm_score,0.5);
		coef_1_pos_norm_score=pow(coef_1_pos_norm_score,0.5);
		coef_1_neg_norm_score=pow(coef_1_neg_norm_score,0.5);
		coef_2_pos_norm_score=pow(coef_2_pos_norm_score,0.5);
		coef_2_neg_norm_score=pow(coef_2_neg_norm_score,0.5);
		score_1=reduce_cost_sum_1/coef_1_pos_cnt/coef_1_max_value;
		score_2=reduce_cost_sum_2/coef_2_pos_cnt/coef_2_max_value;
		score_1=1/(coef_1_pos_norm_score)*violation_complementary_1;//*coef_1_norm_score*coef_1_norm_score
		score_2=1/(coef_2_pos_norm_score)*violation_complementary_2;//*coef_2_norm_score*coef_2_norm_score
		if (coef_1_pos_cnt==0 || coef_2_pos_cnt==0)
		{
			status=-1;
			goto TERMINATE;
		}
		
		*direction=score_1>score_2?1:2;
	}
	if(*direction==0)
	{
		status=-1;
		fprintf(stderr," fixing_direction(): unexpected direction\n");
		goto TERMINATE;
	}
TERMINATE:
	free_and_null((char**)&tableau_1);
	free_and_null((char**)&tableau_2);
	free_and_null((char**)&cstat);
	free_and_null((char**)&rstat);
	free_and_null((char**)&reduce_cost);
	free_and_null((char**)&val);
	return(status);
}

//
int readdata_format_1(char* fileDir, 
			 int *param_n_p,
			 int *param_m_p,
			 int *param_k_p, 
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
	int param_n,param_m,param_k;
	double *param_p=NULL;
	int i,n;
	char filename[100];
	FILE *in=NULL;
	//read c.txt
	strcpy(filename,fileDir);
	strcat(filename,"\\c.txt");
	in=fopen(filename,"r");
	if ( in == NULL ) 
	{
		status = -1;
		fprintf(stderr," unable to open the file c.txt\n");
		goto TERMINATE;
	}else
	{
		if ((status = readarray_row(in, &n,c_coef_p)))
		{
			fprintf(stderr," Input data error, please check c.txt\n");
			goto TERMINATE;
		}
		param_n=n;
		fclose(in);
	}


	//read d.txt
	strcpy(filename,fileDir);
	strcat(filename,"\\d.txt");
	in=fopen(filename,"r");
	if ( in == NULL ) 
	{
		status = -1;
		fprintf(stderr," unable to open the file d.txt\n");
		goto TERMINATE;
	}else
	{
		if ((status = readarray_row(in, &n,d_coef_p)))
		{
			fprintf(stderr," Input data error, please check d.txt\n");
			goto TERMINATE;
		}
		param_m=n;
		fclose(in);
	}


	//read b.txt
	strcpy(filename,fileDir);
	strcat(filename,"\\f.txt");
	in=fopen(filename,"r");
	if ( in == NULL ) 
	{
		status = -1;
		fprintf(stderr," unable to open the file f.txt\n");
		goto TERMINATE;
	}else
	{
		if ((status = readarray_col(in, &n,b_coef_p)))
		{
			fprintf(stderr," Input data error, please check f.txt\n");
			goto TERMINATE;
		}
		param_k=n;
		fclose(in);
	}

	

	//read q.txt
	strcpy(filename,fileDir);
	strcat(filename,"\\q.txt");
	in=fopen(filename,"r");
	if ( in == NULL ) 
	{
		status = -1;
		fprintf(stderr," unable to open the file q.txt\n");
		goto TERMINATE;
	}else
	{
		if ((status = readarray_col(in, &n,q_coef_p)))
		{
			fprintf(stderr," Input data error, please check q.txt\n");
			goto TERMINATE;
		}
		if(n!=param_m)
		{
			fprintf(stderr," Input data error, please check q.txt\n");
			status=-1;
			goto TERMINATE;
		}
		fclose(in);
	}

	*param_n_p=param_n;
	*param_m_p=param_m;
	*param_k_p=param_k;
	*matrix_A_p=(double**)calloc(param_k,sizeof(double*));
	*matrix_B_p=(double**)calloc(param_k,sizeof(double*));
	*matrix_N_p=(double**)calloc(param_m,sizeof(double*));
	*matrix_M_p=(double**)calloc(param_m,sizeof(double*));
	if (*matrix_A_p==NULL || *matrix_B_p==NULL || *matrix_N_p==NULL || *matrix_M_p==NULL)
	{
		status = NO_MEMORY;
		fprintf(stderr," Unable to allocate enough memory for the input data...\n");
		goto TERMINATE;
	}
	//read A.txt
	fprintf(stdout,"A:\n");
	strcpy(filename,fileDir);
	strcat(filename,"\\A.txt");
	in=fopen(filename,"r");
	if ( in == NULL ) 
	{
		status = -1;
		fprintf(stderr," unable to open the file A.txt\n");
		goto TERMINATE;
	}else
	{
		i=0;
		while (!feof(in))
		{
			if ( (status = readarray_row(in, &n, (*matrix_A_p)+i)) ) 
			{
				fprintf(stderr," Input data error, please check A.txt\n");
				goto TERMINATE;
			}
			i++;
			if (n!=param_n)
			{
				status=-1;
				fprintf(stderr," Input data error, please check A.txt\n");
				goto TERMINATE;
			}
		}
		if (i!=param_k)
		{
			status=-1;
			fprintf(stderr," Input data error, please check A.txt\n");
			goto TERMINATE;
		}
		fclose(in);
	}
	

	//read B.txt
	fprintf(stdout,"B:\n");
	strcpy(filename,fileDir);
	strcat(filename,"\\B.txt");
	in=fopen(filename,"r");
	if ( in == NULL ) 
	{
		status = -1;
		fprintf(stderr," unable to open the file B.txt\n");
		goto TERMINATE;
	}else
	{
		i=0;
		while (!feof(in))
		{
			if ( (status = readarray_row(in, &n, (*matrix_B_p)+i)) ) 
			{
				fprintf(stderr," Input data error, please check B.txt\n");
				goto TERMINATE;
			}
			i++;
			if (n!=param_m)
			{
				status=-1;
				fprintf(stderr," Input data error, please check B.txt\n");
				goto TERMINATE;
			}
		}
		if (i!=param_k)
		{
			status=-1;
			fprintf(stderr," Input data error, please check B.txt\n");
			goto TERMINATE;
		}
		fclose(in);
	}

	
	//read N.txt
	fprintf(stdout,"N:\n");
	strcpy(filename,fileDir);
	strcat(filename,"\\N.txt");
	in=fopen(filename,"r");
	if ( in == NULL ) 
	{
		status = -1;
		fprintf(stderr," unable to open the file N.txt\n");
		goto TERMINATE;
	}else
	{
		i=0;
		while (!feof(in))
		{
			if ( (status = readarray_row(in, &n, (*matrix_N_p)+i)) ) 
			{
				fprintf(stderr," Input data error, please check N.txt\n");
				goto TERMINATE;
			}
			i++;
			if ( n != param_n ) 
			{	
				status = -1;	
				fprintf(stderr," Input data error, please check N.txt\n");
				goto TERMINATE;
			}
		}
		if (i!=param_m)
		{
			status = -1;
			fprintf(stderr," Input data error, please check N.txt\n");
			goto TERMINATE;
		}
		fclose(in);
	}
	
	//read M.txt
	fprintf(stdout,"M:\n");
	strcpy(filename,fileDir);
	strcat(filename,"\\M.txt");
	in=fopen(filename,"r");
	if ( in == NULL ) 
	{
		status = -1;
		fprintf(stderr," unable to open the file M.txt\n");
		goto TERMINATE;
	}else
	{
		i=0;
		while (!feof(in))
		{
			if ( (status = readarray_row(in, &n, (*matrix_M_p)+i)) ) 
			{
				fprintf(stderr," Input data error, please check M.txt\n");
				goto TERMINATE;
			}
			i++;
			if ( n != param_m) 
			{	
				status = -1;	
				fprintf(stderr," Input data error, please check M.txt\n");
				goto TERMINATE;
			}
		}
		if (i!=param_m)
		{
			status = -1;
			fprintf(stderr," Input data error, please check M.txt\n");
			goto TERMINATE;
		}
		fclose(in);
	}
	


TERMINATE:
	return (status);

}

int readdata_format_2(char* file, 
					  int *param_n_p,
					  int *param_m_p,
					  int *param_k_p, 
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
	int param_n,param_m,param_k;
	double *param_p=NULL;
	int i,n;
	char ch;
	FILE *in=NULL;
	in=fopen(file,"r");
	if ( in == NULL ) 
	{
		status = -1;
		fprintf(stderr," unable to open the file... \n");
		goto TERMINATE;
	}
	/*read parameter*/
	if ( (status = readarray(in, &n,&param_p)) ) 
	{
		fprintf(stderr," Input data error, please check the input data again...\n");
		goto TERMINATE;
	}
	if ( n != 3 )
	{
		status = -1;
		fprintf(stderr," Input data error, please check the input data again...\n");
		goto TERMINATE;
	}
	param_n=(int) floor(param_p[0]);
	param_m=(int) floor(param_p[1]);
	param_k=(int) floor(param_p[2]);
	*param_n_p=param_n;
	*param_m_p=param_m;
	*param_k_p=param_k;
	*matrix_A_p=(double**)calloc(param_k,sizeof(double*));
	*matrix_B_p=(double**)calloc(param_k,sizeof(double*));
	*matrix_N_p=(double**)calloc(param_m,sizeof(double*));
	*matrix_M_p=(double**)calloc(param_m,sizeof(double*));
	if (*matrix_A_p==NULL || *matrix_B_p==NULL || *matrix_N_p==NULL || *matrix_M_p==NULL)
	{
		status = NO_MEMORY;
		fprintf(stderr," Unable to allocate enough memory for the input data...\n");
		goto TERMINATE;
	}

	//
	/*read c_coef*/
	if ( (status = readarray(in, &n, c_coef_p)) ) 
	{
		fprintf(stderr," Input data error, please check the input data again...\n");
		goto TERMINATE;
	}
	if ( n != param_n) {
		status = -1;
		fprintf(stderr," Input data error, please check the input data again...\n");
		goto TERMINATE;
	}
	/*read d_coef*/
	if ( (status = readarray(in, &n, d_coef_p)) ) 
	{
		fprintf(stderr," Input data error, please check the input data again...\n");
		goto TERMINATE;
	}
	if ( n != param_m) {
		status = -1;
		fprintf(stderr," Input data error, please check the input data again...\n");
		goto TERMINATE;
	}
	/*read b_coef*/
	if ( (status = readarray(in, &n, b_coef_p)) ) 
	{
		fprintf(stderr," Input data error, please check the input data again...\n");
		goto TERMINATE;
	}
	if ( n != param_k) {
		status = -1;
		fprintf(stderr," Input data error, please check the input data again...\n");
		goto TERMINATE;
	}
	/*read q_coef*/
	if ( (status = readarray(in, &n,q_coef_p)) ) 
	{
		fprintf(stderr," Input data error, please check the input data again...\n");
		goto TERMINATE;
	}
	if ( n != param_m) {
		status = -1;
		fprintf(stderr," Input data error, please check the input data again...\n");
		goto TERMINATE;
	}
	/*read matrix_A*/
	for (;;)
	{
		fscanf (in, "%c", &ch);
		if ( ch == '\t' ||
			ch == '\r' ||
			ch == ' '  ||
			ch == '\n'   ) continue;
		if ( ch == '[' ) break;
		status = -1;
		fprintf(stderr," Input data error, please check the input data again...\n");
		goto TERMINATE;
	}

	for ( i = 0; i < param_k; i++ )
	{
		if ( (status = readarray(in, &n, (*matrix_A_p)+i)) ) 
		{
			fprintf(stderr," Input data error, please check the input data again...\n");
			goto TERMINATE;
		}
		if ( n != param_n ) 
		{	
			status = -1;	
			fprintf(stderr," Input data error, please check the input data again...\n");
			goto TERMINATE;
		}
		do {
			fscanf (in, "%c", &ch);
		} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');

		if ( ch == ']' ) 
		{
			if(i!=param_k-1)
			{
				status=-1;
				fprintf(stderr," Input data error, please check the input data again...\n");
				goto TERMINATE;
			}
		}else if ( ch != ',' ) {
			status = -1;
			goto TERMINATE;
		}
	}
	/*read matrix_B*/
	for (;;)
	{
		fscanf (in, "%c", &ch);
		if ( ch == '\t' ||
			ch == '\r' ||
			ch == ' '  ||
			ch == '\n'   ) continue;
		if ( ch == '[' ) break;
		status = -1;
		fprintf(stderr," Input data error, please check the input data again...\n");
		goto TERMINATE;
	}

	for ( i = 0; i < param_k; i++ )
	{
		if ( (status = readarray(in, &n, (*matrix_B_p)+i)) ) 
		{
			fprintf(stderr," Input data error, please check the input data again...\n");
			goto TERMINATE;
		}
		if ( n != param_m) 
		{	
			status = -1;	
			fprintf(stderr," Input data error, please check the input data again...\n");
			goto TERMINATE;
		}
		do {
			fscanf (in, "%c", &ch);
		} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');

		if ( ch == ']' ) 
		{
			if(i!=param_k-1)
			{
				status=-1;
				fprintf(stderr," Input data error, please check the input data again...\n");
				goto TERMINATE;
			}
		}else if ( ch != ',' ) {
			status = -1;
			goto TERMINATE;
		}
	}
	/*read matrix_N*/
	for (;;)
	{
		fscanf (in, "%c", &ch);
		if ( ch == '\t' ||
			ch == '\r' ||
			ch == ' '  ||
			ch == '\n'   ) continue;
		if ( ch == '[' ) break;
		status = -1;
		fprintf(stderr," Input data error, please check the input data again...\n");
		goto TERMINATE;
	}

	for ( i = 0; i < param_m; i++ )
	{
		if ( (status = readarray(in, &n, (*matrix_N_p)+i)) ) 
		{
			fprintf(stderr," Input data error, please check the input data again...\n");
			goto TERMINATE;
		}
		if ( n != param_n ) 
		{	
			status = -1;	
			fprintf(stderr," Input data error, please check the input data again...\n");
			goto TERMINATE;
		}
		do {
			fscanf (in, "%c", &ch);
		} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');

		if ( ch == ']' ) 
		{
			if(i!=param_m-1)
			{
				status=-1;
				fprintf(stderr," Input data error, please check the input data again...\n");
				goto TERMINATE;
			}
		}else if ( ch != ',' ) {
			status = -1;
			goto TERMINATE;
		}
	}

	/*read matrix_M*/
	for (;;)
	{
		fscanf (in, "%c", &ch);
		if ( ch == '\t' ||
			ch == '\r' ||
			ch == ' '  ||
			ch == '\n'   ) continue;
		if ( ch == '[' ) break;
		status = -1;
		fprintf(stderr," Input data error, please check the input data again...\n");
		goto TERMINATE;
	}

	for ( i = 0; i < param_m; i++ )
	{
		if ( (status = readarray(in, &n, (*matrix_M_p)+i)) ) 
		{
			fprintf(stderr," Input data error, please check the input data again...\n");
			goto TERMINATE;
		}
		if ( n != param_m) 
		{	
			status = -1;	
			fprintf(stderr," Input data error, please check the input data again...\n");
			goto TERMINATE;
		}
		do {
			fscanf (in, "%c", &ch);
		} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');

		if ( ch == ']' ) 
		{
			if(i!=param_m-1)
			{
				status=-1;
				fprintf(stderr," Input data error, please check the input data again...\n");
				goto TERMINATE;
			}
		}else if ( ch != ',' ) {
			status = -1;
			goto TERMINATE;
		}
	}

TERMINATE:

	if( in != NULL ) 
	{
		fclose(in);
	}
	free_and_null ((char **)&param_p);
	return (status);

}

int readdata_format_3(char* file, 
					  int *param_n_p,
					  int *param_m_p,
					  int *param_k_p, 
					  double **c_coef_p,
					  double **d_coef_p,
					  double **b_coef_p, 
					  double **q_coef_p,  
					  MATRIX *_matrix_A_p,
					  MATRIX *_matrix_B_p,
					  MATRIX *_matrix_N_p,
					  MATRIX *_matrix_M_p)
{
	int status = 0;
	int param_n,param_m,param_k;
	int n_row,n_col,n_nz;
	double *param_p=NULL;
	double *A_param_p=NULL;
	double *B_param_p=NULL;
	double *M_param_p=NULL;
	double *N_param_p=NULL;
	int n;
	char ch;
	FILE *in=NULL;
	in=fopen(file,"r");
	if ( in == NULL ) 
	{
		status = -1;
		fprintf(stderr," unable to open the file... \n");
		goto TERMINATE;
	}
	/*read parameter*/
	if ( (status = readarray(in, &n,&param_p)) ) 
	{
		fprintf(stderr," Input data error param1, please check the input data again...\n");
		goto TERMINATE;
	}
	if ( n != 3 )
	{
		status = -1;
		fprintf(stderr," Input data error param2, please check the input data again...\n");
		goto TERMINATE;
	}
	param_n=(int) floor(param_p[0]);
	param_m=(int) floor(param_p[1]);
	param_k=(int) floor(param_p[2]);
	*param_n_p=param_n;
	*param_m_p=param_m;
	*param_k_p=param_k;
	
	//
	/*read c_coef*/
	if ( (status = readarray(in, &n, c_coef_p)) ) 
	{
		fprintf(stderr," Input data error c_coef1, please check the input data again...\n");
		goto TERMINATE;
	}
	if ( n != param_n) {
		status = -1;
		fprintf(stderr," Input data error c_coef2, please check the input data again...\n");
		goto TERMINATE;
	}
	/*read d_coef*/
	if ( (status = readarray(in, &n, d_coef_p)) ) 
	{
		fprintf(stderr," Input data error d_coef1, please check the input data again...\n");
		goto TERMINATE;
	}
	if ( n != param_m) {
		status = -1;
		fprintf(stderr," Input data error d_coef2, please check the input data again...\n");
		goto TERMINATE;
	}
	/*read b_coef*/
	if ( (status = readarray(in, &n, b_coef_p)) ) 
	{
		fprintf(stderr," Input data error b_coef1, please check the input data again...\n");
		goto TERMINATE;
	}
	if ( n != param_k) {
		status = -1;
		fprintf(stderr," Input data error b_coef2, please check the input data again...\n");
		goto TERMINATE;
	}
	/*read q_coef*/
	if ( (status = readarray(in, &n,q_coef_p)) ) 
	{
		fprintf(stderr," Input data error q_coef1, please check the input data again...\n");
		goto TERMINATE;
	}
	if ( n != param_m) {
		status = -1;
		fprintf(stderr," Input data error q_coef2, please check the input data again...\n");
		goto TERMINATE;
	}
	/*read matrix_A*/
	free_matrix(_matrix_A_p);
	for (;;)
	{
		fscanf (in, "%c", &ch);
		if ( ch == '\t' ||
			ch == '\r' ||
			ch == ' '  ||
			ch == '\n'   ) continue;
		if ( ch == '[' ) break;
		status = -1;
		fprintf(stderr," Input data error A0, please check the input data again...\n");
		goto TERMINATE;
	}
	//parameter
	if ( (status = readarray(in, &n, &A_param_p)) ) 
	{
		fprintf(stderr," Input data error A1, please check the input data again...\n");
		goto TERMINATE;
	}
	if ( n !=3 ) 
	{	
		status = -1;	
		fprintf(stderr," Input data error A2, please check the input data again...\n");
		goto TERMINATE;
	}else
	{
		n_row=(int) floor(A_param_p[0]);
		n_col=(int) floor(A_param_p[1]);
		n_nz=(int) floor(A_param_p[2]);
		_matrix_A_p->param=1;
		_matrix_A_p->n_row=n_row;
		_matrix_A_p->n_col=n_col;
		_matrix_A_p->nnz=n_nz;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) {
		status = -1;
		goto TERMINATE;
	}
	//matbeg
	if ( (status = readarray_int(in, &n, &(_matrix_A_p->matbeg))) ) 
	{
		fprintf(stderr," Input data error A3, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_row)
	{
		status = -1;	
		fprintf(stderr," Input data error A4, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
    if ( ch != ',' ) {
		status = -1;
		goto TERMINATE;
	}
	//matcnt
	if ( (status = readarray_int(in, &n, &(_matrix_A_p->matcnt))) ) 
	{
		fprintf(stderr," Input data error A5, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_row)
	{
		status = -1;	
		fprintf(stderr," Input data error A6, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) {
		status = -1;
		goto TERMINATE;
	}
	//matind
	if ( (status = readarray_int(in, &n, &(_matrix_A_p->matind))) ) 
	{
		fprintf(stderr," Input data error A7, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_nz)
	{
		status = -1;	
		fprintf(stderr," Input data error A8, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) {
		status = -1;
		goto TERMINATE;
	}
	//matval
	if ( (status = readarray(in, &n, &(_matrix_A_p->matval))) ) 
	{
		fprintf(stderr," Input data error A9, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_nz)
	{
		status = -1;	
		fprintf(stderr," Input data error A10, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ']' ) 
	{
		status=-1;
		fprintf(stderr," Input data error A11, please check the input data again...\n");
		goto TERMINATE;
	}

	/*read matrix_B*/
	free_matrix(_matrix_B_p);
	for (;;)
	{
		fscanf (in, "%c", &ch);
		if ( ch == '\t' ||
			ch == '\r' ||
			ch == ' '  ||
			ch == '\n'   ) continue;
		if ( ch == '[' ) break;
		status = -1;
		fprintf(stderr," Input data error B0, please check the input data again...\n");
		goto TERMINATE;
	}
	//parameter
	if ( (status = readarray(in, &n, &A_param_p)) ) 
	{
		fprintf(stderr," Input data error B1, please check the input data again...\n");
		goto TERMINATE;
	}
	if ( n !=3 ) 
	{	
		status = -1;	
		fprintf(stderr," Input data error B2, please check the input data again...\n");
		goto TERMINATE;
	}else
	{
		n_row=(int) floor(A_param_p[0]);
		n_col=(int) floor(A_param_p[1]);
		n_nz=(int) floor(A_param_p[2]);
		_matrix_B_p->param=1;
		_matrix_B_p->n_row=n_row;
		_matrix_B_p->n_col=n_col;
		_matrix_B_p->nnz=n_nz;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) 
	{
		status=-1;
		fprintf(stderr," Input data error B3, please check the input data again...\n");
		goto TERMINATE;
	}
	//matbeg
	if ( (status = readarray_int(in, &n, &(_matrix_B_p->matbeg))) ) 
	{
		fprintf(stderr," Input data error B4, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_row)
	{
		status = -1;	
		fprintf(stderr," Input data error B5, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) 
	{
		status=-1;
		fprintf(stderr," Input data error B6, please check the input data again...\n");
		goto TERMINATE;
	}
	//matcnt
	if ( (status = readarray_int(in, &n, &(_matrix_B_p->matcnt))) ) 
	{
		fprintf(stderr," Input data error B7, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_row)
	{
		status = -1;	
		fprintf(stderr," Input data error B8, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) 
	{
		status=-1;
		fprintf(stderr," Input data error B9, please check the input data again...\n");
		goto TERMINATE;
	}
	//matind
	if ( (status = readarray_int(in, &n, &(_matrix_B_p->matind))) ) 
	{
		fprintf(stdout," n, n_nz: %i, %i\n",n,n_nz);
		fprintf(stderr," Input data error B10, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_nz)
	{
		status = -1;	
		fprintf(stderr," Input data error B11, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) 
	{
		status=-1;
		fprintf(stderr," Input data error B12, please check the input data again...\n");
		goto TERMINATE;
	}
	//matval
	if ( (status = readarray(in, &n, &(_matrix_B_p->matval))) ) 
	{
		fprintf(stderr," Input data error B13, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_nz)
	{
		status = -1;	
		fprintf(stderr," Input data error B14, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ']' ) 
	{
		status=-1;
		fprintf(stderr," Input data error B15, please check the input data again...\n");
		goto TERMINATE;
	}
	/*read matrix_N*/
	free_matrix(_matrix_N_p);
	for (;;)
	{
		fscanf (in, "%c", &ch);
		if ( ch == '\t' ||
			ch == '\r' ||
			ch == ' '  ||
			ch == '\n'   ) continue;
		if ( ch == '[' ) break;
		status = -1;
		fprintf(stderr," Input data error N0, please check the input data again...\n");
		goto TERMINATE;
	}
	//parameter
	if ( (status = readarray(in, &n, &A_param_p)) ) 
	{
		fprintf(stderr," Input data error N1, please check the input data again...\n");
		goto TERMINATE;
	}
	if ( n !=3 ) 
	{	
		status = -1;	
		fprintf(stderr," Input data error N2, please check the input data again...\n");
		goto TERMINATE;
	}else
	{
		n_row=(int) floor(A_param_p[0]);
		n_col=(int) floor(A_param_p[1]);
		n_nz=(int) floor(A_param_p[2]);
		_matrix_N_p->param=1;
		_matrix_N_p->n_row=n_row;
		_matrix_N_p->n_col=n_col;
		_matrix_N_p->nnz=n_nz;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) 
	{
		status=-1;
		fprintf(stderr," Input data error N3, please check the input data again...\n");
		goto TERMINATE;
	}
	//matbeg
	if ( (status = readarray_int(in, &n, &(_matrix_N_p->matbeg))) ) 
	{
		fprintf(stderr," Input data error N4, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_row)
	{
		status = -1;	
		fprintf(stderr," Input data error N5, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) 
	{
		status=-1;
		fprintf(stderr," Input data error N6, please check the input data again...\n");
		goto TERMINATE;
	}
	//matcnt
	if ( (status = readarray_int(in, &n, &(_matrix_N_p->matcnt))) ) 
	{
		fprintf(stderr," Input data error N7, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_row)
	{
		status = -1;	
		fprintf(stderr," Input data error N8, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) 
	{
		status=-1;
		fprintf(stderr," Input data error N9, please check the input data again...\n");
		goto TERMINATE;
	}
	//matind
	if ( (status = readarray_int(in, &n, &(_matrix_N_p->matind))) ) 
	{
		fprintf(stderr," Input data error N10, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_nz)
	{
		status = -1;	
		fprintf(stderr," Input data error N11, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) 
	{
		status=-1;
		fprintf(stderr," Input data error N12, please check the input data again...\n");
		goto TERMINATE;
	}
	//matval
	if ( (status = readarray(in, &n, &(_matrix_N_p->matval))) ) 
	{
		fprintf(stderr," Input data error N13, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_nz)
	{
		status = -1;	
		fprintf(stderr," Input data error N14, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ']' ) 
	{
		status=-1;
		fprintf(stderr," Input data error N15, please check the input data again...\n");
		goto TERMINATE;
	}
	/*read matrix_M*/
	free_matrix(_matrix_M_p);
	for (;;)
	{
		fscanf (in, "%c", &ch);
		if ( ch == '\t' ||
			ch == '\r' ||
			ch == ' '  ||
			ch == '\n'   ) continue;
		if ( ch == '[' ) break;
		status = -1;
		fprintf(stderr," Input data error M0, please check the input data again...\n");
		goto TERMINATE;
	}
	//parameter
	if ( (status = readarray(in, &n, &A_param_p)) ) 
	{
		fprintf(stderr," Input data error M1, please check the input data again...\n");
		goto TERMINATE;
	}
	if ( n !=3 ) 
	{	
		status = -1;	
		fprintf(stderr," Input data error M2, please check the input data again...\n");
		goto TERMINATE;
	}else
	{
		n_row=(int) floor(A_param_p[0]);
		n_col=(int) floor(A_param_p[1]);
		n_nz=(int) floor(A_param_p[2]);
		_matrix_M_p->param=1;
		_matrix_M_p->n_row=n_row;
		_matrix_M_p->n_col=n_col;
		_matrix_M_p->nnz=n_nz;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) 
	{
		status=-1;
		fprintf(stderr," Input data error M3, please check the input data again...\n");
		goto TERMINATE;
	}
	//matbeg
	if ( (status = readarray_int(in, &n, &(_matrix_M_p->matbeg))) ) 
	{
		fprintf(stderr," Input data error M4, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_row)
	{
		status = -1;	
		fprintf(stderr," Input data error M5, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) 
	{
		status=-1;
		fprintf(stderr," Input data error M6, please check the input data again...\n");
		goto TERMINATE;
	}
	//matcnt
	if ( (status = readarray_int(in, &n, &(_matrix_M_p->matcnt))) ) 
	{
		fprintf(stderr," Input data error M7, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_row)
	{
		status = -1;	
		fprintf(stderr," Input data error M8, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) 
	{
		status=-1;
		fprintf(stderr," Input data error M9, please check the input data again...\n");
		goto TERMINATE;
	}
	//matind
	if ( (status = readarray_int(in, &n, &(_matrix_M_p->matind))) ) 
	{
		fprintf(stderr," Input data error M10, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_nz)
	{
		status = -1;	
		fprintf(stderr," Input data error M11, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ',' ) 
	{
		status=-1;
		fprintf(stderr," Input data error M12, please check the input data again...\n");
		goto TERMINATE;
	}
	//matval
	if ( (status = readarray(in, &n, &(_matrix_M_p->matval))) ) 
	{
		fprintf(stderr," Input data error M13, please check the input data again...\n");
		goto TERMINATE;
	}
	if (n!=n_nz)
	{
		status = -1;	
		fprintf(stderr," Input data error M14, please check the input data again...\n");
		goto TERMINATE;
	}
	do {
		fscanf (in, "%c", &ch);
	} while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
	if ( ch != ']' ) 
	{
		status=-1;
		fprintf(stderr," Input data error M15, please check the input data again...\n");
		goto TERMINATE;
	}

TERMINATE:

	if( in != NULL ) 
	{
		fclose(in);
	}
	free_and_null ((char **)&param_p);
	free_and_null((char **)&A_param_p);
	free_and_null((char **)&B_param_p);
	free_and_null((char **)&M_param_p);
	free_and_null((char **)&N_param_p);
	return (status);

}
//end reading data

//This function is used to convert binary integer program from MPS format to our format
//therefore we can use our code to test the instances from MIPLIB 
int convertMPS(char* file, 
			   int *param_n_p,
			   int *param_m_p,
			   int *param_k_p, 
			   double **c_coef_p,
			   double **d_coef_p,
			   double **b_coef_p, 
			   double **q_coef_p,  
			   double ***matrix_A_p,
			   double ***matrix_B_p, 
			   double ***matrix_N_p,
			   double ***matrix_M_p)
{
	int status=0;
	int lp_num_row;
	int lp_num_col;
	int param_m;
	int param_k;
	int param_n;
	int i,j;
	char *var_type=NULL;
	char *con_sense=NULL;
	double tmp_val;
	double *obj_coef=NULL;
	double *var_lb=NULL;
	double *var_ub=NULL;
	double *con_rhs=NULL;
	CPXENVptr     env = NULL;
	CPXLPptr      lp = NULL;
	env = CPXopenCPLEX (&status);
	if ( env == NULL ) {
		char  errmsg[1024];
		fprintf (stderr, "Could not open CPLEX environment.\n");
		CPXgeterrorstring (env, status, errmsg);
		fprintf (stderr, "%s", errmsg);
		goto TERMINATE;
	}
	lp = CPXcreateprob (env, &status, "mps_lp");
	if ( lp == NULL ) {
		fprintf (stderr, "Failed to create LP.\n");
		goto TERMINATE;
	}
	status=CPXreadcopyprob(env,lp,file,NULL);
	if(status) goto TERMINATE;
	lp_num_col=CPXgetnumcols(env,lp);
	lp_num_row=CPXgetnumrows(env,lp);
	//count the dimension of x and y
	var_type=(char*)malloc(lp_num_col*sizeof(char));
	con_sense=(char*)malloc(lp_num_row*sizeof(char));
	obj_coef=(double*)malloc(lp_num_col*sizeof(double));
	var_lb=(double*)malloc(lp_num_col*sizeof(double));
	var_ub=(double*)malloc(lp_num_col*sizeof(double));
	con_rhs=(double*)malloc(lp_num_row*sizeof(double));
	if(var_type==NULL 
		|| con_sense==NULL
		|| obj_coef==NULL
		|| var_lb==NULL
		|| var_ub==NULL
		|| con_rhs==NULL)
	{
		status=-1;
		fprintf(stderr," convertMPS(): unable to allocate memory\n");
		goto TERMINATE;
	}
	status=CPXgetctype(env,lp,var_type,0,lp_num_col-1);
	if(status)goto TERMINATE;
	status=CPXgetsense(env,lp,con_sense,0,lp_num_row-1);
	if(status)goto TERMINATE;
	status=CPXgetobj(env,lp,obj_coef,0,lp_num_col-1);
	if(status)goto TERMINATE;
	status=CPXgetlb(env,lp,var_lb,0,lp_num_col-1);
	if(status)goto TERMINATE;
	status=CPXgetub(env,lp,var_ub,0,lp_num_col-1);
	if(status)goto TERMINATE;
	status=CPXgetrhs(env,lp,con_rhs,0,lp_num_row-1);
	if(status)goto TERMINATE;
	//first decide param_m, param_n, param_k
	param_k=0;
	param_m=0;
	param_n=0;
	for(i=0;i<lp_num_col;i++)
	{
		if (fabs(var_lb[i]-var_ub[i])<ZERO_TOLERANCE || var_type[i]=='C')
		{
			param_n++;
			if (var_lb[i]>-CPX_INFBOUND)
			{
				param_k++;
			}
			if (var_ub[i]<CPX_INFBOUND)
			{
				param_k++;
			}

		}else if ((var_type[i]=='B' && fabs(var_lb[i])<ZERO_TOLERANCE && fabs(var_ub[i]-1)<ZERO_TOLERANCE)|| (var_type[i]=='I' && fabs(var_lb[i])<ZERO_TOLERANCE && fabs(var_ub[i]-1)<ZERO_TOLERANCE))
		{
			param_m++;
		}else
		{	
			status=-1;
			fprintf(stderr," convertMPS(): unexpected var type: %c, lb: %f, ub: %f\n",var_type[i],var_lb[i],var_ub[i]);
			goto TERMINATE;
		}
	}
	for (i=0;i<lp_num_row;i++)
	{
		if (con_sense[i]=='L' || con_sense[i]=='G')
		{
			param_k++;
		}else if (con_sense[i]=='E')
		{
			param_k++;
			param_k++;
		}else 
		{
			status=-1;
			fprintf(stderr," convertMPS(): unexpected con_sense:%c\n",con_sense[i]);
			goto TERMINATE;
		}
	}
	*param_n_p=param_n;
	*param_m_p=param_m;
	*param_k_p=param_k;
	*c_coef_p=(double*)calloc(param_n,sizeof(double));
	*d_coef_p=(double*)calloc(param_m,sizeof(double));
	*b_coef_p=(double*)calloc(param_k,sizeof(double));
	*q_coef_p=(double*)calloc(param_m,sizeof(double));
	*matrix_A_p=(double**)calloc(param_k,sizeof(double*));
	*matrix_B_p=(double**)calloc(param_k,sizeof(double*));
	*matrix_N_p=(double**)calloc(param_m,sizeof(double*));
	*matrix_M_p=(double**)calloc(param_m,sizeof(double*));
	if (*c_coef_p==NULL
		|| *d_coef_p==NULL
		|| *b_coef_p==NULL
		|| *q_coef_p==NULL
		|| *matrix_A_p==NULL 
		|| *matrix_B_p==NULL 
		|| *matrix_N_p==NULL 
		|| *matrix_M_p==NULL)
	{
		status = -1;
		fprintf(stderr," convertMPS(): unable to allocate memory\n");
		goto TERMINATE;
	}
	for (i=0;i<param_k;i++)
	{
		(*matrix_A_p)[i]=(double*)calloc(param_n,sizeof(double));
		(*matrix_B_p)[i]=(double*)calloc(param_m,sizeof(double));
		if ((*matrix_A_p)[i]==NULL || (*matrix_B_p)[i]==NULL)
		{
			status=-1;
			fprintf(stderr," convertMPS(): unable to allocate memory\n");
			goto TERMINATE;
		}
	}
	for (i=0;i<param_m;i++)
	{
		(*matrix_N_p)[i]=(double*)calloc(param_n,sizeof(double));
		(*matrix_M_p)[i]=(double*)calloc(param_m,sizeof(double));
		if ((*matrix_N_p)[i]==NULL || (*matrix_M_p)[i]==NULL)
		{
			status=-1;
			fprintf(stderr," convertMPS(): unable to allocate memory\n");
			goto TERMINATE;
		}
	}
	//start to initialize 
	param_k=0;
	param_m=0;
	param_n=0;
	for(i=0;i<lp_num_col;i++)
	{
		if (fabs(var_lb[i]-var_ub[i])<ZERO_TOLERANCE || var_type[i]=='C')
		{

			(*c_coef_p)[param_n]=obj_coef[i];
			if (var_lb[i]>-CPX_INFBOUND)
			{
				(*matrix_A_p)[param_k][param_n]=1;
				(*b_coef_p)[param_k]=var_lb[i];
				param_k++;
			}
			if (var_ub[i]<CPX_INFBOUND)
			{
				(*matrix_A_p)[param_k][param_n]=-1;
				(*b_coef_p)[param_k]=-var_ub[i];
				param_k++;
			}
			param_n++;
		}else if ((var_type[i]=='B' && fabs(var_lb[i])<ZERO_TOLERANCE && fabs(var_ub[i]-1)<ZERO_TOLERANCE)|| (var_type[i]=='I' && fabs(var_lb[i])<ZERO_TOLERANCE && fabs(var_ub[i]-1)<ZERO_TOLERANCE))
		{
			(*d_coef_p)[param_m]=obj_coef[i];
			param_m++;	
		}else
		{	
			status=-1;
			fprintf(stderr," convertMPS(): unexpected var type: %c, lb: %f, ub: %f\n",var_type[i],var_lb[i],var_ub[i]);
			goto TERMINATE;
		}
	}
	for (j=0;j<lp_num_row;j++)
	{
		param_m=0;
		param_n=0;
		if (con_sense[j]=='L')
		{
			for(i=0;i<lp_num_col;i++)
			{
				if (fabs(var_lb[i]-var_ub[i])<ZERO_TOLERANCE || var_type[i]=='C')
				{
					status=CPXgetcoef(env,lp,j,i,&tmp_val);
					if(status)goto TERMINATE;
					(*matrix_A_p)[param_k][param_n]=-tmp_val;
					param_n++;
				}else if ((var_type[i]=='B' && fabs(var_lb[i])<ZERO_TOLERANCE && fabs(var_ub[i]-1)<ZERO_TOLERANCE)|| (var_type[i]=='I' && fabs(var_lb[i])<ZERO_TOLERANCE && fabs(var_ub[i]-1)<ZERO_TOLERANCE))
				{
					status=CPXgetcoef(env,lp,j,i,&tmp_val);
					if(status)goto TERMINATE;
					(*matrix_B_p)[param_k][param_m]=-tmp_val;
					param_m++;	
				}else
				{	
					status=-1;
					fprintf(stderr," convertMPS(): unexpected var type: %c, lb: %f, ub: %f\n",var_type[i],var_lb[i],var_ub[i]);
					goto TERMINATE;
				}
			}
			(*b_coef_p)[param_k]=-con_rhs[j];
			param_k++;
		}else if (con_sense[j]=='G')
		{
			for(i=0;i<lp_num_col;i++)
			{
				if (fabs(var_lb[i]-var_ub[i])<ZERO_TOLERANCE || var_type[i]=='C')
				{
					status=CPXgetcoef(env,lp,j,i,&tmp_val);
					if(status)goto TERMINATE;
					(*matrix_A_p)[param_k][param_n]=tmp_val;
					param_n++;
				}else if ((var_type[i]=='B' && fabs(var_lb[i])<ZERO_TOLERANCE && fabs(var_ub[i]-1)<ZERO_TOLERANCE)|| (var_type[i]=='I' && fabs(var_lb[i])<ZERO_TOLERANCE && fabs(var_ub[i]-1)<ZERO_TOLERANCE))
				{
					status=CPXgetcoef(env,lp,j,i,&tmp_val);
					if(status)goto TERMINATE;
					(*matrix_B_p)[param_k][param_m]=tmp_val;
					param_m++;	
				}else
				{	
					status=-1;
					fprintf(stderr," convertMPS(): unexpected var type: %c, lb: %f, ub: %f\n",var_type[i],var_lb[i],var_ub[i]);
					goto TERMINATE;
				}
			}
			(*b_coef_p)[param_k]=con_rhs[j];
			param_k++;
		}else if (con_sense[j]=='E')
		{
			for(i=0;i<lp_num_col;i++)
			{
				if (fabs(var_lb[i]-var_ub[i])<ZERO_TOLERANCE || var_type[i]=='C')
				{
					status=CPXgetcoef(env,lp,j,i,&tmp_val);
					if(status)goto TERMINATE;
					(*matrix_A_p)[param_k][param_n]=-tmp_val;
					(*matrix_A_p)[param_k+1][param_n]=tmp_val;
					param_n++;
				}else if ((var_type[i]=='B' && fabs(var_lb[i])<ZERO_TOLERANCE && fabs(var_ub[i]-1)<ZERO_TOLERANCE)|| (var_type[i]=='I' && fabs(var_lb[i])<ZERO_TOLERANCE && fabs(var_ub[i]-1)<ZERO_TOLERANCE))
				{
					status=CPXgetcoef(env,lp,j,i,&tmp_val);
					if(status)goto TERMINATE;
					(*matrix_B_p)[param_k][param_m]=-tmp_val;
					(*matrix_B_p)[param_k+1][param_m]=tmp_val;
					param_m++;	
				}else
				{	
					status=-1;
					fprintf(stderr," convertMPS(): unexpected var type: %c, lb: %f, ub: %f\n",var_type[i],var_lb[i],var_ub[i]);
					goto TERMINATE;
				}
			}
			(*b_coef_p)[param_k]=-con_rhs[j];
			(*b_coef_p)[param_k+1]=con_rhs[j];
			param_k++;
			param_k++;
		}else 
		{
			status=-1;
			fprintf(stderr," convertMPS(): unexpected con_sense:%c\n",con_sense[j]);
			goto TERMINATE;
		}
	}
	for (i=0;i<(*param_m_p);i++)
	{
		(*q_coef_p)[i]=1;
		(*matrix_M_p)[i][i]=-1;
	}
	//status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
	//if ( status )
	//{
	//	fprintf (stderr, 
	//		"Failure to turn on screen indicator, error %d.\n", status);
	//	goto TERMINATE;
	//}
	//status=CPXmipopt(env,lp);
	//if(status) goto TERMINATE;

TERMINATE:
	if ( lp != NULL ) CPXfreeprob (env, &lp);
	if ( env != NULL ) CPXcloseCPLEX (&env);
	free_and_null((char**)&var_type);
	free_and_null((char**)&con_sense);
	free_and_null((char**)&obj_coef);
	free_and_null((char**)&var_lb);
	free_and_null((char**)&var_ub);
	free_and_null((char**)&con_rhs);
	return(status);
}


int node_pre_process(CPXENVptr env, 
					 CPXLPptr rx_lp, 
					 const int param_n,
					 const int param_m, 
					 const int param_k, 
					 NODE *activenodeptr, 
					 double *ub, 				
					 int pre_rx_lp_condition,
					 int *pre_process_status)
{
	//pre_process_status=1 means this activenode is pruned
	int status=0;
	int i,j;
	int fix_var_cnt;
	int fix_check_1;
	int fix_check_2;
	int row_id;
	const int ncols=CPXgetnumcols(env,rx_lp);
	const int nrows=CPXgetnumrows(env,rx_lp);
	const int num_var= param_n+param_m+param_m+param_k;
	int lp_stat;
	double *lp_soln=NULL;
	double *reduce_cost=NULL;
	double *tableau_1=NULL;
	double *tableau_2=NULL;
	double lp_obj;
	STARTINFO lp_start;
	*pre_process_status=0;
	init_startinfo(&lp_start);
	lp_soln=(double*)malloc(num_var*sizeof(double));
	reduce_cost=(double*)malloc(ncols*sizeof(double));
	tableau_1=(double*)malloc(ncols*sizeof(double));
	tableau_2=(double*)malloc(ncols*sizeof(double));
	if (lp_soln==NULL)
	{
		status=-1;
		fprintf(stderr," node_pre_process(): unable to allocate memory\n");
		goto TERMINATE;
	}
	status=setupnode_rx_lp(env,rx_lp,param_n,param_m,param_k,*activenodeptr);
	if(status) goto TERMINATE;
	status=cplexSolveLP(env,rx_lp,activenodeptr->nodestartinfo,num_var,lp_soln,&lp_obj,&lp_stat,&lp_start,CPX_ALG_DUAL,1,pre_rx_lp_condition);
	if(status) goto TERMINATE;
	if (lp_stat==CPX_STAT_OPTIMAL)
	{
		fix_var_cnt=0;
		for (i=0;i<param_m;i++)
		{
			if (activenodeptr->nodebranchinfo[i]==0)
			{
				if (lp_soln[param_n+i]>ZERO_TOLERANCE && lp_soln[param_n+param_m+i]>ZERO_TOLERANCE)
				{
					fix_check_1=0;
					fix_check_2=0;
					status = CPXgetijrow(env,rx_lp,-1,param_n+i,&row_id);
					if (status) goto TERMINATE;
					status = CPXbinvarow(env,rx_lp,row_id,tableau_1);
					if (status) goto TERMINATE;
					status = CPXgetijrow(env,rx_lp,-1,param_n+param_m+i,&row_id);
					if (status) goto TERMINATE;
					status = CPXbinvarow(env,rx_lp,row_id,tableau_2);
					if(status) goto TERMINATE;
					for (j=0;j<ncols;j++)
					{
						switch(activenodeptr->nodestartinfo.cstat[j])
						{
						case CPX_BASIC:
							break;
						case CPX_AT_LOWER:
							if(tableau_1[j]>ZERO_TOLERANCE)fix_check_1=1;
							if(tableau_2[j]>ZERO_TOLERANCE)fix_check_2=1;
							break;
						case CPX_AT_UPPER:
							if(tableau_1[j]<-ZERO_TOLERANCE)fix_check_1=1;
							if(tableau_2[j]<-ZERO_TOLERANCE)fix_check_2=1;
							break;
						case CPX_FREE_SUPER:
							if(fabs(tableau_1[j])>ZERO_TOLERANCE)fix_check_1=1;
							if(fabs(tableau_2[j])>ZERO_TOLERANCE)fix_check_2=1;
							break;
						default:
							status=-1;
							fprintf(stderr," node_pre_process(): unexpected cstat: %d\n",activenodeptr->nodestartinfo.cstat[j]);
							goto TERMINATE;
						}
						if(fix_check_1 && fix_check_2)break;
					}

					if (fix_check_1==0 || fix_check_2==0)
					{
						fix_var_cnt++;
						if(fix_check_1==0)activenodeptr->nodebranchinfo[i]=2;
						if(fix_check_2==0)activenodeptr->nodebranchinfo[i]=1;
						if (fix_check_1==0 && fix_check_2==0)
						{
							*pre_process_status=1;
							fprintf(stdout," node is pruned during node_pre_process\n");
							goto TERMINATE;
						}
					}
				}else
				{
					if (lp_soln[param_n+i]>ZERO_TOLERANCE)
					{
						fix_check_1=0;
						status = CPXgetijrow(env,rx_lp,-1,param_n+i,&row_id);
						if (status) goto TERMINATE;
						status = CPXbinvarow(env,rx_lp,row_id,tableau_1);
						if (status) goto TERMINATE;
						for (j=0;j<ncols;j++)
						{
							switch(activenodeptr->nodestartinfo.cstat[j])
							{
							case CPX_BASIC:
								break;
							case CPX_AT_LOWER:
								if(tableau_1[j]>ZERO_TOLERANCE)fix_check_1=1;
								break;
							case CPX_AT_UPPER:
								if(tableau_1[j]<-ZERO_TOLERANCE)fix_check_1=1;
								break;
							case CPX_FREE_SUPER:
								if(fabs(tableau_1[j])>ZERO_TOLERANCE)fix_check_1=1;
								break;
							default:
								status=-1;
								fprintf(stderr," node_pre_process(): unexpected cstat: %d\n",activenodeptr->nodestartinfo.cstat[j]);
								goto TERMINATE;
							}
							if(fix_check_1)break;
						}
						if (fix_check_1==0)
						{
							fix_var_cnt++;
							activenodeptr->nodebranchinfo[i]=2;
						}
					}
					if (lp_soln[param_n+param_m+i]>ZERO_TOLERANCE)
					{
						fix_check_2=0;
						status = CPXgetijrow(env,rx_lp,-1,param_n+param_m+i,&row_id);
						if (status) goto TERMINATE;
						status = CPXbinvarow(env,rx_lp,row_id,tableau_2);
						if(status) goto TERMINATE;
						for (j=0;j<ncols;j++)
						{
							switch(activenodeptr->nodestartinfo.cstat[j])
							{
							case CPX_BASIC:
								break;
							case CPX_AT_LOWER:
								if(tableau_2[j]>ZERO_TOLERANCE)fix_check_2=1;
								break;
							case CPX_AT_UPPER:
								if(tableau_2[j]<-ZERO_TOLERANCE)fix_check_2=1;
								break;
							case CPX_FREE_SUPER:
								if(fabs(tableau_2[j])>ZERO_TOLERANCE)fix_check_2=1;
								break;
							default:
								status=-1;
								fprintf(stderr," node_pre_process(): unexpected cstat: %d\n",activenodeptr->nodestartinfo.cstat[j]);
								goto TERMINATE;
							}
							if(fix_check_2)break;
						}
						if (fix_check_2==0)
						{
							fix_var_cnt++;
							activenodeptr->nodebranchinfo[i]=1;
						}
					}
				}
			}
		}
		if(fix_var_cnt>0)fprintf(stdout,"fix %d complementarity during node_pre_process\n",fix_var_cnt);
	}else
	{
		status=-1;
		fprintf(stderr," node_pre_process(): unexpected lp_stat: %d\n",lp_stat);
		goto TERMINATE;
	}
TERMINATE:
	free_and_null((char**)&lp_soln);
	free_and_null((char**)&reduce_cost);
	free_and_null((char**)&tableau_1);
	free_and_null((char**)&tableau_2);
	free_startinfo(&lp_start);
	return(status);
}

int writeSolution(int LPCC_condition,
				  int param_n,	
				  int param_m,	
				  int param_k, 	
				  int solve_node,
				  PARAM preprocess_param,
				  double preprocessing_time,
				  double branching_time,
				  double wall_preproc_time,
				  double wall_solve_time,
				  double lp_relaxation,
				  double improved_lb,
				  double best_lb,
				  double opt_obj,	
				  double *lpcc_soln)
{
	int status=0;
	int i;
	FILE *pFile=NULL;
	fprintf(stdout,"Save solution to file %s\n",preprocess_param.solnfile);
	pFile=fopen(preprocess_param.solnfile,"w");
	fprintf(pFile,"-------------------------------------------------------------------------\n");
	switch (preprocess_param.solve_mode)
	{
	case 'c':
		fprintf(pFile,"Solve LPCC by CPLEX\n");
		fprintf(pFile,"Preprocessing time: %f\n",preprocessing_time);
		fprintf(pFile,"Solving time: %f\n",branching_time);
		fprintf(pFile,"Total time: %f\n",preprocessing_time+branching_time);
                fprintf(pFile,"Preprocessing time (wall): %.0f\n",wall_preproc_time);
                fprintf(pFile,"Solving time (wall): %.0f\n",wall_solve_time);
                fprintf(pFile,"Total time (wall): %.0f\n",wall_preproc_time+wall_solve_time);

		fprintf(pFile,"Total solved node: %d\n",solve_node);
		fprintf(pFile,"Preprocessing parameter:\n");
		fprintf(pFile,"-Feasibility_recovery_breath: %f\n",preprocess_param.feasibility_recovery_breath);
		fprintf(pFile,"-Feasibility_recovery_depth: %d\n",preprocess_param.feasibility_recovery_depth);
		fprintf(pFile,"-BoundCut generate loop: %d\n",preprocess_param.boundcut_loop);
		fprintf(pFile,"-BoundCut split: %d\n",preprocess_param.boundcut_p);
		fprintf(pFile,"-DisjunctiveCut generate loop: %d\n",preprocess_param.disjunctivecut_loop_param);
		fprintf(pFile,"-SimpleCut generate loop: %d\n",preprocess_param.simplecut_loop_param);
		fprintf(pFile,"-MccormickRefineCnt: %d\n",preprocess_param.MccormickRefineCnt);
		fprintf(pFile,"-------------------------------------------------------------------------\n");
		if(lp_relaxation!=MAX_VALUE || lp_relaxation!=MIN_VALUE )fprintf(pFile,"LP relaxation: %f\n",lp_relaxation);
		if(improved_lb!=MAX_VALUE || improved_lb!=MIN_VALUE)fprintf(pFile,"Improved lower bound:%f\n",improved_lb);
		switch(LPCC_condition)
		{
		case LPCC_INFEASIBLE:
			fprintf(pFile,"LPCC is infeasible\n");
			break;
		case LPCC_INForUNBD:
			fprintf(pFile,"LPCC is infeasible or unbounded\n");
			break;
		case LPCC_OPTIMAL:
			fprintf(pFile,"LPCC is feasible and bounded\n");
			fprintf(pFile,"Optimal value: %f\n", opt_obj);
			fprintf(pFile,"Optimal solution:\n");
			for (i=0;i<param_n;i++)
			{
				fprintf(pFile,"x[%d]: %12.12f  ",i,lpcc_soln[i]);
				fprintf(pFile,"\n");
			}
			for (i=0;i<param_m;i++)
			{
				fprintf(pFile,"y[%d]: %12.12f  ",i,lpcc_soln[i+param_n]);
				fprintf(pFile,"\n");
			}
			break;
		case LPCC_TERMINATE_WITH_SOLN:
			fprintf(pFile,"CPLEX time limit reached, and feasible solution exists\n");
			fprintf(pFile,"Best lower bound: %f\n",best_lb);
			fprintf(pFile,"Best objective value; %f",opt_obj);
			fprintf(pFile,"Best solution:\n");
			for (i=0;i<param_n;i++)
			{
				fprintf(pFile,"x[%d]: %12.12f  ",i,lpcc_soln[i]);
				fprintf(pFile,"\n");
			}
			for (i=0;i<param_m;i++)
			{
				fprintf(pFile,"y[%d]: %12.12f  ",i,lpcc_soln[i+param_n]);
				fprintf(pFile,"\n");
			}
			break;
		case LPCC_TERMINATE_WITHOUT_SOLN:
			fprintf(pFile,"CPLEX time limit reached, and no feasible solution\n");
			fprintf(pFile,"Best lower bound: %f\n",best_lb);
			break;
		case LPCC_UNBOUNDED:
			fprintf(pFile,"LPCC is unbounded\n");
			break;
		default:
			fprintf(pFile,"Unexpected MIP CPLEX LPCC solution stauts 3: %d\n",LPCC_condition);
			status=-1;
			goto TERMINATE;
		}
		break;
	case 'b':
		fprintf(pFile,"Solve LPCC by Branch and cut.\n");
		fprintf(pFile,"Preprocessing time: %f\n",preprocessing_time);
		fprintf(pFile,"Solving time: %f\n",branching_time);
		fprintf(pFile,"Total time: %f\n",preprocessing_time+branching_time);
                fprintf(pFile,"Preprocessing time (wall): %.0f\n",wall_preproc_time);
                fprintf(pFile,"Solving time (wall): %.0f\n",wall_solve_time);
                fprintf(pFile,"Total time (wall): %.0f\n",wall_preproc_time+wall_solve_time);
		fprintf(pFile,"Total solved node: %d\n",solve_node);
		fprintf(pFile,"Preprocessing parameter:\n");
		fprintf(pFile,"-Feasibility_recovery_breath: %f\n",preprocess_param.feasibility_recovery_breath);
		fprintf(pFile,"-Feasibility_recovery_depth: %d\n",preprocess_param.feasibility_recovery_depth);
		fprintf(pFile,"-BoundCut generate loop: %d\n",preprocess_param.boundcut_loop);
		fprintf(pFile,"-BoundCut split: %d\n",preprocess_param.boundcut_p);
		fprintf(pFile,"-DisjunctiveCut generate loop: %d\n",preprocess_param.disjunctivecut_loop_param);
		fprintf(pFile,"-SimpleCut generate loop: %d\n",preprocess_param.simplecut_loop_param);
		fprintf(pFile,"-MccormickRefineCnt: %d\n",preprocess_param.MccormickRefineCnt);
		fprintf(pFile,"-------------------------------------------------------------------------\n");
		if(lp_relaxation!=MAX_VALUE || lp_relaxation!=MIN_VALUE )fprintf(pFile,"LP relaxation: %f\n",lp_relaxation);
		if(improved_lb!=MAX_VALUE || improved_lb!=MIN_VALUE)fprintf(pFile,"Improved lower bound:%f\n",improved_lb);
		switch (LPCC_condition)
		{
		case LPCC_UNBOUNDED:
			fprintf(pFile,"LPCC is unbounded\n");
			fprintf(pFile,"Unbouned ray:\n");
			for (i=0;i<param_n;i++)
			{
				fprintf(pFile,"x[%d]: %12.12f  ",i,lpcc_soln[i]);
				fprintf(pFile,"\n");
			}
			for (i=0;i<param_m;i++)
			{
				fprintf(pFile,"y[%d]: %12.12f  ",i,lpcc_soln[i+param_n]);
				fprintf(pFile,"\n");
			}
			break;
		case LPCC_INFEASIBLE:
			fprintf(pFile,"LPCC is infeasible\n");
			break;
		case LPCC_TERMINATE_WITHOUT_SOLN:
			fprintf(pFile,"Time limit reached, and no feasible solution\n");
			fprintf(pFile,"Best lower bound: %f\n",best_lb);
			break;
		case LPCC_TERMINATE_WITH_SOLN:
			fprintf(pFile,"Time limit reached, and feasible solution exists\n");
			fprintf(pFile,"Best lower bound: %f\n",best_lb);
			fprintf(pFile,"Best objective value; %f",opt_obj);
			fprintf(pFile,"Best solution:\n");
			for (i=0;i<param_n;i++)
			{
				fprintf(pFile,"x[%d]: %12.12f  ",i,lpcc_soln[i]);
				fprintf(pFile,"\n");
			}
			for (i=0;i<param_m;i++)
			{
				fprintf(pFile,"y[%d]: %12.12f  ",i,lpcc_soln[i+param_n]);
				fprintf(pFile,"\n");
			}
			break;
		case LPCC_OPTIMAL:
			fprintf(pFile,"LPCC is feasible and bounded\n");
			fprintf(pFile,"Optimal value: %f\n", opt_obj);
			fprintf(pFile,"Optimal solution:\n");
			for (i=0;i<param_n;i++)
			{
				fprintf(pFile,"x[%d]: %12.12f  ",i,lpcc_soln[i]);
				fprintf(pFile,"\n");
			}
			for (i=0;i<param_m;i++)
			{
				fprintf(pFile,"y[%d]: %12.12f  ",i,lpcc_soln[i+param_n]);
				fprintf(pFile,"\n");
			}
			break;
		}
		break;
	}
	fclose(pFile);
TERMINATE:
	return(status);
}

int printResult(  int LPCC_condition,
				  int param_n,	
				  int param_m,	
				  int param_k, 	
				  int solve_node,
				  PARAM preprocess_param,
				  double preprocessing_time,
				  double branching_time,
				  double wall_preproc_time,
				  double wall_solve_time,
				  double lp_relaxation,
				  double improved_lb,
				  double best_lb,
				  double opt_obj)
{

	int status=0;
	fprintf(stdout,"\n-------------------------------------------------------------------------\n");
	switch (preprocess_param.solve_mode)
	{
	case 'c':
		fprintf(stdout,"Solve LPCC by CPLEX.\n");
		fprintf(stdout,"Preprocessing time: %f\n",preprocessing_time);
		fprintf(stdout,"Solving time: %f\n",branching_time);
		fprintf(stdout,"Total time: %f\n",preprocessing_time+branching_time);
		fprintf(stdout,"Preprocessing time (wall): %.0f\n",wall_preproc_time);
		fprintf(stdout,"Solving time (wall): %.0f\n",wall_solve_time);
		fprintf(stdout,"Total time (wall): %.0f\n",wall_preproc_time+wall_solve_time);
		fprintf(stdout,"Total solved node: %d\n",solve_node);
		fprintf(stdout,"-------------------------------------------------------------------------\n");
		if(lp_relaxation!=MAX_VALUE || lp_relaxation!=MIN_VALUE )fprintf(stdout,"LP relaxation: %f\n",lp_relaxation);
		if(improved_lb!=MAX_VALUE || improved_lb!=MIN_VALUE)fprintf(stdout,"Improved lower bound:%f\n",improved_lb);
		switch(LPCC_condition)
		{
		case LPCC_INFEASIBLE:
			fprintf(stdout,"LPCC is infeasible\n");
			break;
		case LPCC_INForUNBD:
			fprintf(stdout,"LPCC is infeasible or unbounded\n");
			break;
		case CPXMIP_OPTIMAL_TOL:
		case LPCC_OPTIMAL:
			fprintf(stdout,"LPCC is feasible and bounded\n");
			fprintf(stdout,"Optimal value: %f\n", opt_obj);
			break;
		case LPCC_TERMINATE_WITH_SOLN:
			fprintf(stdout,"CPLEX time limit reached, and feasible solution exists\n");
			fprintf(stdout,"Best lower bound: %f\n",best_lb);
			fprintf(stdout,"Best objective value; %f",opt_obj);
			break;
		case LPCC_TERMINATE_WITHOUT_SOLN:
			fprintf(stdout,"CPLEX time limit reached, and no feasible solution\n");
			fprintf(stdout,"Best lower bound: %f\n",best_lb);
			break;
		case LPCC_UNBOUNDED:
			fprintf(stdout,"LPCC is unbounded\n");
			break;
        /*  case CPX_STAT_OPTIMAL:
            fprintf(stdout,"LPCC is feasible and bounded\n");
            fprintf(stdout,"Optimal value: %f\n", opt_obj);
            break;  */
		default:
			fprintf(stdout,"Unexpected MIP CPLEX solution stauts 4: %d\n",LPCC_condition);
			status=-1;
			goto TERMINATE;
		}
		break;
	case 'b':
		fprintf(stdout,"Solve LPCC by Branch and cut.\n");
		fprintf(stdout,"Preprocessing time: %f\n",preprocessing_time);
		fprintf(stdout,"Solving time: %f\n",branching_time);
		fprintf(stdout,"Total time: %f\n",preprocessing_time+branching_time);
		fprintf(stdout,"Preprocessing time (wall): %.0f\n",wall_preproc_time);
		fprintf(stdout,"Solving time (wall): %.0f\n",wall_solve_time);
		fprintf(stdout,"Total time (wall): %.0f\n",wall_preproc_time+wall_solve_time);
		fprintf(stdout,"Total solved node: %d\n",solve_node);
		fprintf(stdout,"-------------------------------------------------------------------------\n");
		if(lp_relaxation!=MAX_VALUE || lp_relaxation!=MIN_VALUE )fprintf(stdout,"LP relaxation: %f\n",lp_relaxation);
		if(improved_lb!=MAX_VALUE || improved_lb!=MIN_VALUE)fprintf(stdout,"Improved lower bound:%f\n",improved_lb);
		switch (LPCC_condition)
		{
		case LPCC_UNBOUNDED:
			fprintf(stdout,"LPCC is unbounded\n");
			break;
		case LPCC_INFEASIBLE:
			fprintf(stdout,"LPCC is infeasible\n");
			break;
		case LPCC_TERMINATE_WITHOUT_SOLN:
			fprintf(stdout,"Time limit reached, and no feasible solution\n");
			fprintf(stdout,"Best lower bound: %f\n",best_lb);
			break;
		case LPCC_TERMINATE_WITH_SOLN:
			fprintf(stdout,"Time limit reached, and feasible solution exists\n");
			fprintf(stdout,"Best lower bound: %f\n",best_lb);
			fprintf(stdout,"Best objective value; %f",opt_obj);
			break;
		case LPCC_OPTIMAL:
			fprintf(stdout,"LPCC is feasible and bounded\n");
			fprintf(stdout,"Optimal value: %f\n", opt_obj);
			break;
		}
		break;
	}
TERMINATE:
	return(status);
}
