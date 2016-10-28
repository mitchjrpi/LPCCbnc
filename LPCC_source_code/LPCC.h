#ifndef _LPCC_H
#define _LPCC_H
#include <ilcplex/cplex.h>
//define paremeter
#define ZERO_TOLERANCE								1e-06	
#define WEAKEN_SCALE								1e-08
#define STRENGTHEN_SCALE							1e-08
#define INF_BOUND									CPX_INFBOUND
#define MAX_VALUE									CPX_INFBOUND
#define MIN_VALUE									-CPX_INFBOUND
#define DEFAULT_FEASIBILITY_RECOVERY_BREATH			0.6
#define DEFAULT_FEASIBILITY_RECOVERY_DEPTH			5
#define DEFAULT_FEASIBILITY_RECOVERY_REFINEMENT		1
#define DEFAULT_FEASIBILITY_RECOVERY_MAXLOOP		10
#define DEFAULT_BOUNDCUT_MAXLOOP					3
#define DEFAULT_BOUNDCUT_P							1
#define DEFAULT_DISJUNCTIVECUT_MAXLOOP				0.1
#define DEFAULT_SIMPLECUT_MAXLOOP					0.1
#define DEFAULT_MCCORMICK_REFINEMENT_CNT			0
#define DEFAULT_TIME_LIMIT							3600
#define DEFAULT_OPTIMAL_TOLERANCE					1e-12
#define DEFAULT_PARTITION_X							1
#define DEFAULT_PARTITION_Y							1
//define error message						
#define NO_MEMORY									CPXERR_NO_MEMORY
//define solution status
//LP condition
#define STAT_OPTIMAL								CPX_STAT_OPTIMAL		
#define STAT_NUM_BEST								CPX_STAT_NUM_BEST	
#define STAT_INFEASIBLE								CPX_STAT_INFEASIBLE		
#define STAT_INForUNBD								CPX_STAT_INForUNBD		
#define STAT_UNBOUNDED								CPX_STAT_UNBOUNDED		
#define STAT_OPTIMAL_INFEAS							CPX_STAT_OPTIMAL_INFEAS	
//LP pre-condition
#define LP_PRECONDITION_BOUNDED						0	
#define LP_PRECONDITION_UNBOUNDED					1
#define LP_PRECONDITION_UNKNOWN						2
//LPCC pre-condition
#define LPCC_PRECONDITION_INFEASIBLE				-1
#define LPCC_PRECONDITION_BOUNDED					0
#define LPCC_PRECONDITION_OPT						1
#define LPCC_PRECONDITION_UNBOUNDED					2
//LPCC condition
#define LPCC_UNBOUNDED								-2
#define LPCC_INFEASIBLE								-1
#define LPCC_TERMINATE_WITHOUT_SOLN					0
#define LPCC_OPTIMAL								1
#define LPCC_TERMINATE_WITH_SOLN					2
#define LPCC_INForUNBD								3
//define data structure
struct matrix {
	//param=0: group by column
	//param=1: group by row
	int n_row;
	int n_col;
	int nnz;
	int param;
	int *matbeg;
	int *matcnt;
	int *matind;
	double *matval;
};
typedef struct matrix MATRIX;
struct startinfo {
	int ccnt;
	int *cindices;
	int *cstat;
	int rcnt;
	int *rindices;
	int *rstat;
};
typedef struct startinfo STARTINFO;
struct constraint_set {
	//constraint matrix group by row
	MATRIX cons_matrix;
	double *rhs;
	char *sen;
};
typedef struct constraint_set CONSTRAINT_SET;
struct node {
	STARTINFO nodestartinfo;
	CONSTRAINT_SET nodecuts;
	int n;
	int m;
	int k;
	int *nodebranchinfo;
	int branchIndexofParentNode;
	int branchDirectionofParentNode;
	int x_start_index;
	int y_bar_start_index;
	int mcc_index_row;
	int mcc_index_col;
	double *x_lb;
	double *x_ub;
	double *y_bar_lb;
	double *y_bar_ub;
	double parentNodeValue;
	double parentNodeBranchReduceCostNorm;
	double parentNodeBranchCoefNorm;
	double parentNodeBranchViolation;
	double *nodex;
	double nodevalue;
	double nodelb;
	struct node *link;
};
typedef struct node NODE;
struct param
{
	double feasibility_recovery_breath;
	int feasibility_recovery_depth;
	int feasibility_recovery_refinement;
	int feasibility_recovery_maxloop;
	int boundcut_loop;
	int boundcut_p;
	double disjunctivecut_loop_param;
	double simplecut_loop_param;
	int MccormickRefineCnt;
	char solve_mode;
	int format_mode;
	char profilefile[1024];
	char inputfile[1024];
	char solnfile[1024];
	double time_limit;
	double opt_tolerance;
	int partition_x;
	int partition_y;
};
typedef struct param PARAM;
struct branchHist
{
	int branch_var_cnt;
	long int *branch_1_cnt;
	long int *branch_2_cnt;
	double *average_violation_1_cnt;
	double *average_violation_2_cnt;
	double *average_obj_improve_1;
	double *average_obj_improve_2;
};
typedef struct branchHist BRANCHHIST;
#endif

