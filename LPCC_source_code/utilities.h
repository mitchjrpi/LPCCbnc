#ifndef _UTILITIES_H
#define _UTILITIES_H
int get_matrix(const double **matrix_A,
			   const int n_row,
			   const int n_col, 
			   int param, 
			   MATRIX *matrix_A_p);
int print_matrix ( MATRIX *matrix_A );
int print_full_matrix(MATRIX *matrix_A_p);
int transpose_matrix (MATRIX *matrix_A_p);
int get_matrix_element (MATRIX *matrix_A_p,
						const int index_i, 
						const int index_j, 
						double *val_p);
void init_matrix ( MATRIX *ptr);
void free_matrix ( MATRIX *ptr);
int combine_matrix(MATRIX matrix_A, 
				   MATRIX matrix_B, 
				   MATRIX *matrix_AB, 
				   int method );
int zero_matrix (const int param,
				 const int n_row,
				 const int n_col,
				 MATRIX *matrix_A);
int diag_matrix (const int param,
				 const int matrix_size,
				 const double *val,
				 MATRIX *matrix_A);
int identity_matrix (const int param,
					 const int matrix_size,
					 const double val,
					 MATRIX *matrix_A);
int copy_matrix (MATRIX *matrix_A_p, 
				 MATRIX *matrix_A_cpy_p);
void init_startinfo ( STARTINFO *ptr);
void free_startinfo ( STARTINFO *ptr);
int copy_startinfo (STARTINFO *start_p, 
					STARTINFO *start_cpy_p);
void init_constraint_set ( CONSTRAINT_SET *ptr);
void free_constraint_set ( CONSTRAINT_SET *ptr);
int copy_constraint_set (CONSTRAINT_SET *cons_p, 
						 CONSTRAINT_SET *cons_cpy_p);
int get_constraint(CONSTRAINT_SET *cons_p,
				   const int index, 
				   CONSTRAINT_SET **cons_i_p);
int add_constraint(CONSTRAINT_SET *ptr,
				   int nnz, 
				   int *ind, 
				   double *val,
				   double rhs);
int add_cuts(CONSTRAINT_SET *ptr,
			 int new_cols, 
			 int nnz, 
			 int *ind,
			 double *val,
			 double rhs) ;
int add_constraint_a (CONSTRAINT_SET *ptr,
					  int cnt,
					  double *original_array,
					  double rhs);
int convertArray (int cnt,
				  double *original_array,
				  int *nnz_p, 
				  int **ind_p, 
				  double ** val_p);
int checkComplementary(int index, double *x);
int checkallComplementary(const int param_n,
						  const int param_m,
						  const int param_k, 
						  double *x) ;
int countallComplementary(const int param_n,
						  const int param_m,
						  const int param_k,
						  double *x) ;
double weaken_value(const double value);
double strengthen_value(const double value);
void free_and_null (char **ptr);
int pushnode_ordered(NODE **head_ptr,
					 NODE newnode);
int pushnode_unordered( NODE **head_ptr,
					   NODE newnode);
void popnode(NODE **head_ptr, 
			 NODE **topnode_ptr);
int copy_node(NODE *o_node_p,
			 NODE *cpynode_ptr);
int init_branch_node(NODE activenode,
					 NODE **branchnode_ptr);
int init_branch_nodes(NODE activenode,
					  NODE **branchnode_1_ptr, 
					  NODE **branchnode_2_ptr);
void deletenode(NODE **node_ptr);
void filternode(NODE **head_ptr, 
				double val);
int fixedcnt_node(NODE *node_ptr);
void init_node(NODE *ptr);
void free_node( NODE *node_ptr);
void displaynode( NODE *head_ptr);
double getLastnodevalue(NODE **head_ptr);
double getaveragenodevalue(NODE **head_ptr);
double getmediannodevalue(NODE **head_ptr);
void  usage(char *progname);
int readarray(FILE *in, 
			  int *num_p,
			  double **data_p);
int readarray_int (FILE *in,
				   int *num_p,
				   int **data_p);
int readarray_row (FILE *in,
				 int *num_p,
				 double **data_p);
int readarray_col(FILE *in,
				  int *num_p,
				  double **data_p);
int get_sort_index(const int cnt,
				   double *sort_val,
				   int *sorted_index,
				   const int sort_order);
double sum_max_t(const int cnt,
				 double *val,
				 const int t);
int symmetrizing_matrix(MATRIX *matrix_A_p, 
						MATRIX *sym_matrix_A_p, 
						int *sym_status);
int matrix_addition(MATRIX matrix_A, 
					MATRIX matrix_B,
					MATRIX *matrix_sum);
int test_matrix_addition_symmetrize();
void vector_sum(double *vector_1,
				double *vector_2,
				double cnt,
				double *sum_vector);
void init_branchhist ( BRANCHHIST *ptr);
void free_branchhist ( BRANCHHIST *ptr);
int init_branchhist_info(BRANCHHIST *ptr, 
						 int branch_var_cnt);
int get_branchhist_info(BRANCHHIST *ptr,
						int index, 
						int direction, 
						double *violation_cost, 
						double *improve_cost);
int update_branchhist_info(BRANCHHIST *ptr, 
						   int index, 
						   int direction,
						   double average_violation, 
						   double average_distance);
double vectorproduct(double **Vector_1, 
					 double **Vector_2, 
					 const int dimension, 
					 const int Vector_1_index, 
					 const int Vector_2_index);
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
					int *total_partition_cnt);
int readProfile(PARAM *parameter);
void exit_with_help();
void parse_command_line(int argc, char **argv, PARAM *parameter);
void printParam(PARAM parameter);
#endif

