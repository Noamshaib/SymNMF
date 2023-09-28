# ifndef SYM_NMF
# define SYM_NMF

double **symnmf_c(double **H, double **W, int k, int N, int dim);
double **sym_c(double **X, int N, int dim);
double **ddg_c(double **X, int N, int dim);
double **norm_c(double **X, int N, int dim);
void free_Matrix(double **head);
void print_Matrix(double **head); /*only to check the code, printing in this project is from s-f.py or s-f.c*/


# endif
