# ifndef SYM_NMF
# define SYM_NMF

double **symnmf_c(double **pointsMatrix, int N, int dim, int k);
double **sym_c(double **pointsMatrix, int N, int dim);
double **ddg_c(double **pointsMatrix, int N, int dim);
double **norm_c(double **pointsMatrix, double **pointsMatrix, int N, int dim);
void freeMatrix(double **Matrix);


# endif