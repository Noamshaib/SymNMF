#include <stdio.h>


double **defX(char *filename);
void printMatrix(double **Matrix);
double **calH(double **normMatrix, int k);

double **symnmf_c(double **H0_Matrix,double **WMatrix, int N, int dim, int k);
double **sym_c(double **pointsMatrix, int N, int dim);
double **ddg_c(double **pointsMatrix, int N, int dim);
double **norm_c(double **pointsMatrix, int N, int dim);

double **symnmf_c(double **H0_Matrix,double **WMatrix, int N, int dim, int k){
    double **A, **D, **W;
    int N, dim;

    

    return calH(W, k);    
}



void main(int argc, char* argv[]) {
    double **X, **A, **D, **W;
    int N, dim;

    /* command-line arguments */
    char* goal = argv[1];
    char* inputfile = argv[2];

    X = defX(inputfile);
    if(goal == "sym"){
        printMatrix(sym_c(X, N, dim));
    }
    else{
        if(goal == "ddg"){
            printMatrix(ddg_c(X, N, dim));
        }
        else{   /* goal == "norm" */
            printMatrix(norm_c(X,  N, dim));
        }
    }
    
}
