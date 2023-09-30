#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double **sym_c(double **pointsMatrix, int N, int dim);   /* TESTED UNTIL EXAMPLE 27 */
double **ddg_c(double **pointsMatrix, int N, int dim);   /* TESTED UNTIL EXAMPLE 27 */
double **norm_c(double **pointsMatrix, int N, int dim);   /* TESTED UNTIL EXAMPLE 27 */
double **symnmf_c(double **initialH, double **W, int N, int k);   /* TESTED UNTIL EXAMPLE 27 */

int my_strcmp(char *str1, char *str2);     /* TESTED */
void ptintMatrix(double **mtrx, int rows, int cols);   /* TESTED */
double squaredDistance(double *p1, double *p2, int dim);   /* TESTED */
double calcSqueredFronenius(double **A, int rows, int cols);   /* WASNT TESTED BUT WAS VERIFIED WITH OTHED CODE */
double **mallocMatrix(int rows, int cols);   /* WASNT TESTED BUT WAS VERIFIED WITH OTHED CODE */
void freeMatrix(double **mtrx, int rows);   /* WASNT TESTED BUT WAS VERIFIED WITH OTHED CODE */
int main(int argc, char **argv);

int main(int argc, char **argv){
    int i, rows, cols, row, col;
    double n;
    char c;
    double **dataPointsMatrix, **returnMatrix;
    char *goal;
    FILE *file;

    if (argc != 3)
    {
        printf("An Error Has Occurred");
        exit(1);
    }

    goal = argv[1];
    file = fopen(argv[2], "r");

    if (file == NULL)
    {
        printf("An Error Has Occurred");
        exit(1);
    }

    /* Count the number of lines in the file to determine the matrix size */
    rows = 0;
    cols = 0;
    while (fscanf(file, "%lf%c", &n, &c) == 2)
    {
        if (rows == 0){
            cols++;
        }
        if (c == '\n'){
            rows++;
        }
    }

    /* Rewind the file pointer to read from the beginning again */
    rewind(file);

    /* Allocate memory for the matrix that will contain the input */          /* MAYBE CHANCE THE NUMBER OF COLUMNS */
    dataPointsMatrix = (double **)malloc(rows * sizeof(double *));
    if (dataPointsMatrix == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for (i = 0; i < rows; i++)
    {
        dataPointsMatrix[i] = (double *)malloc(cols * sizeof(double));
        if (dataPointsMatrix[i] == NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }

    /* Rewind the file pointer to read from the beginning again */
    rewind(file);

    /* Read the file and populate the matrix (with the input data) */
    row = 0, col = 0;
    while (fscanf(file, "%lf%c", &n, &c) == 2)
    {   
        dataPointsMatrix[row][col] = n;
        col++;
        if (c == '\n')
        {
            row++;
            col = 0;
        }
    }

    /* Close the file when finished */
    fclose(file);
    
    /* Print a matrix accoarding to the value of "goal" */
    if(my_strcmp(goal, "sym") == 0){
        returnMatrix = sym_c(dataPointsMatrix, rows, cols);
        ptintMatrix(returnMatrix, rows, rows);
        freeMatrix(returnMatrix, rows);
    }
    if(my_strcmp(goal, "ddg") == 0){
        returnMatrix = ddg_c(dataPointsMatrix, rows, cols);
        ptintMatrix(returnMatrix, rows, rows);
        freeMatrix(returnMatrix, rows);
    }
    if(my_strcmp(goal, "norm") == 0){
        returnMatrix = norm_c(dataPointsMatrix, rows, cols);
        ptintMatrix(returnMatrix, rows, rows);
        freeMatrix(returnMatrix, rows);
    }

    /* Free memory used for the dataPointsMatrix */
    freeMatrix(dataPointsMatrix, rows);

    return 0;
}


int my_strcmp(char *str1, char *str2) {
/* returns 0 if srt1=str2 */
/* reutns negative int if str1 comes before str2 in lexicographic order */
/* else returns a positive int. */
    while (*str1 != '\0' && *str2 != '\0') {
        if (*str1 != *str2) {
            return *str1 - *str2;
        }
        str1++;
        str2++;
    }
    return *str1 - *str2;
}

void ptintMatrix(double **mtrx, int rows, int cols){
    int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            if (j < cols-1){
                printf("%.4f,", mtrx[i][j]);
            }
            else{
                printf("%.4f", mtrx[i][j]);
            }

        }
        printf("\n");
    }
}

double squaredDistance( double *p1, double *p2, int dim) {
    /* Returns the squered Euqlidian distance */
    int i;
    double sum, diff;
    sum = 0.0;
    for (i = 0; i < dim; i++) {
        diff = p1[i] - p2[i];
        sum += pow(diff, 2);
    }
    return sum;
}

double calcSqueredFronenius(double **A, int rows, int cols){
    int i, j;
    double norm = 0;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            norm += pow(A[i][j], 2);
            }
        }
    return norm;
}

double **mallocMatrix(int rows, int cols){
    int i, j;
    double **A = (double **)malloc(rows * sizeof(double *));
    if (A == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for (i = 0; i < rows; i++) {
        A[i] = (double *)malloc(cols * sizeof(double));
        if (A[i] == NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
        for(j = 0; j < cols; j++) {
            A[i][j] = 0.0;
        }
    }
    return A;
}

void freeMatrix(double **mtrx, int rows){
    int i;
    for (i = 0; i < rows; i++){
        free(mtrx[i]);
    }
    free(mtrx);
}

double **sym_c(double **pointsMatrix, int N, int dim) {
    int i, j;

    /* Allocate memory for similarityMatrix */
    double **similarityMatrix = mallocMatrix(N, N);

    /* Populate similarityMatrix */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (i != j) {
                double dist = squaredDistance(pointsMatrix[i], pointsMatrix[j], dim);
                similarityMatrix[i][j] = exp((-dist) / 2.0);
            } else {
                similarityMatrix[i][j] = 0.0;
            }
        }
    }

    return similarityMatrix;
}

double **ddg_c(double **pointsMatrix, int N, int dim) {
    int i, j;
    double degree;
    double **similarityMatrix, **degreeMatrix;

    /* Calculate the similarity matrix A using sym_c function */
    similarityMatrix = sym_c(pointsMatrix, N, dim);

    /* Allocate memory for the degree matrix D */
    degreeMatrix = mallocMatrix(N, N);

    /* Compute the degrees d1,...,dn and form the degree matrix */
    for (i = 0; i < N; i++) {
        degree = 0.0;
        for (j = 0; j < N; j++) {
            degree += similarityMatrix[i][j];
        }
        degreeMatrix[i][i] = degree;
    }

    /* Free memory used for the similarity matrix A */
    freeMatrix(similarityMatrix, N);

    return degreeMatrix;
}

double **norm_c(double **pointsMatrix, int N, int dim){
    int i, j;
    double d;

    /* Calculate the similarity matrix A using sym_c function */
    double **similarityMatrix = sym_c(pointsMatrix, N, dim);

    /* Calculate the diagonal degree matrix D using ddg_c function */
    double **degreeMatrix = ddg_c(pointsMatrix, N, dim);

    /* Allocate memory for the normalized matrix W */
    double **normalizedMatrix = mallocMatrix(N, N);

    /* Populate normalizedMatrix W */

    /* left multiplation ( D^(-0.5) * A ) */
    for (i = 0; i < N; i++) {
        d = pow(degreeMatrix[i][i], -0.5);
        for (j = 0; j < N; j++) {
            normalizedMatrix[i][j] = similarityMatrix[i][j]*d;
        }
    }

    /* right multiplation ( (D^(-0.5)*A) * D^(-0.5) ) */
    for (j = 0; j < N; j++) {
        d = pow(degreeMatrix[j][j], -0.5);
        for (i = 0; i < N; i++) {
            normalizedMatrix[i][j] = normalizedMatrix[i][j]*d;
        }
    }

    freeMatrix(similarityMatrix, N);
    freeMatrix(degreeMatrix, N);

    return normalizedMatrix;
}


double **symnmf_c(double **initialH, double **W, int N, int k){
    double **prevH, **nextH, **WH, **HHTH, **temp, **conditor;
    int i, j, m, iter;
    double beta, c, squeredFronenius, eps;

    beta = 0.5;
    eps = 0.0001;

    /* Allocate memory for prevH, nextH, WH, HHTH */
    prevH = mallocMatrix(N, k);
    nextH = mallocMatrix(N, k);
    WH = mallocMatrix(N, k);
    HHTH = mallocMatrix(N, k);
    temp = mallocMatrix(N, N);
    conditor = mallocMatrix(N, k);


    /* Initialize first prevH to be initialH */
    for (i = 0; i < N; i++) {
        for (j = 0; j < k; j++) {
            prevH[i][j] = initialH[i][j];
        }
    }

    /* Iterate to compute H(i) */
    for (iter = 0; iter < 300; iter++){
        /* Compmute WH (H is the prevH) */
        for (i = 0; i < N; i++) {        /* CHECKED */
            for (j = 0; j < k; j++) {
                WH[i][j] = 0.0;
                for(m = 0; m < N; m++){
                    WH[i][j] += W[i][m] * prevH[m][j];
                }
            }
        }
        
        /* Copmute HHTH (H is the prevH) */
        /* temp = (H)(HT) */
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                temp[i][j] = 0.0;
                for(m = 0; m < k; m++){
                    temp[i][j] += prevH[i][m] * prevH[j][m];
                }
            }
        }

        /* (HHT)H (computes temp*H) */
        for (i = 0; i < N; i++) {
            for (j = 0; j < k; j++) {
                HHTH[i][j]=0.0;
                for(m = 0; m < N; m++){
                    HHTH[i][j] += temp[i][m] * prevH[m][j];
                }
            }
        }

        /* Populate nextH */
        for (i = 0; i < N; i++) {
            for (j = 0; j < k; j++) {
                nextH[i][j] = 0.0;
                c = (WH[i][j])/(HHTH[i][j]);
                nextH[i][j] = prevH[i][j] * (1.0 - beta + beta*c);
            }
        }

        /* Check stopping condition */
        for (i = 0; i < N; i++) {
            for (j = 0; j < k; j++) {
                conditor[i][j] = nextH[i][j]-prevH[i][j];
            }
        }

        squeredFronenius = calcSqueredFronenius(conditor, N, k);
        if(squeredFronenius < eps){
            break;
        }

        /* update prevH */
        for (i = 0; i < N; i++) {
            for (j = 0; j < k; j++) {
                prevH[i][j] = nextH[i][j];
            }
        }

    }
    /* Free memory */
    freeMatrix(prevH, N);
    freeMatrix(WH, N);
    freeMatrix(HHTH, N);
    freeMatrix(temp, N);
    freeMatrix(conditor, N);
    

    return nextH;
}



/*from now on the code is for kmeans (for the analysis part)*/

/*this fucntion returns the index of the most close centroid to the vector x_i*/
int CLASSIFY_CENTROID(double* X_i,int K, int dimension,double* CENTROIDS){
    int RESULT = -1;
    double MINIMUM = __DBL_MAX__;
    double MINIMUM_CONTENDER;
    int j = 0;
    double* centroids_number_j_ptr;

    for( j = 0; j < K ; j++) {

    centroids_number_j_ptr = CENTROIDS + j*(dimension);
        
        MINIMUM_CONTENDER = squaredDistance(X_i , centroids_number_j_ptr , dimension);

        if (MINIMUM_CONTENDER < MINIMUM) {
            MINIMUM = MINIMUM_CONTENDER;
            RESULT = j;
        }
        
    }

    return RESULT;

}

/*this function will update the new Centroids, we will use the funcion CLASSIFY_CENTROID so we know which vector we want to add to each centroid
the function will return 0 if for each centroid : LENGTH(new_centroid[i],centroid[i]) < epsilon , else it will return 1
as a result , we will know we arrived to the wanted precision of centroids and we can finifsh the program*/

int UPDATE_CENTROIDS(double** matrix, int k, int dimension,int num_rows ,double* centroids ,double* new_centroids, int* num_of_elem_in_cluster,double eps){
    int enough_precision = 0;
 
    long double EPSILON = eps;

    int i = 0;
    int j = 0;
    int a = 0;
    int b = 0;
    int index = 0;



    for ( i = 0; i < k ; i++){
        num_of_elem_in_cluster[i] = 0;
    }

    
    for( a = 0 ; a < k ; a++){
            for ( b = 0; b < dimension ; b ++){
                new_centroids[a*(dimension) + b] = 0;
            }
    }


    for( i = 0 ; i < num_rows ; i++ ){

        index = CLASSIFY_CENTROID(matrix[i],k,dimension,centroids); /*identifying the index of the closest centroid and adding it to the new centroid*/
        num_of_elem_in_cluster[index] +=1;

        for ( j = 0; j < dimension ; j++){
            new_centroids[index*(dimension) + j] += matrix[i][j];

        }

    }



    for ( i = 0 ; i < k ; i++){
        for (  j =0 ; j < dimension ; j++){
            new_centroids[i*(dimension) + j] = ( (new_centroids[i*(dimension) + j])*((double) 1 / (double) ( num_of_elem_in_cluster[i])) );

        }
    }

        for( i = 0; i < k ; i++){
            if (squaredDistance(centroids + i*dimension , new_centroids + i*dimension , dimension) >= EPSILON) {
                enough_precision = 1;

            }
        }
    
    for ( i = 0 ; i < k ; i++){
        for (  j =0 ; j < dimension ; j++){
            centroids[i*(dimension) + j] = ((new_centroids[i*(dimension) + j]));

        }

    }

return enough_precision;

}

/*this function will return the k centroids that will satisfy the parameters
it will use the tool functions : LENGTH ,UPDATE_CENTROIDS classify_CENTOIRDS*/
double* kmeans_c(int Pyk, int Pyiter,int Pynumofrows,int Pydim ,double Pyepsilon,double** C_Matrix, double* C_Centroids){

    int k = Pyk;
    int iteration = Pyiter;
    int num_of_rows = Pynumofrows;
    int dimension = Pydim;
    double epsilon = Pyepsilon*Pyepsilon;
    int i = 0;

    int precision = 1;
    int* num_of_rows_ptr = & num_of_rows;
    int* dimension_ptr = & dimension;

    double** MATRIX =  C_Matrix;
    double* centroids = C_Centroids;
    double* new_centroids;
    int* num_of_elem_in_kluster;
    int kluster_times_dimensions;

    kluster_times_dimensions = *(dimension_ptr)*k;

    new_centroids = malloc(kluster_times_dimensions*sizeof(double));

      if (new_centroids ==NULL){
        printf("failed to allocate memory new centroids\n");
        exit(1);
    }

    num_of_elem_in_kluster = malloc(k*sizeof(int));

      if (num_of_elem_in_kluster ==NULL){
        printf("failed to allocate memory klusters\n");
        exit(1);
    }

    i  = 0;
  
    while ( i < iteration && precision == 1){ // the algorithm
        precision = UPDATE_CENTROIDS(MATRIX,k,*(dimension_ptr),*(num_of_rows_ptr),centroids,new_centroids,num_of_elem_in_kluster,epsilon);
        i +=1;
    
    }
   
    for( i = 0; i < *(num_of_rows_ptr) ; i++) {
    free(MATRIX[i]);
    }
    
    free(num_of_elem_in_kluster);
    free(new_centroids);
    free(MATRIX);

    return centroids;
}
