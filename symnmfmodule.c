# define PY_SSIZE_T_CLEAN
#include<stdlib.h>
#include<Python.h>
#include "symnmf.h"

static PyObject *symnmf(PyObject *self, PyObject* args);
static PyObject *sym(PyObject *self, PyObject* args);
static PyObject *ddg(PyObject *self, PyObject* args);
static PyObject *norm(PyObject *self, PyObject* args);

PyMODINIT_FUNC PyInit_mysymnmfsp(void);
static double **convert_PyMatrix_To_CMatrix(PyObject *pyMatrix, int N, int dim);
static PyObject *convert_CMatrix_To_PyMatrix(double **ret_matrix, int N, int dim);


/*symnmf called from symnmf.py with the arguments - H0, W, dim, k, N */
static PyObject *symnmf(PyObject *self, PyObject *args) {
    PyObject *pyW;                          /* PyObject* - normal Matrix */
    PyObject *pyH;                          /* PyObject* Matrix H */
    int k, N, dim;                          /* C Objects - int arguments */
    double **W, **H;                    /* C Objects - W, H Matrix*/

    
    /* Parse the arguments to retrieve the Python objects and additional int arguments */
    if (!PyArg_ParseTuple(args, "OOiii", &pyW, &pyH, &k, &N, &dim)) {
        return NULL;  
    }

    /* Convert the PyObject* (list of lists) to Cmatrix* (double**) using appropriate conversion functions */
    W = convert_PyMatrix_To_CMatrix(pyW, N, N);
    H = convert_PyMatrix_To_CMatrix(pyH, N, k);

    /* Calculates the H matrix from symnmf algo using func symnmf_c from symnmf.c */
    H = symnmf_c(H, W, N, k);
    
    /* Convert the Cmatrix* (double**) to PyObject* (list of lists) using appropriate conversion functions*/
    pyH = convert_CMatrix_To_PyMatrix(H, N, k);
    
    /* Cleanup memory */
    freeMatrix(W, N);
    freeMatrix(H, N);
    
    /* Return the H matrix as a PyObject* (list of lists) */
    return pyH;
}

/*sym called from symnmf.py with the arguments - data point, dim, N */
static PyObject *sym(PyObject *self, PyObject *args) {
    PyObject *pyX, *pyA;                /* PyObject* - similarity Matrix */
    int N, dim;                         /* C Objects - int arguments */
    double **X, **A;                    /* C Objects - X, A Matrix*/
    
    
    /* Parse the arguments to retrieve the Python objects and additional int arguments */
    if (!PyArg_ParseTuple(args, "Oii", &pyX, &N, &dim)) {
        return NULL;  
    }

    /* Convert the PyObject* (list of lists) to Cmatrix* (double**) using appropriate conversion functions */
    X = convert_PyMatrix_To_CMatrix(pyX, N, dim);

    /* Calculates the similarity matrix using func sym_c from symnmf.c */
    A = sym_c(X, N, dim);
    
    /* Convert the Cmatrix* (double**) to PyObject* (list of lists) using appropriate conversion functions*/
    pyA = convert_CMatrix_To_PyMatrix(A, N, N);

    /* Cleanup memory */
    freeMatrix(X, N);
    freeMatrix(A, N);
    
    /* Return the A matrix as a PyObject* (list of lists) */
    return pyA;
}

/*ddg called from symnmf.py with the arguments - data point, dim, N */
static PyObject *ddg(PyObject *self, PyObject *args) {
    PyObject *pyX, *pyD;                /* PyObject* - diagonal degree Matrix */
    int N, dim;                         /* C Objects - int arguments */
    double **X, **D;                    /* C Objects - X, D Matrix*/
    
    
    /* Parse the arguments to retrieve the Python objects and additional int arguments */
    if (!PyArg_ParseTuple(args, "Oii", &pyX, &N, &dim)) {
        return NULL;  
    }

    /* Convert the PyObject* (list of lists) to Cmatrix* (double**) using appropriate conversion functions */
    X = convert_PyMatrix_To_CMatrix(pyX, N, dim);

    /* Calculates the diagonal degree matrix using func ddg_c from symnmf.c */
    D = ddg_c(X, N, dim);
    
    /* Convert the Cmatrix* (double**) to PyObject* (list of lists) using appropriate conversion functions*/
    pyD = convert_CMatrix_To_PyMatrix(D, N, N);
    
    /* Cleanup memory */
    freeMatrix(X, N);
    freeMatrix(D, N);
    
    /* Return the D matrix as a PyObject* (list of lists) */
    return pyD;
}

/*norm called from symnmf.py with the arguments - data point, dim, N */
static PyObject *norm(PyObject *self, PyObject *args) {
    PyObject *pyX, *pyW;                /* PyObject* - normal Matrix */
    int N, dim;                         /* C Objects - int arguments */
    double **X, **W;                    /* C Objects - X, W Matrix*/
    
    
    /* Parse the arguments to retrieve the Python objects and additional int arguments */
    if (!PyArg_ParseTuple(args, "Oii", &pyX, &N, &dim)) {
        return NULL;  
    }

    /* Convert the PyObject* (list of lists) to Cmatrix* (double**) using appropriate conversion functions */
    X = convert_PyMatrix_To_CMatrix(pyX, N, dim);

    /* Calculates the normal matrix using func norm_c from symnmf.c */
    W = norm_c(X, N, dim);
    
    /* Convert the Cmatrix* (double**) to PyObject* (list of lists) using appropriate conversion functions*/
    pyW = convert_CMatrix_To_PyMatrix(W, N, N);
    
    /* Cleanup memory */
    freeMatrix(X, N);
    freeMatrix(W, N);
    
    /* Return the D matrix as a PyObject* (list of lists) */
    return pyW;
}

static PyMethodDef module_methods[] = {
    { "symnmf", (PyCFunction) symnmf, METH_VARARGS, "Perform full the symNMF and output H" },
    { "sym", (PyCFunction) sym, METH_VARARGS, "Calculate and output the similarity matrix" },
    { "ddg", (PyCFunction) ddg, METH_VARARGS, "Calculate and output the Diagonal Degree Matrix" },
    { "norm", (PyCFunction) norm, METH_VARARGS, "Calculate and output the normalized similarity matrix" },
    { NULL, NULL, 0, NULL }  /* Sentinel indicating the end of the array */
};


static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "mysymnmfsp",       /* name of module */
    NULL,               /* module documentation */
    -1,                 /* the module keeps state in global variables */
    module_methods      /* the PyMethodDef array */
};

PyMODINIT_FUNC PyInit_mysymnmfsp(void){
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m) {
        return NULL;
    }
    return m;
}

static double **convert_PyMatrix_To_CMatrix(PyObject *pyMatrix, int m, int n){
    PyObject *curr_lst, *curr_cord;
    double **CMatrix;
    int i,j;
    
    /* providing space for the C Matrix - mXn*/
    CMatrix = (double**) malloc(sizeof(double **)*(m));
    if (CMatrix == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    
    /* fill the C matrix with values of py matrix */
    for (i = 0; i < m; i++){

        curr_lst = PyList_GET_ITEM(pyMatrix,i);     
        CMatrix[i] = (double*) malloc(sizeof(double *)*(n));
        if (CMatrix[i] == NULL){
            printf("An Error Has Occurred\n");
            exit(1);
        }

        for(j = 0; j < n; j++){
            curr_cord = (PyObject *)PyList_GET_ITEM(curr_lst,j);
            CMatrix[i][j] = PyFloat_AsDouble(curr_cord);
        }

    }

    return CMatrix;
}

static PyObject *convert_CMatrix_To_PyMatrix(double **CMatrix, int m, int n){
    PyObject *pyMatrix, *curr_lst, *pycord;
    int i,j;
    
    /* Create the Python list to hold the pyMatrix (list of lists) */
    pyMatrix = PyList_New(m);  

    for(i = 0 ; i < m ; i++){

        curr_lst = PyList_New(n);
        for (j = 0 ; j < n ; j++){
            pycord = Py_BuildValue("d", CMatrix[i][j]);
            PyList_SetItem(curr_lst,j,pycord);
        }

        PyList_SetItem(pyMatrix,i,curr_lst);        
    }
    return pyMatrix;
}
