# define PY_SSIZE_T_CLEAN
#include<stdlib.h>
#include<Python.h>
#include "symnmf.h"

/*initiallize new classes - used as linked list*/

static PyObject *symnmf(PyObject *self, PyObject* args);
static PyObject *sym(PyObject *self, PyObject* args);
static PyObject *ddg(PyObject *self, PyObject* args);
static PyObject *norm(PyObject *self, PyObject* args);

PyMODINIT_FUNC PyInit_mykmeanssp(void);
static void free_Matrix(double **head);
static double **convert_PyMatrix_To_CMatrix(PyObject *pyMatrix, int N, int dim);
static PyObject *convert_CMatrix_To_PyObject(double **ret_matrix, int N, int dim);
static void print_Matrix(double **head);


/*symnmf called from symnmf.py with the arguments - data point, k, dim, N */
static PyObject *symnmf(PyObject *self, PyObject *args) {
    PyObject *pyPointsMatrix;             /* PyObject* - Matrix of points */
    PyObject *pyH;                        /* PyObject* Matrix H */
    int k;                                /* C Objects - int arguments */
    double **pointsMatrix, **ret_H;       /* C Objects - Matrix pointer*/

    
    /*providing space for pointsMatrix head*/
    pointsMatrix = malloc(sizeof(double));
    if (pointsMatrix == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    
    /* Parse the arguments to retrieve the Python objects and additional int arguments */
    if (!PyArg_ParseTuple(args, "Osi", &pyPointsMatrix, &k)) {
        return NULL;  
    }

    /* Convert the PyObject* to struct vector* using appropriate conversion functions */
    pointsMatrix = convert_PyMatrix_To_CMatrix(pyPointsMatrix, M, dim);

    /* Perform the fitting algorithm using the extracted data and additional arguments*/
    ret_H = symnmf_c(pointsMatrix, k);
    

    /* Convert the struct vector* to PyObject* using appropriate conversion functions*/
    pyH = convert_CMatrix_To_PyObject(ret_H, k, dim);
    
    /* Cleanup memory */
    free_Matrix(pointsMatrix);
    

    /* Return the new centroids list as a PyObject* */
    return pyH;

}
/*sym called from symnmf.py with the arguments - data point, k, dim, N */
static PyObject *symnmf(PyObject *self, PyObject *args) {
    PyObject *pyPointsMatrix;             /* PyObject* - Matrix of points */
    PyObject *pyH;                        /* PyObject* Matrix H */
    int k;                                /* C Objects - int arguments */
    double **pointsMatrix, **ret_H;       /* C Objects - Matrix pointer*/

    
    /*providing space for pointsMatrix head*/
    pointsMatrix = malloc(sizeof(double));
    if (pointsMatrix == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    
    /* Parse the arguments to retrieve the Python objects and additional int arguments */
    if (!PyArg_ParseTuple(args, "Osi", &pyPointsMatrix, &k)) {
        return NULL;  
    }

    /* Convert the PyObject* to struct vector* using appropriate conversion functions */
    pointsMatrix = convert_PyMatrix_To_CMatrix(pyPointsMatrix, M, dim);

    /* Perform the fitting algorithm using the extracted data and additional arguments*/
    ret_H = sym_c(pointsMatrix, k);
    

    /* Convert the struct vector* to PyObject* using appropriate conversion functions*/
    pyH = convert_CMatrix_To_PyObject(ret_H, k, dim);
    
    /* Cleanup memory */
    free_Matrix(pointsMatrix);
    

    /* Return the new centroids list as a PyObject* */
    return pyH;

}
static PyMethodDef module_methods[] = {
    { "symnmf", (PyCFunction) symnmf, METH_VARARGS, "Perform full the symNMF and output H" },
    { "sym", (PyCFunction) sym, METH_VARARGS, "Calculate and output the similarity matrix" },
    { "ddg", (PyCFunction) ddg, METH_VARARGS, "Calculate and output the Diagonal Degree Matrix" },
    { "norm", (PyCFunction) norm, METH_VARARGS, "Calculate and output the normalized similarity matrix" },
    { NULL, NULL, 0, NULL }  /* Sentinel indicating the end of the array */
};

static struct PyModuleDef fitmodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",       /* name of module */
    NULL,               /* module documentation */
    -1,                 /* the module keeps state in global variables */
    module_methods      /* the PyMethodDef array */
};

PyMODINIT_FUNC PyInit_mykmeanssp(void){
    PyObject *m;
    m = PyModule_Create(&fitmodule);
    if (!m) {
        return NULL;
    }
    return m;
}

static struct vector* convert_PyPointsList_To_VectorPointer(PyObject *pylist, int N, int dim){
    PyObject *item;
    PyObject *curr_lst;
    double num;
    int i,j;
    struct vector *head_vec, *curr_vec, *next_curr_vec;
    struct cord *curr_cord, *next_curr_cord;

    /* providing space for the head_vec - first vec*/
    curr_vec = malloc(sizeof(struct vector));
    if (curr_vec == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    curr_vec->next = NULL;

    /*the first vec is the head of the linked list*/
    head_vec = curr_vec;

    for (i = 0; i < N; i++) {

        /*gets the i'th arg in pylist - a list of doubles*/
        curr_lst = PyList_GetItem(pylist, i);  
        if (!PyList_Check(curr_lst)) {
            printf("Invalid list object encountered. Exiting.\n");
            freeVectors(head_vec);
            return NULL;
        }
        
        /* providing space for the first cord*/
        curr_cord = malloc(sizeof(struct cord));
        if (curr_cord == NULL) {
            printf("An Error Has Occurred\n");
            freeVectors(head_vec);
            exit(1);
        }
        curr_cord->next = NULL;


        /*iterativly add cord to the vec */
        for (j = 0; j < dim; j++){
            
            if (j == 0){    /*the first cord is the cords of the vec*/
            curr_vec->cords = curr_cord;
            }
            
            /*converting the py-float to c-double*/
            item = PyList_GetItem(curr_lst, j); /*gets the j'th double in the list*/
            if (!PyFloat_Check(item)) {
                printf("Invalid float object encountered. Exiting.\n");
                freeVectors(head_vec);
                return NULL;
            }
            num = PyFloat_AsDouble(item);    /*converts the py-float to c-double */
            
            curr_cord->value = num;         /*add the cord value*/

            if (j != dim-1){
                /* providing space for the next cord */
                next_curr_cord = malloc(sizeof(struct cord));
                if (next_curr_cord == NULL) {
                    printf("An Error Has Occurred\n");
                    freeVectors(head_vec);
                    exit(1);
                }
                next_curr_cord->next = NULL;
                curr_cord->next = next_curr_cord;
                curr_cord = next_curr_cord;
            
            }
        }

        /* providing space for the next vec */
        if (i != N-1){
        next_curr_vec = malloc(sizeof(struct vector));
        if (next_curr_vec == NULL) {
            printf("An Error Has Occurred\n");
            freeVectors(head_vec);
            exit(1);
        }
        next_curr_vec->next = NULL;
        curr_vec->next = next_curr_vec;
        curr_vec = next_curr_vec;
        }
        
    }
    return head_vec;
}

static PyObject *convert_VectorPointer_To_PyObject(struct vector *ret_centroidsList, int k, int dim){
    
    PyObject *pylist, *cords_list, *cord_obj;
    struct vector *curr_vec;
    struct cord *curr_cord;
    int i,j;
    double val;

    pylist = PyList_New(k);  /* Create the Python list to hold the vectors lists */
    
    if (pylist == NULL) {
        printf("Failed to create Python list. Exiting.\n");
        return NULL;
    }

    /* pointer to run over the vectors, points to the first vector */
    curr_vec = ret_centroidsList;

    for (i = 0; i < k; i++) {

        /* Create a new Python list to hold the cordinates of the i'th vector */
        cords_list = PyList_New(dim);

        if (cords_list == NULL) {
            printf("Failed to create Python list for coordinates. Exiting.\n");
            return NULL;
        }   

         /* pointer to run over the cordinates of the curr vector */
        curr_cord = curr_vec->cords;

        /*iterativly add cord to the vec */
        for (j = 0; j < dim; j++){

            /* Convert the cordinate to a Python float object */
            val = curr_cord->value;
            cord_obj = Py_BuildValue("d", val);  
            
            if (cord_obj == NULL) {
                printf("Failed to create Python float object. Exiting.\n");
                return NULL;
            }
            
            /* Add the j'th cordinate to the cordinates list */
            PyList_SetItem(cords_list, j, cord_obj);
            /* steping to the next cord of the vec */  
            curr_cord = curr_cord->next;
        }

        /* Add the cords list to the i'th vec */
        PyList_SetItem(pylist, i, cords_list);

        /* steping to the next vec of the vecs */  
        curr_vec = curr_vec->next;
    } 
    return pylist;
}

