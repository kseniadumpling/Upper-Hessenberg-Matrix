#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 * Func creates new matrix with random elements
 */
void create_random_matrix(int size, double *matrix){
    srand(time(NULL));
    for (int i = 0; i < size*size; i++){
        matrix[i] = rand()%137;
    }
}

/*
 * Func prints the matrix
 */
void print_matrix(int size, double *matrix){
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            printf("%12.6lf ", matrix[size*i+j]);
        }
        printf("\n");
    }
}

/*
 * Func that calculates vector V. It's preparation
 * for future building of Householder matrix (which will
 * be based on vector V)
 * 
 * Input:
 * {int} size - size of main matrix
 * {double*} matrix - the main matrix
 * {int} sizeH - size of future Householder matrix
 * 
 * Output:
 * {double*} vectorV 
 */
void get_vector_v(const int size, double *matrix, const int sizeH, double *vectorV){
    int order = size - sizeH;
    int offset = order-1;
    double sigma = 0; 

    double vectorX[sizeH];
    for (int i = 0; i < sizeH; i++) {
        vectorX[i] = matrix[(i+order)*size + offset];
        sigma += pow(vectorX[i], 2);
    }
    sigma = sqrt(sigma);

    double vectorE[sizeH];
    vectorE[0] = sigma;
    for (int i = 1; i < sizeH; i++){
        vectorE[i] = 0;
    }

    for (int i = 0; i < sizeH; i++){
        vectorV[i] = vectorX[i] + vectorE[i];
    }
}

/*
 * Func builds Householder matrix which is 
 * based on vector V (see get_vector_v())
 * 
 * Input:
 * {int} sizeH - size of Householder matrix
 * {double*} vectorV - vector v
 * 
 * Output:
 * {double*} matrixH - Householder Matrix 
 */
void build_housheholder_matrix(const int sizeH, double *vectorV, double *matrixH){
    double matrixE[sizeH*sizeH];
    for (int i = 0; i < sizeH; i++){
        for (int j = 0; j < sizeH; j++){
            if (i != j){
                matrixE[i*sizeH + j] = 0;
            }
            else {
                matrixE[i*sizeH + j] = 1;
            }
        }
    }

    double normV = 0;
    for (int i = 0; i < sizeH; i++){
        normV += pow(vectorV[i], 2);
    }
    
    double matrixTemp[sizeH*sizeH];
    for (int i = 0; i < sizeH; i++){
        for (int j = 0; j < sizeH; j++){
            matrixTemp[i*sizeH +j] = vectorV[i] * vectorV[j];
            matrixTemp[i*sizeH +j] = matrixTemp[i*sizeH +j]*2/normV;
        }
    }

    for (int i = 0; i < sizeH*sizeH; i++){
        matrixH[i] = matrixE[i] - matrixTemp[i]; 
    }
}


/*
 * Func builds transformation matrix U.
 * It will be multiplied with the main matrix A
 * in order to build Hessenberg matrix (iteratively)
 * 
 * Input:
 * {int} sizeH - size of Householder matrix
 * {double*} matrixH - Householder Matrix
 * {int} size - size of matrix U (and also of matrix A)
 * 
 * Output:
 * {double*} matrixU - transformation matrix
 */
void build_u_matrix(const int sizeH, double *matrixH, const int size, double *matrixU){
    int offset = size - sizeH;

    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            if (i < offset){
                if (i == j){
                    matrixU[i*size + j] = 1;
                }
                else {
                    matrixU[i*size + j] = 0;
                }
            }
            else {
                if (j < offset) {
                    matrixU[i*size + j] = 0;
                }
                else {
                    matrixU[i*size + j] = matrixH[(i-offset)*sizeH + (j-offset)];
                }
            }
        }
    }
}


/*
 * Func multiplies transformation matrix U and
 * the main matrix A by order U*A*U.
 * As a result, transformed matrix A will become
 * more and more similar with "almost" triangular 
 * matrix 
 * 
 * Input:
 * {int} size - size of matrix A
 * {double*} matrixA - the main matrix
 * {double*} matrixU - transformation matrix
 * 
 * Output:
 * {double*} matrixA - transformated matrix A
 */
void multiply_matrices_UAU(const int size, double *matrixA, double *matrixU){
    double matrixTemp[size*size];
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            matrixTemp[i*size+j] = 0;
            for (int k = 0; k < size; k++){
                matrixTemp[i*size+j] += matrixU[i*size+k]*matrixA[k*size+j];
            }
        }
    }

    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            matrixA[i*size+j] = 0;
            for (int k = 0; k < size; k++){
                matrixA[i*size+j] += matrixTemp[i*size+k]*matrixU[k*size+j];
            }
        }
    }
}

int main(){

    /*
     * This part olny for inputs from file
     */ 
    /*
    FILE *fp;
    char *fname = "input.txt";
    fp = fopen(fname, "r");
    if (fp == NULL){
        printf("\nError. File '%s' wasn't opened\n", fname);
        return -1;
    }
    int size;
    fscanf(fp, "%d", &size);
    
    double matrixA[size*size];
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            fscanf(fp, "%lf", &matrixA[size*i+j]); 
        }
    }
    printf("\nInput matrix A:\n");
    print_matrix(size, matrixA);
    */

    int size = 5;
    double matrixRand[size*size];
    create_random_matrix(size, matrixRand);
    printf("\nInput: matrix with random elements\n");
    print_matrix(size, matrixRand);

    
    for (int i = 1; i < size-1; i++){
        double vectorV[size-i];
        get_vector_v(size, matrixRand, size-i, vectorV);
        double matrixH[(size-i)*(size-i)];
        build_housheholder_matrix(size-i, vectorV, matrixH);
        double matrixU[size*size];
        build_u_matrix(size-i, matrixH, size, matrixU);
        multiply_matrices_UAU(size, matrixRand, matrixU);
    }
    printf("\n\nOutput: Hessenberg matrix\n");
    print_matrix(size, matrixRand);
}
