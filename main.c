#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Util.c"
#include "ReadMatrixFromFile.c"

int main(){
    
    double **M;
    double *b;
    int size_order;
    int size_order2;

    read_matrix_from_file("input_files/matrix.in", &M, &b, &size_order);
    read_vector_from_file("input_files/x_vector.in", &b, &size_order2);
    
    print_matrix(M, size_order, size_order);
    print_vector(b, size_order);
    double *p = prod_matrix_vector(M, b, size_order, size_order);
    double *q = prod_vector_matrix(b, M, size_order, size_order);
    
    
    print_vector(p, size_order);
    
    print_vector(q, size_order);
    
    double **C = prod_matrix(M, M, size_order, size_order, size_order);
    printf("\n\n M before \n");
    print_matrix(M, size_order, size_order);
    free_matrix_double(M,size_order);
    printf("\n\n M after \n");
    print_matrix(M, size_order, size_order);
    M = C;
    printf("\n\n C  \n");
    print_matrix(C, size_order, size_order);
    printf("\n\n M  \n");
    print_matrix(M, size_order, size_order);
    
    // printf("\n\n %lf \n", normm(0, M, size_order));
    // printf("\n\n %lf \n", normm(1, M, size_order));
    // printf("\n\n %lf \n", normm(2, M, size_order));
    // printf("\n\n 1 %lf \n", normv(1, b, size_order));
    // printf("\n\n 2 %lf \n", normv(2, b, size_order));
    // printf("\n\n 3 %lf \n", normv(3, b, size_order));
    // printf("\n\n 4 %lf \n", normv(4, b, size_order));
    return 0;
}