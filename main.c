#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Util.c"
#include "Cholesky.c"
#include "ReadMatrixFromFile.c"


int main(){
    
    double **M;
    double *b;
    int size_order;
    int size_order2;

    read_matrix_from_file("input_files/matrix.in", &M, &b, &size_order);
    read_vector_from_file("input_files/x_vector.in", &b, &size_order2);
    
    // double **Mt = transpose(M, size_order);
    
    bool x = is_pos_definite_chol(M, size_order);
    printf("%s \n\n", x ? "true" : "false");

    double **L = chol_fatoration(M, size_order);
    print_matrix(L, size_order, size_order);
    
    // printf("\n\n %lf \n", normm(0, M, size_order));
    // printf("\n\n %lf \n", normm(1, M, size_order));
    // printf("\n\n %lf \n", normm(2, M, size_order));
    // printf("\n\n 1 %lf \n", normv(1, b, size_order));
    // printf("\n\n 2 %lf \n", normv(2, b, size_order));
    // printf("\n\n 3 %lf \n", normv(3, b, size_order));
    // printf("\n\n 4 %lf \n", normv(4, b, size_order));
    return 0;
}