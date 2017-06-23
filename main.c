#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Util.c"
#include "Cholesky.c"
#include "LU.c"
#include "SOR.c"
#include "QR.c"
#include "ReadMatrixFromFile.c"


int main(){
    
    double **M;
    double *b;
    int size_order;
    int size_order2;

    read_matrix_from_file("input_files/matrix.in", &M, &b, &size_order);
    read_vector_from_file("input_files/x_vector.in", &b, &size_order2);
    
    // double **Mt = transpose(M, size_order);
    
    double **invert = invert_matrix_n(M, size_order);
    print_matrix(invert, size_order, size_order);
    double c_num0 = condition_number(0,M,size_order);
    double c_num1 = condition_number(1,M,size_order);
    double c_num2 = condition_number(2,M,size_order);
    double c_num3 = condition_number_est(M ,size_order);
    printf("\n0: %lf\n 1: %lf\n 2: %lf\n est: %lf\n ",c_num0, c_num1, c_num2, c_num3 );
    
    bool x = criterio_norma(M, size_order);
    printf("\n %s \n\n", x ? "true" : "false");

    double *c = SOR_method(M,b,size_order, 1.058823529, 0.000001);
    print_vector(c, size_order);
    
    double **Q = allocate_matrix_double(2,2);
    double **R = allocate_matrix_double(2,2);
    QR_decompose_2(M,b,Q,R,size_order);
    print_matrix(Q,size_order,size_order);
    print_matrix(R,size_order,size_order);
    double *qr = QR_decomposition(M,b,size_order);
    print_vector(c, size_order);
    // printf("\n\n %lf \n", normm(0, M, size_order));
    // printf("\n\n %lf \n", normm(1, M, size_order));
    // printf("\n\n %lf \n", normm(2, M, size_order));
    // printf("\n\n 1 %lf \n", normv(1, b, size_order));
    // printf("\n\n 2 %lf \n", normv(2, b, size_order));
    // printf("\n\n 3 %lf \n", normv(3, b, size_order));
    // printf("\n\n 4 %lf \n", normv(4, b, size_order));
    return 0;
}