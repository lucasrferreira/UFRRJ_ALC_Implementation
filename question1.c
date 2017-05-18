#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Util.c"
#include "ReadMatrixFromFile.c"

int main(){
    
    double **M;
    double *b;
    double *x;
    int size_order;
    // int size_order2;

    read_matrix_from_file("input_files/matrix.in", &M, &b, &size_order);
    read_vector_from_file("input_files/x_vector.in", &x, &size_order);
    
   print_matrix(M,size_order,size_order);
   print_vector(b, size_order);
   print_vector(x, size_order);
   
    return 0;
}