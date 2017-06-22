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
    double erro_maximo = 0.1;

    read_matrix_from_file("input_files/matrix.in", &M, &b, &size_order);
    read_vector_from_file("input_files/x_vector.in", &x, &size_order);

	
	if( is_correct(M, b, x, erro_maximo, size_order) ){
		printf("Correto");
	}else{
		printf("Incorreto");
	}
	
	free_matrix_double(M,size_order);
	free(b);
	free(x);
	
    return 0;
}
