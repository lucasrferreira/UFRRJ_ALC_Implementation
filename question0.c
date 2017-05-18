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
    double *Mx;
    double erro;
    double erro_maximo = 0.1;
    // int size_order2;

    read_matrix_from_file("input_files/matrix.in", &M, &b, &size_order);
    read_vector_from_file("input_files/x_vector.in", &x, &size_order);

    
    Mx = prod_matrix_vector(M,x,size_order,size_order); 

	// erro = || b-Mx || norma 2
	erro = normv(2, sub_v(b,Mx,size_order) ,size_order);

	if(erro <= erro_maximo){
		printf("Correto");
	}else{
		printf("Incorreto");
	}
	
	
	free_matrix_double(M,size_order);
	free(b);
	free(x);
	free(Mx);
	
    return 0;
}
