#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void decompose(double **M, double*b, double **L, double **U, int order){
	int i;
    int j;
    int k;
    double m;
    generate_identity(L,order);
	clone_matrix(U,M,order);
    
    for(j = 0; j < order-1; j++){
	
    	for(i = j+1; i < order; i++){
    		m = U[i][j]/U[j][j];
    		L[i][j] = m;

    		for(k = 0; k < order; k++){
    			U[i][k] = U[i][k] - (U[j][k] * m);
			}
		}
	}
}

double *lu_decomposition(double **M, double*b, int order){
	
    double *x;
    double *y;
	double **U = allocate_matrix_double(order,order);
    double **L = allocate_matrix_double(order,order);
    
    if(isSingular(M, order)){
    	printf("A matriz deve ser nao singular.");
		exit(1);
	}
    decompose(M,b,L,U,order);
    
	// Ly = b
	y = forward_substitution(L,b,order);
	
	//Utiliza o y encontrado pra encontrar x ( Ux = y )
	x = back_substitution(U,y,order);
		
    return x;
}
