#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void decompose(double **M, double*b, double **L, double **U, int order){
	int i;
    int j;
    int k;
    double m;
    // generate_identity(L,order);
// 	clone_matrix(U,M,order);
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            L[i][j] = 0;
            U[i][j] = 0;
        }
    }
    
    for(i = 0; i < order; i++){
    	for(j = 0; j < order; j++){
    		if(i <= j)
    		{
    		    double sum = 0;
    		    for (int k = 1; k < i; k++) {
    		        sum += U[k][j] * L[i][k];
    		    }
    		    U[i][j] = M[i][j] - sum;
    		}
    		if(i >= j) 
    		{
    		    double sum = 0;
    		    for (int k = 1; k < i; k++) {
    		        sum += U[k][j] * L[i][k];
    		    }
    		    L[i][j] = (1/U[j][j]) * (M[i][j] - sum);
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
	
	free_matrix_double(L,order);
	free(y);
	free_matrix_double(U, order);
    
    return x;
}

double condition_number_est(double **M, int order){
    double *w = allocate_vet_double(order);
    fill_double_vector_1(w, order);
    double *c = lu_decomposition(M,w,order);
    
    double c_norm = normv(1,c,order);
    double M_norm = normm(1,M,order);
    double w_norm = normv(1,w,order);
    free(c);
    free(w);
    double condition_number_est = (M_norm*c_norm)/w_norm;
    return condition_number_est;
}