#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void QR_decompose_2(double **M, double*b, double **Q, double **R, int order){
	if(order != 2){
	    printf("A matriz deve ser de ordem 2.");
		exit(1);
	}
	
	double a11 = M[0][0];
	double a21 = M[1][0];
	double sqrt_a = sqrt( a11*a11 + a21*a21  );
	double q_sin = a21 / sqrt_a;
	double q_cos = a11 / sqrt_a;
	Q[0][0] = q_cos;
	Q[0][1] = q_sin;
	Q[1][0] = q_sin * -1;
	Q[1][1] = q_cos;

    double **auxR = prod_matrix(Q, M, order, order, order);
    clone_matrix(R,auxR,order);
    free_matrix_double(auxR, order);
}

double *QR_decomposition(double **M, double*b, int order){
	
    double *x;
    double *y;
	double **Q = allocate_matrix_double(order,order);
    double **R = allocate_matrix_double(order,order);
    
    if(order != 2){
	    printf("A matriz deve ser de ordem 2.");
		exit(1);
	}
	
    decompose(M,b,Q,R,order);
    
    //QRx=b => y=QTb => Rx = y
	// Qy = b
	double **Qt = transpose(Q,order);
	y = prod_matrix_vector(Qt, b, order, order);
	
	//Utiliza o y encontrado pra encontrar x ( Rx = y )
	x = back_substitution(R,y,order);
	
	free_matrix_double(Qt,order);
	free_matrix_double(Q,order);
	free(y);
	free_matrix_double(R, order);
    
    return x;
}

