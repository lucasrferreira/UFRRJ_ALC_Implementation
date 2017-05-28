#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *gauss_seidel_method(double **M, double*b, int order, double erro_maximo){

    double *x;
    double *x_anterior;
    int i;
    int j;
    bool convergiu = false;
    double soma;
    
    x = allocate_vet_double(order);
    
    //TODO: precisa disso aqui?
	/*if(main_diagonal_has_zeros(M,order)){
		printf("Os elementos da diagonal principal devem ser diferentes de zero");
		exit(1);
	}*/
	
	//valores iniciais para x = (0, ... ,0)
	zerar_vetor_double(x, order);
	
	while(convergiu == false){
		
		for(i = 0; i < order; i++){
			soma = 0.0;
			
			for(j = 0; j < order; j++){
				if(j != i){
					soma += M[i][j] * x[j];
				}
			}
			
			x[i] = (b[i] - soma)/M[i][i];
		}
		
		convergiu = is_correct(M, b, x, erro_maximo, order);
	}
	
    return x;
}
