#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double *jacobi_method(double **M, double*b, int order, double erro_maximo){

    double *x;
    double *x_anterior;
    int i;
    int j;
    bool convergiu = false;
    double soma;
    
    clock_t tempo_inicial;
	clock_t tempo_atual;
	float tempo_corrido = 0;
    
    x = allocate_vet_double(order);
    x_anterior = allocate_vet_double(order);
    
	if(main_diagonal_has_zeros(M,order)){
		printf("Os elementos da diagonal principal devem ser diferentes de zero");
		exit(1);
	}
	
	//valores iniciais para x = (0, ... ,0)
    zerar_vetor_double(x_anterior, order);
    
   	tempo_inicial = clock();
	
	while(convergiu == false && tempo_corrido < 15){
		
		for(i = 0; i < order; i++){
			soma = 0.0;
			
			for(j = 0; j < order; j++){
				if(j != i){
					soma += M[i][j] * x_anterior[j];
				}
			}
			
			x[i] = (b[i] - soma)/M[i][i];
		}
		
		tempo_atual = clock();
		
		//Calcula o tempo de execu��o do programa em segundos
		tempo_corrido = (tempo_atual - tempo_inicial)/CLOCKS_PER_SEC;
		convergiu = is_correct(M, b, x, erro_maximo, order);
		clone_vector(x_anterior, x, order);
	}
	free(x_anterior);
    return x;
}
