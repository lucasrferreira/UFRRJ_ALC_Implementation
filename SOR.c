#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<time.h>

double *SOR_method(double **M, double*b, int order, double w, double erro_maximo){

    double *x;
    int i;
    int j;
    bool convergiu = false;
    double soma;
    
	clock_t tempo_inicial;
	clock_t tempo_atual;
	float tempo_corrido = 0;

    x = allocate_vet_double(order);
    
	if(main_diagonal_has_zeros(M,order)){
		printf("Os elementos da diagonal principal devem ser diferentes de zero");
		exit(1);
	}

	zerar_vetor_double(x, order);
	tempo_inicial = clock();
	
	//tempo_corrido eh pra interromper a execucao em 15s, caso nao convirja
	while(convergiu == false && tempo_corrido < 15){
		
		for(i = 0; i < order; i++){
			soma = 0.0;
			
			for(j = 0; j < order; j++){
				if(j != i){
					soma += M[i][j] * x[j];
				}
			}
			
			x[i] = (1-w) * x[i] + w * ((b[i] - soma)/M[i][i]);
		}
		tempo_atual = clock();
		//Calcula o tempo de execucao do programa em segundos
		tempo_corrido = (tempo_atual - tempo_inicial)/CLOCKS_PER_SEC;

		convergiu = is_correct(M, b, x, erro_maximo, order);
	}
	
    return x;
}
