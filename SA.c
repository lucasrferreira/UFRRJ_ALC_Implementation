#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


void shaking(double *x, int i){
	x[i] += 0.1;
}

double custo(double **M, double *x, double *b, int order){
	double *Mx;
	double erro;
	Mx = prod_matrix_vector(M,x,order,order); 
		
	// erro = || b-Mx || norma 2
	erro = normv(2, sub_v(b,Mx,order) ,order);
	free(Mx);
	return erro;
}

bool aceitar(double *x1, double *x2, double t, double **M, double *b, int order){
	double c1,c2;
	double prob;
	double random;
	c1 = custo(M,x1,b,order);
	c2 = custo(M,x2,b,order);
	
	if(c2 <= c1){
		prob = 1;
	}else{
		prob = exp( -(c2-c1)/t);
	}
	
	random = (rand()%100) / 100.0;
	
	return random <= prob;
}

void update_best(double **M, double *x1, double *x_best, double *b,int order){
	if(custo(M,x1,b,order) <  custo(M,x_best,b,order)){
		clone_vector(x_best,x1,order);
		//printf("aceitou best:\n");
		print_vector(x_best,order);
	}
}

double *simulated_annealing(double **M, double*b, int order, double erro_maximo){

    double *x1;
    double *x2;
    double *x_best;
    int i;
    int j;
    bool convergiu = false;
    double t = 100;
    double taxa_t = 0.5;
    
    
    clock_t tempo_inicial;
	clock_t tempo_atual;
	float tempo_corrido = 0;
    
    x1 = allocate_vet_double(order);
    x2 = allocate_vet_double(order);
    x_best = allocate_vet_double(order);
    
//	if(main_diagonal_has_zeros(M,order)){
//		printf("Os elementos da diagonal principal devem ser diferentes de zero");
//		exit(1);
//	}
	
	//valores iniciais para x = (0, ... ,0)
    zerar_vetor_double(x1, order);
    zerar_vetor_double(x_best, order);
    
   	tempo_inicial = clock();
	
	while(convergiu == false && tempo_corrido < 20){
		
		for(i = 0; i < order; i++){
			clone_vector(x2,x1,order);
			shaking(x2,i);
				
			if(aceitar(x1,x2,t,M,b,order) ){
				clone_vector(x1,x2,order);
				//se o custo de x1 < x_best, x_best = x1
				update_best(M,x1,x_best,b, order);
			}
			
		}
		//so pra nao ficar negativo
		t -= fabs(taxa_t);
		
		tempo_atual = clock();
		
		//Calcula o tempo de execucao do programa em segundos
		tempo_corrido = (tempo_atual - tempo_inicial)/CLOCKS_PER_SEC;
		convergiu = is_correct(M, b, x_best, erro_maximo, order);
	}
	
	free(x1);
	free(x2);
    return x_best;
}
