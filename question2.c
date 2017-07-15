#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Util.c"
#include "ReadMatrixFromFile.c"

void ler_arquivo(char *path, double **P, double *xAnterior,int *t){
	int i;
	int j;
	double valor;
	FILE *myFile;
	
    myFile = fopen(path, "r");
    generate_identity(P,4);
	fscanf(myFile, "%d", t);
    
	for(i = 0; i< 4; i++){
    	fscanf(myFile, "%lf", &valor);
    	xAnterior[i] = valor;
	}
    
	for(i=0; i<4; i++){
		for(j=0; j<4; j++){
			fscanf(myFile, "%lf", &P[i][j]);
		}
	}

    fclose(myFile);
}

int main(){
    
    double *xAnterior;
    double *xAtual;
    char estado[5] = "SNYA";    
    int t;
    int i;
	
	xAnterior = allocate_vet_double(4);
	double **P = allocate_matrix_double(4,4);
	ler_arquivo("input_files/probabilities.in",P,xAnterior,&t);
	
	printf("Numero de individuos em t = 0:\n");
	for(i = 0; i < 4; i++){
		printf("%c[0]= %lf\n",estado[i], xAnterior[i]);
	}
	
	for(i = 0; i < t; i++){
		xAtual = prod_matrix_vector(P,xAnterior,4,4);
		clone_vector(xAnterior,xAtual,4);
	}

	printf("\nNumero de individuos em t = %d:\n",t);
	for(i = 0; i < 4; i++){
		printf("%c[%d]= %lf\n",estado[i], t, xAtual[i]);
	}
	
	free(xAnterior);
	free(xAtual);
	free_matrix_double(P,4);
    return 0;
}
