#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Util.c"
#include "ReadMatrixFromFile.c"

double *normalizar(double *v, int order){
	double norma;
	norma = normv(2,v,order);
	return prod_vector_scalar(v,(1/norma),order);
}

double *orthogonal_proj(double **V, double **U, int row, int order){
	int i;
	double scalar;
	double *aux;
	double *result = allocate_vet_double(order);
	
	zerar_vetor_double(result, order);
	//proj = (V[row] . U[i-1]) . U[i-1] + (V[row] . U[i-2]) . U[i-2] + ...
	for(i = row; i > 0; i--){	
		scalar = scalar_prod(V[row],U[i-1],order);
		aux = prod_vector_scalar(U[i-1], scalar,order);
		result = sum_v(result,aux,order);
	}
	free(aux);
	return result;
}

double **orthonormal_vectors(double **V, int order){
	int i;
	double *y = allocate_vet_double(order);
	double *projAux;
	double **U = allocate_matrix_double(order,order);
	
	clone_vector(U[0], normalizar(V[0],order), order);

	for(i = 1; i < order; i++){
		projAux = orthogonal_proj(V,U,i,order);
		//y = V[i] - projecao
		y = sub_v(V[i], projAux, order);
		//U[i] = y normalizado. Ou seja: ( y/norma(y) )
		clone_vector(U[i],  normalizar(y, order)   , order);
	}

	free(y);
	free(projAux);
	
	return U;
}

int main(){
    double **M;
    double *b;
    int size_order;
    double **X;
    double **U;

    read_matrix_from_file("input_files/matrix.in", &M, &b, &size_order);

	//transformar todos os vetores em uma matriz, onde cada vetor eh uma coluna
	//se o det for != 0, sao linearmente independentes
	X = transpose(M,size_order);
	if( is_linearmente_independentes(X,size_order)){
		printf("Os vetores sao LI.\n\n");
		
		U = orthonormal_vectors(M,size_order);
		printf("Vetores ortonormais: (cada linha corresponde a um vetor)\n\n");
		print_matrix(U,size_order,size_order);
	}else{
		printf("Os vetores sao LD\n");
	}
	
	free_matrix_double(M,size_order);
	free_matrix_double(X,size_order);
	free_matrix_double(U,size_order);
	free(b);
	
    return 0;
}
