#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Util.c"
#include "ReadMatrixFromFile.c"
#include "LU.c"
#include "jacobi.c"
#include "gauss_seidel.c"
#include "SOR.c"
#include "Cholesky.c"
#include "QR.c"

#define MAX_ERROR 0.0001

int main(){
    
    double **M;
    double *b;
    int size_order;

    read_matrix_from_file("input_files/matrix.in", &M, &b, &size_order);
    
    printf("Propriedades da matriz: (a)\n");
	printf("Positiva Definida: ");
	if(is_pos_definite(M,size_order)){
    	printf("True\n");
	}else{
		printf("False\n");
	}
	
	printf("Ortogonal: ");
	if(isOrthogonal(M,size_order)){
    	printf("True\n");
	}else{
		printf("False\n");
	}
	
//	printf("Matriz Banda: ");
//	if(bandMatrixBandwidth(M,size_order)){
//    	printf("True \n");
//	}else{
//		printf("False\n");
//	}

	printf("Tridiagonal: ");
	if(isTridiagonal(M,size_order)){
    	printf("True\n");
	}else{
		printf("False\n");
	}
	
	printf("Matriz de Hilbert: ");
	if(isHilbert(M,size_order)){
    	printf("True\n");
	}else{
		printf("False\n");
	}
	
	printf("Singular: ");
	if(isSingular(M,size_order)){
    	printf("True\n");
	}else{
		printf("False\n");
	}
	
	//LU
	printf("\nFatoracao LU: (b)\n");
	double **L = allocate_matrix_double(size_order,size_order);
	double **U = allocate_matrix_double(size_order,size_order);
	double *xLU;
	decompose(M,b,L,U,size_order);
	printf("L:\n");
	print_matrix(L,size_order,size_order);
	printf("\nU:\n");
	print_matrix(U,size_order,size_order);
	printf("\nResultado: (h2)\n");
	xLU = lu_decomposition(M,b,size_order);
	print_vector(xLU,size_order);
	printf("\nNumero condicao estimado(LU): (c)\n");
	printf("%lf \n", condition_number_est(M,size_order));
	//liberar LU da memoria
	free_matrix_double(L,size_order);
	free_matrix_double(U,size_order);
	free(xLU);
	
	//Numero Condicao
	printf("\nNumero Condicao: (d)\n");
	printf("Utilizando Norma 1: %lf\n", condition_number(1,M,size_order));
	printf("Utilizando Norma Infinito: %lf\n", condition_number(2,M,size_order));
	printf("Utilizando Norma Frobenius: %lf\n", condition_number(0,M,size_order));
	
	printf("\nCriterios: (e)\n");
	printf("Criterio das linhas: ");
	if(criterio_linhas(M,size_order) ){
    	printf("True\n");
	}else{
		printf("False\n");
	}
	
	printf("Criterio das colunas: ");
	if(criterio_colunas(M,size_order) ){
    	printf("True\n");
	}else{
		printf("False\n");
	}
	
	printf("Criterio de Sassenfeld: ");
	if(sassenfeld(M,size_order) ){
    	printf("True\n");
	}else{
		printf("False\n");
	}
	
	
	printf("\nMetodos Iterativos: (f)\n");
	double *xGauss;
	xGauss = gauss_seidel_method(M,b,size_order,MAX_ERROR);
	printf("Por Gauss-Seidel:\n");
	print_vector(xGauss,size_order);
	printf("Distancia entre sol. aproximada e sol. do sistema: %lf\n", solution_distance(M,b,xGauss,size_order));
	free(xGauss);
	
	double *xJacobi;
	xJacobi = jacobi_method(M,b,size_order,MAX_ERROR);
	printf("\nPor Jacobi:\n");
	print_vector(xJacobi,size_order);
	printf("Distancia entre sol. aproximada e sol. do sistema: %lf\n", solution_distance(M,b,xJacobi,size_order));
	free(xJacobi);
	
	double *xSOR;
	xSOR = SOR_method(M,b,size_order,1.5,MAX_ERROR);
	printf("\nPor SOR:\n");
	print_vector(xSOR,size_order);
	printf("Distancia entre sol. aproximada e sol. do sistema: %lf\n", solution_distance(M,b,xSOR,size_order));
	free(xSOR);
	
	printf("\nCholesky: (h)\n");
	if(is_pos_definite(M,size_order)){
    	double **RT = chol_fatoration(M,size_order);
    	double **R = transpose(RT,size_order);
    	
    	//Rty = b
    	double *choleskyY = forward_substitution(RT,b,size_order);
    	//Rx = y
    	double *choleskyX = back_substitution(R,choleskyY,size_order);
    	
    	printf("Fator de Cholesky:\n");
    	print_matrix(R,size_order,size_order);
    	
		printf("Resultado:\n");
		print_vector(choleskyX,size_order);
		
		free_matrix_double(RT,size_order);
		free_matrix_double(R,size_order);
		free(choleskyX);
		free(choleskyY);
    	
	}else{
		printf("Nao eh possivel utilizar a Fatoracao de Cholesky, pois a matriz nao eh Positiva definida.\n");
	}
	
	printf("\nMaior angulo entre vetores: (i)\n");
	max_angle_between_vec(M,size_order);
	
	printf("\n\nFatoracao QR: (j)\n");
	double *xQR = QR_decomposition(M,b,size_order);
	if(xQR != NULL){
		printf("Resultado: \n");
		print_vector(xQR,size_order);
		free(xQR);
	}
	
//   print_matrix(M,size_order,size_order);
//   print_vector(b, size_order);
//   print_vector(x, size_order);
	free_matrix_double(M,size_order);
	free(b);
    return 0;
}
