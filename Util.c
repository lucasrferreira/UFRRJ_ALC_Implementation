#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print_matrix( double **M, int size_column, int size_line){
    
    for (int i = 0; i < size_line; i++) {
        for (int j = 0; j < size_column; j++) {
            printf("%.20f ", M[i][j]);
        }
        printf("\n");
    }
}

void print_vector( double *v , int size) {
    for (int j = 0; j < size; j++) {
            printf("%.20f \n", v[j]);
    }
}

double **allocate_matrix_double(int Rows, int Cols)
{    
    // allocate Rows rows, each row is a pointer to int
    double **board = (double **)malloc(Rows * sizeof(double *)); 
    int row;

    // for each row allocate Cols ints
    for (row = 0; row < Rows; row++) {
        board[row] = (double *)malloc(Cols * sizeof(double));
    }

    return board;
}
double *allocate_vet_double(int Rows)
{    
    double *board = (double *)malloc(Rows * sizeof(double));
    
    return board;
}

double **prod_matrix( double **A, double **B, int A_column_size, int B_line_size, int n_size) {
    
    double **M = allocate_matrix_double(A_column_size, B_line_size);
    
    double value = 0;
    for (int i = 0; i < A_column_size; i++) {
        for (int j = 0; j < B_line_size; j++) {
            value = 0;
            for (int k = 0; k < n_size; k++) {
                value += A[i][k] * B[k][j];   
            }
            M[i][j] = value;
        }
    }
    return M;
}

void free_matrix_double(double **board, int Rows) 
{
    int row;

    // first free each row
    for (row = 0; row < Rows; row++) {
         free(board[row]);
    }

    // Eventually free the memory of the pointers to the rows
    free(board);
 }

double *prod_matrix_vector(double **A, double *b, int A_column_size,  int line_size){
    double *M = allocate_vet_double(line_size);
    
    double value = 0;
    for (int i = 0; i < A_column_size; i++) {
        value = 0;
        for (int j = 0; j < line_size; j++) {
            value += A[i][j] * b[j];   
        }
        M[i] = value;
    }
    return M;
}
double *prod_vector_matrix(double *b, double **A, int column_size,  int A_line_size){
    double *M = allocate_vet_double(column_size);
    
    double value = 0;
    for (int i = 0; i < A_line_size; i++) {
        value = 0;
        for (int j = 0; j < column_size; j++) {
            value += A[j][i] * b[j];   
        }
        M[i] = value;
    }
    return M;
}

double normm(int norm_id, double **M, int matrix_order)
{
    if(norm_id == 0) //frobenius
    {
        double sum_res = 0;
        for (int i = 0; i < matrix_order; i++) {
            for (int j = 0; j < matrix_order; j++) {
                sum_res += M[i][j]*M[i][j];
            }
        }
        return sqrt(sum_res);
    }
    else if(norm_id == 1) //line
    {
        double atual_max_norm = 0;
        for (int i = 0; i < matrix_order; i++) {
            double atual_norm = 0;
            for (int j = 0; j < matrix_order; j++) {
                atual_norm += abs(M[i][j]);
            }
            
            if(atual_max_norm < atual_norm)
            {
                atual_max_norm = atual_norm;
            }
        }
        return atual_max_norm;
    }
    else if(norm_id == 2) //column
    {
        double atual_max_norm = 0;
        for (int i = 0; i < matrix_order; i++) {
            double atual_norm = 0;
            for (int j = 0; j < matrix_order; j++) {
                atual_norm += abs(M[j][i]);
            }
            
            if(atual_max_norm < atual_norm)
            {
                atual_max_norm = atual_norm;
            }
        }
        return atual_max_norm;
    }
}

double normv(int norm_order, double *v, int vet_order)
{
    double sum_res = 0;
    for (int i = 0; i < vet_order; i++) {
        sum_res += pow(abs(v[i]), norm_order);
    }
    return pow(sum_res, 1.0/norm_order);

}
double normv_max(double *v, int vet_order)
{
    double sum_res = 0;
    for (int i = 0; i < vet_order; i++) {
        double temp = abs(v[i]);
        if(temp > sum_res)
            sum_res = temp;
        
    }
    return sum_res;

}

double *sum_v(double *v1, double *v2, int size){
    double *result = allocate_vet_double(size);
    
    for (int i = 0; i < size; i++) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}
double **sum_m(double **v1, double **v2, int size){
    double **result = allocate_matrix_double(size, size);
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            result[i][j] = v1[i][j] + v2[i][j];
        }
    }
    return result;
}


double *sub_v(double *v1, double *v2, int size){
    double *result = allocate_vet_double(size);
    
    for (int i = 0; i < size; i++) {
        result[i] = v1[i] - v2[i];
    }
    return result;
}

void zerar_vetor_double(double *v, int n){
	int i;
	for(i=0; i < n; i++){
		v[i] = 0.0;
	}
}

double det_ordem_inferior_a_4(double **m, int n){
	int i;
	int j;
	double resultado = 0.0;
	//a quantidade de diagonais nao eh proporcional
	//entao nao deu pra arrumar uma formula que englobasse n=2 e n=3
	if(n == 2){
		double dp;
		double ds;
		dp = m[0][0]*m[1][1];
		ds = m[0][1]*m[1][0];
		return dp - ds;
	}else{
		double dp[n];
		double ds[n];
		
		//inicializando os vetores
		for(i = 0; i < n; i++){
			dp[i] = 1.0;
			ds[i] = 1.0;
		}
		
		for(i = 0; i < n; i++){
			
			for(j = 0; j < n; j++){
				//(i+j)%n eh pra dar a volta na matriz (sarrus)
				dp[j] *= m[i][(i+j)%n];
				ds[j] *= m[(n-1)-i][(i+j+1)%n];
			}
		}

		for(i = 0; i < n; i++){
			resultado += dp[i] - ds[i];
		}
		return resultado;
	}
}



int isSingular(double **M){
	//calcular o determinante
	//se for 0, eh singular
}
