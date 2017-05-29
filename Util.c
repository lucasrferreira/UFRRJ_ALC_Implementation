#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>


bool approximately(double v1, double v2, double precision){
    return fabs((double) v1 -  (double) v2) < (double) precision;
}
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
                atual_norm += fabs(M[i][j]);
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
                atual_norm += fabs(M[j][i]);
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
        sum_res += pow(fabs(v[i]), norm_order);
    }
    
    return pow(sum_res, 1.0/norm_order);
}

double normv_max(double *v, int vet_order)
{
    double sum_res = 0;
    for (int i = 0; i < vet_order; i++) {
        double temp = fabs(v[i]);
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

double **transpose(double **M, int order){
    double **Mt = allocate_matrix_double(order, order);
    
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            Mt[j][i]=M[i][j];
        }
    }
    return Mt;
}

bool isIdentity(double **M, int order){
    bool _isIdentity = true;
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            if(i == j){
                if(M[i][j] != 1){
                    _isIdentity = false;
                }
            }
            else if(i != j){
                if(M[i][j] != 0) {
                    _isIdentity = false;
                }
            }
        }
    }
    return _isIdentity;
}

bool isOrthogonal(double **M, int order){
    double **Mt = transpose(M, order);
    double **I = prod_matrix(M, Mt, order, order, order);
    return isIdentity(I, order);
}

int bandMatrixBandwidth(double **M, int order){
    int n = 0;
    int j = 0;
    int i = 0;
    bool should_break = false;
    for (n = 0; n < order; n++) {
        for (i = 0; i <= n; i++) {
            j = (order - 1) - n + i;
            if (M[i][j] != 0 || M[j][i] != 0){
                should_break = true;
            }
            if(should_break)
                break;
        }
        if(should_break)
            break;
    }
    return order - 1 - n;
}

bool isTridiagonal(double **M,int order ){
    int k = bandMatrixBandwidth(M, order);
    if (k == 1){
        return true;
    }
    else{
        return false;
    }
}

bool isHilbert(double **M, int order){
    double hilbertValue = 0;
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            double denom = (double)i + 1 + (double)j + 1 - 1;
            hilbertValue = (double)1/denom;
            if(!approximately(M[i][j], hilbertValue , 0.001))
                return false;
        }
    }
    return true;
}


bool isSingular(double **M, int order){
    double det = det_ordem_inferior_a_4(M, order);
    return det == 0;
}




double *forward_substitution(double **L, double *b, int order){
	double *x = allocate_vet_double(order);
	int i;
	int j;
	double soma;
	
	x[0] = b[0]/L[0][0];
	for(i=1; i < order; i++){
		soma = 0.0;
		for(j=0; j < i; j++){
			soma += L[i][j] * x[j];
		}
		x[i] = (b[i] - soma)/L[i][i];
	}
	
	return x;
}


double *back_substitution(double **U, double *b, int order){
	double *x = allocate_vet_double(order);
	int i;
	int j;
	int n = order;
	double soma;
	
	x[n-1] = b[n-1]/U[n-1][n-1];
	for(i=n-2; i >= 0; i--){
		soma = 0.0;
		for(j=i+1; j < n; j++){
			soma += U[i][j] * x[j];
		}
		x[i] = (b[i] - soma)/U[i][i];
	}
	
	return x;
}

bool is_correct(double **M, double *b, double *x, double erro_maximo, int order){
	double *Mx;
	double erro;
	Mx = prod_matrix_vector(M,x,order,order); 
	
	// erro = || b-Mx || norma 2
	erro = normv(2, sub_v(b,Mx,order) ,order);
	
	if(erro <= erro_maximo){
		return true;
	}else{
		return false;
	}
}

bool main_diagonal_has_zeros(double **A, int order){
	int i;
	for(i = 0; i < order; i++){
		if(A[i][i] == 0.0)
			return true;
	}
	return false;
}

void *clone_vector(double *v1, double *v2, int order){
	int i;
	for(i = 0; i < order; i++){
		v1[i] = v2[i];
	}
}

bool criterio_linhas(double **A, int order){
	int i;
	int j;
	double soma_linha;
	for(i=0; i < order; i++){
		soma_linha = 0;
		for(j=0; j < order; j++){
			if(i != j){
				soma_linha += fabs(A[i][j]);
			}
		}
		if(fabs(A[i][i]) <= soma_linha){
			return false;
		}
	}
	return true;
}

bool criterio_colunas(double **A, int order){
	int i;
	int j;
	double soma_coluna;
	for(j=0; j < order; j++){
		soma_coluna = 0;
		for(i=0; i < order; i++){
			if(i != j){
				soma_coluna += fabs(A[i][j]);
			}
		}
		if(fabs(A[j][j]) <= soma_coluna){
			return false;
		}
	}
	return true;
}
