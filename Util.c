#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define PI 3.14159265358979323846
#define MIN 0.0001

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

void copy_matrix(double **M_to, double **M_from, int order){
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            M_to[i][j] = M_from[i][j];
        }
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

double **prod_matrix_scalar(double **A, double b, int A_column_size,  int line_size){
    double **M = allocate_matrix_double(line_size, A_column_size);
    double value = 0;
    for (int i = 0; i < A_column_size; i++) {
        for (int j = 0; j < line_size; j++) {
            M[i][j] = A[i][j] * b;   
        }
    }
    return M;
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

//norm_id = 0 for frobenius, 1 for norm-1 and 2 for infinite, 
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
    else if(norm_id == 2) //infinite
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
    else if(norm_id == 1) //norm 1
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

double determinant_order_n(double **a,int n)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = allocate_matrix_double(n,n);
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,j1+2.0) * a[0][j1] * determinant_order_n(m,n-1);
         for (i=0;i<n;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}

double **coFactor(double **a,int n)
{
    int i,j,ii,jj,i1,j1;
    double det;
    double **c;
    double **b = allocate_matrix_double(n,n);
    c = allocate_matrix_double(n,n);
    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) {
            /* Form the adjoint a_ij */
            i1 = 0;
            for (ii=0;ii<n;ii++) {
                if (ii == i)
                    continue;
                j1 = 0;
                for (jj=0;jj<n;jj++) {
                    if (jj == j)
                        continue;
                    c[i1][j1] = a[ii][jj];
                    j1++;
                }
                i1++;
            }   
            /* Calculate the determinate */
            det = determinant_order_n(c,n-1);
            /* Fill in the elements of the cofactor */
            b[i][j] = pow(-1.0,i+j+2.0) * det;
        }
    }
    free_matrix_double(c,n);
    return b;
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

//TODO: fazer para ordem > 3
bool isSingular(double **M, int order){
    double det = determinant_order_n(M, order);
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
	double *diff;
	double erro;
	Mx = prod_matrix_vector(M,x,order,order); 
	
	// erro = || b-Mx || norma 2
	diff = sub_v(b,Mx,order);
	erro = normv(2, diff ,order);
	
	free(Mx);
	free(diff);
	
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

//v1 = v2
void clone_vector(double *v1, double *v2, int order){
	int i;
	for(i = 0; i < order; i++){
		v1[i] = v2[i];
	}
}

//M1 = M2
void clone_matrix(double **M1, double **M2, int order){
	int i;
	int j;
	for(i = 0; i < order; i++){
		for(j = 0; j < order; j++){
			M1[i][j] = M2[i][j];
		}
	}
}

bool criterio_norma(double **M, int order){
    double **auxB = allocate_matrix_double(order, order);
    double **B = prod_matrix_scalar(auxB,0,order,order);
    free_matrix_double(auxB, order);
    
    for(int i=0; i < order; i++){
        for(int j=0; j < order; j++){
            if(i!=j)  
                B[i][j] = -M[i][j]/M[i][i];
        }
    }
    
    print_matrix(B, order, order);
    //satisfaz se alguma norma for menor que 1 
    bool satisfy = false;
    for (int i = 0; i < 3; i++) {
        satisfy = (normm(i, B, order) < 1) || satisfy;
    }
    return satisfy;
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

//preenche o vetor com 1
void fill_double_vector_1(double *v, int n){
	int i;
	for(i=0; i < n; i++){
		v[i] = 1.0;
	}
}

double max_value(double *v, int order){
	double max = v[0];
	int i;
	
	//ta redundante, mas pelo menos nao tem indexOutOfBound (caso order == 1)
	for(i = 0; i < order; i++){
		if(v[i] > max){
			max = v[i];
		}
	}
	return max;
}

bool sassenfeld(double **A, int order){
	int i;
	int j;
	double soma_linha;
	double *beta;
	beta = allocate_vet_double(order);
	
	//preenche o vetor dos betas com 1
	fill_double_vector_1(beta, order);
	
	for(i=0; i < order; i++){
		soma_linha = 0;
		for(j=0; j < order; j++){
			if(i != j){
				soma_linha += fabs( A[i][j] * beta[j]);
			}
		}
		beta [i] = soma_linha/fabs(A[i][i]);
	}
	
	return max_value(beta, order) < 1;
}

void generate_identity(double **M, int order){
	int i;
	int j;
	
	for(i = 0; i < order; i++){
		for(j = 0; j < order; j++){
			if(i == j){
				M[i][j] = 1;
			}else{
				M[i][j] = 0;
			}
		}
	}
}


double **invert_matrix_n(double **m, int n){
    double idet = 1/determinant_order_n(m,n);
    double **cofac = coFactor(m,n);
    double **adj = transpose(cofac,n);
    free_matrix_double(cofac,n);
    
    double **inverse = prod_matrix_scalar(adj, idet, n, n);
    free_matrix_double(adj,n);
    return inverse;
}

double condition_number(int norm_id, double **M, int matrix_order)
{
    double **iM = invert_matrix_n(M,matrix_order);
    double c_num = normm(norm_id,M,matrix_order) * normm(norm_id,iM,matrix_order);
    free_matrix_double(iM,matrix_order);
    return c_num;
}

bool is_linearmente_independentes(double **M, int order){
	double det = determinant_order_n(M,order);
	//MIN eh aproximadamente 0
	return fabs(det)  > MIN;
}

double *prod_vector_scalar(double *v, double scalar, int order){
	int i;
	double *result = allocate_vet_double(order);
	for(i=0; i < order; i++){
		result[i] = v[i]*scalar;
	}
	return result;
}

double scalar_prod(double *v1, double *v2, int order){
	int i;
	double result = 0;
	for(i = 0; i < order; i++){
		result += v1[i] * v2[i];
	}
	return result;
}

double angle_between_vec(double *v1, double *v2, int order){
	double cos_theta;
	cos_theta = scalar_prod(v1,v2,order) / ( normv(2,v1,order) * normv(2,v2,order) );
	// *180/PI -> conversao de rad para graus
	return (acos(cos_theta)*180)/PI;
}

void max_angle_between_vec(double **M, int order){
	double cos_phi;
	int i,j;
	int *vi;
	int *vj;
	int cont = 0;
	int tam_vect = order;
	double angle;
	double max_angle = -1.0;
	double **X = transpose(M,order);
	
	//alocacao inicial dos pares de colunas com maior angulo
	vi = (int *)malloc(tam_vect * sizeof(int));
	vj = (int *)malloc(tam_vect * sizeof(int));
	
	for(i=0; i<order-1; i++){
		
		for(j=i+1; j<order; j++){
			angle = angle_between_vec(X[i],X[j],order);
			if(angle > max_angle){
				max_angle = angle;
				//se encontrou um angulo maior, recomecar a lista de pares
				vi[0] = i;
				vj[0] = j;
				cont = 1;
			}else if(angle == max_angle){
				//adicionar par a lista
				if(cont >= tam_vect){
					//alocar mais memoria
					tam_vect++;
					vi = (int *)realloc(vi, tam_vect * sizeof(int));
					vj = (int *)realloc(vj, tam_vect * sizeof(int));
				}
				
				vi[cont] = i;
				vj[cont] = j;
				cont++;
			}
		}
	}
	
	printf("Pares de colunas que formam o maior angulo:");
	for(i = 0; i < cont; i++){
		printf("\nColunas %d e %d\n", vi[i], vj[i]);
		print_vector(X[vi[i]], order);
		printf("\n");
		print_vector(X[vj[i]], order);
	}
	
	printf("\nMaior angulo: %lf", max_angle);
	free(X);
	free(vi);
	free(vj);
}

bool is_pos_definite(double **A, int order){
    
    double **L = allocate_matrix_double(order, order);
    copy_matrix(L,A,order);
    int n = order;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < (i+1); j++) {
            double s = 0;
            for (int k = 0; k < j; k++)
                s += L[i][k] * L[j][k];
                
            if(i == j){
                double sum = A[j][j] - s;
                if (sum <= 0){ 
                    return false;
                }
                L[j][j] = sqrt(sum); 
            }else{
                L[i][j] = (A[i][j] - s) / L[j][j];
                L[j][i] = 0;
                
            }
            
        }
    }
    return true;
}

double solution_distance(double **M, double *b, double *x, int order){
	double *Mx;
	double *diff;
	double erro;
	Mx = prod_matrix_vector(M,x,order,order); 
	
	// erro = || b-Mx || norma 2
	//erro = normv(2, sub_v(b,Mx,order) ,order);
	diff = sub_v(b,Mx,order);
	erro = normv(2, diff ,order);
	
	free(Mx);
	return erro;
}
