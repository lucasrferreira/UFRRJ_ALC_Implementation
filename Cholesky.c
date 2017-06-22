#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// #include "Util.c"

bool is_pos_definite_chol(double **A, int order);
double **chol_fatoration(double **A, int order);

bool is_pos_definite_chol(double **A, int order){
    
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

double **chol_fatoration(double **A, int order){
    
    double **L = allocate_matrix_double(order, order);
    copy_matrix(L,A,order);
    int n = order;
    for (int i = 0; i < n; i++){
        for (int j = 0; j <= i ; j++) {
            double s = 0;
            for (int k = 0; k < j; k++)
                s += L[i][k] * L[j][k];
                
            if(i == j){
                double sum = A[j][j] - s;
                L[j][j] = sqrt(sum); 
            }else{
                L[i][j] = (A[i][j] - s) / L[j][j];
                L[j][i] = 0;
            }
            
        }
    }
    return L;
}