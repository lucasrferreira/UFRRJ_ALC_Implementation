#include <stdio.h>
#include <stdlib.h>


void read_matrix_from_file(char *path_To_file, double ***M, double **b, int *matrix_order)
{
    
    // double **M;

    FILE *myFile;
    myFile = fopen(path_To_file, "r");

    fscanf(myFile, "%d", matrix_order);

    *M = malloc( *matrix_order * sizeof(double**));
    for (int i = 0; i < *matrix_order; i++) {
        (*M)[i] = malloc( *matrix_order * sizeof(double*));
        for (int j = 0; j < *matrix_order; j++) {
            // M[i][j] = malloc(sizeof(double));
            // (*M)[i][j] = 8465341;
            fscanf(myFile, "%lf", &(*M)[i][j]) ;
        }
    }
    
    
    *b = malloc(*matrix_order * sizeof(double*));
    for (int j = 0; j < *matrix_order; j++) {
        // (*b)[j] = 0;
        fscanf(myFile, "%lf", &(*b)[j]);
    }
    
    // print_matrix(M, matrix_order, matrix_order);
    // print_vector(b, matrix_order);
    
    
    fclose(myFile);

}

void read_vector_from_file(char *path_To_file, double **b, int *v_order)
{
    
    // double **M;

    FILE *myFile;
    myFile = fopen(path_To_file, "r");

    fscanf(myFile, "%d", v_order);
    
    *b = malloc(*v_order * sizeof(double*));
    for (int j = 0; j < *v_order; j++) {
        // (*b)[j] = 0;
        fscanf(myFile, "%lf", &(*b)[j]);
    }
    
    fclose(myFile);

}
