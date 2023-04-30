/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: fespunbe7.alumnes
 *
 * Created on 25 / d’octubre / 2018, 09:00
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int get_number_rows(double*** Matrix){
    for(int i = 0; i !=-1; i++){
        if(Matrix[i] == 0){
            return i;
        }
    }
}

int get_number_columns(double*** Matrix){
    for(int i = 0; i !=-1; i++){
        if(Matrix[0][i] == 0){
            return i;
        }
    }
}

void printM(double*** Matrix){
    
    int exists = 1;
    
    if(Matrix == 0){
        exists = 0;
    }else if(Matrix[0] == 0){
        exists = 0;
    }
    
    if(exists == 0){
        printf("La matriu no existeix");
            return;
    }
    
    int i = 0;
    int j = 0;
    
    while(Matrix[i] != 0){
        while(Matrix[0][j] != 0){
            printf("%7.3f   ", *Matrix[i][j]);
            j++;
        }
        printf("\n");
        i++;
        j = 0;
    }
    printf("\n");
}

void randomFill(double*** Matrix){
    for(int i = 0; Matrix[i] != 0; i ++){
        for(int j = 0; Matrix[i][j] != 0; j++){
            double f = (double)(rand() % 10);//(((double)(rand()%10000000))/100000.0);
            *Matrix[i][j] = f;
        }
    }
    
}

void gauss(double*** Matrix, double*** Matrix2, int n){
    
    double zero = 0.01;
    
    
    
    if(Matrix[n] == 0){
        return;
    }
    
    double biggest;
    int biggestrow;
    int column = n;
    double* first = Matrix[0][0];
    int done = 0;
    
    
    while(fabs(biggest) < zero && Matrix[0][column] != 0){
        done = 1;
        biggestrow = n;
        biggest = *Matrix[n][column];
        for (int i = n+1; Matrix[i] != 0; i++){
            if (fabs(*Matrix[i][column]) > fabs(biggest)){
                double k = *Matrix[i][column];
                biggestrow = i;
                biggest = (*Matrix[i][column]);
            }
        }
        column++;
    }
    
    if(done == 1){
        column--;
    }
    
    double* test = Matrix[0][column];
    
    if (test == 0){
        return;
    }
    
    
    double** temp = Matrix[n];
    Matrix[n] = Matrix[biggestrow];
    Matrix[biggestrow] = temp;
    
    if(Matrix2 != 0){
        temp = Matrix2[n];
        Matrix2[n] = Matrix2[biggestrow];
        Matrix2[biggestrow] = temp;
    }
    
    
    
    if(fabs(biggest) < zero){
        return;
    }
    
    for(int i = n + 1; Matrix[i] != 0; i++){
        if(fabs(*Matrix[i][column]) >= zero){
            double multiplier = *Matrix[i][column] / *Matrix[n][column];
            for(int j = 0; Matrix[i][j] != 0; j++){
                *Matrix[i][j] = *Matrix[i][j] - (*Matrix[n][j])*multiplier;
                if(Matrix2 != 0){
                    *Matrix2[i][j] = *Matrix2[i][j] - (*Matrix2[n][j])*multiplier;
                }
            }
        }
    }
    
    
    gauss(Matrix, Matrix2, n+1);
    
}

void jordan(double*** Matrix){
    
    gauss(Matrix, 0, 0);
    
    
    int n = get_number_rows(Matrix);
    int m = get_number_columns(Matrix);
    double zero = 0.01;
    
    double lastrow = -1;
    double lastrowAmpliada = -1;
    
    if(m < 2){
        return;
    }
    
    for(int i = n - 1; i >= 0; i--){
        if(lastrow == -1 && fabs(*Matrix[i][m-1]) > zero){
            lastrow = i;
        }

        if(lastrowAmpliada == -1 && fabs(*Matrix[i][m-2]) > zero){
            lastrowAmpliada = i;
        }
    }
    
    if(lastrow != lastrowAmpliada){
        printf("\nsistema incompatible");
        return;
    }
    
    if(lastrow < m-2){
        printf("\nsistema compatible indeterminat");
        return;
    }
    
    printf("\nsistema compatible determinat");
    
    
    for(int i = lastrow; i > 0; i--){
        for(int k = i-1; k >= 0; k--){
            
            double multiplier = (*Matrix[k][i])/(*Matrix[i][i]);
            
            for(int j = i; j < m; j++){
                *Matrix[k][j] = *Matrix[k][j] - (*Matrix[i][j])*multiplier;
            }
        }
        
    }
    
    for(int i = 0; i <= lastrow; i++){
        *Matrix[i][m-1] = (*Matrix[i][m-1]) / (*Matrix[i][i]);
        *Matrix[i][i] = 1.0;
    }
}

double*** getidentity(int n){
    double*** Matrix;
    Matrix = (double***)malloc((n+1)*sizeof(double**));
    for(int e = 0; e<n;e++){
        Matrix[e] = (double**)malloc((n+1)*sizeof(double*));
        for(int i = 0; i < n; i++){
            (Matrix[e][i]) = (double*)malloc(sizeof(double));
            if(e == i){
                *(Matrix[e][i]) = 1;
            }else{
                *(Matrix[e][i]) = 0; 
            }
        }
        Matrix[e][n] = 0;
    }
    (Matrix[n]) = 0;
    
    return Matrix;
}

double*** get0(int n, int m){
    double*** Matrix;
    Matrix = (double***)malloc((n+1)*sizeof(double**));
    for(int e = 0; e<n;e++){
        Matrix[e] = (double**)malloc((m+1)*sizeof(double*));
        for(int i = 0; i < m; i++){
            (Matrix[e][i]) = (double*)malloc(sizeof(double));
            *Matrix[e][i] = 0;
        }
        Matrix[e][m] = 0;
    }
    (Matrix[n]) = 0;
    
    return Matrix;
}

double*** copy(double*** Matrix){
    
    double*** MatrixCopy = get0(get_number_rows(Matrix), get_number_columns(Matrix));
    
    for(int i = 0; Matrix[i] != 0; i++){
        for(int j = 0; Matrix[i][j] != 0; j++){
            *MatrixCopy[i][j] = *Matrix[i][j];
        }
    }
    
    return MatrixCopy;
}

void squarejordan(double*** Matrix, double*** Matrix2){
    
    int m = get_number_rows(Matrix);
    
    for(int i = m-1; i > 0; i--){
        for(int j = i-1; j >= 0; j--){
            double multiplier = (*Matrix[j][i]) / (*Matrix[i][i]);
            for(int k = 0; k < m; k++){
                *Matrix[j][k] = (*Matrix[j][k]) - (*Matrix[i][k])*multiplier;
                *Matrix2[j][k] = (*Matrix2[j][k]) - (*Matrix2[i][k])*multiplier;
            }
        }
    }
    
    for(int i = 0; i < m; i++){
        for(int j = 0; j < m; j++){
            *Matrix2[i][j] = (*Matrix2[i][j])/(*Matrix[i][i]);   
        }
        *Matrix[i][i] = (*Matrix[i][i])/(*Matrix[i][i]);
    }
}

double*** product(double*** Matrix1, double*** Matrix2){
    /*NO SE SI VA*/
    
    int n1 = get_number_rows(Matrix1);
    int m1 = get_number_columns(Matrix1);
    
    int n2 = get_number_rows(Matrix2);
    int m2 = get_number_columns(Matrix2);
    
    if(m1 != n2){
        return 0;
    }
    
    double*** product = get0(n1, m2);
    
    for(int i = 0; product[i] != 0; i++){
        for(int j = 0; product[i][j] != 0; j++){
            for(int k = 0; Matrix1[k] != 0; k++){
                *product[i][j] = *product[i][j] + (*Matrix1[i][k])*(*Matrix2[k][j]);
            }
        }
    }
    
    return product;
    
}

double*** transpose(double*** Matrix){
    
    int n = get_number_rows(Matrix);
    int m = get_number_columns(Matrix);
    double*** transposed = get0(m, n);
    
    for(int i = 0; Matrix[i] != 0; i++){
        for(int j = 0; Matrix[i][j] != 0; j++){
            transposed[j][i] = Matrix[i][j];
        }
    }
    
    return transposed;
}

double*** alu (double*** Matrix, int n){
    
    double zero = 0.01;
    
    if(Matrix[n] == 0){
        return Matrix;
    }
    
    if(fabs(*Matrix[n][n]) < zero){
        return 0;
    }
    
    for(int i = n + 1; Matrix[i] != 0; i++){
        double multiplier = (*Matrix[i][n]) / (*Matrix[n][n]);
        for(int j = n + 1; Matrix[i][j] != 0; j++){
            *Matrix[i][j] = (*Matrix[i][j]) - (*Matrix[n][j]) * multiplier;
        }
        *Matrix[i][n] = multiplier;
    }
    
    alu(Matrix, n+1);
}

double*** choleski(double*** Matrix){
    
    double zero = 0.01;
    
    double*** MatrixAlu = alu(Matrix, 0);
    if (MatrixAlu == 0){
        return 0;
    }
    
    double multiplier;
    
    for(int i = 0; Matrix[i] !=  0; i++){
        if(*Matrix[i][i] < zero){
            return 0;
        }
        
        multiplier = sqrt(*MatrixAlu[i][i]);
        
        for(int j = i; Matrix[i][j] != 0; j++){
            *MatrixAlu[j][i] = (*MatrixAlu[j][i])*multiplier;
            *MatrixAlu[i][j] = (*MatrixAlu[j][i]);
            
        }
        
        *MatrixAlu[i][i] = multiplier;
    }
    
    return MatrixAlu;
    
}

double*** choleski2(double*** Matrix){
    
    double zero = 0.01;
    
    for(int i = 0; Matrix[i] != 0; i++){
        for(int k = 0; k < i; k++){
            
            if (Matrix[i][k] == 0){
                return 0;
            }
            
            for (int j = 0; j < k; j++){
                *Matrix[i][k] = *Matrix[i][k] - (*Matrix[i][j])*(*Matrix[k][j]);
            }
            *Matrix[i][k] = (*Matrix[i][k])/(*Matrix[k][k]);
            
        }
        
        for(int j = 0; j < i; j++){
            *Matrix[i][i] = *Matrix[i][i] - (*Matrix[i][j])*(*Matrix[i][j]);
        }
        *Matrix[i][i] = sqrt(*Matrix[i][i]);
        
        
    }
    return Matrix;
    
    
}

int main(int argc, char** argv) {
    
    int n = 3;
    
    double*** Matrix1 = getidentity(n);
    
    *Matrix1[0][0] = 4;
    *Matrix1[0][1] = -1;
    *Matrix1[0][2] = 1;
    *Matrix1[1][0] = -1;
    *Matrix1[1][1] = 4.25;
    *Matrix1[1][2] = 2.75;
    *Matrix1[2][0] = 1;
    *Matrix1[2][1] = 2.75;
    *Matrix1[2][2] = 3.5;
    
    //randomFill(Matrix1);
    
    printM(Matrix1);
    
    double*** choleski1 = choleski2(Matrix1);
    printM(choleski1);
    
    


    return (EXIT_SUCCESS);
}