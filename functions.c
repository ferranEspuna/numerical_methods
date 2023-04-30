#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double PI = 3.14159265359;

// get n by m zero matrix
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

//get n by n triangular zero matrix
double*** getTriangular0(int n){
    
    double*** Matrix;

    Matrix = (double***)malloc((n+1)*sizeof(double**));
    
    for(int i = 0; i < n; i++){
        
        Matrix[i] = (double**)malloc((i+2)*sizeof(double*));
        
        for(int j = 0; j <= i; j++){
            
            (Matrix[i][j]) = (double*)malloc(sizeof(double));
            *Matrix[i][j] = 0;
            
        }
        
        Matrix[i][i+1] = 0;
    }
    Matrix[n] = 0;
    
    return Matrix;
}

//get 1 by n zero matrix, or "zero vector", n-long array
double** get0array(int n){
    return get0(1, n)[0];
}

//number of rows in a matrix
int get_number_rows(double*** Matrix){
    for(int i = 0; i !=-1; i++){
        if(Matrix[i] == 0){
            return i;
        }
    }
    return 0;
}

//number of columns in a matrix
int get_number_columns(double*** Matrix){
    for(int i = 0; i !=-1; i++){
        if(Matrix[0][i] == 0){
            return i;
        }
    }
    return 0;
}

//make a copy of a matrix
double*** copy(double*** Matrix){
    
    double*** MatrixCopy = get0(get_number_rows(Matrix), get_number_columns(Matrix));
    
    for(int i = 0; Matrix[i] != 0; i++){
        for(int j = 0; Matrix[i][j] != 0; j++){
            *MatrixCopy[i][j] = *Matrix[i][j];
        }
    }
    
    return MatrixCopy;
}

//print a matrix
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
        while(Matrix[i][j] != 0){
            printf("%9.6f   ", *Matrix[i][j]);
            j++;
        }
        printf("\n");
        i++;
        j = 0;
    }
    
    printf("\n");
}

//print an array
void printarray(double** Matrix){
    
    int exists = 1;
    
    if(Matrix == 0){
        exists = 0;
    }
    
    if(exists == 0){
        printf("L'array no existeix");
            return;
    }
    
    int i = 0;
    
    while(Matrix[i] != 0){
        printf("%7.15f   ", *Matrix[i]);
        i++;
    }
    printf("\n");
    
}

//get an array of n+1 equidistant points (with n spaces between them)
double** equidistants(int grau, double start, double end){
    
    double** abscisses = get0array(grau + 1);
    
    int numpoints = grau + 1;
    double range = end - start;
    double step = range/grau;
    
    
    for(int k = 0; k < numpoints; k++){
        double value = start + (double) k * step;
        *abscisses[k] = value;
    }
    
    return abscisses;
    
}

//get the roots of the n+1st txebitxev polynomial, remapping the usual interval [-1, 1] to the desired interval [start, end]
//this is useful to get the nth degree interpolating polynomial to a function with minimal error
double** Txebitxev(int grau, double start, double end){
    
    double** abscisses = get0array(grau + 1);
    
    
    for(int k = 0; k <= grau; k++){
        double value = cos((PI/2.0 + k*PI)/(grau+1)) ;
        *abscisses[grau - k] = value;
    }
    
    double multiply = (end - start)/2.0;
    double add = (end + start)/2.0;
    
    for(int i = 0; abscisses[i] != 0; i++){
        *abscisses[i] *= multiply;
        *abscisses[i] += add;
    }
    
    return abscisses;
    
    
}

//evaluate a function fn at a set of n+1 points between start and end, selected according to one of the previous two funcitons
//the purpose of this is to then interpolate the function with an nth degree polynomial that coincides with it at the desired points
//returns both x and y values
double*** evalfunc(double (*fn) (double), int degree, double start, double end, int txebitxev){
    
    int numpoints = degree + 1;
    
    double*** vectors = get0(2, degree + 1);
    
    if(txebitxev == 0){
        vectors[0] = equidistants(degree, start, end);
    }else{
        vectors[0] = Txebitxev(degree, start, end);
    }
    
    for(int i = 0; i < numpoints; i++){
        *vectors[1][i] = fn(*vectors[0][i]);
    }
    
    return vectors;
}

//the same as before, but with two functions, which are supposed to represent a function and its derivative
double*** evalfunciderivada(double (*fn) (double), double (*der) (double), int degree, double start, double end, int txebitxev){
    
    double*** uno = evalfunc(fn, degree, start, end, txebitxev);
    double*** dos = evalfunc(der, degree, start, end, txebitxev);
    
    double*** resposta = get0(3, degree + 1);
    
    resposta[0] = uno[0];
    resposta[1] = uno[1];
    resposta[2] = dos[1];
    
    return resposta;
}

//evaluate the lagrange interpolating polynomial of a set of x and y values at a desired x value. if x, y is in dades, evalLagrange(dades, x) = y
double evalLagrange(double*** dades, double x){
    
    double** abscisses = dades[0];
    double** ys = dades[1];
    
    double total = 0;
    double current;
    
    for(int i = 0; abscisses[i] != 0; i++){
        
        current = 1;
        
        for(int j = 0; abscisses[j] != 0; j++){
            
            if(i != j){
                current *= (x - *abscisses[j]);
                current /= (*abscisses[i] - *abscisses[j]);
            }
        }
        
        total += current * (*ys[i]);
        
    }
    
    return total;
    
}

//a particular function and its derivative
double func(double x){
    return pow(x, 2.0) - 2.0;
}

double derivada(double x){
    return 2.0*x;
}

//calculate newton's divided differences for particular x and y values
double*** Newton (double*** dades){
    
    int degree = 0;
    while(dades[0][degree + 1] != 0){
        degree++;
    }
    
    for(int pas =  1; pas <= degree; pas++){
        for(int k = degree; k >= pas; k--){
            int i = k - pas; 
            
            *dades[1][k] = (*dades[1][k] - *dades[1][k-1])/(*dades[0][k] - *dades[0][i]);
        }
    }
    
    return dades;
}

//given the above computation, evaluate newton's interpolating polynomial at a point x
double evalNewton(double*** difdiv, double x){
    
    int i = get_number_columns(difdiv) - 1;
    
    double total = 0;
    
    for(i = i; i >= 0; i--){
        
        total *= (x - *difdiv[0][i]);
        total += *difdiv[1][i];
    }
    return total;
    
}

//same for hermite's divided diferences, which interpolate both values for y(x_i) and y'(x_i)
double *** Hermite(double*** dades){
    
    double degree = get_number_columns(dades) - 1;
    double*** difdiv = get0(2, 2*degree + 2);
    
    for(int i = 0; dades[0][i] != 0; i++){
        
        *difdiv[0][2*i] = *dades[0][i];
        *difdiv[0][2*i + 1] = *dades[0][i];
        *difdiv[1][2*i] = *dades[1][i];
        *difdiv[1][2*i + 1] = *dades[1][i];
    }
    
    
    for(int pas = 1; pas <= 2*degree + 1; pas++){
        
        for(int k = 2*degree + 1; k >= pas; k--){
            int i = k - pas;
            
            if(pas == 1 && i%2 == 0){
                *difdiv[1][k] = *dades[2][i/2];
            }else{
                *difdiv[1][k] = (*difdiv[1][k] - *difdiv[1][k-1])/(*difdiv[0][k] - *difdiv[0][i]);
            }
            
        }
        
        
    }
    
    return difdiv;
}

//write a matrix to a file
void tofile(char* str, double*** dades){
    
    FILE * fp;
    fp = fopen(str ,"w");
    
    for(int i = 0; dades[0][i] != 0; i++){
        fprintf (fp, "%f %f \n", *dades[0][i], *dades[1][i]);
    }
}

//symmetrical numerical derivative (f(x+dx)-f(x-dx))/(2dx) of a function f at a point x, with dx = h
double derivada_simetrica(double (*fn) (double), double x, double h){
    
    return (fn(x+h) - fn(x-h))/(2*h);
    
}

//calculate richardson's increasingly better approximations of the derivative, which are better than the symmetrical option before, by cleverly solving approximate systems of equations
double Richardson(double (*fn) (double), double x, double h0, double q, double precisio, int limit){

    double*** Matrix = getTriangular0(limit);
    double h = h0;
    
    for(int i = 0; Matrix[i] != 0; i++){
        *Matrix[i][0] = derivada_simetrica(fn, x, h);
        
        for(int j = 1; j <= i; j++){
            
            *Matrix[i][j] = *Matrix[i][j-1] + (1.0/(pow(q, (double) j) - 1)) * (*Matrix[i][j-1] - *Matrix[i-1][j-1]);
            
            if(j == i){
                
                double diferencia = fabs(*Matrix[i][j] - *Matrix[i][j-1]);
                if(diferencia < precisio){
                    
                    printM(Matrix);
                    return *Matrix[i][j];
                }
            }
        }
        
        h /= q;
    }
    
    printM(Matrix);
    return *Matrix[limit - 1][limit - 1];
}

//calculate a better approximation than trapezis for the intergral of a function by doubling the "resolution" of a first trapezoidal aproximation, without recalculating everything
double mestrapezis(double (*fn) (double), double trapezis, int n, double start, double end){
    
     
    double** punts = equidistants(2*n, start, end);
    double total = 0;
    
    for(int i = 0; i < n; i++){
        
            total += fn(*punts[2*i + 1]);
            
        }
    
    total *= ((end - start)/n);
    total += trapezis;
    total /= 2.0;
    
    
    return total;
}

//calculate the area of a slice under a function via a trapezoidal approximation
double trapezi(double (*fn) (double), double start, double end){
 
    return (fn(start) + fn(end)) * (end - start) / 2.0;
}

//similar to richardson's trick, but for integrals. Uses the above functions to avoid evaluating the function too many times
double Romberg(double (*fn) (double), double start, double end, double precisio, int limit){

    double*** Matrix = getTriangular0(limit);
    int n = 1;
    double q = 2; //DO NOT CHANGE. I HAVE NOT GENERALIZED THIS. (mestrapezis would have to be more general, which I have not done)
    double trapezis = 0;
    
    for(int i = 0; Matrix[i] != 0; i++){
        
        if(n == 1){
            trapezis = trapezi(fn, start, end);
        }else{
            trapezis = mestrapezis(fn, trapezis, n, start, end);
        }
        
        *Matrix[i][0] = trapezis;
        
        for(int j = 1; j <= i; j++){
            
            *Matrix[i][j] = *Matrix[i][j-1] + (1.0/(pow(q, (double) j) - 1)) * (*Matrix[i][j-1] - *Matrix[i-1][j-1]);
            
            if(j == i){
                
                double diferencia = fabs(*Matrix[i][j] - *Matrix[i][j-1]);
                if(diferencia < precisio){
                
                    printM(Matrix);
                    return *Matrix[i][j];
                }
            }
        }
        
        n *= q;
    }
    
    printM(Matrix);
    return *Matrix[limit - 1][limit - 1];
}

//Newton-Rhapson's root-finding algorithm
double**** NR(double (*fn) (double), double (*der) (double), double start, int imax, double prec1, double prec2){
    
    
    double**** ret;
    
    
    ret = (double****)malloc((3)*sizeof(double***));
    
    
    double convergencia = 0;
    
    double*** llista = get0(2, imax);
    *llista[0][0] = start;
    *llista[1][0] = fn(start);
    
    ret[2] = llista;
    
    double x = start;
    
    

    for (int i = 1; llista[0][i] != 0; i++){
        
        
        double prevx = *llista[0][i - 1];
        double prevy = *llista[1][i - 1];
        
        
        
        double newx = prevx - prevy / der(prevx);
        
        double newy = fn(newx);
        
        *llista[0][i] = newx;
        *llista[1][i] = newy;
        
        
        
        x = newx;
        
        
        if( (fabs(prevx - newx) < prec1) || (fabs(newy) < prec2) ){
            
            
            
            llista[0][i + 1] = 0;
            llista[1][i + 1] = 0;
            
           convergencia = 1;
           break;
            
        }
    }
    
    double*** mat0 = get0(1,1);
    double*** mat1 = get0(1,1);
    
    *mat0[0][0] = convergencia;
    *mat1[0][0] = x;
    
    ret[0] = mat0;
    ret[1] = mat1;
    
        
    return ret;
    
}

//secant root-finding algorithm, useful when we can't calculate the derivative
double**** secant(double (*fn) (double), double start1, double start2, int imax, double prec1, double prec2){
    
    
    double**** ret;
    
    
    ret = (double****)malloc((3)*sizeof(double***));
    
    
    double convergencia = 0;
    
    double*** llista = get0(2, imax);
    *llista[0][0] = start1;
    *llista[1][0] = fn(start1);
    *llista[0][1] = start2;
    *llista[1][1] = fn(start2);
    
    ret[2] = llista;
    
    double x = start2;
    
    

    for (int i = 2; llista[0][i] != 0; i++){
        
        
        double prevx = *llista[0][i - 1];
        double prevy = *llista[1][i - 1];
        
        double prevprevx = *llista[0][i - 2];
        double prevprevy = *llista[1][i - 2];
        
        double der = (prevy - prevprevy)/(prevx - prevprevx);
        
        
        
        double newx = prevx - prevy / der;
        
        double newy = fn(newx);
        
        *llista[0][i] = newx;
        *llista[1][i] = newy;
        
        
        
        x = newx;
        
        
        if( (fabs(prevx - newx) < prec1) || (fabs(newy) < prec2) ){
            
            
            
            llista[0][i + 1] = 0;
            llista[1][i + 1] = 0;
            
           convergencia = 1;
           break;
            
        }
    }
    
    double*** mat0 = get0(1,1);
    double*** mat1 = get0(1,1);
    
    *mat0[0][0] = convergencia;
    *mat1[0][0] = x;
    
    ret[0] = mat0;
    ret[1] = mat1;
    
        
    return ret;
    
}

//try to find approximate locations of roots of a function by 'sweeping' an interval and then looking at sign changes and halving the interval each time
//we assume the initial sweep is fine enough as to always capture roots as sign changes. This may not be the case, especially for "even" roots.
double** escombrada(double (*fn) (double), int n, double minx, double maxx, double pas){
    
    double** ret = get0array(n);
    
    int i = 0;
    
    double x = minx;
    double fx = fn(x);
    
    double newx;
    double newfx;
    
    double zero;
    
    while(x < maxx){
        
        
        if (ret[i] == 0){
            
            return ret;
            
        }
        
        newx = x + pas;
        newfx = fn(newx);
        
        if(fx * newfx <= 0){
            
            zero = (x + newx)/2.0;
            *ret[i] = zero;
            
            i += 1;
            
        }
        
        x = newx;
        fx = newfx;
        
        
    }
    
    
    ret[i] = 0;
    return ret;
    
}

//perform newton-rhapson at the points found with the previous funciton
double** manyNR(double (*fn) (double), double(*der) (double), int n, double minx, double maxx, double pas){
    
    int imax = 10;
    double prec1 = 0.0000001;
    double prec2 = 0.0000001;
    
    double** startpoints = escombrada(fn, n, minx, maxx, pas);
    
    for(int i = 0; startpoints[i] != 0; i++){
        
        double startx = *startpoints[i];
        *startpoints[i] = *NR(fn, der, startx, imax, prec1, prec2)[1][0][0];
        
    }
    
    return startpoints;
    
    
}

//same thing with the secant method
double** manySecant(double (*fn) (double), int n, double minx, double maxx, double pas){
    
    int imax = 10;
    double prec1 = 0.0000001;
    double prec2 = 0.0000001;
    
    double** startpoints = escombrada(fn, n, minx, maxx, pas);
    
    for(int i = 0; startpoints[i] != 0; i++){
        
        double startx1 = *startpoints[i] - pas/2.0;
        double startx2 = *startpoints[i] + pas/2.0;
        *startpoints[i] = *secant(fn, startx1, startx2, imax, prec1, prec2)[1][0][0];
        
    }
    
    return startpoints;
    
    
}

//evaluate a polynomial at a single point efficiently (O(n))
double hornernatural(double** polinomi, double x){
    
    int i;
    
    for (i = 0; polinomi[i] != 0; i++){
    }
    
    i--;
    
    double total = 0;
    
    for(i = i; i >= 0; i--){
        
        total *= x;
        total += *polinomi[i];
        
    }
    
    return total;
    
}

//same, but we calculate the derivative while we're at it
double hornerderivada(double** polinomi, double x){
    
    int i;
    
    for (i = 0; polinomi[i] != 0; i++){
    }
    
    i -= 2;
    
    double total = 0;
    
    for(i = i; i >= 0; i--){
        
        total *= x;
        total += *polinomi[i + 1] * (i + 1);
        
    }
    
    return total;
    
}

//bound for the roots of a polynomial (I don't know how to prove it or where it comes from)
double rootbound(double** polinomi){
    
    double cota = 0;
    
    int i;
    for (i = 0; polinomi[i] != 0; i++){
    }
    i--;
    
    for(int j = i - 1; j >= 0; j--){
        
        cota += fabs(*polinomi[j] / *polinomi[i]);
    }
    
    if(cota < 1){
        return 1;
    }
    return cota;
    
    
}


//(try to) find all roots of a polynomial with newton's method
double** findrootsNR(double** polinomi, double pas){
    
    int i;
    //find the degree of the polynomial
    for (i = 0; polinomi[i] != 0; i++){
    }
    i--;
    
    double bound = rootbound(polinomi);
    
    double fn(double x){
        return hornernatural(polinomi, x);
    }
    
    double der(double x){
        return hornerderivada(polinomi, x);
    }
    
    return manyNR(fn, der, i, - bound, bound, pas);
    
    
}

//(try to) find all roots of a polynomial with the secant method
double** findrootsSecant(double** polinomi, double pas){
    
    int i;
    for (i = 0; polinomi[i] != 0; i++){
    }
    i--;
    
    double bound = rootbound(polinomi);
    
    double fn(double x){
        return hornernatural(polinomi, x);
    }    
    return manySecant(fn, i, - bound, bound, pas);
    
    
}

//defining some example functions for testing
/*  
double composicio(double x, int n){
    
    if (n == 0){
        
        return x;
        
    }else{
        
        double prev = composicio(x, n - 1);
        return 4.0*prev*(1.0 - prev);
        
    }
    
}


double polinomin(double x, int n){
    
    return composicio(x, n) - x;
    
}

double derivada1(double x){
    
    return 4.0 - 8.0*x;
    
}

double derivadacomp(double x, int n){
    
    double der;
    
    
    if (n == 0){
        
        der = 1;
        
    }else{
        
        double a = composicio(x, n - 1);
        double der1 = derivada1(a);
        double der2 = derivadacomp(x, n - 1);
        
        der = der1 * der2;
        
    }
    
    return der;
    
}

double derivadan(double x, int n){
    
    return derivadacomp(x, n) - 1;
    
}



int composicions = 2;
double polinomiconcret(double x){
    
    return polinomin(x, composicions);
    
}

double derivadaconcreta(double x){
    
    return derivadan(x, composicions);
    
    
}
*/

/*
double polinomitest(double x){
    
    double** polinomi = get0array(4);
    *polinomi[0] = 0.9;
    *polinomi[1] = -3;
    *polinomi[2] = 0;
    *polinomi[3] = 1;
    
    return hornernatural(polinomi, x);
    
}

double derivadatest(double x){
    
    double** polinomi = get0array(4);
    *polinomi[0] = 0.9;
    *polinomi[1] = -3;
    *polinomi[2] = 0;
    *polinomi[3] = 1;
    
    return hornerderivada(polinomi, x);
    
}
*/

//example use cases
int main(int argc, char** argv){

    /*
    double**** ret1 = NR(func, derivada, 5, 10, 0.0001, 0.000001);
    double**** ret2 = secant(func, 3, 4, 10, 0.00001, 0.0000001);
    
    printM(ret2[2]);
    */
    
    
    /*
    composicions = 6;
    printarray(manyNR(polinomiconcret, derivadaconcreta, escombrada(polinomiconcret, 100, -1, 2, 0.0001)));
    */
    
    
    /*
    printf("%f\n", hornernatural(polinomi, 0.5));
    printf("%f\n", hornerderivada(polinomi, 0.5));
    */
    
    
    /*
    printarray(manyNR(polinomitest, derivadatest, 10, -2, 2, 0.1));
    printarray(manySecant(polinomitest, 10, -2, 2, 0.1));
    */
    
    
    /*
    double** polinomi = get0array(4);
    *polinomi[0] = 0.9;
    *polinomi[1] = -3;
    *polinomi[2] = 0;
    *polinomi[3] = 1;
    
    printarray(findrootsNR(polinomi, 0.1));
    printarray(findrootsSecant(polinomi, 0.1));
    */
    
    
}

