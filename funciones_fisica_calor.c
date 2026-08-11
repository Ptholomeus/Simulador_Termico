/*
	p66 del pdf
	Sea el dominio L. Luego 0 <= x <= L.
	L se divide en M regiones iguales, de ancho Dx = L/M.
	Se tienen entonces M+1 nodos, de i=0 a i=M.
	
	La formula discretizada es T_i-1 - 2*T_i + T_i+1 + (1/k)*Dx*g_i = 0.
	
	Si T_0 = 300, T_M = 400, la formula es valida para i entre 1 y M-1.
	
	usar k = 1.0 y g = 0.
	
	Se busca un vector v solucion (T_0, T_1, ···, T_M).
	Este vector se sabe que es (300, T_1, ···, 400).
	
	Cada T_i conseguido con la fórmula tiene un coeficiente. los coeficientes
	de T_0 y T_M son 1.0. Se tiene entonces una matriz A de coeficientes.
	
	Cada fila de la matriz corresponde a un i, entonces hay tantas filas como i's.
	Cada colunba de la matriz corresponde a un T_i. Entonces hay tantas columnas
	como i's.
	La matriz es cuadrada. dim (M+1)x(M+1)
	
	Este problema es de la forma Av = b, donde b es un vector constante.
	La matriz A se construye con los coeficientes de T_i's.

*/
#include <math.h>

#include "funciones_matrices.h"
#include "funciones_fisica_calor.h"


// matriz "hardcoded"
Matriz crear_matriz_coefs_1Dsteadystate (double L, int M){
	
	double Dx = L / (double) M;
	
	Matriz coefs;
	coefs = crear_matriz(M+1, M+1);
	
	//Inicializar matriz
	
	coefs.matriz[0][0] = 1.0;
	coefs.matriz[M][M] = 1.0;
	
	for (int i = 1; i < M; i++){
		coefs.matriz[i][i-1] =  1.0;
		coefs.matriz[i][i+0] = -2.0;
		coefs.matriz[i][i+1] =  1.0;
	}
		
	return coefs;
	
}

// matriz "hardcoded"

Matriz crear_matriz_coefs_1Dsteadystate_2 (double L, int M){
	
	double Dx = L / (double) M;
	
	Matriz coefs;
	coefs = crear_matriz(M+1, M+1);
	
	//Inicializar matriz
	
	// Fila 0
	coefs.matriz[0][0] =  0.0;
	coefs.matriz[0][1] = -2.0;
	coefs.matriz[0][2] =  1.0;
	
	// Fila M
	coefs.matriz[M][M-2] =  1.0;
	coefs.matriz[M][M-1] = -2.0;
	coefs.matriz[M][M-0] =  0.0;
	
	for (int i = 1; i < M; i++){
		coefs.matriz[i][i-1] =  1.0;
		coefs.matriz[i][i+0] = -2.0;
		coefs.matriz[i][i+1] =  1.0;
	}
		
	return coefs;
	
}


// aaaaaaaaaaaaaaaaaaa
Matriz crear_matriz_const_1Dsteadystate (int M){
	
	Matriz constante;
	constante = crear_matriz(M+1, 1);
	//iniciar_matriz_cero (constante);
	
	constante.matriz[0][0] = 300.0; // T_0
	constante.matriz[M][0] = 400.0; // T_M
	
	return constante;
	
}

// aaaaaaaaaaaaaaaaaaa

Matriz crear_matriz_const_1Dsteadystate_2 (int M, double T_0, double T_M){
	
	Matriz constante;
	constante = crear_matriz(M+1, 1);
	//iniciar_matriz_cero (constante);
	
	constante.matriz[0][0] = -T_0; // T_0
	constante.matriz[M][0] = -T_M; // T_M
	
	return constante;
	
}
 


/*
	p131 libro Osizik
	Resolver transitorio, dado L, M, N, T(0,n) y T(M,n). Es decir, con potencial
	prescrito en la frontera. Nota: Es similar entonces a un estacionario local
	en el tiempo?
*/
/*
	Argumentos:
	L : Longitud de la barra
	n : numero del instante
	Dt : longitud del instante
	M : cantidad de nodos espaciales en L
	N : cantidad de nodos temporales en tan
	T_0 : Temperatura de frontera izquierda
	T_M : Temperatura de frontera derecha
	acumulado: Matriz que acumula soluciones
	alfa : numero alfa 
	solucion_fila : fila que almacena la solucion
	
	1 Se calcula Dx 
	2 Se calcula r, se valida, se almacena
	3 Se asigna el valor del nodo 0 a la solucion
	4 Se valida el numero del instante, luego se calculan los nodos intermedios
	5 Se asigna el valor al nodo M de la solucion
	
	Si n = 1 se tiene el caso donde se puede usar el acumulado como memoria del paso previo.
	
*/
void matriz_1DSimpleTransient(double L, int n, double Dt, int M, int N, double T_0, double T_M, Matriz *acumulado, double alfa, Matriz *solucion_fila) {
	
	double Dx = L / (double) M;
	double r  = alfa * Dt / (Dx * Dx);
	
	// resolver 
	solucion_fila->matriz[0][0] = T_0;	
	if (n < 0){
		printf("Error! Tiempo negativo!\n");
		solucion_fila->datos = NULL;
		solucion_fila->matriz = NULL;
		}
	if (n == 0) {
		for(int i = 0; i <= M; i++)
			solucion_fila->matriz[0][i] = acumulado->matriz[0][i];
		}	
	if (n >= 1) {
		// printf("entre al bucle\n");
		//solucion_fila->matriz[0][1] = r*T_0 + (1.0-2.0*r)*acumulado->matriz[n-1][1] + r*acumulado->matriz[n-1][2];
		for(int idx = 1; idx < M; idx++) {			
			solucion_fila->matriz[0][idx] = r*acumulado->matriz[n-1][idx-1] + (1.0-2.0*r)*acumulado->matriz[n-1][idx+0] + r*acumulado->matriz[n-1][idx+1];
		}
	}
	solucion_fila->matriz[0][M] = T_M;
	//Inicializar matriz
	
}

/*
	A definir...
*/

/*
	Orden correcto de la malla:
	i es y
	j es x
*/
int es_frontera(int i, int j, double malla_x, double malla_y){
	if (i == 0 || i == malla_y || j == 0 || j == malla_x)
		return 1;
	else 
		return 0;
}

double T_pared(int i,int j, double malla_x, double malla_y, double T_N, double T_S, double T_E, double T_O){
	
	if (i == 0 && (j > 0 && j < malla_x))
		return T_N; 
	if (i == malla_y && (j > 0 && j < malla_x))
		return T_S;
	if (j == 0 && (i > 0 && i < malla_y))
		return T_O;
	if (j == malla_x && (i > 0 && i < malla_y))
		return T_E;
	
	if (j == 0 && i == 0)
		return 0.5 * (T_N + T_O);
	if (j == 0 && i == malla_y)
		return 0.5 * (T_S + T_O);
	if (j == malla_x && i == 0)
		return 0.5 * (T_N + T_E);
	if (j == malla_x && i == malla_y)
		return 0.5 * (T_S + T_E);
	
	return -1.0;
}

Matriz crear_coefs_multiDimSteadyState_2D(double x_len, double y_len, int malla_x, int malla_y){

	/*
		Ec: div(grad(T)) = 0;
	*/
	
	/*
		La grilla es una malla de x_len * y_len, separado en malla_x+1 nodos en x y malla_y+1 nodos en y.
		
		Hay 4 paredes, las 4 tendran temperaturas distintas.
		Serán entonces, pared Norte, Sur, Este, Oeste.
		
		En la malla, estas son:
		
		Norte:
			desde i=0,j=1 hasta i=0, j=malla_x-1
		Sur:
			desde i=malla_y,j=1 hasta i=malla_y, j=malla_x-1
		Oeste:
			desde i=1,j=0 hasta i=malla_y-1, j=0
		Este: 
			desde i=1,j=malla_x hasta i=malla_y-1, j=malla_x
			
		Las esquinas serán:
		
		NO:
			0.5*(i=0,j=1 + i=1,j=0)
		NE: 
			0.5*(i=0,j=malla_x-1 + i=1, j=malla_x)
		SO:
			0.5*(i=malla_y-1 j=0 + i=malla_y, j=1)
		SE:
			0.5*(i=malla_y-1 j=malla_x + i=malla_y, j=malla_x-1)
	*/
	
	double Dx;
	double Dy;
	
	int k;
	
	double beta;
	
	Matriz coeficientes;
	int dim_matriz_coef = (malla_x + 1)* (malla_y + 1);
	coeficientes = crear_matriz(dim_matriz_coef, dim_matriz_coef);
	
	Dx = x_len / (double) malla_x;
	Dy = y_len / (double) malla_y;
	
	beta = Dx / Dy;
	
	/*
		Llenar la matriz de coeficientes 
	*/
	
	for (int i = 0; i <= malla_y ; i++){
		for (int j = 0; j <= malla_x; j++){
			// Ver bien por qué
			k = i * (malla_x + 1) + j;
			if (es_frontera(i, j, malla_x, malla_y)){
				coeficientes.matriz[k][k] = 1.0;
			} else {
				// Por qué hay que cambiar de lugar los coeficientes???
				coeficientes.matriz[k][k] = -2.0 * (1.0 + beta * beta);
				coeficientes.matriz[k][k-1] = 1.0;
				coeficientes.matriz[k][k+1] = 1.0;
				coeficientes.matriz[k][k-(malla_x+1)] = beta * beta;
				coeficientes.matriz[k][k+(malla_x+1)] = beta * beta;
			}
		}
	}

	return coeficientes;
}

Matriz crear_const_multiDimSteadyState_2D(int malla_x, int malla_y, double T_N, double T_S, double T_E, double T_O){

	/*
		Ec: div(grad(T)) = 0;
	*/
	
	/*
		La grilla es una malla de x_len * y_len, separado en malla_x+1 nodos en x y malla_y+1 nodos en y.
		
		Hay 4 paredes, las 4 tendran temperaturas distintas.
		Serán entonces, pared Norte, Sur, Este, Oeste.
		
		En la malla, estas son:
		
		Norte:
			desde i=0,j=1 hasta i=0, j=malla_x-1
		Sur:
			desde i=malla_y,j=1 hasta i=malla_y, j=malla_x-1
		Oeste:
			desde i=1,j=0 hasta i=malla_y-1, j=0
		Este: 
			desde i=1,j=malla_x hasta i=malla_y-1, j=malla_x
			
		Las esquinas serán:
		
		NO:
			0.5*(i=0,j=1 + i=1,j=0)
		NE: 
			0.5*(i=0,j=malla_x-1 + i=1, j=malla_x)
		SO:
			0.5*(i=malla_y-1 j=0 + i=malla_y, j=1)
		SE:
			0.5*(i=malla_y-1 j=malla_x + i=malla_y, j=malla_x-1)
	*/
	
	int k;
	
	Matriz constantes;
	int dim_matriz_coef = (malla_y + 1)*(malla_x + 1);
	constantes = crear_matriz(dim_matriz_coef, 1);

	/*
		Llenar la matriz de constantes 
	*/
	
	for (int i = 0; i <= malla_y ; i++){
		for (int j = 0; j <= malla_x; j++){
			k = i * (malla_x + 1) + j;
			if (es_frontera(i, j, malla_x, malla_y)){
				constantes.matriz[k][0] = T_pared(i, j, malla_x, malla_y, T_N, T_S, T_E, T_O);
			} else {
				constantes.matriz[k][0] = 0.0;
			}
		}
	}

	return constantes;
}

void matriz_2DSimpleTransient(double W, double H, int n, MallaSim *Malla, double T_N, double T_S, double T_E, double T_O, Matriz *previo, double alfa, Matriz *actual, double k, double (*g)(double x, double y, double t), int usar_frontera) {
	
	int M     = Malla->M;
	int N     = Malla->N;
	double Dx = Malla->Dx;
	double Dy = Malla->Dy;
	double Dt = Malla->Dt;
	
	double rx  = alfa * Dt / (Dx * Dx);
	double ry  = alfa * Dt / (Dy * Dy);
	
	double coef_g = alfa * Dt / k;
	
	// She's got the devil in her heart
	// I can't believe that she ever, ever go
	
	// resolver 
	if (n < 0){
		printf("Error! Tiempo negativo!\n");
		actual->datos = NULL;
		actual->matriz = NULL;
		}
	
	if (n == 0) {
		for(int i = 0; i <= M; i++)
			copiar_matriz(previo, actual);
		}
	
	// En el programa, j es x, i es y
	if (n >= 1) {
		for (int i = 0; i < actual->fil; i++){
			for(int j = 0; j < actual->col; j++) {
				if (es_frontera(i, j, M, N)) {
					if (usar_frontera > 0) {
						actual->matriz[i][j] = T_pared(i, j, M, N, T_N, T_S, T_E, T_O);
					} else {
						// No se actualiza si la frontera es constante en el tiempo
						actual->matriz[i][j] = previo->matriz[i][j]; 
					}
				} else {
					actual->matriz[i][j] = previo->matriz[i][j] +
										   rx * (previo->matriz[i][j-1] - 2.0 * previo->matriz[i][j] + previo->matriz[i][j+1]) + 
										   ry * (previo->matriz[i-1][j] - 2.0 * previo->matriz[i][j] + previo->matriz[i+1][j]) + coef_g * g(i*Dx, j*Dy, n*Dt);
					if(actual->matriz[i][j] < 0.0)	actual->matriz[i][j] = 0.0;
				}
			}
		}
	}
	//Inicializar matriz
	
}

