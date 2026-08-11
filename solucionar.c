#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include <string.h>
#include <math.h>
#include "funciones_matrices.h"
#include "funciones_fisica_calor.h"
#include "testear.h"
#include "solucionar.h"
#include "exportar_csv.h"
#include "gnuplot_i.h"
#include "ipc_calor.h"

double cero_C_en_K = 273.15;

typedef struct {
	
	double kp, ki, kd;
	double dt;
	double acumulado;
	double e_prev;
	
} PID_parametros;

double k_func(double x, double y);
double g(double x, double y, double t);
void testing_Dt(double *Dt, double test_Dt);
double calcular_pid(PID_parametros *pid, double setpoint, double val_actual, double out_min, double out_max);



void resolver() {
		
	double L = 5.0;
	int M;
	double k = 1.0;
	
	printf("Introducir tamano de discretizacion. \nM = ");
	scanf("%d", &M);
	
	// Defino matrices
	Matriz solucion;
	Matriz coeficientes;
	Matriz inv_coefs;
	Matriz constantes;
	
	Matriz test;
	Vector_N vector;
	
	Matriz temp_1;
	
	// crear matriz de soluciones	
		// 13.12.25 Despues de cambiar la definicion de mult escalar matriz, no necesito crearla.
	//solucion = crear_matriz(M+1, 1); // vector columna
	
	// crear matriz de coeficientes
	coeficientes = crear_matriz_coefs_1Dsteadystate(L, M);
	
	/*
	// conseguir inversa de la matriz de coeficientes
	inv_coefs = crear_inversa(&coeficientes);
	*/
	
	// crear matriz constante
	constantes = crear_matriz_const_1Dsteadystate(M);
	
	// solucionar
	//imprimir_matriz(coeficientes);
	//imprimir_matriz(constantes);
	

	/*
	//imprimir y probar la inversa
	test = copiar_matriz(coeficientes);
	
	printf("Triangular coefs: \n");
	triangular_matriz(&test, &vector);
	imprimir_matriz(test);
	
	printf("Inversa\n");
	imprimir_matriz(inv_coefs);
	
	printf("Inversa de coeficientes: \n");
	imprimir_matriz(inv_coefs);
	
	printf("Producto de coeficientes e inversa: \n");
	test = crear_mult_matriz(inv_coefs, coeficientes);
	imprimir_matriz(test);
	*/
	
	/*
	// crear matriz temporal, que almacena (1/k)*A^-1
	// 13.12.25 Acá también podría haber problema, por eso cambio la def de esc mult matr
	temp_1 = crear_escalar_mult_matriz( (1/k) , inv_coefs );
	*/
	
	/*
	// solucion es (1/k)*A^-1*b
	solucion = crear_mult_matriz(temp_1, constantes);
	
	*/
	
	temp_1 = crear_escalar_mult_matriz( (1/k) , constantes);
	crear_solucion_sistema(&solucion, &coeficientes, &temp_1);
	
	// solucion horizontal
	Matriz sol_h;
	sol_h = crear_traspuesta(solucion);
	
	//imprimir_matriz(sol_h);
	
	exportar_solucion(solucion, "texto.csv");
	printf("\nListo!");
	
	getchar();
	getchar();
	
	// destruirrr
	destruir_matriz(&solucion);
	destruir_matriz(&coeficientes);
	destruir_matriz(&constantes);
	// destruir_matriz(&inv_coefs);
	destruir_matriz(&temp_1);
	destruir_matriz(&sol_h);
	// destruir_matriz(&test);
	// free(vector.vector);
	
	return;
	
}

void resolver_2() {
		
	double L = 5.0;
	int M;
	double k = 1.0;
	
	printf("Introducir tamano de discretizacion. \nM = ");
	scanf("%d", &M);
	
	// Defino matrices
	Matriz solucion;
	Matriz coeficientes;
	Matriz inv_coefs;
	Matriz constantes;
	
	Matriz test;
	Vector_N vector;
	
	Matriz temp_1;
	
	// crear matriz de soluciones	
		// 13.12.25 Despues de cambiar la definicion de mult escalar matriz, no necesito crearla.
	//solucion = crear_matriz(M+1, 1); // vector columna
	
	// crear matriz de coeficientes
	coeficientes = crear_matriz_coefs_1Dsteadystate_2(L, M);
	
	
	// crear matriz constante
	constantes = crear_matriz_const_1Dsteadystate_2(M, 100.0, 300.0);
	
	// solucionar
	//imprimir_matriz(coeficientes);
	//imprimir_matriz(constantes);
	
	
	temp_1 = crear_escalar_mult_matriz( (1/k) , constantes);
	crear_solucion_sistema(&solucion, &coeficientes, &temp_1);
	
	// solucion horizontal
	Matriz sol_h;
	sol_h = crear_traspuesta(solucion);
	
	//imprimir_matriz(sol_h);
	
	exportar_solucion(solucion, "texto.csv");
	printf("\nListo!");
	
	getchar();
	getchar();
	
	// destruirrr
	destruir_matriz(&solucion);
	destruir_matriz(&coeficientes);
	destruir_matriz(&constantes);
	// destruir_matriz(&inv_coefs);
	destruir_matriz(&temp_1);
	destruir_matriz(&sol_h);
	// destruir_matriz(&test);
	// free(vector.vector);
	
	return;
	
}

// Para testear transitorio
void resolver_3() {
		
	double L = 5.0;
	double alfa = 1.0;
	double Dt;
	int M;
	double k = 1.0;
	double T_A = 100.0;
	double T_B = 200.0;
	int t_idx;
	int N; // t_max
	
	printf("Introducir tamano de discretizacion. \nM = ");
	scanf("%d", &M);
	printf("Introducir tiempo maximo. \nN = ");
	scanf("%d", &N);
	printf("Introducir dt. \ndt = ");
	scanf("%lf", &Dt);
	
	double Dx = L / (double) M;
	double r_temp = alfa * Dt / (Dx * Dx);

	if(r_temp > 0.5){
		printf("Dt no es util, r vale %lf\n", r_temp);
		
		Dt = 0.5*(Dx * Dx)/alfa - 0.01;
		
		printf("Se usa valor para dt por defecto: dt = %lf\n", Dt);
	}
	
	Matriz acumulado;
	acumulado = crear_matriz(N+1, M+1);
	
	// crear, inicializar y acumular la matriz de condiciones iniciales
	Matriz condicion_inicial;
	condicion_inicial = crear_matriz(1, M+1);
	condicion_inicial.matriz[0][0] = T_A;
	for(int i = 1; i < M; i++)
		condicion_inicial.matriz[0][i] = 0.0;
	condicion_inicial.matriz[0][M] = T_B;
	
	// Acumulador
	acumular_soluciones(&acumulado, &condicion_inicial, 0);
	
	Matriz temp_1;
	temp_1 = crear_matriz(1, M+1);
	
	int n = 1;
		while(true){
			
			matriz_1DSimpleTransient(L, n, Dt, M, N, T_A, T_B, &acumulado, alfa, &temp_1);
			acumular_soluciones(&acumulado, &temp_1, n);
			
			n++;
			if(n == N+1)
				break;
		}
	
	exportar_matriz(acumulado, "transitorio.csv");
	printf("Exportado!\n");
	
	// solucion
	
	getchar();
	getchar();
	
	// destruirrr
	destruir_matriz(&temp_1);
	destruir_matriz(&condicion_inicial);
	destruir_matriz(&acumulado);
	
	return;
	
}

void resolver_4() {
		
	double L = 5.0;
	double alfa = 1.0;
	double Dt;
	int M;
	double k = 1.0;
	double T_A = 100.0;
	double T_B = 200.0;
	int t_idx;
	int N; // t_max
	
	printf("Introducir tamano de discretizacion. \nM = ");
	scanf("%d", &M);
	printf("Introducir tiempo maximo. \nN = ");
	scanf("%d", &N);
	printf("Introducir dt. \ndt = ");
	scanf("%lf", &Dt);
	
	double Dx = L / (double) M;
	double r_temp = alfa * Dt / (Dx * Dx);

	if(r_temp > 0.5){
		
		printf("Dt no es util, r vale %lf\n", r_temp);
		
		Dt = 0.5*(Dx * Dx)/alfa - 0.01;
		
		printf("Se usa valor para dt por defecto: dt = %lf\n", Dt);
	}
	
	// Tiene dos filas, en la 1 la solucion anterior, en la 2 la actual
	// Acumulado con previo y actual
	Matriz solucion_previa;
	solucion_previa = crear_matriz(1, M+1);
	
	// crear, inicializar y acumular la matriz de condiciones iniciales
	Matriz condicion_inicial;					// Declarar matriz 
	condicion_inicial = crear_matriz(1, M+1);   // Crear matriz 
	condicion_inicial.matriz[0][0] = T_A;       // Nodo inicial
	for(int i = 1; i < M; i++)					// Nodos intermedio
		condicion_inicial.matriz[0][i] = 0.0;
	condicion_inicial.matriz[0][M] = T_B;	    // Nodo final 
	
	// No acumular, sino, conservar el previo
	copiar_matriz(&solucion_previa, &condicion_inicial);
	
	Matriz temp_1;
	temp_1 = crear_matriz(1, M+1);
	
	double tiempo = 0;
	int n = 1;
	if(N == -1){
		getchar();
		printf("Presiona tecla cualquiera!\n");
		getchar();
	}
		while(true || N == -1){
			
			if(N == -1){
				printf("Modo infinito!\n\n");
				printf("Tiempo transcurrido: %6.2lf segundos \n\n", tiempo += Dt);
			}
			matriz_1DSimpleTransient(L, 1, Dt, M, N, T_A, T_B, &solucion_previa, alfa, &temp_1);
			acumular_soluciones(&solucion_previa, &temp_1, 0);
			imprimir_matriz(solucion_previa);
			system("cls");
			
			n++;
			if(n == N+1)
				break;
		}
	
	printf("Terminado!\n");
	
	// solucion
	
	getchar();
	getchar();
	
	// destruirrr
	destruir_matriz(&temp_1);
	destruir_matriz(&condicion_inicial);
	destruir_matriz(&solucion_previa);
	
	return;
	
}

void resolver_4_grafico() {
		
	double L = 5.0;
	double alfa = 1.0;
	double Dt;
	int M;
	double k = 1.0;
	double T_A = 100.0;
	double T_B = 200.0;
	int t_idx;
	int N; // t_max
	
	////// probar pipes a la gemini
	HANDLE hPipe;
	char buffer[128];
	DWORD dwRead;
	
	int contador_ipc = 0;
	int frecuencia_ipc = 10; // Revisar solo cada x iteraciones
	
	hPipe = CreateNamedPipe(
		"\\\\.\\pipe\\ControlCalor",
		PIPE_ACCESS_INBOUND, //solo entrada
		PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_NOWAIT, // no bloqueante
		1,    // 1 instancia
		1024, // buffer salida
		1024, // buffer entrada
		0,    // timeout
		NULL  // seguridad por defecto
	);
	
	if (hPipe == INVALID_HANDLE_VALUE) {
        printf("Error al crear la tuberia: %ld\n", GetLastError());
    }	
	
	///////////////////////////////
	
	// Configurar gnuplot_i
	gnuplot_ctrl *handle = gnuplot_init();
	gnuplot_setstyle(handle, "set term windows noraise");
	gnuplot_setstyle(handle, "lines");
	
	gnuplot_cmd(handle, "set xrange [%lf:%lf]", - L * 0.1, L * 1.1);
	//gnuplot_cmd(handle, "set yrange [%lf:%lf]", -T_B * 0.1, T_B * 1.1);
	gnuplot_cmd(handle, "set yrange [%lf:%lf]", 0.0, 1000.0);
	gnuplot_cmd(handle, "set grid");
	gnuplot_cmd(handle, "set mxtics 5");
	gnuplot_cmd(handle, "set mxtics 5");
	gnuplot_cmd(handle, "set grid mxtics mytics");
	
	printf("Introducir tamano de discretizacion. \nM = ");
	scanf("%d", &M);
	printf("Introducir tiempo maximo. \nN = ");
	scanf("%d", &N);
	printf("Introducir dt. \ndt = ");
	scanf("%lf", &Dt);
	
	double Dx = L / (double) M;
	double r_temp = alfa * Dt / (Dx * Dx);

	testing_Dt(&Dt, Dx*Dx*0.5 / alfa);
		
	// numero de iteraciones:
	int freq = 1 / Dt;
	if (freq == 0) freq = 2;
	
	printf("Dt usado: %lf\n", Dt);
	printf("freq usado: %i\n", freq);
	getchar();
	
	// Tiene dos filas, en la 1 la solucion anterior, en la 2 la actual
	// Acumulado con previo y actual
	Matriz solucion_previa;
	solucion_previa = crear_matriz(1, M+1);
	
	// crear, inicializar y acumular la matriz de condiciones iniciales
	Matriz condicion_inicial;					// Declarar matriz 
	condicion_inicial = crear_matriz(1, M+1);   // Crear matriz 
	condicion_inicial.matriz[0][0] = T_A;       // Nodo inicial
	for(int i = 1; i < M; i++)					// Nodos intermedio
		condicion_inicial.matriz[0][i] = 0.0;
	condicion_inicial.matriz[0][M] = T_B;	    // Nodo final 
	
	// No acumular, sino, conservar el previo
	copiar_matriz(&solucion_previa, &condicion_inicial);
	
	Matriz temp_1;
	temp_1 = crear_matriz(1, M+1);
	
	double x_coord[M+1];
	for (int i = 0; i < M+1; i++)
		x_coord[i] = i * Dx;
	
	double tiempo = 0;
	int n = 1;
	if(N == -1){
		getchar();
		printf("Presiona tecla cualquiera!\n");
		getchar();
	}
	
	while((true || N == -1) && T_A >= 0.0){
		
		////// probar pipes a la gemini
		ConnectNamedPipe(hPipe, NULL);
		
		DWORD bytesAvail = 0;
		if (contador_ipc++ >= frecuencia_ipc) {
			contador_ipc = 0;
			ConnectNamedPipe(hPipe, NULL); 
			DWORD bytesAvail = 0;
			  
			if (PeekNamedPipe(hPipe, NULL, 0, NULL, &bytesAvail, NULL) && bytesAvail > 0){
				if(ReadFile(hPipe, buffer, sizeof(buffer)-1, &dwRead, NULL)){
					buffer[dwRead] = '\0';
					T_A = atof(buffer); // Actualiza T_A
					printf("\a"); //avisa el cambio
				}
				DisconnectNamedPipe(hPipe);
			}
		}
		///////////////////////////////
		
		if(N == -1){
			printf("Modo infinito!\n\n");
			printf("Tiempo transcurrido: %6.2lf segundos \n\n", tiempo += Dt);
		}
		
		matriz_1DSimpleTransient(L, 1, Dt, M, N, T_A, T_B, &solucion_previa, alfa, &temp_1);
		acumular_soluciones(&solucion_previa, &temp_1, 0);
		
		if ( n % 5 == 1){
			// imprimir_matriz(solucion_previa);
			gnuplot_resetplot(handle);
			gnuplot_plot_coordinates(handle, x_coord, solucion_previa.matriz[0], M+1, "Calor");
		}
		
		// system("cls");
		n++;
		if(n == N+1)
			break;
	}
	////// probar pipes a la gemini
	CloseHandle(hPipe);
	///////////////////////////////

	printf("Terminado!\n");
	
	// solucion
	
	getchar();
	getchar();
	
	// destruirrr
	destruir_matriz(&temp_1);
	destruir_matriz(&condicion_inicial);
	destruir_matriz(&solucion_previa);
	
	getchar();
	gnuplot_close (handle);
	
	return;
	
}

// Solver de periodo estacionario de una placa 2D con T fronteras constante
void resolver_5() {
		
	double H = 5.0;
	double W = 6.0;
	
	int M;
	int N;
	
	double T_N = 400.0;
	double T_S = 300.0;
	double T_E = 200.0;
	double T_O = 100.0;
	
	double k = 1.0;
	
	printf("Introducir tamano de malla en X. \nM = ");
	scanf("%d", &M);
	printf("Introducir tamano de malla en Y. \nN = ");
	scanf("%d", &N);
	
	gnuplot_ctrl *handle = gnuplot_init();
	
	// Defino matrices
	Matriz solucion;
	Matriz coeficientes;
	Matriz constantes;
	
	Matriz temp_1;

	// crear matriz de coeficientes
	coeficientes = crear_coefs_multiDimSteadyState_2D(H, W, M, N);

	// crear matriz constante
	constantes = crear_const_multiDimSteadyState_2D(M, N, T_N, T_S, T_E, T_O);
	
	// solucionar
	// imprimir_matriz(coeficientes);
	// imprimir_matriz(constantes);

	temp_1 = crear_escalar_mult_matriz( (1/k) , constantes);
	crear_solucion_sistema(&solucion, &coeficientes, &temp_1);
	
	// solucion horizontal
	Matriz sol_h;
	sol_h = crear_traspuesta(solucion);
	
	// Recrear matriz
	Matriz replegada;
	replegada = crear_matriz_desde_fila(N+1, M+1, &sol_h);
	
	//imprimir_matriz(sol_h);
	
	exportar_matriz(replegada, "2D_estacionario.txt");
	printf("\nListo!");
	
	////////////////////////////////////////////////////////////////////////////////
	// 1. Asignar un bloque de memoria estrictamente contiguo
    double *grid_plana = (double *)malloc((N+1) * (M+1) * sizeof(double));
    
    // 2. Extraer los datos reales iterando sobre los punteros seguros de la matriz
    for (int i = 0; i < N + 1; i++) {
        for (int j = 0; j < M + 1; j++) {
            grid_plana[i * (M + 1) + j] = replegada.matriz[i][j];
        }
    }
	////////////////////////////////////////////////////////////////////////////////
	
	gnuplot_splot_grid(handle, grid_plana, N+1, M+1, "A");
	
	imprimir_matriz(replegada);
	
	getchar();
	getchar();
	
	// destruirrr
	destruir_matriz(&solucion);
	destruir_matriz(&coeficientes);
	destruir_matriz(&constantes);
	destruir_matriz(&temp_1);
	destruir_matriz(&sol_h);
	destruir_matriz(&replegada);
	
	gnuplot_close (handle);
	free(grid_plana);
	
	return;
}


double g(double x, double y, double t){
	
	double z;
	
	// z = 150.0*(sin(x+y+t))*(1.5*t+3.0*y+2.0*x)/sqrt(t*t+x*x+y*y);
	
	z = 1.0 + 0.1 * x + 0.02 * y;
	
	// z = 0;
	
	return z;
}


void testing_Dt(double *Dt, double test_Dt){
	
	if( *Dt > test_Dt){
	
		printf("Dt = %lf no es util, debe ser menor a %lf\n", *Dt, test_Dt);
		
		// Algoritmo kk
		int temp_time;
		temp_time = (int)test_Dt;
		int test_time_idx = 0;
		
		while(true){
			if (temp_time != 0){
				break;
			} else {
				test_time_idx++;
				temp_time = (int) (test_Dt * pow(10.0, (double) test_time_idx));
			}
		}
		
		*Dt = test_Dt - 1.0 / pow(10.0, (double) (test_time_idx + 1)); // No me gusta mucho!
		
		printf("Se usa valor para dt por defecto: dt = %lf\n", *Dt);
		
		getchar();
		
	}
}

void resolver_6_grafico() {
		
	double W = 10.0;
	double H = 12.0;
	
	double k    = 0.05;
	double alfa = 0.05;

	int M;
	int N;
	
	double Dx;
	double Dy;
	double Dt;
	
	MallaSim Malla;
	
	double T_N = 400.0;
	double T_S = 300.0;
	double T_E = 450.0;
	double T_O = 250.0;
	
	int t_idx;
	int Tmax;
	
	// Configurar gnuplot_i
	gnuplot_ctrl *handle = gnuplot_init();
	gnuplot_setstyle(handle, "set term windows noraise");
	gnuplot_setstyle(handle, "lines");
	
		
	printf("Introducir tamano de malla en X. \nM = ");
	scanf("%d", &M);
	printf("Introducir tamano de malla en Y. \nN = ");
	scanf("%d", &N);
	
	Dx = W / (double) M;
	Dy = H / (double) N; // It's only love and that is all
	
	printf("Introducir tiempo maximo. \nN = ");
	scanf("%d", &Tmax);
	
	printf("Introducir dt. \ndt = ");
	scanf("%lf", &Dt);
	
	double test_Dt = (0.5 / alfa) * (Dx * Dx * Dy * Dy) / (Dx * Dx + Dy * Dy);

	testing_Dt(&Dt, test_Dt);
		
	// numero de iteraciones:
	int freq = 1 / Dt;
	if (freq == 0) freq = 2;
	
	printf("Dt usado: %lf\n", Dt);
	printf("freq usado: %i\n", freq);
	getchar();
	
	// Tiene dos filas, en la 1 la solucion anterior, en la 2 la actual
	// Acumulado con previo y actual
	Matriz solucion_previa;
	solucion_previa = crear_matriz(N+1, M+1);
	
	// crear, inicializar y acumular la matriz de condiciones iniciales
	Matriz condicion_inicial;					// Declarar matriz 
	condicion_inicial = crear_matriz(N+1, M+1);   // Crear matriz
	for(int i = 0; i <= M; i++){
		for (int j = 0; j <= N; j++){
			condicion_inicial.matriz[i][j] = 200.0;
		}
	}
	/////////////////////////////////////////////////////////////////
	llenar_fil(&condicion_inicial, 0, T_N);
	llenar_fil(&condicion_inicial, N, T_S);
	llenar_col(&condicion_inicial, 0, T_O);
	llenar_col(&condicion_inicial, M, T_E);
	condicion_inicial.matriz[0][0] = 0.5*(T_N+T_O);
	condicion_inicial.matriz[0][M] = 0.5*(T_N+T_E);
	condicion_inicial.matriz[N][0] = 0.5*(T_S+T_O);
	condicion_inicial.matriz[N][M] = 0.5*(T_S+T_E);
	/////////////////////////////////////////////////////////////////
	
	// No acumular, sino, conservar el previo
	copiar_matriz(&solucion_previa, &condicion_inicial);
	
	Matriz temp_1;
	temp_1 = crear_matriz(N+1, M+1);
	
	double x_coord[M+1];
	for (int i = 0; i < M+1; i++)
		x_coord[i] = i * Dx;
	double y_coord[N+1];
	for (int i = 0; i < N+1; i++)
		y_coord[i] = i * Dy;
	
	double tiempo = 0;
	int n = 1;
	if(N == -1){
		printf("Presiona tecla cualquiera!\n");
		getchar();
	}
	
	Malla.M  = M;
	Malla.N  = N;
	Malla.Dx = Dx;
	Malla.Dy = Dy;
	Malla.Dt = Dt;
	
	if(Tmax == -1)
		printf("Modo infinito!\n\n");
	
	while(true || Tmax == -1){
		
		if(Tmax == -1){
			printf("Tiempo transcurrido: %6.2lf segundos \n\n", tiempo += Dt);
		}
		
		matriz_2DSimpleTransient(W, H, n, &Malla, T_N, T_S, T_E, T_O, &solucion_previa, alfa, &temp_1, k, g, 1);
		copiar_matriz(&solucion_previa, &temp_1);
		
		// condicionnnnn
		if ( 1 == 1){
			gnuplot_resetplot(handle);
			// gnuplot_splot_grid(handle, temp_1.datos, N+1, M+1, "A");
			//gnuplot_splot_grid_fixed(handle, temp_1.datos, N+1, M+1, Dx, Dy, -10.0, 600.0, "A");
			gnuplot_termografia_2D(handle, temp_1.datos, N+1, M+1, Dx, Dy, -200.0, 2000.0, "AAAA");
			// gnuplot_dual_view_2D(handle, temp_1.datos, N+1, M+1, Dx, Dy, 0.0, 1200.0, "uuuuh");
			Sleep(100);
		}
		
		//system("cls");
		n++;
		//if(n == Tmax + 1 || n == 5000)
			// break;
	}

	printf("Terminado!\n");
	
	// solucion
	
	getchar();
	getchar();
	
	// destruirrr
	destruir_matriz(&temp_1);
	destruir_matriz(&condicion_inicial);
	destruir_matriz(&solucion_previa);
	
	getchar();
	gnuplot_close (handle);
	
	system("pause");
	return;
	
}

double calcular_pid(PID_parametros *pid, double setpoint, double val_actual, double out_min, double out_max){
	
	double P, I, D;
	double error;
	
	error = setpoint - val_actual;
	P = pid->kp * error;
	D = pid->kd * (error - pid->e_prev) / pid->dt;
	
	//pid->acumulado += error * pid->dt;
	//I = pid->ki * pid->acumulado;
	
	double I_tentativa = pid->ki * (pid->acumulado + error * pid->dt);
    double out_tentativa = P + I_tentativa + D;
	
	// pid->e_prev = error;
	
	//double out = P+I+D;
	
	/*
	if (out > out_max) {
		out = out_max;
		// Back-calculation: deshacer la integracion de este paso
		pid->acumulado -= error * pid->dt; 
	} else if (out < out_min) {
		out = out_min;
		pid->acumulado -= error * pid->dt; 
	}
	*/ 
	
	// Clamping (Anti-Windup)
    if (out_tentativa > out_max) {
        // Solo acumulamos si el error ayuda a desaturar (error < 0)
        if (error < 0) pid->acumulado += error * pid->dt;
        return out_max;
    } else if (out_tentativa < out_min) {
        // Solo acumulamos si el error ayuda a desaturar (error > 0)
        if (error > 0) pid->acumulado += error * pid->dt;
        return out_min;
    } else {
        // Operación lineal normal
        pid->acumulado += error * pid->dt;
        pid->e_prev = error;
        return out_tentativa;
    }	
	
	//return out;
	
}

/*
	Conozco:-k*D_x·T = q_0 y k*D_x·T = q_L 
	[q_i] = W/m^2
	[T]   = K
	[k]   = W/mK
*/

void resolver_7_PID() {
		
	double L = 5.0;
	double alfa = 1.0;
	double Dt;
	int M;
	double k = 1.0;
	double T_A = 25.0 + cero_C_en_K;
	double q_0 = 100.0;
	double T_B = 35.0 + cero_C_en_K;
	int t_idx;
	int N; // t_max
	
	double q0_min = 0.0;
	double q0_max = 50000.0;
	
	PID_parametros pid;
	
	////// probar pipes a la gemini
	HANDLE hPipe;
	char buffer[128];
	DWORD dwRead;
	
	int contador_ipc = 0;
	int frecuencia_ipc = 10; // Revisar solo cada x iteraciones
	
	hPipe = CreateNamedPipe(
		"\\\\.\\pipe\\ControlCalor",
		PIPE_ACCESS_INBOUND, //solo entrada
		PIPE_TYPE_MESSAGE | PIPE_READMODE_MESSAGE | PIPE_NOWAIT, // no bloqueante
		1,    // 1 instancia
		sizeof(ComandoControlTransitorio2D), // buffer salida
		sizeof(ComandoControlTransitorio2D), // buffer entrada
		0,    // timeout
		NULL  // seguridad por defecto
	);
	
	if (hPipe == INVALID_HANDLE_VALUE) {
        printf("Error al crear la tuberia: %ld\n", GetLastError());
    }	
	
	///////////////////////////////
	
	// Configurar gnuplot_i
	gnuplot_ctrl *handle = gnuplot_init();
	gnuplot_setstyle(handle, "set term windows noraise");
	gnuplot_setstyle(handle, "lines");
	
	gnuplot_cmd(handle, "set xrange [%lf:%lf]", - L * 0.1, L * 1.1);
	//gnuplot_cmd(handle, "set yrange [%lf:%lf]", -T_B * 0.1, T_B * 1.1);
	gnuplot_cmd(handle, "set yrange [%lf:%lf]", 270.0, 500.0);
	gnuplot_cmd(handle, "set grid");
	gnuplot_cmd(handle, "set mxtics 5");
	gnuplot_cmd(handle, "set mxtics 5");
	gnuplot_cmd(handle, "set grid mxtics mytics");
	
	printf("Introducir tamano de discretizacion. \nM = ");
	scanf("%d", &M);
	printf("Introducir tiempo maximo. \nN = ");
	scanf("%d", &N);
	printf("Introducir dt. \ndt = ");
	scanf("%lf", &Dt);
	
	double Dx = L / (double) M;
	
	testing_Dt(&Dt, Dx*Dx*0.5 / alfa);
	double r_temp = alfa * Dt / (Dx * Dx);

	
	// numero de iteraciones:
	int freq = 1 / Dt;
	if (freq == 0) freq = 2;
	
	printf("Dt usado: %lf\n", Dt);
	printf("freq usado: %i\n", freq);
	getchar();
	
	// Tiene dos filas, en la 1 la solucion anterior, en la 2 la actual
	// Acumulado con previo y actual
	Matriz solucion_previa;
	solucion_previa = crear_matriz(1, M+1);
	
	// crear, inicializar y acumular la matriz de condiciones iniciales
	Matriz condicion_inicial;					// Declarar matriz 
	condicion_inicial = crear_matriz(1, M+1);   // Crear matriz 
	condicion_inicial.matriz[0][0] = T_A;       // Nodo inicial, con temperatura inicial T_A
	for(int i = 1; i < M; i++)					// Nodos intermedio, 25 ºC cond inicial
		condicion_inicial.matriz[0][i] = cero_C_en_K + 25.0;
	condicion_inicial.matriz[0][M] = T_B;	    // Nodo final 
	
	// No acumular, sino, conservar el previo
	copiar_matriz(&solucion_previa, &condicion_inicial);
	
	Matriz temp_1;
	temp_1 = crear_matriz(1, M+1);
	
	double x_coord[M+1];
	for (int i = 0; i < M+1; i++)
		x_coord[i] = i * Dx;
	
	double tiempo = 0;
	int n = 1;
	if(N == -1){
		printf("Presiona tecla cualquiera!\n");
		
		getchar();
		printf("\rModo infinito!\n\n");
	}
	
	pid.dt = Dt;
	pid.acumulado = 0;
	pid.e_prev = 0;
	
	pid.kp = 0.00000000001;
	pid.ki = 0.0;
	pid.kd = 0.0;
	
	double T_SP = T_B;
	
	while((true || N == -1) && T_SP >= 0.0){
		
		////// probar pipes a la gemini
		ConnectNamedPipe(hPipe, NULL);
		
		DWORD bytesAvail = 0;
		if (contador_ipc++ >= frecuencia_ipc) {
			contador_ipc = 0;
			ConnectNamedPipe(hPipe, NULL); 
			DWORD bytesAvail = 0;
			  
			if (PeekNamedPipe(hPipe, 
							  NULL, 
							  0, 
							  NULL, 
							  &bytesAvail, 
							  NULL) && bytesAvail >= sizeof(ComandoControlTransitorio2D)){
				
				ComandoControlTransitorio2D struct_datos;
				
				if(ReadFile(hPipe, &struct_datos, sizeof(ComandoControlTransitorio2D), &dwRead, NULL)){
					if(dwRead == sizeof(ComandoControlTransitorio2D)){
						pid.kp = struct_datos.kp;
						pid.ki = struct_datos.ki;
						pid.kd = struct_datos.kd;
						T_B  = struct_datos.TB  + cero_C_en_K;
						T_SP = struct_datos.TSP + cero_C_en_K;
						printf("Cambio aplicado\n");
					}
				}
				DisconnectNamedPipe(hPipe);
			}
		}
		///////////////////////////////
		
		if(N == -1){
			
			printf("\rkp: %2.2lf ki: %2.2lf kd: %2.2lf T_SP: %2.2lf, T_A: %2.2lf, T_mid: %2.2lf, T_B: %2.2lf q_0: %2.2lf, t: %2.2lf", 
					pid.kp,
					pid.ki,
					pid.kd,
					T_SP - cero_C_en_K, 
					T_A - cero_C_en_K, 
					solucion_previa.matriz[0][(int) (M / 2.0)] - cero_C_en_K, 
					T_B - cero_C_en_K,
					q_0, 
					tiempo += Dt);
			fflush(stdout);
			// system("cls");
		}
		
		q_0 = calcular_pid(&pid, T_SP, solucion_previa.matriz[0][(int) (M / 2.0)], q0_min, q0_max);
		
		T_A = (1-2.0*r_temp)*solucion_previa.matriz[0][0] +
							  2.0*r_temp*solucion_previa.matriz[0][1] +
							  2.0*r_temp*Dx * q_0 / k;
		
		matriz_1DSimpleTransient(L, 1, Dt, M, N, T_A, T_B, &solucion_previa, alfa, &temp_1);
		acumular_soluciones(&solucion_previa, &temp_1, 0);
		
		if ( n % 5 == 1){
			// imprimir_matriz(solucion_previa);
			gnuplot_resetplot(handle);
			gnuplot_plot_coordinates(handle, x_coord, solucion_previa.matriz[0], M+1, "Calor");
			// Sleep(50);
		}
		
		// system("cls");
		n++;
		if(n == N+1)
			break;
	}
	////// probar pipes a la gemini
	CloseHandle(hPipe);
	///////////////////////////////

	printf("Terminado!\n");
	
	// solucion
	
	getchar();
	getchar();
	
	// destruirrr
	destruir_matriz(&temp_1);
	destruir_matriz(&condicion_inicial);
	destruir_matriz(&solucion_previa);
	
	getchar();
	gnuplot_close (handle);
	
	return;
	
}


/*

	Ec de difusion:
		
		rho*C*dT/dt = div(k * grad ( T ) ) + g 
		
		rho*C*dT/dt = d(k*dT/dx)/dx + d(k*dT/dy)/dy + g
	
	Ley de Fourier:
		
		q' = -k*dT/dx 
		
		(Usado en condición de Neumann - 2da especie)
		
	Discretizada la dif.:
		
		A = Dt / ( rho * C ) = alfa * Dt / k
		
		T_prox[i][j] = T[i][j] +
					   A*k(i, j)*( T[i-1][j] - 2 * T[i][j] + T[i+1][j] )/(Dx*Dx) + 
					   A*k(i, j)*( T[i][j-1] - 2 * T[i][j] + T[i][j+1] )/(Dy*Dy) +
					   A*(k(i+1, j)-k(i-1,j))*(T[i+1][j] - T[i-1][j])/(4*Dx*Dx)  +
					   A*(k(i, j+1)-k(i,j-1))*(T[i][j+1] - T[i][j-1])/(4*Dy*Dy)  +
					   A*g(i,j);
		
		Estudiemos la pared izquierda.
		Al evaluar esto en i=0, aparece un término i=-1, luego
		discretizar la ec. de Fourier centrada en i=0
		
		-k(0,j)*( T[1][j] - T[-1][j] )/(2*Dx) = q'
		
		ó
		
		T[-1][j] = T[1][j] + 2 * Dx * q' / k[0][j]
		
		(q' es el que da el PID)
		
		De aquí se despeja T[-1][j] y se reemplaza, lo que da (usar derivada de k en x hacia delante):
		
		T_prox[0][j] = T[0][j] +
					   A*k(0, j)*( T[1][j] + 2 * Dx * q' / k[0][j] - 2 * T[0][j] + T[1][j] )/(Dx*Dx) + 
					   A*k(0, j)*( T[0][j-1] - 2 * T[0][j] + T[0][j+1] )/(Dy*Dy) +
					   A*(k(1, j)-k(0,j))*q'/(k(0, j)*Dx)  +
					   A*(k(0, j+1)-k(0,j-1))*(T[0][j+1] - T[0][j-1])/(4*Dy*Dy)  +
					   A*g(0,j);

*/

// Solo oeste
void frontera_2da_especie(Matriz *previo, Matriz *actual, Matriz *K_map,
						  double (*g)(double, double, double), 
						  double Dx, double Dy, double Dt, 
						  double C, double q, double t){
	
	int N = actual->fil-1;
	double A_const = Dt / C;
		
	for (int i = 1; i < N; i++){
		
		double y = i * Dy;
		
		double k_local = K_map->matriz[i][0];
        double k_der_x = K_map->matriz[i][1];
        double k_sup_y = K_map->matriz[i+1][0];
        double k_inf_y = K_map->matriz[i-1][0];
		
		actual->matriz[i][0] = previo->matriz[i][0] +
            A_const * k_local * (2.0 * previo->matriz[i][1] - 2.0 * previo->matriz[i][0] + (2.0 * Dx * q) / k_local) / (Dx * Dx) +
            A_const * k_local * (previo->matriz[i+1][0] - 2.0 * previo->matriz[i][0] + previo->matriz[i-1][0]) / (Dy * Dy) -
            A_const * ((k_der_x - k_local) / Dx) * (q / k_local) +
            A_const * ((k_sup_y - k_inf_y) / (2.0 * Dy)) * (previo->matriz[i+1][0] - previo->matriz[i-1][0]) / (2.0 * Dy) +
            A_const * g(0.0, y, t);
		
	}
	
	return;
}

double k_func(double x, double y){
	
	
	return 1.0;
}

void resolver_8_grafico_PID_2D() {
	
	/////////////////////////
	// Dimensiones fisicas //
	/////////////////////////
	
	// Chapa
	double W = 0.30;
	double H = 0.20;
	
	// Coeficientes de calor
	double k    = k_func(0.0, 0.0);
	double alfa = 0.000003;
	
	// Espacio
	double Dx;
	double Dy;
	double Dt;
	
	// Temperaturas de frontera (fijas)
	double T_N = 40.0 + cero_C_en_K;
	double T_S = 35.0 + cero_C_en_K;
	double T_E = 45.0 + cero_C_en_K;
	double T_O = 25.0 + cero_C_en_K;
	
	// Control de tiempo
	int t_idx;
	int Tmax;
	
	////////////////////////////
	// Dimensiones simulacion //
	////////////////////////////
	
	// Cantidad de nodos
	int M;
	int N;
	
	double q;
	double q0_min = 0.0;
	double q0_max = 50000.0;
	
	PID_parametros pid;
	
	//////////////////////////
	// Named PIPES de Win32 //
	//////////////////////////
	
	HANDLE hPipe;
	char buffer[128];
	DWORD dwRead;
	
	int contador_ipc = 0;
	int frecuencia_ipc = 10; // Revisar solo cada x iteraciones
	
	hPipe = CreateNamedPipe(
		"\\\\.\\pipe\\ControlCalor",
		PIPE_ACCESS_INBOUND, //solo entrada
		PIPE_TYPE_MESSAGE | PIPE_READMODE_MESSAGE | PIPE_NOWAIT, // no bloqueante
		1,    // 1 instancia
		sizeof(ComandoControlTransitorio2D), // buffer salida
		sizeof(ComandoControlTransitorio2D), // buffer entrada
		0,    // timeout
		NULL  // seguridad por defecto
	);
	
	if (hPipe == INVALID_HANDLE_VALUE) {
        printf("Error al crear la tuberia: %ld\n", GetLastError());
    }	
	
	//////////////////////
	// Seteo de GNUPLOT //
	//////////////////////
	
	// Configurar gnuplot_i
	gnuplot_ctrl *handle = gnuplot_init();
	gnuplot_setstyle(handle, "set term windows noraise");
	gnuplot_setstyle(handle, "lines");
	
	//////////////////////////
	// Sección de prompting //
	//////////////////////////
	
	printf("Introducir tamano de malla en X. \nM = ");
	scanf("%d", &M);
	printf("Introducir tamano de malla en Y. \nN = ");
	scanf("%d", &N);
	
	printf("Introducir tiempo maximo. \nt_max = ");
	scanf("%d", &Tmax);
	
	printf("Introducir dt. \ndt = ");
	scanf("%lf", &Dt);
	
	///////////////////////////////////
	// Establecimiento de parametros //
	//        con los prompts        //
	///////////////////////////////////
	
	Dx = W / (double) M;
	Dy = H / (double) N; // It's only love and that is all
	
	// Testeo de Dt 
	
	double test_Dt = (0.5 / alfa) * (Dx * Dx * Dy * Dy) / (Dx * Dx + Dy * Dy);
	testing_Dt(&Dt, test_Dt);
		
	// COSA RARA CON FRECUENCIAS PARA NO PLOTEAR TODOOOOO
	int freq = 1 / Dt;
	if (freq == 0) freq = 2;
	
	////////////////////////////
	// Print en consola de Dt //
	////////////////////////////
	
	printf("Dt usado: %lf\n", Dt);
	printf("freq usado: %i\n", freq);
	getchar();
	
	
	//////////////////////////////////////////
	// Declaración y definicion de matrices //
	//////////////////////////////////////////
	
	// Tiene dos filas, en la 1 la solucion anterior, en la 2 la actual
	// Acumulado con previo y actual
	Matriz solucion_previa;
	solucion_previa = crear_matriz(N+1, M+1);
	
	// crear, inicializar y acumular la matriz de condiciones iniciales
	Matriz condicion_inicial;					// Declarar matriz 
	condicion_inicial = crear_matriz(N+1, M+1);   // Crear matriz
	for(int i = 0; i <= M; i++){
		for (int j = 0; j <= N; j++){
			condicion_inicial.matriz[i][j] = 200.0;
		}
	}
	
	/////////////////////////////////////////////////////////////////
	///////// Acá se llenan las fronteras de 1ra especie ////////////
	/////////////////////////////////////////////////////////////////
	llenar_fil(&condicion_inicial, 0, T_N);
	llenar_fil(&condicion_inicial, N, T_S);
	llenar_col(&condicion_inicial, 0, T_O);
	llenar_col(&condicion_inicial, M, T_E);
	condicion_inicial.matriz[0][0] = 0.5*(T_N+T_O);
	condicion_inicial.matriz[0][M] = 0.5*(T_N+T_E);
	condicion_inicial.matriz[N][0] = 0.5*(T_S+T_O);
	condicion_inicial.matriz[N][M] = 0.5*(T_S+T_E);
	/////////////////////////////////////////////////////////////////
	
	Matriz temp_1;
	temp_1 = crear_matriz(N+1, M+1);
	
	
	// No acumular, sino, conservar el previo
	copiar_matriz(&solucion_previa, &condicion_inicial);
	copiar_matriz(&temp_1, &condicion_inicial);

	///////////////////////
	//// COSAS  / / / // //
	///////////////////////
	
	double tiempo = 0;   // Tiempo a mostrar en consola
	int n = 1;           // Contador de iteraciones
	
	
	//////////////////////////////////////////
	/////// Estructura para los pipes ////////
	//////////////////////////////////////////
	
	/*
	struct {
	
	int M, N;
	double Dx, Dy, Dt;
	
	} MallaSim;
	
	*/
	
	MallaSim Malla;
	
	Malla.M  = M;
	Malla.N  = N;
	Malla.Dx = Dx;
	Malla.Dy = Dy;
	Malla.Dt = Dt;
	
	/////////////////////////////////
	//////  Matriz K          ///////
	/////////////////////////////////
	
	Matriz matriz_k = crear_matriz(N+1, M+1);
	for(int i = 0; i <= N; i++){
		for(int j = 0; j <= M; j++){
			matriz_k.matriz[i][j] = k_func(j * Dx, i * Dy);
		}
	}
	
	pid.dt = Dt;
	pid.acumulado = 0;
	pid.e_prev = 0;
	
	pid.kp = 1.0;
	pid.ki = 0.0;
	pid.kd = 0.0;
	
	double T_B = cero_C_en_K + 20.0;
	double T_SP = T_B;
	
	/////////////////////////////////
	//////  Modo infinito!!   ///////
	/////////////////////////////////
	
	if(Tmax == -1)
		printf("Modo infinito!\n\n");
	
	/////////////////////////////////
	///// Bucle de solución     /////
	/////////////////////////////////
	t_idx = 0;
	while(true){
		
		if(Tmax == -1){
			printf("Tiempo transcurrido: %6.2lf segundos \n\n", tiempo += Dt);
		}
		
	
		///////////////////////////////
		////// PIPING ! !!!! ! !   ////
		///////////////////////////////
				
		ConnectNamedPipe(hPipe, NULL);
		
		DWORD bytesAvail = 0;
		if (contador_ipc++ >= frecuencia_ipc) {
			contador_ipc = 0;
			ConnectNamedPipe(hPipe, NULL); 
			DWORD bytesAvail = 0;
			  
			if (PeekNamedPipe(hPipe, 
							  NULL, 
							  0, 
							  NULL, 
							  &bytesAvail, 
							  NULL) && bytesAvail >= sizeof(ComandoControlTransitorio2D)){
				
				ComandoControlTransitorio2D struct_datos;
				
				if(ReadFile(hPipe, &struct_datos, sizeof(ComandoControlTransitorio2D), &dwRead, NULL)){
					if(dwRead == sizeof(ComandoControlTransitorio2D)){
						pid.kp = struct_datos.kp;
						pid.ki = struct_datos.ki;
						pid.kd = struct_datos.kd;
						T_B  = struct_datos.TB  + cero_C_en_K;
						T_SP = struct_datos.TSP + cero_C_en_K;
						printf("Cambio aplicado\n");
					}
				}
				DisconnectNamedPipe(hPipe);
			}
		}

		///////////////////////////////
		////// SOLVING! !!!! ! !   ////
		///////////////////////////////
		
		
		q = calcular_pid(&pid, T_SP, solucion_previa.matriz[(int) (N / 2.0)][(int)(M / 2.0)], q0_min, q0_max);
		matriz_2DSimpleTransient(W, H, t_idx, &Malla, T_N, T_S, T_E, T_O, &solucion_previa, alfa, &temp_1, k, g, -1);
		frontera_2da_especie(&solucion_previa, &temp_1, &matriz_k, g, Dx, Dy, Dt, k / (1.0 * alfa), q, t_idx*Dt);
		swap_matrices(&solucion_previa, &temp_1);
		
		///////////////////////////////
		//// PLOT!!!!!!!!!!!    ///////
		///////////////////////////////
		int plot_freq = (int)(0.1 / Dt); 
		if (plot_freq < 1) plot_freq = 1;
		if (t_idx % plot_freq == 0){
			gnuplot_resetplot(handle);
			//gnuplot_splot_grid(handle, temp_1.datos, N+1, M+1, "A");
			gnuplot_splot_grid_fixed(handle, temp_1.datos, N+1, M+1, Dx, Dy, 0.0, 1500.0, "A");
			//gnuplot_termografia_2D(handle, temp_1.datos, N+1, M+1, Dx, Dy, cero_C_en_K, 1000, "AAAA");
			//gnuplot_dual_view_2D(handle, temp_1.datos, N+1, M+1, Dx, Dy, cero_C_en_K, 1200.0, "uuuuh");
			Sleep(10);
		}
		
		t_idx++;
		if (t_idx > Tmax && Tmax != -1) break;
	}

	CloseHandle(hPipe);

	printf("Terminado!\n");
	
	// solucion
	
	getchar();
	getchar();
	
	// destruirrr
	destruir_matriz(&temp_1);
	destruir_matriz(&condicion_inicial);
	destruir_matriz(&solucion_previa);
	destruir_matriz(&matriz_k);
	
	getchar();
	gnuplot_close (handle);
	
	system("pause");
	return;
	
}