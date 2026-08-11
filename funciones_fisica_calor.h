typedef struct {
	
	double W, H;
	double alfa, k;
	
} Chapa2D;

typedef struct {
	
	int M, N;
	double Dx, Dy, Dt;
	
} MallaSim;

typedef struct {
	
	double N, S, E, O;
	
} CondicionesFronteras;

Matriz crear_matriz_coefs_1Dsteadystate (double L, int M);
Matriz crear_matriz_coefs_1Dsteadystate_2 (double L, int M);
Matriz crear_matriz_const_1Dsteadystate (int M);
Matriz crear_matriz_const_1Dsteadystate_2 (int M, double T_0, double T_M);
void matriz_1DSimpleTransient(double L, int n, double Dt, int M, int N, double T_0, double T_M, Matriz *acumulado, double alfa, Matriz *solucion_fila);
int es_frontera(int i, int j, double malla_x, double malla_y);
Matriz crear_coefs_multiDimSteadyState_2D(double x_len, double y_len, int malla_x, int malla_y);
Matriz crear_const_multiDimSteadyState_2D(int malla_x, int malla_y, double T_N, double T_S, double T_E, double T_O);
double T_pared(int i, int j, double malla_x, double malla_y, double T_N, double T_S, double T_E, double T_O);
void matriz_2DSimpleTransient(double W, double H, int n, MallaSim *Malla, double T_N, double T_S, double T_E, double T_O, Matriz *previo, double alfa, Matriz *actual, double k, double (*g)(double x, double y, double t), int usar_frontera);
