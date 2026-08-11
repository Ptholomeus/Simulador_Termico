void resolver();
void resolver_2();
void resolver_3();
void resolver_4();
void resolver_4_grafico();
void resolver_5();
void resolver_6_grafico();
void resolver_7_PID();
void resolver_8_grafico_PID_2D();
void frontera_2da_especie(Matriz *previo, Matriz *actual, Matriz *K_map,
						  double (*g)(double, double, double), 
						  double Dx, double Dy, double Dt, 
						  double C, double q, double t);