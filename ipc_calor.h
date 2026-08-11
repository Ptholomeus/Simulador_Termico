#ifndef ESTRUCTURA_IPC_CALOR
#define ESTRUCTURA_IPC_CALOR

#pragma pack(push, 1)
typedef struct {
	
	double TB;
	double TSP;
	double kp, ki, kd;
	
} ComandoControlTransitorio2D;
#pragma pack(pop)

#endif
