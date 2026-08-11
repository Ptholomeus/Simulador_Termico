#include <windows.h>
#include <stdio.h>
#include "ipc_calor.h"

int main() {
    HANDLE hPipe;
    DWORD dwWritten;
    ComandoControlTransitorio2D struct_datos_enviar;
	
	
	
    printf("--- Controlador de Temperatura ---\n");
	
	/*
    while (1) {
        hPipe = CreateFile("\\\\.\\pipe\\ControlCalor", GENERIC_WRITE, 0, NULL, OPEN_EXISTING, 0, NULL);
        
        if (hPipe != INVALID_HANDLE_VALUE) {
            printf("Ingrese nueva T_A: ");
            scanf("%s", buffer);

            WriteFile(hPipe, buffer, strlen(buffer), &dwWritten, NULL);
            CloseHandle(hPipe);
            printf("Enviado!\n");
        } else {
            printf("Esperando a que el simulador se inicie...\n");
            Sleep(2000);
        }
    }
    return 0;
	*/
	
	while (1) {
        // Primero pedimos el dato al usuario
        printf("\nIngrese nueva kp ki kd T_B y T_SP: ");
        if (scanf("%lf %lf %lf %lf %lf",
			&struct_datos_enviar.kp,
			&struct_datos_enviar.ki,
			&struct_datos_enviar.kd,
			&struct_datos_enviar.TB, 
			&struct_datos_enviar.TSP) == 5) {
            
            // Intentamos abrir la tubería solo para este envío
            hPipe = CreateFile("\\\\.\\pipe\\ControlCalor", 
				GENERIC_WRITE, 
				0, 
				NULL, 
				OPEN_EXISTING, 
				0, 
				NULL);
            
            if (hPipe != INVALID_HANDLE_VALUE) {
                if (WriteFile(hPipe, 
							  &struct_datos_enviar, 
							  sizeof(ComandoControlTransitorio2D), 
							  &dwWritten, 
							  NULL)
					) {
                    printf("Enviado correctamente: %lu bites enviados\n", dwWritten);
                }
                CloseHandle(hPipe); // Cerramos para liberar la instancia
            } else {
                printf("Error: No se pudo conectar. ¿Esta corriendo el simulador?\n");
            }
        } else {
			printf("Error en el comando de entrada. \n");
			int c;
			while ((c = getchar()) != '\n' && c != EOF);
		}
    }
    return 0;
}