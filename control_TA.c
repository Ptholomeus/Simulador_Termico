#include <windows.h>
#include <stdio.h>

int main() {
    HANDLE hPipe;
    DWORD dwWritten;
    char buffer[128];
	
	
	
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
        printf("\nIngrese nueva T_A: ");
        if (scanf("%127s", buffer) == 1) {
            
            // Intentamos abrir la tubería solo para este envío
            hPipe = CreateFile("\\\\.\\pipe\\ControlCalor", GENERIC_WRITE, 0, NULL, OPEN_EXISTING, 0, NULL);
            
            if (hPipe != INVALID_HANDLE_VALUE) {
                if (WriteFile(hPipe, buffer, (DWORD)strlen(buffer), &dwWritten, NULL)) {
                    printf("Enviado correctamente: %s\n", buffer);
                }
                CloseHandle(hPipe); // Cerramos para liberar la instancia
            } else {
                printf("Error: No se pudo conectar. ¿Esta corriendo el simulador?\n");
            }
        }
    }
    return 0;
}