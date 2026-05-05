# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Run

```sh
make          # Compila todo -> programa.exe
make clean    # Limpia .exe, .o, .d
```

Compilación manual (sin make):
```sh
gcc main.c funciones_matrices.c funciones_fisica_calor.c testear.c solucionar.c exportar_csv.c funciones_vectores.c gnuplot_i.c -o programa.exe
```

Flags usados: `-O3 -Wall -Wpedantic -Wextra -Wconversion`

Programas auxiliares (compilar por separado):
```sh
gcc ipc_transitorio_2d.c -o control_T_A.exe    # Controlador con PID via named pipe
gcc control_TA.c -o control_TA.exe              # Controlador simple via named pipe
```

## Arquitectura

El programa es una simulación numérica de transferencia de calor por diferencias finitas. `main.c` presenta un menú interactivo que despacha a las funciones `resolver()` … `resolver_7_PID()` definidas en `solucionar.c`.

### Módulos principales

| Archivo | Rol |
|---|---|
| `funciones_matrices.c/.h` | Álgebra lineal: creación de matrices, descomposición PLU, eliminación gaussiana, sustitución hacia adelante/atrás, resolución de sistemas `Ax=b` |
| `funciones_vectores.c/.h` | Tipo `Vector_N` (vector dinámico de doubles), usado para registrar permutaciones de filas en PLU |
| `funciones_fisica_calor.c/.h` | Construcción de matrices de coeficientes y constantes para calor 1D estacionario, 1D transitorio, 2D estacionario y 2D transitorio. Define `Chapa2D`, `MallaSim`, `CondicionesFronteras` |
| `solucionar.c/.h` | Flujo de cada problema: discretización → ensamble → resolución → exportación/graficado. `resolver_7_PID()` incluye un controlador PID en tiempo real |
| `testear.c/.h` | Pruebas de operaciones matriciales (`test_matrices_0/1/2`) |
| `exportar_csv.c/.h` | Exporta `Matriz` a archivos CSV |
| `gnuplot_i.c/.h` | Interfaz con gnuplot vía pipes. Biblioteca de terceros con extensiones propias (`gnuplot_plot_coords_smooth`, `gnuplot_termografia_2D`, `gnuplot_dual_view_2D`, etc.) |
| `ipc_calor.h` | Struct compartida `ComandoControlTransitorio2D` para comunicación entre procesos |
| `ipc_transitorio_2d.c` | Programa independiente: envía parámetros PID y condiciones de frontera al simulador 2D transitorio vía Windows named pipe `\\.\pipe\ControlCalor` |

### Flujo de resolución de un sistema lineal

1. `crear_matriz_coefs_*()` → matriz de coeficientes A  
2. `crear_matriz_const_*()` → vector de constantes b  
3. `crear_solucion_sistema(&x, &A, &b)` → internamente hace PLU (`obtener_P_L_U`) + sustitución (`crear_fwd_subs` / `crear_bwd_subs`)

### IPC para el problema transitorio 2D con PID (opción 11)

`resolver_7_PID()` abre un Windows named pipe como servidor. El programa `ipc_transitorio_2d.exe` (cliente) envía structs `ComandoControlTransitorio2D` en caliente para modificar kp/ki/kd y condiciones de frontera sin reiniciar la simulación.

## Plataforma

El código usa `windows.h` (named pipes, `Sleep`, `system("cls")`). Solo compila en Windows. Gnuplot debe estar instalado y en el PATH para las opciones de graficado.
