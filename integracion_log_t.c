/* INTEGRACION NUMERICA POR EL METODO DE LOS TRAPECIOS
 * ENTRADA: NINGUNA
 * SALIDA: ESTIMACION DE LA INTEGRAL DESDE a HASTA b DE f(x)
 * USANDO EL METODO DE LOS TRAPECIOS CON n TRAPECIOS */
#include <stdio.h>
#include <math.h>
#include <mpi.h>

main(int argc, char** argv) {
    float integral; // RESULTADO DE LA INTEGRAL
    float a = 1.0; // EXTREMO IZQUIERDO
    float b = 100000.0; // EXTREMO DERECHO
    int n = 7927920; // NUMERO DE TRAPECIOS
    float h; // LONGITUD DE LA BASE DEL TRAPECIO
    float x;

    int mi_rango;
    int p;
    int i;
    float a_local;
    float b_local;
    int n_local;
    float integral_local;
    float x_local;
    float trapecio;
    int fuente;
    int dest = 0;
    int etiqueta = 0;

    double inicio;
    double fin;

    MPI_Status  status;
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &mi_rango);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if(n % p != 0)  {
        if(mi_rango == 0) {
            printf("El numero de trapecios debe ser multiplo del numero de\t"
                    "procesos.\n");
        }
        MPI_Finalize();
    }

    h = (b - a) / n;
    a_local = a + mi_rango * h * n / p;
    b_local = a + (mi_rango + 1) * h * n / p;
    n_local = n / p;

    float f(float x); /* FUNCION QUE ESTAMOS INTEGRANDO */

    /* 
    integral_local = (f(a_local) + f(b_local)) / 2.0;
    x_local = a_local;
    for(i = 1; i <= n_local - 1; i++) {
        x_local = x_local + h;
        integral_local = integral_local + f(x_local);
    }
    integral_local = integral_local * h;
    */

    if(mi_rango == 0) {
        inicio = MPI_Wtime();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    integral_local = 0;
    x_local = a_local;
    for(i = 0; i <= n_local - 1; i++) {
        trapecio = (f(x_local) + f(x_local + h)) / 2.0 * h;
        integral_local = integral_local + trapecio;
        x_local = x_local + h;
    }
    
    if(mi_rango == 0) {
        integral = integral_local;
        for(fuente = 1; fuente < p; fuente++) {
            MPI_Recv(&integral_local, 1, MPI_FLOAT, fuente, etiqueta,
                    MPI_COMM_WORLD, &status);
	        integral = integral + integral_local;
        }
    } else {
        MPI_Send(&integral_local, 1, MPI_FLOAT, dest, etiqueta, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(mi_rango == 0) {
        fin = MPI_Wtime() - inicio;
        printf("NUMERO DE PROCESOS: %d. TIEMPO TOTAL: %f\n", p, fin);
        printf("ESTIMACION USANDO n=%d TRAPECIOS,\n", n);
        printf("DE LA INTEGRAL DESDE %f HASTA %f = %f\n", a, b, integral);
        printf("\nESTIMACION DE PI: %f\n", 2 * integral);
    }

    MPI_Finalize();
} /* MAIN */

/* FUNCION QUE ESTAMOS INTEGRANDO */
float f(float x) {
    float return_val;
    /* CALCULA f(x) Y DEVUELVE SU VALOR */
    return_val = log(x);
    return return_val;
} /* f */
