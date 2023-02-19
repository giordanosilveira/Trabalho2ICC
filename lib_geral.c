#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "lib_geral.h"
#include "lib_sislin.h"

double timestamp(void)
{

  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
  return((double)(tp.tv_sec + tp.tv_nsec*1.0e-9));
}

int verificarArgumentos(int *flags, char *argumentosFaltantes)
{

    // Percorre a lista de flags e verifica se alguns dos argumentos obrigatórios está faltando
    for (int i = 0; i < N_FLAGS; ++i)
    {
        if (i != 4 && flags[i] == 0)
        {
            switch (i)
            {
            case 0:
                strcat(argumentosFaltantes, "n ");
                break;
            case 1:
                strcat(argumentosFaltantes, "k ");
                break;
            case 2:
                strcat(argumentosFaltantes, "p ");
                break;
            case 3:
                strcat(argumentosFaltantes, "i ");
                break;
            case 5:
                strcat(argumentosFaltantes, "o ");
                break;
            }
        }
    }

    if (strlen(argumentosFaltantes) > 0)
    {
        argumentosFaltantes[strlen(argumentosFaltantes) - 1] = '\0';
        return -1;
    }

    return 0;
}


int validarArgumentos(char **argumentos, int *tamanhoSL, int *k_diagonais, int *nInteracoes, int *pre_condicionador, real_t *erro)
{

    char *lixo;

    *(tamanhoSL) = strtol(argumentos[0], &lixo, 10);
    *(k_diagonais) = strtol(argumentos[1], &lixo, 10);
    *(pre_condicionador) = strtol(argumentos[2], &lixo, 10);
    *(nInteracoes) = strtol(argumentos[3], &lixo, 10);

    if (argumentos[4] != NULL)
        sscanf(argumentos[4], "%lg", erro);
    else {
        *(erro) = -1.0;
    }

    if (*(tamanhoSL) < TAM_SL_MIN)
    {
        fprintf(stderr, "A dimensão do Sistema Linear deve ser maior que 10\n");
        return -1;
    }

    if ((*(k_diagonais) < 1) || (*(k_diagonais) % 2 == 0))
    {
        fprintf(stderr, "O número de diagonais da matriz A deve ser maior que 1 e deve ser ímpar\n");
        return -1;
    }

    if (*(pre_condicionador) < 0)
    {
        fprintf(stderr, "O pré-condicionador deve ser 0 (sem pré-condicionador) ou maior que 0 (pré-condicionador de Jacobi)\n");
        return -1;
    }

    if (*(nInteracoes) <= 0)
    {
        fprintf(stderr, "O número de interações do método deve ser maior ou igual a 0\n");
        return -1;
    }

    if (*erro > -1.0) { 
        if (((*(erro) < DBL_EPSILON) && (*(erro) > -DBL_EPSILON)) || (*(erro) < ERRO_MIN) || (*(erro) >= ERRO_MAX))
        {
            fprintf(stderr, "O erro deve ser um número 0 < x < 1\n");
            return -1;
        }
    }

    return 0;
}


void *alocarVetor(int tamanho, int size)
{

    void *v;

    // Aloca espaço na memória e verifica se foi alocado.
    v = calloc(tamanho, size);
    if (!v)
    {
        fprintf(stderr, "Não foi possível alocadr espaço para o vetor\n");
        exit(1);
    }

    return v;
}

// Coeficientes_t *calcularTransposta(SistLinear_t *orig){

//     Coeficientes_t *orig_transp = alocarCoeficiente(orig->n, orig->k);
//     if (! orig_transp) {
//         perror("Não foi possível alocar espaço para a matriz transposta");
//         return NULL;
//     }

//     for (int i = 0; i < (orig->k / 2); ++i){
//         for (int j = 0; j < orig->n; ++j) {
//             if (j + i + 1 < orig->n) {
//                 orig_transp->diagonais_superiores[i][j] = orig->A->diagonais_inferiores[i][j + i + 1];
//                 orig_transp->diagonais_inferiores[i][j + i + 1] = orig->A->diagonais_superiores[i][j];
//             }
//         } 
//     }
   
//     for (int j = 0; j < orig->n; ++j) {
//         orig_transp->diagonal_princial[j] = orig->A->diagonal_princial[j];
//     } 
    

//     return orig_transp;
// }

// void liberarVetor(void *v){
//     if (v)
//         free(v);
// }


// void liberarMatriz(real_t **matriz) {
//     if (matriz[0])
//         free(matriz[0]);
//     if (matriz)
//         free(matriz);
// }


real_t **alocarMatriz(unsigned int n, unsigned int k, unsigned int tam_ptr, unsigned int tam_ele) {

    real_t **A;
    
    A = malloc(n*tam_ptr);
    if (A) {
        
        A[0] = (real_t*)calloc(n*k, tam_ele);
        if (!A[0]) {
            fprintf(stderr, "Não foi possível alocar espaço para o vetor de elementos\n");
            return NULL;
        }
        
        for (int i = 1; i < n; ++i) {
            A[i] = A[0] + i * k;
        }

        
    }
    return A;
    
}


// void cpyMatriz(real_t**dest, real_t**ori, unsigned int n) {

//     for (int i = 0; i < n; ++i){
//         for (int j = 0; j < n; ++j) {
//             dest[i][j] = ori[i][j];
//         }
//     }
    
// }


void cpyVetor(real_t *dest, real_t *orig, unsigned int *tam)
{

    for (int i = 0; i < *(tam); ++i)
    {
        dest[i] = orig[i];
    }
}




// real_t multiplicarVtxV(real_t *vt, real_t* v, unsigned int *tam)
// {

//     real_t soma = 0.0;

//     if (!v) {
        
//         // o produto de uma matriz 1xN por sua tranposta é o somatório dos elementos ao quadrado
//         for (int i = 0; i < *(tam); ++i)
//         {
//             soma = soma + vt[i] * vt[i];          //overflow
//             if (isnan(soma) || isinf(soma))
//             {
//                 fprintf(stderr, "Erro soma(calcularNumeradorEscalarA): %g é NaN ou +/-Infinito\n", soma);
//                 exit(1);
//             }
//         }
//         return soma;
//     }

//     else {

//         //Somatório do produto de dois vetores
//         for (int i = 0; i < *(tam); ++i)
//         {
//             soma = soma + vt[i] * v[i];                //overflow
//             if (isnan(soma) || isinf(soma))
//             {
//                 fprintf(stderr, "Erro soma(calcularNumeradorEscalarA): %g é NaN ou +/-Infinito\n", soma);
//                 exit(1);
//             }
//         }
//         return soma;
//     }


// }


// void calcularMxb(real_t *b, real_t *M, unsigned int n) {

//     for (int i = 0; i < n; ++i) {
//         b[i] = M[i] * b[i];
//         if (isnan(b[i]) || isinf(b[i]))
//         {
//             fprintf(stderr, "Erro SL->b[i](calcularMxb): %g é NaN ou +/-Infinito, Linha: %i \n", b[i], i);
//             exit(1);
//         }
//     }

// }


// void calcularMxA(real_t **A, real_t *M, unsigned int n) {

//     for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < n; ++j) {
//             A[i][j] = A[i][j]*M[i];
//             if (isnan(A[i][j]) || isinf(A[i][j]))
//             {
//                 fprintf(stderr, "Erro A[i][j](calcularMxA): %g é NaN ou +/-Infinito, Linha: %i, Coluna: %i \n", A[i][j], i, j);
//                 exit(1);
//             }
//         }
//     }

// }


// void prnVetor (FILE* arq_saida, real_t *v, unsigned int n)
// {
//     int i;

//     fprintf (arq_saida, "\n");
//     for(i=0; i < n; ++i)
//         fprintf (arq_saida, "%.15g ", v[i]);
//     fprintf (arq_saida, "\n\n");

// }


