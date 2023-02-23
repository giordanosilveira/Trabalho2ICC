#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "lib_geral_v1.h"
#include "lib_sislin_v1.h"

double timestamp(void)
{
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
    return ((double)(tp.tv_sec + tp.tv_nsec * 1.0e-9));
}

void liberarVetor(void *v)
{
    if (v)
    {
        free(v);
    }
}

void liberarMatriz(real_t **matriz)
{
    if (matriz[0])
        free(matriz[0]);
    if (matriz)
        free(matriz);
}

real_t **alocarMatriz(unsigned int n, unsigned int tam_ptr, unsigned int tam_ele)
{

    real_t **A;

    A = malloc(n * tam_ptr);
    if (!A)
    {
        fprintf(stderr, "Não foi possível alocar espaço para o vetor de ponteiros\n");
        exit(1);
    }

    A[0] = (real_t *)calloc(n * n, tam_ele);
    if (!A[0])
    {
        fprintf(stderr, "Não foi possível alocar espaço para o vetor de elementos\n");
        free(A);
        exit(1);
    }

    // Ajusta os ponteiro da matriz
    for (int i = 1; i < n; ++i)
    {
        A[i] = A[i - 1] + n;
    }

    return A;
}

void cpyMatriz(real_t **dest, real_t **ori, unsigned int n)
{

    // Percorre a matrizes e copia o valor da origem para o destino
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            dest[i][j] = ori[i][j];
        }
    }
}

void cpyVetor(real_t *dest, real_t *orig, unsigned int *tam)
{

    // Percorre os vetores e copia os elementos da origem para o destino.
    for (int i = 0; i < *(tam); ++i)
    {
        dest[i] = orig[i];
    }
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

real_t multiplicarVtxV(real_t *vt, real_t *v, unsigned int *tam)
{

    real_t soma = 0.0;

    if (!v)
    {

        // o produto de um vetor Nx1 por sua tranposta 1xN gera uma matriz 1x1, ou seja,
        // um número real.
        for (int i = 0; i < *(tam); ++i)
        {
            soma = soma + vt[i] * vt[i];

            // Teste para ver se não foi gerado um NaN ou um número infinito.
            if (isnan(soma) || isinf(soma))
            {
                fprintf(stderr, "Erro soma(calcularNumeradorEscalarA): %g é NaN ou +/-Infinito\n", soma);
                exit(1);
            }
        }
        return soma;
    }

    else
    {

        // o produto de um vetor Nx1 por sua tranposta 1xN gera uma matriz 1x1, ou seja,
        // um número real.
        for (int i = 0; i < *(tam); ++i)
        {
            soma = soma + vt[i] * v[i];
            // Teste para ver se não foi gerado um NaN ou um número infinito.
            if (isnan(soma) || isinf(soma))
            {
                fprintf(stderr, "Erro soma(calcularNumeradorEscalarA): %g é NaN ou +/-Infinito\n", soma);
                exit(1);
            }
        }
        return soma;
    }
}

void calcularMxb(real_t *b, real_t *M, unsigned int n)
{

    // Percorre o vetor b. Bi = Mi * Bi
    for (int i = 0; i < n; ++i)
    {
        b[i] = M[i] * b[i];

        // Teste para ver se não foi gerado um NaN ou um número infinito.
        if (isnan(b[i]) || isinf(b[i]))
        {
            fprintf(stderr, "Erro SL->b[i](calcularMxb): %g é NaN ou +/-Infinito, Linha: %i \n", b[i], i);
            exit(1);
        }
    }
}

void calcularMxA(real_t **A, real_t *M, unsigned int n)
{

    // Percorre a matriz A. Aij = Aij * Mi
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            A[i][j] = A[i][j] * M[i];

            // Teste para ver se não foi gerado um NaN ou um número infinito.
            if (isnan(A[i][j]) || isinf(A[i][j]))
            {
                fprintf(stderr, "Erro A[i][j](calcularMxA): %g é NaN ou +/-Infinito, Linha: %i, Coluna: %i \n", A[i][j], i, j);
                exit(1);
            }
        }
    }
}

void prnVetor(FILE *arq_saida, real_t *v, unsigned int n)
{
    int i;

    fprintf(arq_saida, "\n");
    for (i = 0; i < n; ++i)
        fprintf(arq_saida, "%.15g ", v[i]);
    fprintf(arq_saida, "\n\n");
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
    else
    {
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

    if (*erro > -1.0)
    {
        if (((*(erro) < DBL_EPSILON) && (*(erro) > -DBL_EPSILON)) || (*(erro) < ERRO_MIN) || (*(erro) >= ERRO_MAX))
        {
            fprintf(stderr, "O erro deve ser um número 0 < x < 1\n");
            return -1;
        }
    }

    return 0;
}