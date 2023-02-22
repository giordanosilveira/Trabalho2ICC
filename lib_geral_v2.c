#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "lib_geral_v2.h"
#include "lib_sislin_v2.h"

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


void liberarVetor(void *v){
    if (v)
        free(v);
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


void cpyVetor(real_t *dest, real_t *orig, unsigned int *tam)
{

    int i = 0;
    for (; i < *(tam)-*(tam)%UNROLL; i += UNROLL)
    {
        dest[i] = orig[i];
        dest[i+1] = orig[i+1];
        dest[i+2] = orig[i+2];
        dest[i+3] = orig[i+3];
    }
    for (; i < *(tam); ++i)
        dest[i] = orig[i];
}


real_t multiplicarMesmoVtxV(real_t *vt, unsigned int *tam) {

    real_t soma = 0.0;
    real_t soma_v[UNROLL] = {0.0,0.0,0.0,0.0};
    int i;

    for (i = 0; i < *(tam)-*(tam) % UNROLL; i += UNROLL){
        soma_v[0] += vt[i] * vt[i];
        soma_v[1] += vt[i + 1] * vt[i + 1];
        soma_v[2] += vt[i + 2] * vt[i + 2];
        soma_v[UNROLL - 1] += vt[i + 3] * vt[i + 3];
    }

    soma = soma_v[0] + soma_v[1] + soma_v[2] + soma_v[3];
    for (; i < (*tam); ++i)
        soma += vt[i]*vt[i];

    return soma;
}

real_t multiplicarVtxV(real_t* restrict vt, real_t* restrict v, unsigned int *tam) {
    
    int i = 0;
    real_t soma = 0.0;
    real_t soma_v[UNROLL] = {0.0,0.0,0.0,0.0};
    for (; i < *(tam)-*(tam) % UNROLL; i += UNROLL){
        soma_v[0] += vt[i] * v[i];
        soma_v[1] += vt[i + 1] * v[i + 1];
        soma_v[2] += vt[i + 2] * v[i + 2];
        soma_v[UNROLL - 1] += vt[i + 3] * v[i + 3];
    }
    soma = soma_v[0] + soma_v[1] + soma_v[2] + soma_v[3];
    for (; i < (*tam); ++i)
        soma += vt[i]*v[i];

    return soma;
}


void multiplicarMatrizPorVetor(SistLinear_t *SL, real_t *v, real_t *resultado) {

    int i = 0;
    int j;
    int k = 0;
    real_t soma;
    for (; i < (SL->k/2); ++i) {
        soma = 0.0;
        for (j = (SL->k/2) - i; j < SL->k; ++j) {
            soma += SL->A[i][j] * v[k];
            ++k;
        }
        resultado[i] = soma;
        k = 0;
    }

    int deslocamento = 0;
    for (; i < SL->n - (SL->k / 2); ++i) {
        soma = 0.0;
        for (j = 0; j < SL->k; ++j) {
            soma += SL->A[i][j] * v[j + deslocamento];
        }
        resultado[i] = soma;
        ++deslocamento;
    }

    k = 1;
    for (; i < SL->n; ++i){
        soma = 0.0;
        for (j = 0; j < SL->k - k; ++j) {
           soma += SL->A[i][j] * v[j + deslocamento]; 
        }
        ++deslocamento;
        resultado[i] = soma;
        ++k;
    }

}


void multiplicarVetorPorVetor(real_t* restrict v1, real_t* restrict v2, unsigned int *n){

    int i = 0;
    for (; i < *(n)-*(n)%UNROLL; i += UNROLL){
        v1[i] = v2[i] * v1[i];
        v1[i + 1] = v2[i + 1] * v1[i + 1];
        v1[i + 2] = v2[i + 2] * v1[i + 2];
        v1[i + 3] = v2[i + 3] * v1[i + 3];
    }
    for (; i < *(n); ++i)
        v1[i] = v2[i] * v1[i];
    
    // if (isnan(b[i]) || isinf(b[i]))
    //     {
    //         fprintf(stderr, "Erro SL->b[i](calcularMxb): %g é NaN ou +/-Infinito, Linha: %i \n", b[i], i);
    //         exit(1);
    //     }

}

void multiplicarVetorPorMatriz(real_t *vetor, SistLinear_t *SL){

    int i;
    int j;
    for (i = 0; i < SL->n-(SL->n%UNROLL); i += UNROLL) {
        for (j = 0; j < SL->k; ++j) {
            SL->A[i][j] = vetor[i] * SL->A[i][j];
            SL->A[i+1][j] = vetor[i+1] * SL->A[i+1][j];
            SL->A[i+2][j] = vetor[i+2] * SL->A[i+2][j];
            SL->A[i+3][j] = vetor[i+3] * SL->A[i+3][j];
        }
    }

    for (; i < SL->n; ++i) {
        for (j = 0; j < SL->k; ++j)
            SL->A[i][j] = vetor[i] * SL->A[i][j];
    }

}

void prnVetor (FILE* arq_saida, real_t *v, unsigned int n)
{
    int i;

    fprintf (arq_saida, "\n");
    for(i=0; i < n; ++i)
        fprintf (arq_saida, "%.15g ", v[i]);
    fprintf (arq_saida, "\n\n");

}


