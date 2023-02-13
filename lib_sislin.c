#include "lib_sislin.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "lib_geral.h"


/***********************
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * i,j: coordenadas do elemento a ser calculado (0<=i,j<n)
 * k: numero de diagonais da matriz A
 ***********************/
 static inline double generateRandomA(unsigned int i, unsigned int j, unsigned int k)
{
    static double invRandMax = 1.0 / (double)RAND_MAX;
    return ((i == j) ? (double)(k << 1) : 1.0) * (double)rand() * invRandMax;
}

/***********************
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * k: numero de diagonais da matriz A
 ***********************/
 static inline double generateRandomB(unsigned int k)
{
    static double invRandMax = 1.0 / (double)RAND_MAX;
    return (double)(k << 2) * (double)rand() * invRandMax;
}


void liberarGradientes(Gradiente_t *grad){

    if (grad) {
        if (grad->aux_A[0])
            free(grad->aux_A[0]);
        if (grad->aux_A)
            free(grad->aux_A);
        if (grad->aux_b)
            free(grad->aux_b);
        if (grad->vetor_d)
            free(grad->vetor_d);
        if (grad->vetor_r)
            free(grad->vetor_r);
        if (grad->vetor_r0)
            free(grad->vetor_r0);
        if (grad->vetor_x0)
            free(grad->vetor_x0);
    }
    free(grad);

}


Gradiente_t* alocarGradiente(unsigned int n) {

    Gradiente_t *grad = malloc(sizeof(Gradiente_t));
    if (grad) {
        
        grad->aux_A = alocarMatriz(n, sizeof(real_t*), sizeof(real_t));

        grad->aux_b = (real_t*)alocarVetor(n, sizeof(real_t));
        if (! grad->aux_b) {
            fprintf(stderr, "Não foi possível alocar espaço na memória para o vetor\n");
            free(grad->aux_b);
            exit(1);
        }
        grad->vetor_d = (real_t*)alocarVetor(n, sizeof(real_t));
        if (! grad->vetor_d) {
            fprintf(stderr, "Não foi possível alocar espaço na memória para o vetor\n");
            free(grad->aux_b);
            free(grad->vetor_d);
            exit(1);
        }
        grad->vetor_r0 = (real_t*)alocarVetor(n, sizeof(real_t));
        if (! grad->vetor_r0) {
            fprintf(stderr, "Não foi possível alocar espaço na memória para o vetor\n");
            free(grad->aux_b);
            free(grad->vetor_d);
            free(grad->vetor_r0);
            exit(1);
        }
        grad->vetor_r = (real_t*)alocarVetor(n, sizeof(real_t));
        if (! grad->vetor_r) {
            fprintf(stderr, "Não foi possível alocar espaço na memória para o vetor\n");
            free(grad->aux_b);
            free(grad->vetor_d);
            free(grad->vetor_r0);
            free(grad->vetor_r);
            exit(1);
        }
        grad->vetor_x0 = (real_t*)alocarVetor(n, sizeof(real_t));
        if (! grad->vetor_x0) {
            fprintf(stderr, "Não foi possível alocar espaço na memória para o vetor\n");
            free(grad->aux_b);
            free(grad->vetor_d);
            free(grad->vetor_r0);
            free(grad->vetor_r);
            free(grad->vetor_x0);
            exit(1);
        }
    }
    return grad;

}


void tornarDiagonalDominante(SistLinear_t *SL) {

    real_t soma;

    int i = 0;
    while (i < (SL->k / 2) + 1) {
        soma = 0.0;
        for (int j = 0; j <= (SL->k / 2) + i; ++j) {
            soma = soma + SL->A[i][j];
        }
        SL->A[i][i] = soma + 1.0;
        ++i;
    }

    int m = 1;
    for (i = (SL->k / 2) + 1; i < SL->n; ++i){
        soma = 0.0;
        for (int j = m; j <= SL->k + m && j < SL->n; ++j){
            soma = soma + SL->A[i][j];
        }
        ++m;
        SL->A[i][i] = soma + 1.0;
    }
}


void calcularResiduo(real_t**coef, real_t *residuo, real_t*b, real_t *x, int n)
{

    // Percorre a matriz, realiza a soma da linha e tira do resíduo.
    for (int i = 0; i < n; ++i)
    {
        real_t soma = 0.0;
        for (int j = 0; j < n; ++j)
        {
            soma = soma + coef[i][j] * x[j];               //overflow
            if (isnan(soma) || isinf(soma))
            {
                fprintf(stderr, "Erro soma(calcularResiduo): %g é NaN ou +/-Infinito\n", soma);
                exit(1);
            } 
        }
        residuo[i] = b[i] - soma;
    }
}


real_t *calcularAtxB(SistLinear_t *SL) {

    real_t *b = (real_t *)alocarVetor(SL->n, sizeof(real_t));

    // Percorre a matriz de coeficientes e os termos independentes.
    for (int i = 0; i < SL->n; ++i) {
        for (int j = 0; j < SL->n; ++j) {
            // Realiza os cálculos.
            b[i] += SL->A[j][i] * SL->b[j];
            if (isnan(b[i]) || isinf(b[i]))
            {
                fprintf(stderr, "Erro b[i](calcularBxAt): %g é NaN ou +/-Infinito, Linha: %i, Coluna: %i\n", b[i], i, j);
                exit(1);
            }
        }
    }

    return b;

}


real_t **calcularMatrizAtxA(SistLinear_t* SL) {

    real_t **A = alocarMatriz(SL->n, sizeof(real_t*), sizeof(real_t));

    //Percorre a matriz A
    for (int i = 0; i < SL->n; ++i){
        for (int j = 0; j < SL->n; ++j){
            A[i][j] = 0.0;
            // Efetua o cálculo.
            for (int k = 0; k < SL->n; ++k) {
                A[i][j] = A[i][j] + SL->A[k][i] * SL->A[k][j];
                if (isnan(A[i][j]) || isinf(A[i][j]))
                {
                    fprintf(stderr, "Erro A[i][j](calcularMatrizAxAt): %g é NaN ou +/-Infinito, Linha: %i, Coluna: %i\n", A[i][j], i, j);
                    exit(1);
                }
            }
        }
    }
    
    return A;
}


void prnSisLin(FILE* arq_saida, SistLinear_t *SL)
{
    int n = SL->n;

    for (int i = 0; i < n; ++i)
    {
        fprintf(arq_saida, "\n  ");
        for (int j = 0; j < n; ++j)
            fprintf(arq_saida, "%10g", SL->A[i][j]);
        fprintf(arq_saida, "   |   %g", SL->b[i]);
    }
    fprintf(arq_saida, "\n\n");
}


void liberarSisLin(SistLinear_t *SL)
{

    if (SL)
    {

        if (SL->A)
        {
            if (SL->A[0])
                free(SL->A[0]);
            free(SL->A);
        }

        if (SL->b)
            free(SL->b);

        free(SL);
    }
}


SistLinear_t *alocarSisLin(unsigned int n, unsigned int k, unsigned int p)
{

    // Alocação do ponteiro do Sistema Linear
    SistLinear_t *SL = (SistLinear_t *)malloc(sizeof(SistLinear_t));

    if (SL)
    {

        SL->n = n;
        SL->k = k;
        SL->p = p;

        // Alocação das estruturas internar do Sistema Linear.
        SL->A = (real_t **)malloc(n * sizeof(real_t *));
        SL->b = (real_t *)malloc(n * sizeof(real_t));

        // Verifica se não foi alocado.
        if (!(SL->A) || !(SL->b))
        {
            liberarSisLin(SL);
            return NULL;
        }

        SL->A[0] = (real_t *)malloc(n * n * sizeof(real_t));
        // Verifica se não foi alocado.
        if (!(SL->A[0]))
        {
            liberarSisLin(SL);
            return NULL;
        }

        for (int i = 1; i < n; ++i)
        {
            SL->A[i] = SL->A[i - 1] + n;
        }
    }
    return SL;
}


void initSisLin(SistLinear_t *SL)
{

    // Percorre a matriz de coeficiente e o vetor de termos independentes e inicializa as estruturas
    for (int i = 0; i < SL->n; ++i)
    {
        for (int j = 0; j < SL->n; ++j)
        {
            // SL->A[i][j] = generateRandomA(i, j, SL->k);
            if (i == j)
                SL->A[i][j] = generateRandomA(i, j, SL->k);
            else if (i > j)
            {

                int resp = j + (SL->k / 2);
                if (resp >= i)
                {
                    SL->A[i][j] = generateRandomA(i, j, SL->k);
                }
                else
                    SL->A[i][j] = 0.0;
            }
            else
            {

                int resp = j - (SL->k / 2);
                if (resp <= i)
                {
                    SL->A[i][j] = generateRandomA(i, j, SL->k);
                }
                else
                    SL->A[i][j] = 0.0;
            }
        }
        SL->b[i] = generateRandomB(SL->k);
    }
}
