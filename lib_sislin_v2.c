#include "lib_sislin_v2.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "lib_geral_v2.h"


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



void liberarSisLin(SistLinear_t *sistl) {

    if (sistl) {

        if (sistl->b)
            free(sistl->b);

        if (sistl->A) {
            if (sistl->A[0])
                free(sistl->A[0]);
            free(sistl->A);
        }

        free(sistl);
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

        SL->A = alocarMatriz(n, k, sizeof(real_t*), sizeof(real_t));
        SL->b = (real_t *)malloc(n * sizeof(real_t));

        // Verifica se não foi alocado.
        if (!(SL->A) || !(SL->b)) {
            liberarSisLin(SL);
            return NULL;
        }

        
    }
    return SL;
}

void initSistLinear(SistLinear_t *sistl) {

    int j = 0;
    
    int k = (sistl->k) / 2;
    for (; j < (sistl->k / 2); ++j){
        int i = 0;

        for (; i < sistl->n; ++i) {
            if (i < k)
                sistl->A[i][j] = 0.0;
            else
                sistl->A[i][j] = generateRandomA(i, j, sistl->k);
        }
        --k; 
    }

    for (int i = 0; i < sistl->n; ++i){
        sistl->A[i][j] = generateRandomA(i, i, sistl->k);
    }
    j++;

    k = 1;
    for (; j < sistl->k; ++j){
        int i = 0;
        for (; i < sistl->n; ++i) {
            if (i < sistl->n - k)
                sistl->A[i][j] = generateRandomA(i, j, sistl->k);
            else
                sistl->A[i][j] = 0.0;
        }
        ++k; 
    }

    for (int i = 0; i < sistl->n; ++i)
        sistl->b[i] = generateRandomB(sistl->k);

}

/**
 * @brief 
 * 
 * @param orig 
 * @param dest 
 */
void copiarSistLinear(SistLinear_t*orig, SistLinear_t *dest) {

    for (int i = 0; i < orig->n; ++i){
        for (int j = 0; j < orig->k; ++j){
            dest->A[i][j] = orig->A[i][j];
        }
    }

}


SistLinear_t* calcularTransposta(SistLinear_t *original) {

    SistLinear_t* transposto =  alocarSisLin(original->n, original->k, original->p);

    if (transposto) {

        copiarSistLinear(original, transposto);
        for (int i = 0; i < original->n - (original->k/2) ; ++i){
            
            int deslocamento = 1;
            double aux;
            for (int j = (original->k/2) + 1; j < original->k; ++j){

                aux = transposto->A[i + deslocamento][(original->k/2) - deslocamento]; 
                transposto->A[i + deslocamento][(original->k/2) - deslocamento] = transposto->A[i][j];
                transposto->A[i][j] = aux;
                
                deslocamento++; 
            }
            deslocamento = 1;
        }

        return transposto;

    }
    return transposto;

}


void liberarGradientes(Gradiente_t *grad){

    if (grad) {

        if (grad->vetor_d)
            free(grad->vetor_d);
        if (grad->vetor_r)
            free(grad->vetor_r);
        if (grad->vetor_r0)
            free(grad->vetor_r0);
        if (grad->vetor_x0)
            free(grad->vetor_x0);
        
        free(grad);

    }

}


Gradiente_t* alocarGradiente(unsigned int n) {

    Gradiente_t *grad = malloc(sizeof(Gradiente_t));
    if (grad) {

        grad->vetor_d = (real_t*)alocarVetor(n, sizeof(real_t));
        grad->vetor_r = (real_t*)alocarVetor(n, sizeof(real_t));
        grad->vetor_r0 = (real_t*)alocarVetor(n, sizeof(real_t));
        grad->vetor_x0 = (real_t*)alocarVetor(n, sizeof(real_t));

        if (!grad->vetor_d || !grad->vetor_r || !grad->vetor_r0 || !grad->vetor_x0) {
            fprintf(stderr, "Erro ao alocar as estruturas do gradiente\n");
            liberarGradientes(grad);
        }

        return grad;
    }

    return grad;
}


void tornarDiagonalDominante(SistLinear_t *SL) {

    int i = 0;
    int j;
    real_t soma = 0.0;
    for (; i < SL->n; ++i){
        for (j = 0; j < SL->k; ++j)
            if (j != SL->k/2)
                soma += SL->A[i][j];
        SL->A[i][SL->k / 2] += soma;
    }
    
}


void calcularResiduo(SistLinear_t *SL, real_t *residuo, real_t *x, int n)
{
    int i = 0;
    int j;
    int k = 0;
    real_t soma;
    for (; i < (SL->k/2); ++i) {
        soma = 0.0;
        for (j = (SL->k/2) - i; j < SL->k; ++j) {
            soma += SL->A[i][j] * x[k];
            ++k;
        }
        residuo[i] = SL->b[i] - soma;
        k = 0;
    }

    int deslocamento = 0;
    for (; i < SL->n - (SL->k / 2); ++i) {
        soma = 0.0;
        for (j = 0; j < SL->k; ++j) {
            soma += SL->A[i][j] * x[j + deslocamento];
        }
        residuo[i] = SL->b[i] - soma;
        ++deslocamento;
    }

    k = 1;
    for (; i < SL->n; ++i){
        soma = 0.0;
        for (j = 0; j < SL->k - k; ++j) {
           soma += SL->A[i][j] * x[j + deslocamento]; 
        }
        ++deslocamento;
        residuo[i] = SL->b[i] - soma;
        ++k;
    }

// if (isnan(soma) || isinf(soma))
//             {
//                 fprintf(stderr, "Erro soma(calcularResiduo): %g é NaN ou +/-Infinito\n", soma);
//                 exit(1);
//             } 
}


real_t calcularNormaL2Residuo(SistLinear_t *SL, real_t *x, real_t *tempo)
{

    real_t soma_v[UNROLL] = {0.0, 0.0, 0.0, 0.0};
    real_t soma = 0.0;
    real_t raiz;
    int i = 0;

    // Pecorre o vetor de soluções
    *tempo = timestamp();
    for (; i < SL->n-(SL->n%UNROLL); i += UNROLL) {

        // Soma o quadrado dos elementos das soluções
        soma_v[0] = soma_v[0] + x[i]*x[i];
        soma_v[1] = soma_v[1] + x[i+1]*x[i+1];
        soma_v[2] = soma_v[2] + x[i+2]*x[i+2];
        soma_v[3] = soma_v[3] + x[i+3]*x[i+3];
        // if (isnan(soma) || isinf(soma))
        // {
        //     fprintf(stderr, "Erro soma(calcularNormaL2Residuo): %g é NaN ou +/-Infinito\n", soma);
        //     exit(1);
        // }

    }
    soma = soma_v[0] + soma_v[1] + soma_v[2] + soma_v[3]; 
    for(; i < SL->n; ++i)
        soma = soma + x[i]*x[i];
        
    raiz = sqrt(soma);
    *tempo = timestamp() - *tempo;

    // Retorna a raíz quadrada da soma.
    return raiz;

}


// real_t *calcularAtxB(SistLinear_t *SL) {

//     real_t *b = (real_t *)alocarVetor(SL->n, sizeof(real_t));

//     // Percorre a matriz de coeficientes e os termos independentes.
//     for (int i = 0; i < SL->n; ++i) {
//         for (int j = 0; j < SL->n; ++j) {
//             // Realiza os cálculos.
//             b[i] += SL->A[j][i] * SL->b[j];
//             if (isnan(b[i]) || isinf(b[i]))
//             {
//                 fprintf(stderr, "Erro b[i](calcularBxAt): %g é NaN ou +/-Infinito, Linha: %i, Coluna: %i\n", b[i], i, j);
//                 exit(1);
//             }
//         }
//     }

//     return b;

// }


SistLinear_t *calcularMatrizAtxA(SistLinear_t* SL) {

    SistLinear_t* AtxA = alocarSisLin(SL->n, SL->k + (SL->k - 1), SL->p);
    real_t *linha = (real_t*)alocarVetor(SL->k, sizeof(real_t));
    real_t soma = 0.0;

    if (AtxA && linha) {
        
        int i = 0;
        int j;
        int m;
        
        int deslocamento;
        int start = (AtxA->k / 2);
        for (; i < SL->k/2; ++i){
            deslocamento = 0;
            for (m = i; m < (SL->k + i) && (m < SL->n); ++m){
                soma = 0.0;
                for (j = SL->k/2 - i; j < SL->k; ++j){
                    if (j - deslocamento > -1)
                        soma += SL->A[i][j]*SL->A[m][j - deslocamento];
                }
                AtxA->A[i][start] = soma;
                if (i + deslocamento < SL->n)
                    AtxA->A[i + deslocamento][(AtxA->k / 2) - deslocamento] = soma;

                deslocamento++;
                start++;
            }
            start = (AtxA->k / 2);
        }
        
        
        for (; i < SL->n; ++i){
            deslocamento = 0;
            for(m = i; m < (SL->k + i) && (m < SL->n); ++m){
                soma = 0.0;
                for (j = 0; j < SL->k; ++j) {
                    if (j - deslocamento > -1)
                        soma += SL->A[i][j]*SL->A[m][j - deslocamento];
                }
                AtxA->A[i][start] = soma;
                if (i + deslocamento < SL->n)
                    AtxA->A[i + deslocamento][(AtxA->k / 2) - deslocamento] = soma;

                deslocamento++;
                start++;
            }
            start = (AtxA->k / 2);
        }
        //    if (isnan(soma) || isinf(soma))
        //     {
        //         //fprintf(stderr, "Erro A[i][j](calcularMatrizAxAt): %g é NaN ou +/-Infinito, Linha: %i, Coluna: %i\n", A[i][j], i, j);
        //         liberarSisLin(AtxA);
        //         exit(1);
        //     }
        
    }
    
    return AtxA;
}

// void prnCoef(FILE* arq_saida, Coeficientes_t *A, unsigned int n, unsigned int k)
// {

//     for (int i = 0; i < k / 2; ++i)
//     {
//         fprintf(arq_saida, "\n  Diagonal inferior\n");
//         for (int j = 0; j < n; ++j) {
//             fprintf(arq_saida, "%10g ", A->diagonais_inferiores[i][j]);
//         }
//         fprintf(arq_saida, "\n  Diagonal superiores\n");
//         for (int j = 0; j < n; ++j) {
//             fprintf(arq_saida, "%10g ", A->diagonais_superiores[i][j]);
//         }
//     }

//     fprintf(arq_saida, "\n\n\n");
//     for (int i = 0; i < n; ++i){
//         fprintf(arq_saida, "%10g ", A->diagonal_princial[i]);
//     }
    
// }

void prnSisLin(FILE* arq_saida, SistLinear_t *SL)
{

    for (int i = 0; i < SL->n; ++i)
    {
        for (int j = 0; j < SL->k; ++j)
            fprintf(arq_saida, "%10g ", SL->A[i][j]);
        fprintf(arq_saida, "|%10g ", SL->b[i]);
        fprintf(arq_saida, "\n");
    }     
    fprintf(arq_saida, "\n\n");

}

