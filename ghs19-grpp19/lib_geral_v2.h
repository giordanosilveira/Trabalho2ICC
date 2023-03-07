#ifndef _LIB_GERAL_
#define _LIB_GERAL_

#define N_FLAGS 6       // Número max de flags permitidas
#define TAM_SL_MIN 3    // Tamanho mínimo do sistema linear
#define ERRO_MIN 0.0    // erro mínimo permitido (exclusive)
#define ERRO_MAX 1.0    // erro máximo permitido (exclusive)
#define ERRO_IT -1.0    // erro para testar o número de iterações
#define UNROLL 4        // Unroll

// Módulo de um número
#define ABS(num) ((num) < 0.0 ? -(num) : (num))


#include "lib_sislin_v2.h"


/**
 * @brief Função que retorna 
 * 
 * @return double 
 */
double timestamp(void);

/**
 * @brief Verifica quais argumentos foram passados na linha de comando.
 *
 * @param flags (int*) : Vetor com as flags representando os argumentos que foram passados na linha de comando.
 * @param argumentosFaltantes (char *) : Vetor que indicará, se precisar, quais argumentos obrigatórios estão faltando.
 * @return (int) : 0 se todos os argumentos foram passados; -1 se faltou algum obrigatório.
 */
int verificarArgumentos(int *flags, char *argumentosFaltantes);

/**
 * @brief Valida os argumentos que foram passados na linha de comando.
 *
 * @param argumentos (char**) : Vetor de ponteiros com os argumentos que foram passados.
 * @param tamanhoSL (int*) : Ponteiro onde será colocado o tamanho do sistema linear.
 * @param k_diagonais (int*) : Ponteiro onde será colocado o número de k-diagonais.
 * @param nInteracoes (int*) : Ponteiro onde será colocado o número máximo de interações.
 * @param pre_condicionador (int*) :  Ponteiro onde será colocado os pré-condicionadores.
 * @param erro (real_t*) : Ponteiro onde será coloca o erro máximo.
 * @return (int) : 0 se todos os argumentos são válidos; -1 se algum não está
 */
int validarArgumentos(char **argumentos, int *tamanhoSL, int *k_diagonais, int *nInteracoes, int *pre_condicionador, real_t *erro);

/**
 * @brief Se for possível, libera o vetor 'v'.
 * 
 * @param v (void *) : vetor para ser liberado.
 */
void liberarVetor(void *v);


/**
 * @brief Aloca uma matriz genérica na memória
 * 
 * @param n (unsigned int) : Número de linhas da matriz.
 * @param k (unsigned int) : Número de colunas da matriz.
 * @param tam_ptr (unsigned int) : tamanho dos elementos.
 * @param tam_ele (unsigned int) : tamanho dos elementos.
 * @return (real_t**) : Ponteiro para a matriz 
 */
real_t **alocarMatriz(unsigned int n, unsigned int k, unsigned int tam_ptr, unsigned int tam_ele);


/**
 * @brief Aloca um vetor genério na memória.
 *
 * @param tamanho (int) : tamanho do vetor.
 * @param size (int) : tamanho dos elementos do vetor.
 * @return (void*) : Ponteiro genérico para esse vetor.
 */
void *alocarVetor(int tamanho, int size);



/**
 * @brief Copia um vetor para outro.
 *
 * @param dest (real_t*) : vetor destino.
 * @param orig (real_t*) : vetor origem.
 * @param tam (int) : tamanho dos vetores.
 */
void cpyVetor(real_t *orig, real_t *dest, unsigned int *tam);


/**
 * @brief Calcula o vetor transposto 'vt' pelo vetor 'v'. Isso resulta em uma
 * matriz 1x1, ou seja, um número real.
 * 
 * @param vt (real_t*) : Vetor que será multiplicado.
 * @param tam (unsigned int*) : Tamanho do vetor.
 * @return real_t : o resultado da multiplicação dos dos vetores.
 */
real_t multiplicarMesmoVtxV(real_t *vt, unsigned int *tam);


/**
 * @brief Calcula o vetor transposto 'vt' pelo vetor 'v'. Isso resulta em uma
 * matriz 1x1, ou seja, um número real.
 *
 * @param vt (real_t*) : Vetor com o vt.
 * @param v (real_t*) : Vetor v (pode não ser usado)
 * @param tam (int*) : Ponteiro que indica o tamanho do vetor do resíduo.
 * @return (real_t) : Numerador para o cálculo do alfa.
 */
real_t multiplicarVtxV(real_t* restrict vt, real_t* restrict v, unsigned int *tam);


/**
 * @brief Multiplica uma matriz por um vetor. Levando em conta que não
 * é uma matriz completa, mas sim, só suas diagonais.
 * 
 * @param SL (SistLinear_t*) : O Sistema Linear cuja matriz de coeficientes será usada.
 * @param v (real_t*) : Vetor.
 * @param resultado (real_t*) : Vetor onde será colocado o resultado da multiplicação.
 */
void multiplicarMatrizPorVetor(SistLinear_t *SL, real_t *v, real_t *resultado);

/**
 * @brief Multiplica um vetor por uma matriz. Levando em conta que não
 * é uma matriz completa, mas sim, só suas diagonais.
 * 
 * @param vetor (real_t*) : Vetor.
 * @param SL (SistLinear_t*) : O Sistema Linear cuja matriz de coeficientes será usada.
 */
void multiplicarVetorPorMatriz(real_t *vetor, SistLinear_t *SL);

/**
 * @brief Multiplica um vetor v1 por um vetor v2.
 * 
 * @param v1 (real_t*) : Vetor v1.
 * @param v2 (real_t*) : Vetor v2.
 * @param n (unsigned int *) : Tamanho dos vetores.
 */
void multiplicarVetorPorVetor(real_t* restrict v1, real_t* restrict v2, unsigned int *n);



/**
 * @brief Printa o vetor específicado
 * 
 * @param arq_saida (FILE*) : Onde será printado a saída formatada pedida pelo professor
 * @param v (real_t*) : Vetor.
 * @param n (unsigned int n) : Tamanho do vetor.
 */
void prnVetor (FILE* arq_saida, real_t *v, unsigned int n);




#endif