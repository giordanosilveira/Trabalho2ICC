#ifndef _LIB_GERAL_
#define _LIB_GERAL_

#define N_FLAGS 6     // Número max de flags permitidas
#define TAM_SL_MIN 10 // Tamanho mínimo do sistema linear
#define ERRO_MIN 0.0  // erro mínimo permitido (exclusive)
#define ERRO_MAX 1.0  // erro máximo permitido (exclusive)

// Módulo de um número
#define ABS(num) ((num) < 0.0 ? -(num) : (num))


#include "lib_sislin.h"


/**
 * @brief Função que retorna 
 * 
 * @return double 
 */
double timestamp(void);

/**
 * @brief Se for possível, libera o vetor 'v'
 * 
 * @param v (void *) : vetor para ser liberado.
 */
void liberarVetor(void *v);

/**
 * @brief Libera a matriz 'matriz' se possível.
 * 
 * @param matriz (real_t*) : Matriz para ser liberada.
 */
void liberarMatriz(real_t **matriz);


/**
 * @brief Aloca uma matriz genérica na memória
 * 
 * @param n (unsigned int) : tamanho da matriz.
 * @param tam_ptr (unsigned int) : tamanho dos ponteiros.
 * @param tam_ele (unsigned int) : tamanho dos elementos.
 * @return (void**) : Ponteiro para a matriz 
 */
real_t **alocarMatriz(unsigned int n, unsigned int tam_ptr, unsigned int tam_ele);


/**
 * @brief Copia os elementos da matriz origem (ori)
 * para a matriz destino (dest)
 * 
 * @param dest (real_t**) : Matriz destino.
 * @param ori (real_t**) : Matriz origem.
 * @param n (unsigned int) : Tamanho da matriz.
 */
void cpyMatriz(real_t**dest, real_t**ori, unsigned int n);
 
/**
 * @brief Aloca um vetor genério na memória.
 *
 * @param tamanho (int) : tamanho do vetor.
 * @param size (int) : tamanho dos elementos do vetor.
 * @return (void*) : Ponteiro genérico para esse vetor.
 */
void *alocarVetor(int tamanho, int size);


/**
 * @brief Copia um vetor para outro
 *
 * @param dest (real_t*) : vetor destino
 * @param orig (real_t*) : vetor origem
 * @param tam (int) : tamanho dos vetores
 */
void cpyVetor(real_t *dest, real_t *orig, unsigned int *tam);

/**
 * @brief Calcula o vetor transposto 'vt' pelo vetor 'v'. Isso resulta em uma
 * matriz 1x1, ou seja, um número real.
 *
 * @param vt (real_t*) : Vetor com o vt.
 * @param v (real_t*) : Vetor v (pode não ser usado)
 * @param tam (int*) : Ponteiro que indica o tamanho do vetor do resíduo.
 * @return (real_t) : Numerador para o cálculo do alfa.
 */
real_t multiplicarVtxV(real_t *vt, real_t *v, unsigned int *tam);

/**
 * @brief multiplica a matriz M (vetor com a diagonal principal) com um vetor 'b'.
 * 
 * @param b (real_t*) : Vetor b (geralmente os termos independentes.)
 * @param M (real_t*) : Vetor M (geralmente o vetor contendo os elementos da diagonal principal da matriz).
 * @param n (unsigned int) : Tamanho dos vetores
 */
void calcularMxb(real_t *b, real_t *M, unsigned int n);

/**
 * @brief multiplica a matriz A pelo vetor M.
 * 
 * @param A (real_t**) : Matriz A (geralmente os coeficientes do sistema)
 * @param M (real_t*) : matriz M (geralmente o vetor contendo a diagonal principal da matriz A)
 * @param n (unsigned int) : Dimensão da matriz e do vetor.
 */
void calcularMxA(real_t **A, real_t *M, unsigned int n);


/**
 * @brief Printa o vetor específicado
 * 
 * @param arq_saida (FILE*) : Onde será printado a saída formatada pedida pelo professor
 * @param v (real_t*) : Vetor.
 * @param n (unsigned int n) : Tamanho do vetor.
 */
void prnVetor (FILE* arq_saida, real_t *v, unsigned int n);


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

#endif