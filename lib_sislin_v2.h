#ifndef _LIB_SISLIN_
#define _LIB_SISLIN_
#include <stdio.h>

typedef double real_t;


typedef struct
{
    real_t **A;                     //Coeficientes   
    real_t *b;                      // termos independentes
    unsigned int n;                 // tamanho do sistema linear
    unsigned int k;                 // k diagonais
    unsigned int p;                 // pré condicionador de jacobi

} SistLinear_t;


//Representa as estruturas internas (não todas) usadas dentro do cálculo do método
typedef struct {

    real_t *vetor_r;    //Resíduo
    real_t *vetor_r0;   //Resíduo anterior
    real_t *vetor_d;    //Próxima direção (chamado de 'v' no livro da M.Cristina C. Cunha)
    real_t *vetor_x0;   //Solução anterior
    real_t escalarM;    //Escalar m
    real_t escalarA;    //Escalar a (chamado de 's' no livro da M.Cristina C. Cunha)                      

} Gradiente_t;          //Representa as estruturas internas (não todas) usadas dentro do cálculo do método


/**
 * @brief Libera a estrutura apontada por 'grad'
 * 
 * @param grad (Gradiente_t*) : Ponteiro para a estrutura gradiente.
 */
void liberarGradientes(Gradiente_t *grad);


/**
 * @brief Aloca espaço na memória para o tipo Grandiente e sua estruturas internas
 * 
 * @param n (unsigned int) : Tamanho do dos vetores.
 * @return (Gradiente_t*) : Ponteiro para essa estrutura. 
 */
Gradiente_t* alocarGradiente(unsigned int n);


/**
 * @brief Modifica os coeficientes da matriz A do sistema linear
 * para que ela fique diagonal dominante.
 * 
 * @param SL (SistLinear_t*) : O Sistema Linear.
 */
void tornarDiagonalDominante(SistLinear_t *SL);


/**
 * @brief Calcula o resíduo para um sistema linear.
 * 
 * @param SL (SistLinear_t*) : Sistema Linear.
 * @param residuo (real_t*) : Vetor onde será colocado o resíduo.
 * @param x (real_t*) : Vetor com as soluções.
 * @param n (int) : Tamanho dos vetores. 
 */
void calcularResiduo(SistLinear_t *SL, real_t *residuo, real_t *x, int n);


/**
 * @brief Calcula a norma L2 do vetor com as soluções 'x'
 * 
 * @param tempo (real_t*) : Ponteiro onde será guardado o tempo do cálculo da norma.
 * @param SL (SistLinear_t*) : O sistema linear.
 * @param x (real_t*) : O Vetor com as soluções
 * @return (real_t) : Norma L2 do vetor. 
 */
real_t calcularNormaL2Residuo(SistLinear_t *SL, real_t *x, real_t *tempo);


/**
 * @brief Efetua uma multiplcação de matrizes entre a matriz transposta 
 * A e A. Essa transformação é necessária para ser possível efetuar o métodos dos
 * Gradientes conjugados.
 * 
 * @param SL (SistLinear_t*) : Ponteiro para o sistema linear.
 * @return (SistLinear_t*) : Matriz resultante da multiplicação. 
 */
SistLinear_t *calcularMatrizAtxA(SistLinear_t* SL);


/**
 * @brief Printa o sistema linear.
 * 
 * @param arq_saida (FILE*) : Arquivo onde será printado o Sistema Linear
 * @param SL (SistLinear_t*) : O Sistema Linear
 */
void prnSisLin(FILE * arq_saida, SistLinear_t *SL);


/**
 * @brief Libera as estruturas do sistema linear
 *
 * @param SL (SistLinear_t*) : O sistema Linear
 */
void liberarSisLin(SistLinear_t *SL);

/**
 * @brief Aloca espaço na memória para o Sistema Linear.
 *
 * @param n (unsigned int) : Tamanho do sistema linear.
 * @param k (unsigned int) : Número de k-diagonais.
 * @param p (unsigned int) : Pré-condicionadores.
 * @return (SistLinear_t*) : Um ponteiro para o Sistema linear alocado.
 */
SistLinear_t *alocarSisLin(unsigned int n, unsigned int k, unsigned int p);

/**
 * @brief Inicializa o sistema linear.
 *
 * @param SL (SistLinear_t *) : O sistema linear
 */
void initSistLinear(SistLinear_t *SL);

/**
 * @brief Calcula os sistema linear com os coeficientes transpostos.
 * 
 * @param original (SistLinear_t*) : Sistema Linear original
 * @return SistLinear_t* : Sistema Linear com os coeficientes transpostos. 
 */
SistLinear_t* calcularTransposta(SistLinear_t *original);

/**
 * @brief Copia um sistema linear inteiro.
 * 
 * @param orig (SistLinear_t*) : Sistema Linear originário. 
 * @param dest (SistLinear_t*) : Sistema Linear copiado.
 */
void copiarSistLinear(SistLinear_t*orig, SistLinear_t *dest);


#endif