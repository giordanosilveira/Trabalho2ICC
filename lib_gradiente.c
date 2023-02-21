#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib_geral.h"
#include "lib_gradiente.h"
#include "lib_sislin.h"


/**
 * @brief Calcula a nova solução 'x'. Ou seja, Xk' = Xk + s*Dk.
 * Também guarda a solução antiga em 'X0'.
 *  
 * @param x0 (real_t*) : Guarda a solução anterior.
 * @param x (real_t*) : Guarda a solução atual.
 * @param d (real_t*) : Próxima direção (chamado de 'v' no livro M.Cristina C. Cunha). 
 * @param escalarA (real_t) : Escalar A (chamado de 's' no livro M.Cristina C. Cunha).
 * @param n (unsigned int) : Tamanho dos vetores.
 */
static inline void calcularNovaSolucao(real_t *x0, real_t *x, real_t *d, real_t escalarA, unsigned int *n) {

    // Percorre o vetor de soluções e realiza o calculo Xk' = Xk + s*Dk.
    int i = 0;
    for (; i < (*n)-(*n)%UNROLL; i += UNROLL) {
        x0[i] = x[i];
        x0[i+1] = x[i+1];
        x0[i+2] = x[i+2];
        x0[i+3] = x[i+3];

        x[i] = x[i] + escalarA*d[i];
        x[i+1] = x[i+1] + escalarA*d[i+1];
        x[i+2] = x[i+2] + escalarA*d[i+2];
        x[i+3] = x[i+3] + escalarA*d[i+3];
    }

    for (; i < (*n); ++i) {
        x0[i] = x[i];
        x[i] = x[i] + escalarA*d[i];
    }

}

/**
 * @brief Calcula o novo resíduo. Ou seja, Rk' = Rk - s*A*D.
 * Também guarda o resíduo antigo em 'R0'.
 * 
 * @param r0 (real_t*) : Guarda o resíduo antigo.
 * @param r (real_t*) : Guarda o resíduo atual.
 * @param y (real_t*) : Vetor que é resultado da multiplicação de A*dk.
 * @param escalarA (real_t) : Escalar A (chamado de 's' no livro da M.Cristina C. Cunha).
 * @param n (unsigned int*) : Tamanho dos vetores
 */
static inline void calcularNovoResiduo(real_t *r0, real_t *r, real_t *y, real_t escalarA, unsigned int *n) {

    //Percorre a matriz dos coeficientes
    int i = 0;
    for (; i < (*n)-(*n)%UNROLL; i += UNROLL) {
        r0[i] = r[i];
        r0[i+1] = r[i+1];
        r0[i+2] = r[i+2];
        r0[i+3] = r[i+3];

        r[i] = r[i] - escalarA*y[i];
        r[i+1] = r[i+1] - escalarA*y[i+1];
        r[i+2] = r[i+2] - escalarA*y[i+2];
        r[i+3] = r[i+3] - escalarA*y[i+3];
    }
    for (; i < (*n); ++i) {
        r0[i] = r[i];
        r[i] = r[i] - escalarA*y[i];
    }

}

/**
 * @brief Calcula o valor de z, isto é Zĸ = C⁻¹ * Rĸ
 * 
 * @link: https://math-linux.com/mathematics/linear-systems/article/preconditioned-conjugate-gradient-method
 * 
 * @param z (real_t* restrict) : Vetor onde ser será colocado o resultado de 'Zĸ = C⁻¹ * Rĸ'.
 * @param inverse_c (real_t* restrict) : A matriz inversa do pré-condicionador de Jacobi.
 * @param residue (real_t* restrict) : O vetor de resíduo.
 * @param size (unsigned int) : O tamanho do vetor de Z.
 * @return (real_t*) : O vetor de Z calculado, após Zĸ = C⁻¹. Rĸ.
 */
static inline void calcula_z(real_t* restrict z, real_t* restrict inverse_c, real_t* restrict residue, unsigned int *size){

    int i = 0;
    for (; i < (*size)-(*size)%UNROLL; i += UNROLL){
        z[i] = inverse_c[i] * residue[i];
        z[i+1] = inverse_c[i+1] * residue[i+1];
        z[i+2] = inverse_c[i+2] * residue[i+2];
        z[i+3] = inverse_c[i+3] * residue[i+3];
    }
    for (; i < *(size); ++i)
        z[i] = inverse_c[i] * residue[i];

}


/**
 * @brief Calcula a próxima direção. Ou seja, Dk' = escalarM*Dk + Vk.
 * Lembrando que o vetor 'v', dependendo do método, pode mudar.
 * 
 * @param d (real_t*) : Guarda a direção após o cálculo.
 * @param v (real_t*) : Vetor que pode ser o resíduo 'r' (S/ pré-condicionador) ou vetor 'z' (C/ pré-condicionador).
 * @param escalarM (real_t) : Escalar M.
 * @param n (unsigned int*) : Tamanho dos vetores.
 */
static inline void calcularProxDirecao(real_t * restrict d, real_t * restrict z, real_t escalarM, unsigned int *n) {

    int i = 0;
    for (; i < (*n)-(*n)%UNROLL; i += UNROLL){
        d[i] = escalarM*d[i] + z[i];
        d[i+1] = escalarM*d[i+1] + z[i+1];
        d[i+2] = escalarM*d[i+2] + z[i+2];
        d[i+3] = escalarM*d[i+3] + z[i+3];
    }
    for(; i < (*n); ++i)
        d[i] = escalarM*d[i] + z[i];

}

/**
 * @brief Calcula a inversa do pré-condicionador de Jacobi, isto é: M⁻¹
 * 
 * @link https://math-linux.com/mathematics/linear-systems/article/preconditioned-conjugate-gradient-method
 * 
 * @param SL (SistLinear_t*) : O sistema linear.
 * @param M (real_t*) : Vetor contendo a diagonal principal do Sistema Linear original
 * 
 * @return (real_t) : Matriz inversa de C, C⁻¹.
 */
static inline void inverse_jacobi_preconditioner(SistLinear_t *SL, real_t *M)
{

    int i;
    for (i = 0; i < SL->n-(SL->n%UNROLL); i += UNROLL) {
        M[i] = (real_t)1 / SL->A[i][SL->k/2];
        M[i+1] = (real_t)1 / SL->A[i+1][SL->k/2];
        M[i+2] = (real_t)1 / SL->A[i+2][SL->k/2];
        M[i+3] = (real_t)1 / SL->A[i+3][SL->k/2];
    }
    for (; i < SL->n; ++i)
        M[i] = (real_t)1 / SL->A[i][SL->k/2];

}

/**
 * @brief 
 * 
 * @param x 
 * @param x0 
 * @param maior_erro_max_abs 
 * @param tam 
 */
static inline void normaMaxErroAbsoluto(real_t *x, real_t *x0, real_t *maior_erro_max_abs, unsigned int *tam)
{

    *maior_erro_max_abs = ABS(x[0] - x0[0]);

    // Percorre o vetor de soluções,
    for (int i = 1; i < *(tam); ++i){

        // o maior erro absoluto.
        if (ABS(x[i] - x0[i]) > *maior_erro_max_abs)
            *maior_erro_max_abs = ABS(x[i] - x0[i]);

    }

}


int gradienteConjugadosCPreCondicionadores(FILE*arq_saida, SistLinear_t *SL, SistLinear_t *SLTransposto, real_t *x, real_t* tempo_metodo, real_t *tempo_preparacao, real_t erro, int nInteracoes) {
 
    Gradiente_t *grad = alocarGradiente(SL->n);
    if (! grad) {
        fprintf(stderr, "Não foi possível alocar espaço para as estruturas do gradiente\n");
        exit(1);
    }

    /*Declaração das variáveis utilizadas internamente pela função*/
    real_t *matrix_M = (real_t*)alocarVetor(SL->n, sizeof(real_t));         //Pré-condicionador de Jacobi M = D
    real_t *vetor_z = (real_t*)alocarVetor(SL->n, sizeof(real_t));          //Chamado de 'y' no livro M.Cristina C. Cunha                
    real_t *vetor_z0 = (real_t*)alocarVetor(SL->n, sizeof(real_t));         //Z anterior.
    real_t *vetor_y = (real_t*)alocarVetor(SL->n, sizeof(real_t));
    real_t max_erro_abs_aprox;
    real_t tempo_inicial, tempo_final;

    
    //Tempo inicial para o cálculo do tempo para a preparação do método
    (*tempo_preparacao) = timestamp();
    

    /*Inicializando as variável que serão utilizadas pelo método*/
    inverse_jacobi_preconditioner(SL, matrix_M);                            //Inicializa a "matriz" 'M'
    SistLinear_t *SLTranspXSL = calcularMatrizAtxA(SLTransposto);
    multiplicarMatrizPorVetor(SLTransposto, SL->b, SLTranspXSL->b);
    multiplicarVetorPorMatriz(matrix_M, SLTranspXSL);                       //Calcular M^-1 * A
    multiplicarVetorPorVetor(SLTranspXSL->b,matrix_M, &SLTranspXSL->n);     //Calcular M^-1 * b                               

    //Tempo final da preparação do método
    (*tempo_preparacao) = timestamp() - (*tempo_preparacao);

    cpyVetor(grad->vetor_r, SL->b, &SL->n);                                 //Inicializa o resíduo 'r0' com SL->b (O X inicial é 0)
    calcula_z(vetor_z, matrix_M, grad->vetor_r, &SL->n);                     //Inicializa 'z0'
    cpyVetor(grad->vetor_d, vetor_z, &SL->n);                               //Inicializa 'd0' 
    
    int i = 0;
    (*tempo_metodo) = 0;
    while (i < nInteracoes) {

        tempo_inicial = timestamp();
        multiplicarMatrizPorVetor(SLTranspXSL, grad->vetor_d, vetor_y);

        grad->escalarA = multiplicarVtxV(vetor_z, grad->vetor_r, &SL->n)/multiplicarVtxV(grad->vetor_d, vetor_y, &SL->n);


        // Calcula o novo valor de 'x' e guarda o antigo em 'x0'
        calcularNovaSolucao(grad->vetor_x0, x, grad->vetor_d, grad->escalarA, &SL->n);

        normaMaxErroAbsoluto(x, grad->vetor_x0, &max_erro_abs_aprox, &SL->n);
        fprintf(arq_saida, "# ||%10g||\n", max_erro_abs_aprox);

        // Calcula o resíduo atual e copia o resíduo anterior
        calcularNovoResiduo(grad->vetor_r0, grad->vetor_r, vetor_y, grad->escalarA, &SLTranspXSL->n);

        //Z0 anterior
        cpyVetor(vetor_z0, grad->vetor_r, &SL->n);

        //Calcula novo z
        calcula_z(vetor_z, matrix_M, grad->vetor_r, &SL->n);

        //Calcula escalar M
        grad->escalarM = multiplicarVtxV(vetor_z, grad->vetor_r, &SL->n)/multiplicarVtxV(vetor_z0, grad->vetor_r0, &SL->n);

        calcularProxDirecao(grad->vetor_d, vetor_z, grad->escalarM, &SL->n);

        tempo_final = timestamp() - tempo_inicial;
        *tempo_metodo += tempo_final;
        ++i;
    }
    *tempo_metodo = *tempo_metodo/i;
    liberarGradientes(grad);
    liberarVetor(vetor_z);
    liberarVetor(vetor_z0);
    liberarVetor(matrix_M);
    

    return 0;
}