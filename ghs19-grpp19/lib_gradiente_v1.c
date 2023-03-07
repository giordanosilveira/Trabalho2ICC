#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib_geral_v1.h"
#include "lib_gradiente_v1.h"
#include "lib_sislin_v1.h"



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
static void calcularNovaSolucao(real_t *x0, real_t *x, real_t *d, real_t escalarA, unsigned int n)
{

    // Percorre o vetor de soluções e realiza o calculo Xk' = Xk + s*Dk.
    for (int i = 0; i < n; ++i)
    {
        x0[i] = x[i];
        x[i] = x[i] + escalarA * d[i];
    }
}

/**
 * @brief Calcula o novo resíduo. Ou seja, Rk' = Rk - s*A*D.
 * Também guarda o resíduo antigo em 'R0'.
 *
 * @param r0 (real_t*) : Guarda o resíduo antigo.
 * @param r (real_t*) : Guarda o resíduo atual.
 * @param d (real_t*) : Próxima direção (chamado de 'v' no livro da M.Cristina C. Cunha).
 * @param escalarA (real_t) : Escalar A (chamado de 's' no livro da M.Cristina C. Cunha).
 * @param SL (SistLinear_t *) : O Sistema Linear.
 */
static void calcularNovoResiduo(real_t *r0, real_t *r, real_t *d, real_t escalarA, SistLinear_t *SL)
{

    // Percorre a matriz dos coeficientes
    for (int i = 0; i < SL->n; ++i)
    {
        r0[i] = r[i];
        real_t soma = 0.0;
        for (int j = 0; j < SL->n; ++j)
        {

            // calcula o i-ésimo elemento do vetor apontado como 'z' no livro da M.Cristina C. Cunha.
            soma = soma + SL->A[i][j] * d[j];
        }

        // o i-ésimo elemento de R recebe o i-ésimo elemento de 'z' vezes o escalar 'a'.
        r[i] = r[i] - escalarA * soma;
    }
}

/**
 * @brief Calcula a próxima direção. Ou seja, Dk' = escalarM*Dk + Vk.
 * Lembrando que o vetor 'v', dependendo do método, pode mudar.
 *
 * @param d (real_t*) : Guarda a direção após o cálculo.
 * @param v (real_t*) : Vetor que pode ser o resíduo 'r' (S/ pré-condicionador) ou vetor 'z' (C/ pré-condicionador).
 * @param escalarM (real_t) : Escalar M.
 * @param n (unsigned int) : Tamanho dos vetores.
 */
static void calcularProxDirecao(real_t *d, real_t *v, real_t escalarM, unsigned int n)
{

    for (int i = 0; i < n; ++i)
    {
        d[i] = escalarM * d[i] + v[i];
    }
}

real_t calcularNormaL2Residuo(SistLinear_t *SL, real_t *x, real_t *tempo)
{

    real_t soma = 0.0;
    real_t raiz;

    // Pecorre o vetor de soluções
    *tempo = timestamp();
    for (int i = 0; i < SL->n; ++i)
    {

        // Soma o quadrado dos elementos das soluções
        soma = soma + x[i] * x[i];

        // Teste para ver se não foi gerado um NaN ou um número infinito
        if (isnan(soma) || isinf(soma))
        {
            fprintf(stderr, "Erro soma(calcularNormaL2Residuo): %g é NaN ou +/-Infinito\n", soma);
            exit(1);
        }
    }
    raiz = sqrt(soma);
    *tempo = timestamp() - *tempo;

    // Retorna a raíz quadrada da soma.
    return raiz;
}

real_t normaMaxErroRelativo(real_t *x, real_t *x0, real_t *maior_erro_max_abs, unsigned int *tam)
{

    real_t maior = ABS(x[0] - x0[0]) / ABS(x[0]);
    *maior_erro_max_abs = ABS(x[0] - x0[0]);

    // Percorre o vetor de soluções,
    for (int i = 1; i < *(tam); ++i)
    {
        // acha o maior erro relativo e
        if (ABS(x[i] - x0[i]) / ABS(x[i]) > maior)
        {
            maior = ABS(x[i] - x0[i]) / ABS(x[i]);
            // Teste para ver se não foi gerado um NaN ou um número infinito
            if (isnan(maior) || isinf(maior))
            {
                fprintf(stderr, "Erro maior(normaMaxErroRelativo): %g é NaN ou +/-Infinito\n", maior);
                exit(1);
            }
        }
        // o maior erro absoluto.
        if (ABS(x[i] - x0[i]) > *maior_erro_max_abs)
        {
            *maior_erro_max_abs = ABS(x[i] - x0[i]);
        }
    }

    // Retorna ao final o maior erro absoluto.
    return maior;
}

real_t calcularDenominadorEscalarA(real_t *D, SistLinear_t *SL)
{

    real_t *vetor_resultante = (real_t *)alocarVetor(SL->n, sizeof(real_t));

    // Percorre a matriz SL->A e o vetor D
    for (int i = 0; i < SL->n; ++i)
    {
        real_t soma = 0.0;
        for (int j = 0; j < SL->n; ++j)
        {
            // calcula o iésimo elemento do vetor_resultante (esse vetor é chamado de 'z' no livro M.Cristina C. Cunha)
            soma = soma + SL->A[i][j] * D[j];

            // Teste para ver se não foi gerado um NaN ou um número infinito
            if (isnan(soma) || isinf(soma))
            {
                fprintf(stderr, "Erro soma(calcularDenominadorEscalarA): %g é NaN ou +/-Infinito\n", soma);
                exit(1);
            }
        }
        // O iésimo elemento do vetor resultante recebe soma.
        vetor_resultante[i] = soma;
    }

    // Tendo o vetor resultante, agora é só multiplicar D^t * vetor_resultante.
    // O cálculo entre um vetor transposto e um vetor normal gera uma matriz de 1x1, ou seja, um número real.
    real_t soma = 0.0;
    for (int i = 0; i < SL->n; ++i)
    {
        soma = soma + D[i] * vetor_resultante[i];
        // Teste para ver se não foi gerado um NaN ou um número infinito
        if (isnan(soma) || isinf(soma))
        {
            fprintf(stderr, "Erro soma(calcularDenominadorEscalarA): %g é NaN ou +/-Infinito\n", soma);
            exit(1);
        }
    }
    liberarVetor(vetor_resultante);
    return soma;
}

void inverse_jacobi_preconditioner(SistLinear_t *SL, real_t *M)
{

    // Percorre a diagonal da matriz SL->A
    for (int i = 0; i < SL->n; ++i)
    {

        // O i-ésimo elemento de M recebe o i-ésimo elemento da diagonal princial de SL->A
        M[i] = (real_t)1 / SL->A[i][i];
        if (isnan(M[i]) || isinf(M[i]))
        {
            fprintf(stderr, "Erro M[i](inverse_jacobi_preconditioner): %g é NaN ou +/-Infinito, Linha: %i\n", M[i], i);
            exit(1);
        }
    }
}

void calcula_z(real_t *z, real_t *inverse_c, real_t *residue, unsigned int size)
{

    // Percorre o vetor z
    for (int i = 0; i < size; ++i)
    {

        // o Zi = recebe M-¹ * Ri
        z[i] = inverse_c[i] * residue[i];

        // Teste para ver se não foi gerado um NaN ou um número infinito
        if (isnan(z[i]) || isinf(z[i]))
        {
            fprintf(stderr, "Erro z[i](calcula_z): %g é NaN ou +/-Infinito, Linha: %i\n", z[i], i);
            exit(1);
        }
    }
}

int gradienteConjugadoSPreCondicionadores(FILE *arq_saida, real_t *tempo, SistLinear_t *SL, real_t *x, real_t erro, real_t nInteracoes)
{
    Gradiente_t *grad = alocarGradiente(SL->n);
    if (!grad)
    {
        fprintf(stderr, "Não for possível alocar espaço para a estrutura do gradiente\n");
        exit(1);
    }

    /*Declaração das variáveis utilizadas internamente pela função*/
    real_t max_erro_abs_aprox;
    real_t tempo_inicial, tempo_final;

    /*Transformação do Sistema Ax = b em (A^t)Ax = (A^t)b */
    cpyMatriz(grad->aux_A, SL->A, (unsigned int)SL->n);   // Copia os Coeficientes Originais
    cpyVetor(grad->aux_b, SL->b, (unsigned int *)&SL->n); // Copia os Termos independentes originais
    real_t **coef_transposto = calcularMatrizAtxA(SL);    //(A^t)*A.
    real_t *novo_b = calcularAtxB(SL);                    //(A^t)*b.
    cpyVetor(SL->b, novo_b, &SL->n);                      // Copia os novos termos independentes.
    cpyMatriz(SL->A, coef_transposto, SL->n);             // Copia os novos coeficientes.

    /*Inicialização de variáveis usadas pelo programa*/
    cpyVetor(grad->vetor_r, SL->b, &SL->n);         // Inicializa resíduo com SL->b (O X inicial é 0).
    cpyVetor(grad->vetor_d, grad->vetor_r, &SL->n); // Inicializa o vetor da próxima direção.

    (*tempo) = 0;
    int i = 0;
    while (i < nInteracoes)
    {

        // Tempo inicial do começo do método
        tempo_inicial = timestamp();

        // A variável A é calculada
        grad->escalarA = multiplicarVtxV(grad->vetor_r, NULL, &SL->n) / calcularDenominadorEscalarA(grad->vetor_d, SL); // Div 0

        calcularNovaSolucao(grad->vetor_x0, x, grad->vetor_d, grad->escalarA, SL->n);

        calcularNovoResiduo(grad->vetor_r0, grad->vetor_r, grad->vetor_d, grad->escalarA, SL);

        // Verifica se o erro relativo é menor que o informado, se for, destroi as estruturas internas
        //  e termina o método.
        double max_erro_rel = normaMaxErroRelativo(x, grad->vetor_x0, &max_erro_abs_aprox, &SL->n);
        fprintf(arq_saida, "# iter %d: ||%.15g||\n", i, max_erro_abs_aprox);
        if (max_erro_rel < erro)
        {
            // Copia os valores originais dos coeficientes e dos termos independentes do Sistema Linear.
            cpyMatriz(SL->A, grad->aux_A, (unsigned int)SL->n);
            cpyVetor(SL->b, grad->aux_b, (unsigned int *)&SL->n);

            // Libera as estruturas
            liberarGradientes(grad);
            liberarMatriz(coef_transposto);
            liberarVetor(novo_b);

            // Cálcula o tempo final.
            tempo_final = timestamp() - tempo_inicial;
            (*tempo) = ((*tempo) + tempo_final) / i;
            return 0;
        }

        // A variável M é calculada
        grad->escalarM = multiplicarVtxV(grad->vetor_r, NULL, &SL->n) / multiplicarVtxV(grad->vetor_r0, NULL, &SL->n);

        calcularProxDirecao(grad->vetor_d, grad->vetor_r, grad->escalarM, SL->n);

        // calcula o tempo final do método
        tempo_final = timestamp() - tempo_inicial;
        *tempo += tempo_final;
        ++i;
    }

    // Calcula o tempo final do método
    (*tempo) /= i;

    // Copia os valores originais dos coeficientes e dos termos independentes do Sistema Linear.
    cpyMatriz(SL->A, grad->aux_A, (unsigned int)SL->n);
    cpyVetor(SL->b, grad->aux_b, (unsigned int *)&SL->n);

    // Libera as estruturas
    liberarGradientes(grad);
    liberarMatriz(coef_transposto);
    liberarVetor(novo_b);
    
    // Se o método não convergiu, verifica se o erro foi informado, caso sim, o método
    // realmente não convergiu, ou seja, retorna um erro.
    if (erro != -1)
        return -1;

    // Se o método não convergiu, porém o erro não foi informado, o método foi testado usando o
    // número de iterações como parada. Isso não é um erro, portanto retorna 0.
    return 0;
}

int gradienteConjugadosCPreCondicionadores(FILE *arq_saida, SistLinear_t *SL, real_t *x, real_t *tempo, real_t *tempo_pre_cond, real_t erro, real_t nInteracoes)
{
    Gradiente_t *grad = alocarGradiente(SL->n);
    if (!grad)
    {
        fprintf(stderr, "Não foi possível alocar espaço para as estruturas do gradiente\n");
        exit(1);
    }

    /*Declaração das variáveis utilizadas internamente pela função*/
    real_t *matrix_M = (real_t *)alocarVetor(SL->n, sizeof(real_t)); // Pré-condicionador de Jacobi M = D
    real_t *vetor_z = (real_t *)alocarVetor(SL->n, sizeof(real_t));  // Chamado de 'y' no livro M.Cristina C. Cunha
    real_t *vetor_z0 = (real_t *)alocarVetor(SL->n, sizeof(real_t)); // Z anterior.
    real_t max_erro_abs_aprox;
    real_t tempo_inicial, tempo_final;

    /*Inicialização das variáveis*/
    cpyMatriz(grad->aux_A, SL->A, SL->n); // Inicializa com os coeficientes originais
    cpyVetor(grad->aux_b, SL->b, &SL->n); // Inicializa com os termos independentes originais

    // Tempo inicial para o cálculo do tempo para a preparação do método
    (*tempo_pre_cond) = timestamp();

    /*Inicializando as variável que serão utilizadas pelo método*/
    inverse_jacobi_preconditioner(SL, matrix_M);   // Inicializa a "matriz" 'M'
    real_t **novo_coef_A = calcularMatrizAtxA(SL); //(A^t)*A.
    real_t *novo_b = calcularAtxB(SL);             //(A^t)*b.
    cpyVetor(SL->b, novo_b, &SL->n);               // Copia os novos termos independentes.
    cpyMatriz(SL->A, novo_coef_A, SL->n);          // Copia os novos coeficientes.
    calcularMxA(SL->A, matrix_M, SL->n);           // Calcular M^-1 * A
    calcularMxb(SL->b, matrix_M, SL->n);           // Calcular M^-1 * b

    // Tempo final da preparação do método
    (*tempo_pre_cond) = timestamp() - (*tempo_pre_cond);

    cpyVetor(grad->vetor_r, SL->b, &SL->n);             // Inicializa o resíduo 'r0' com SL->b (O X inicial é 0)
    calcula_z(vetor_z, matrix_M, grad->vetor_r, SL->n); // Inicializa 'z0'
    cpyVetor(grad->vetor_d, vetor_z, &SL->n);           // Inicializa 'd0'

    int i = 0;
    (*tempo) = 0;
    while (i < nInteracoes)
    {

        // Tempo inicial do começo do método
        tempo_inicial = timestamp();

        // A variável A é calculada
        grad->escalarA = multiplicarVtxV(vetor_z, grad->vetor_r, &SL->n) / calcularDenominadorEscalarA(grad->vetor_d, SL);

        calcularNovaSolucao(grad->vetor_x0, x, grad->vetor_d, grad->escalarA, SL->n);

        calcularNovoResiduo(grad->vetor_r0, grad->vetor_r, grad->vetor_d, grad->escalarA, SL);

        // Verifica se o erro relativo é menor que o informado, se for, destroi as estruturas internas
        //  e termina o método.
        double max_erro_rel = normaMaxErroRelativo(x, grad->vetor_x0, &max_erro_abs_aprox, &SL->n);
        fprintf(arq_saida, "# iter %d: ||%.15g||\n", i, max_erro_abs_aprox);
        if (max_erro_rel < erro)
        {

            // Copia os valores originais dos coeficientes e dos termos independentes do Sistema Linear.
            cpyMatriz(SL->A, grad->aux_A, (unsigned int)SL->n);
            cpyVetor(SL->b, grad->aux_b, (unsigned int *)&SL->n);

            // Libera as estruturas
            liberarGradientes(grad);
            liberarMatriz(novo_coef_A);
            liberarVetor(novo_b);
            liberarVetor(vetor_z);
            liberarVetor(vetor_z0);
            liberarVetor(matrix_M);

            // Cálcula o tempo final.
            tempo_final = timestamp() - tempo_inicial;
            (*tempo) = ((*tempo) + tempo_final) / i;

            return 0;
        }

        // Salva o vetor z anterior em z0 e calcula o próximo z.
        cpyVetor(vetor_z0, grad->vetor_r, &SL->n);
        calcula_z(vetor_z, matrix_M, grad->vetor_r, SL->n);

        // A variável M é calculada
        grad->escalarM = multiplicarVtxV(vetor_z, grad->vetor_r, &SL->n) / multiplicarVtxV(vetor_z0, grad->vetor_r0, &SL->n);

        calcularProxDirecao(grad->vetor_d, vetor_z, grad->escalarM, SL->n);

        // calcula o tempo final do método
        tempo_final = timestamp() - tempo_inicial;
        *tempo += tempo_final;
        ++i;
    }
    *tempo = *tempo / i;

    // Copia os valores originais dos coeficientes e dos termos independentes do Sistema Linear.
    cpyMatriz(SL->A, grad->aux_A, (unsigned int)SL->n);
    cpyVetor(SL->b, grad->aux_b, (unsigned int *)&SL->n);

    // Libera as estruturas
    liberarGradientes(grad);
    liberarMatriz(novo_coef_A);
    liberarVetor(novo_b);
    liberarVetor(vetor_z);
    liberarVetor(vetor_z0);
    liberarVetor(matrix_M);

    // Se o método não convergiu, verifica se o erro foi informado, caso sim, o método
    // realmente não convergiu, ou seja, retorna um erro.
    if (erro != -1)
        return -1;

    // Se o método não convergiu, porém o erro não foi informado, o método foi testado usando o
    // número de iterações como parada. Isso não é um erro, portanto retorna 0.
    return 0;
}
