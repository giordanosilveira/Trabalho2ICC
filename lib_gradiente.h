#ifndef _LIB_GRADIENTE_
#define _LIB_GRADIENTE_

// /**
//  * @brief Calcula a norma L2 do vetor com as soluções 'x'
//  * 
//  * @param tempo (real_t*) : Ponteiro onde será guardado o tempo do cálculo da norma.
//  * @param SL (SistLinear_t*) : O sistema linear.
//  * @param x (real_t*) : O Vetor com as soluções
//  * @return (real_t) : Norma L2 do vetor. 
//  */
// real_t calcularNormaL2Residuo(SistLinear_t *SL, real_t *x, real_t *tempo);


// /**
//  * @brief Calcula a norma máxima do erro relativo.
//  *
//  * @param x (real_t*) : Vetor de incógnitas com a solução atual.
//  * @param x0 (real_t*) : Vetor de incógnitas com a solução anterior.
//  * @param tam (int*) : Ponteiro que indica o tamanho do vetor.
//  * @return (real_t) : Norma máxima do erro absoluto.
//  */
// real_t normaMaxErroRelativo(real_t *x, real_t *x0, real_t *maior_erro_max_abs, unsigned int *tam);


// /**
//  * @brief Calcula o denominador para o cálcula do alfa
//  *
//  * @param proxDirecao (real_t*) : Vetor que índica a proxima direção do método.
//  * @param SL (SistLinear_t*) : O sistema linear.
//  * @return (real_t) : Denominador para o cálculo do alfa.
//  */
// real_t calcularDenominadorEscalarA(real_t *proxDirecao, SistLinear_t *SL);


// /**
//  * @brief Calcula a inversa do pré-condicionador de Jacobi, isto é: M⁻¹
//  * 
//  * @link https://math-linux.com/mathematics/linear-systems/article/preconditioned-conjugate-gradient-method
//  * 
//  * @param SL (SistLinear_t*) : O sistema linear.
//  * @param M (real_t*) : Vetor contendo a diagonal principal do Sistema Linear original
//  * 
//  * @return (real_t) : Matriz inversa de C, C⁻¹.
//  */
// void inverse_jacobi_preconditioner(SistLinear_t *SL, real_t *M);



// void calcula_z(real_t *z, real_t *inverse_c, real_t *residue, unsigned int size);


// /**
//  * @brief Resolve o sistema linear pelo método dos Gradientes Conjugados.
//  *
//  * @param arq_saida (FILE*) : Arquivo onde conterá a saída formatada pedida pelos professores.
//  * @param tempo (real_t*) : Ponteiro onde será guardado o tempo de uma iteração.
//  * @param SL (SistLinear_t*) : O sistema linear.
//  * @param x (real_t*) : Vetor onde ficará a solução do sistema linear.
//  * @param erro (real_t*) : Ponteiro que representa o erro máximo.
//  * @param nInteracoes (int) : Número de interações máximas.
//  * @return (int) : 0 se ocorreu tudo certo; -1 se não ocorreu.
//  */
// int gradienteConjugadoSPreCondicionadores(FILE* arq_saida, real_t *tempo, SistLinear_t *SL, real_t *x, real_t erro, real_t nInteracoes);


/**
 * @brief Calcula o método com pré-condicionadores.
 * 
 * @param arq_saida (FILE*) : Arquivo para printar a saída.
 * @param SL (SistLinear_t*) : Sistema Linear original.
 * @param SLTransposto (SistLinear_t*) : Sistema Linear com os coeficientes transposto.
 * @param x (real_t*) : Vetor para ser colocado a solução.
 * @param tempo_metodo (real_t*) : Aonde vai ser colocado o tempo médio total do método.
 * @param tempo_preparacao (real_t*) : Aonde vai ser colocado o tempo total para a preparação do método
 * @param erro (real_t) : Erro.
 * @param nInteracoes (int) : Número de iterações.
 * @return int 
 */
int gradienteConjugadosCPreCondicionadores(FILE*arq_saida, SistLinear_t *SL, SistLinear_t *SLTransposto, real_t *x, real_t* tempo_metodo, real_t *tempo_preparacao, real_t erro, int nInteracoes);

#endif