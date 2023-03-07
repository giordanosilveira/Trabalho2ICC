#ifndef _LIB_GRADIENTE_
#define _LIB_GRADIENTE_


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