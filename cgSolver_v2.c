#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <float.h>
#include "likwid.h"
#include "lib_geral_v2.h"
#include "lib_sislin_v2.h"
#include "lib_gradiente_v2.h"
    //gradienteConjugadosCPreCondicionadores(arq_saida, SL, SLTransp, x, residuo, nInteracoes);
    
    // prnSisLin(stdout, SL);
    // prnSisLin(stdout, SLTransp);

    // SistLinear_t *SLTranspoxSL = calcularMatrizAtxA(SLTransp);
    // multiplicarMatrizPorVetor(SLTransp, SL->b, SLTranspoxSL->b);
    // prnSisLin(stdout, SLTranspoxSL);

    // tornarDiagonalDominante(SLTranspoxSL);
    // prnSisLin(stdout, SLTranspoxSL);

    // if (flagD)
        
        
    // #ifdef SISLIN
    //     FILE *arq = fopen("matriz_inicial_c_pre.txt", "w+");
    //     fprintf(arq, "Sistema antes da aplicação do método\n");
    //     prnSisLin(arq, SL);
    //     fclose(arq);
    // #endif

    

    // #ifdef SISLIN
    //     FILE *arq2 = fopen("matriz_final_c_pre.txt", "w+");
    //     fprintf(arq2, "Sistema depois da aplicação do método\n");
    //     prnSisLin(arq2, SL);
    //     fclose(arq2);
    // #endif
int main(int argc, char *argv[])
{

    real_t erro;                // erro tolerado
    int tamanhoSL;              // Tamanho do Sistema Linear
    int nInteracoes;            // Número de interações
    int k_diagonais;            // Número de diagonais do Sistema
    int pre_condicionador;      // Pré-condicionador
    int *flags;                 // Vetor que indica quais argumentos foram passados;
    char *argumentos[N_FLAGS];  // Valor dos argumentos que foram passados;
    argumentos[4] = NULL;
    char argumentosFaltantes[11] = {'\0'};

    srand(20222);

    flags = (int *)alocarVetor(N_FLAGS, sizeof(int));

    // Pegando as opções que o usuário passa pela linha de comando.
    opterr = 0;
    int opcao;
    while ((opcao = getopt(argc, argv, "n:k:p:i:e:o:d")) != -1)
    {

        //Liga as flags correspondentes para cada opção colocada
        switch (opcao)
        {
        case 'n':
            flags[0] = 1;
            argumentos[0] = optarg;
            break;

        case 'k':
            flags[1] = 1;
            argumentos[1] = optarg;
            break;

        case 'p':
            flags[2] = 1;
            argumentos[2] = optarg;
            break;

        case 'i':
            flags[3] = 1;
            argumentos[3] = optarg;
            break;

        case 'e':
            flags[4] = 1;
            argumentos[4] = optarg;
            break;

        case 'o':
            flags[5] = 1;
            argumentos[5] = optarg;
            break;
        case 'd':
            break;

        default:
            fprintf(stderr, "Está opção não existe\n");
            exit(1);
            break;
        }
    }

    //Verifica se todos os agumentos obrigatórios foram colocados
    if (verificarArgumentos(flags, argumentosFaltantes) < 0)
    {
        fprintf(stderr, "Alguns argumentos obrigatórios ( %s ) não foram especificados\n", argumentosFaltantes);
        exit(1);
    }

    //Valida os argumentos colocados
    if (validarArgumentos(argumentos, &tamanhoSL, &k_diagonais, &nInteracoes, &pre_condicionador, &erro) < 0)
    {
        fprintf(stderr, "Alguns elementos passados não são válidos, encerrando execução.");
        exit(1);
    }

    FILE *arq_saida = fopen(argumentos[5],"w+");
    SistLinear_t *SL;           // Sistema Linear
    SistLinear_t *SLTransp;     // Sistema Linear com os coeficientes transpostos
    real_t *x;                  // Vetor de incógnitas
    real_t *residuo;            // Resíduo
    real_t tempo_metodo;        // Tempo médio para calcular o método
    real_t tempo_preparacao;    // Tempo médio para calcular a preparação para aplicar o método
    real_t tempo_residuo;       // Tempo médio para calcular o resíduo

    fprintf(arq_saida, "# ghs19 Giordano Henrique Silveira\n");
    fprintf(arq_saida, "# grpp19 Gabriel Razzolini Pires de Paula\n\n");

    SL = alocarSisLin(tamanhoSL, k_diagonais, pre_condicionador);
    x = (real_t*)alocarVetor(SL->n, sizeof(real_t));
    residuo = (real_t*)alocarVetor(SL->n, sizeof(real_t));

    initSistLinear(SL);
    tornarDiagonalDominante(SL);
    
    SLTransp = calcularTransposta(SL);


    LIKWID_MARKER_INIT; // Init Likwid markers

    LIKWID_MARKER_START("conj_grad_pre_OTIM"); // Likwid starter conj_grad_pre
    
    gradienteConjugadosCPreCondicionadores(arq_saida, SL, SLTransp, x, &tempo_metodo, &tempo_preparacao, ERRO_IT, nInteracoes);

    LIKWID_MARKER_STOP("conj_grad_pre_OTIM"); // Likwid stopper conj_grad_pre

    LIKWID_MARKER_START("calc_residue_OTIM"); // Likwid starter calc_residue

    calcularResiduo(SL, residuo, x, SL->n);
    real_t norma_residuo = calcularNormaL2Residuo(SL, residuo, &tempo_residuo);

    LIKWID_MARKER_STOP("calc_residue_OTIM"); // Likwid stopper calc_residue

    LIKWID_MARKER_CLOSE; // Shutdown Likwid Markers
    

    fprintf(arq_saida, "# resíduo: %.15g\n", norma_residuo);
    fprintf(arq_saida, "# Tempo PC: %.15g\n", tempo_preparacao);
    fprintf(arq_saida, "# Tempo iter: %.15g\n", tempo_metodo);
    fprintf(arq_saida, "# Tempo norma: %.15g\n", tempo_residuo);
    prnVetor(arq_saida, x, (unsigned int)SL->n);

    //liberar estruturaras
    fclose(arq_saida);
    liberarSisLin(SL);
    liberarSisLin(SLTransp);
    liberarVetor(x);
    liberarVetor(residuo);
    liberarVetor(flags);
    return 0;
}