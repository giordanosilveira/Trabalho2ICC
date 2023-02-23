#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <float.h>
#include <likwid.h>
#include "lib_geral_v1.h"
#include "lib_sislin_v1.h"
#include "lib_gradiente_v1.h"

int main(int argc, char *argv[])
{
    real_t erro;               // erro tolerado
    int flagD = 0;             // Flag opcional D
    int tamanhoSL;             // Tamanho do Sistema Linear
    int nInteracoes;           // Número de interações
    int k_diagonais;           // Número de diagonais do Sistema
    int pre_condicionador;     // Pré-condicionador
    int *flags;                // Vetor que indica quais argumentos foram passados;
    char *argumentos[N_FLAGS]; // Valor dos argumentos que foram passados;
    argumentos[4] = NULL;
    char argumentosFaltantes[11] = {'\0'};
    double time_passed;        // variable that counts time passed in function
 
    srand(20222);

    flags = (int *)alocarVetor(N_FLAGS, sizeof(int));

    // Pegando as opções que o usuário passa pela linha de comando.
    opterr = 0;
    int opcao;
    while ((opcao = getopt(argc, argv, "n:k:p:i:e:o:d")) != -1)
    {

        // Liga as flags correspondentes para cada opção colocada
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
            flagD = 1;
            break;

        default:
            fprintf(stderr, "Está opção não existe\n");
            exit(1);
            break;
        }
    }

    // Verifica se todos os agumentos obrigatórios foram colocados
    if (verificarArgumentos(flags, argumentosFaltantes) < 0)
    {
        fprintf(stderr, "Alguns argumentos obrigatórios ( %s ) não foram especificados\n", argumentosFaltantes);
        exit(1);
    }

    // Valida os argumentos colocados
    if (validarArgumentos(argumentos, &tamanhoSL, &k_diagonais, &nInteracoes, &pre_condicionador, &erro) < 0)
    {
        fprintf(stderr, "Alguns elementos passados não são válidos, encerrando execução.");
        exit(1);
    }

    /*Variaveis utilizadas pelo método*/
    FILE *arq_saida = fopen(argumentos[5], "w+");
    SistLinear_t *SL;   // Sistema Linear
    real_t *x;          // Vetor de incógnitas
    real_t *residuo;    // Resíduo
    real_t tempo;       // Tempo médio do método
    real_t tempo_norma; // Tempo para o cálculo da norma

    /*Inicialização das variáveis*/
    residuo = (real_t *)alocarVetor(tamanhoSL, sizeof(real_t));
    SL = alocarSisLin(tamanhoSL, k_diagonais, pre_condicionador);
    x = (real_t *)alocarVetor(tamanhoSL, sizeof(real_t));
    initSisLin(SL);

    fprintf(arq_saida, "# ghs19 Giordano Henrique Silveira\n");
    fprintf(arq_saida, "# grpp19 Gabriel Razzolini Pires de Paula\n\n");

    // Verifica se o usuário quer o método sem ou com pré-condicionador
    if (pre_condicionador <= 0)
    {

        if (flagD)
            tornarDiagonalDominante(SL);

#ifdef SISLIN
        FILE *arq = fopen("matriz_inicial_s_pre.txt", "w+");
        fprintf(arq, "Sistema antes da aplicação do método\n");
        prnSisLin(arq, SL);
        fclose(arq);
#endif

        if (gradienteConjugadoSPreCondicionadores(arq_saida, &tempo, SL, x, erro, nInteracoes) < 0)
        {
            fprintf(stderr, "O método dos Gradientes Conjugados sem pré-condicionadores não convergiu para uma solução\n");

            fclose(arq_saida);
            liberarVetor(residuo);
            liberarVetor(x);
            liberarSisLin(SL);
            liberarVetor(flags);

            exit(1);
        }

#ifdef SISLIN
        FILE *arq2 = fopen("matriz_final_s_pre.txt", "w+");
        fprintf(arq2, "Sistema depois da aplicação do método\n");
        prnSisLin(arq2, SL);
        fclose(arq2);
#endif

        calcularResiduo(SL->A, residuo, SL->b, x, SL->n);
        real_t norma_residuo = calcularNormaL2Residuo(SL, residuo, &tempo_norma);
        fprintf(arq_saida, "# resíduo: %.15g\n", norma_residuo);
    }
    else
    {
        real_t tempo_pre_cond; // Tempo da preparação das matrizes para o cálculo do método com pré-condicionadores

        if (flagD)
            tornarDiagonalDominante(SL);

#ifdef SISLIN
        FILE *arq = fopen("matriz_inicial_c_pre.txt", "w+");
        fprintf(arq, "Sistema antes da aplicação do método\n");
        prnSisLin(arq, SL);
        fclose(arq);
#endif

        LIKWID_MARKER_INIT; // Init Likwid markers

        LIKWID_MARKER_START("conj_grad_pre_SEM_OTIM"); // Likwid starter conj_grad_pre

        time_passed = timestamp(); // Init of variable that counts time passed in function

        if (gradienteConjugadosCPreCondicionadores(arq_saida, SL, x, &tempo, &tempo_pre_cond, erro, nInteracoes) < 0)
        {
            fprintf(stderr, "O método dos Gradientes Conjugados com pré-condicionadores não convergiu para uma solução\n");

            fclose(arq_saida);
            liberarVetor(residuo);
            liberarVetor(x);
            liberarSisLin(SL);
            liberarVetor(flags);

            exit(1);
        }

        time_passed += timestamp(); // // Shutdown of variable that counts time passed in function

#ifdef SISLIN
        FILE *arq2 = fopen("matriz_final_c_pre.txt", "w+");
        fprintf(arq2, "Sistema depois da aplicação do método\n");
        prnSisLin(arq2, SL);
        fclose(arq2);
#endif

        LIKWID_MARKER_STOP("conj_grad_pre_SEM_OTIM"); // Likwid stopper conj_grad_pre

        LIKWID_MARKER_START("calc_residue_SEM_OTIM"); // Likwid starter calc_residue

        calcularResiduo(SL->A, residuo, SL->b, x, SL->n);
        real_t norma_residuo = calcularNormaL2Residuo(SL, residuo, &tempo_norma);

        LIKWID_MARKER_STOP("calc_residue_SEM_OTIM"); // Likwid stopper calc_residue

        LIKWID_MARKER_CLOSE; // Shutdown Likwid Markers

        fprintf(arq_saida, "# resíduo: %.15g\n", norma_residuo);
        fprintf(arq_saida, "# Tempo PC: %.15g\n", tempo_pre_cond);
    }

    fprintf(arq_saida, "# Tempo iter: %.15g\n", tempo);
    fprintf(arq_saida, "# Tempo norma: %.15g\n", tempo_norma);
    fprintf(arq_saida, "%d", SL->n);
    prnVetor(arq_saida, x, (unsigned int)SL->n);

    fclose(arq_saida);
    liberarVetor(residuo);
    liberarVetor(x);
    liberarSisLin(SL);
    liberarVetor(flags);

    return 0;
}