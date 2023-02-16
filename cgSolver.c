#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <float.h>
#include "lib_geral.h"
#include "lib_sislin.h"
#include "lib_gradiente.h"

int main(int argc, char *argv[])
{

    real_t erro;                // erro tolerado
    int flagD = 0;              // Flag opcional D                     
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
            flagD = 1;
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

    SistLinear_t *SL;           // Sistema Linear 
    Coeficientes_t *A_trasp;
    SL = alocarSisLin(tamanhoSL, k_diagonais, pre_condicionador);
    initSistLinear(SL);
    A_trasp = calcularTransposta(SL);
    prnCoef(stderr, SL->A, SL->n, SL->k);
    prnCoef(stderr, A_trasp, SL->n, SL->k);

    /*Variaveis*/
    //FILE *arq_saida = fopen(argumentos[5],"w+");
    // real_t *x;                  // Vetor de incógnitas
    // real_t *residuo;            // Resíduo
    // real_t tempo;               // Tempo médio do método
    // real_t tempo_norma;         // Tempo para o cálculo da norma
    

    // /*Inicialização das variáveis*/
    // residuo = (real_t*)alocarVetor(tamanhoSL, sizeof(real_t));
    // SL = alocarSisLin(tamanhoSL, k_diagonais, pre_condicionador);
    // x = (real_t *)alocarVetor(tamanhoSL, sizeof(real_t));
    // initSisLin(SL);

    // fprintf(arq_saida, "# ghs19 Giordano Henrique Silveira\n");
    // fprintf(arq_saida, "# grpp19 Gabriel Razzolini Pires de Paula\n\n");
    
    // // Verifica se o usuário quer o método sem ou com pré-condicionador
    // if (pre_condicionador <= 0) {

    //     if (flagD)
    //         tornarDiagonalDominante(SL);

    //     #ifdef SISLIN
    //         FILE *arq = fopen("matriz_inicial_s_pre.txt", "w+");
    //         fprintf(arq, "Sistema antes da aplicação do método\n");
    //         prnSisLin(arq, SL);
    //         fclose(arq);
    //     #endif

    //     gradienteConjugadoSPreCondicionadores(arq_saida, &tempo, SL, x, erro, nInteracoes);

    //     #ifdef SISLIN
    //         FILE *arq2 = fopen("matriz_final_s_pre.txt", "w+");
    //         fprintf(arq2, "Sistema depois da aplicação do método\n");
    //         prnSisLin(arq2, SL);
    //         fclose(arq2);
    //     #endif

    //     calcularResiduo(SL->A, residuo, SL->b, x, SL->n);
    //     real_t norma_residuo = calcularNormaL2Residuo(SL, residuo, &tempo_norma);
    //     fprintf(arq_saida, "# resíduo: %.15g\n", norma_residuo);
    // }
    // else {
    //     real_t tempo_pre_cond;  //Tempo da preparação das matrizes para o cálculo do método com pré-condicionadores

    //     if (flagD)
    //         tornarDiagonalDominante(SL);
        
    //     #ifdef SISLIN
    //         FILE *arq = fopen("matriz_inicial_c_pre.txt", "w+");
    //         fprintf(arq, "Sistema antes da aplicação do método\n");
    //         prnSisLin(arq, SL);
    //         fclose(arq);
    //     #endif

    //     gradienteConjugadosCPreCondicionadores(arq_saida, SL, x, &tempo, &tempo_pre_cond, erro, nInteracoes);

    //     #ifdef SISLIN
    //         FILE *arq2 = fopen("matriz_final_c_pre.txt", "w+");
    //         fprintf(arq2, "Sistema depois da aplicação do método\n");
    //         prnSisLin(arq2, SL);
    //         fclose(arq2);
    //     #endif

    //     calcularResiduo(SL->A, residuo, SL->b, x, SL->n);
    //     real_t norma_residuo = calcularNormaL2Residuo(SL, residuo, &tempo_norma);
    //     fprintf(arq_saida, "# resíduo: %.15g\n", norma_residuo);
    //     fprintf(arq_saida, "# Tempo PC: %.15g\n", tempo_pre_cond);
    // }

    // fprintf(arq_saida, "# Tempo iter: %.15g\n", tempo);
    // fprintf(arq_saida, "# Tempo norma: %.15g\n", tempo_norma);
    // fprintf(arq_saida, "%d", SL->n);
    // prnVetor(arq_saida, x, (unsigned int)SL->n);


    // fclose(arq_saida);
    return 0;
}