Makefile:

    No makefile existe a regra "sislin". Se você compilar utilizando "make sislin", ele vai compilar usando
usando o parâmetro -DSISLIN. Esse parâmetro ativa a macro SISLIN dentro da main. Esta macro serve para
escrever o Sistema Linear antes da aplicação do método e depois da aplicação do método em arquivos diferentes.
    Os arquivos que podem ser criados são:
        * matriz_inicial_s_pre.txt : Matriz antes da aplicação do método sem pré-condicionadores.
        * matriz_final_s_pre.txt : Matriz depois da aplicação do método sem pré-condicionadores.
        * matriz_inicial_c_pre.txt : Matriz antes da aplicação do método com pré-condicionadores.
        * matriz_final_c_pre.txt : Matriz depois da aplicação do método com pré-condicionadores.
O objetivo é que, com o diff, esses arquivos contenham sempres os mesmo sistemas lineares. Tanto porque a matriz
antes do método tem que ser igual a do final do método, tanto porque a matriz gerada tem que ser a mesma, visto que
a seed para a geração dos números aleatórios é a mesma sempre.
    Utilizando "make" você compila somente com as flags "-Wall" e "-g". Uma compilação padrão. Existe as regras
"clean" e "purge" para remover os arquivos objetos e o arquivo executável "cgSolver".


cgSolver:

    Para executar o programa é preciso informar as flags:
        * -n <int> (obrigatório) : Dimensão do sistema linear.
        * -i <int> (obrigatório) : Número max de iterações.
        * -o <arquivo_saida> (obrigatório) : Arquivo onde será escrito a saída formatada.
        * -p <int> (obrigatório) : Se o método é com pré-condicionador (-p <int maior que 1> ), ou sem pré-condicionador (-p 0). 
        * -k <int> (obrigatório) : Número de k diagonais. Somente números ímpares e menor que o número informado em "-n".
        * -e <double> (opcional) : Erro máximo.
        * -d <any> (opcional) : Esse parâmetro não está na especificação original, ele aceita basicamente qualquer coisa, ele só
        precisa estar ligado. Este parâmetro força que os coeficientes do Sistema Linear sejam diagonal dominante. Com isto é possível
        testar o método dos gradientes conjugados com pré-condicionadores. Visto que o pré-condicionador de Jacobe é bom para quando
        a matriz de coeficientes é diagonal dominante.
    
    Exemplo: ./cgSolver -n 500 -k 37 -i 1000 -p 1 -o c_pre_cond.txt -e 1.0e-08 -d 1

Erros:
    
    Caso de erros, número NaN ou Inf, ou caso o método não convirja, aparecera uma mensagem na tela informando o erro, e o programa
    terminara com o código 1.

<arquivo_saida>:

    Este aquivo, que foi informado pela flag "-o", conterá a saída que foi específicada pelos professores.

Arquivos:

    * lib_sislin.h : Contém as definições das estruturas SistLinear_t e Gradiente_t utilizada pelos métodos além dos protótipos das funções. Essas
    funções tem haver mais com o Sistema Linear em si. Funções como alocação, liberação e inicialização do sistema linear, alocação, liberação e inicialização
    da estrutura do gradiente, cálculo de resíduo, entre outras.
    * lib_sislin.c : Contém a implemetação das funções definidas em lib_sislin.h.

    * lib_geral.h : Contém os protótipos das funções que tem um propósito mais geral: alocação e liberação de matriz, alocação e liberação de vetores,
    cálculos de matrizes e vetores. Ou seja, métodos que não necessariamente tem haver com o método.
    * lib_geral.c : Contém a implemetação das funções definidas em lib_geral.h.

    * lib_gradiente.h : Contém os protótipos das funções que tem haver com o cálculo do método.
    * lib_gradiente.c : Contém a implemetação das funções definidas em lib_gradiente.h.

Clareza:

    Tentamos utilizar funções com nomes claros e tentamos comentar o código aonde era preciso para tentar deixar o código mais claro.

