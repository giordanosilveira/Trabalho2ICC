matriz_resultado = [[0,0,0],[0,0,0],[0,0,0]]
matriz_resultado_2 = [[0,0,0],[0,0,0],[0,0,0]]
matriz_resultado_3 = [[0,0,0],[0,0,0],[0,0,0]]
matriz_resultado_4 = [[0,0,0],[0,0,0],[0,0,0]]

matriz = [[1, 2, 3],[4,5,6],[7,8,9]]
matriz_transposta = [[1, 4, 7],[2,5,8],[3,6,9]]

print(matriz)

print(matriz_transposta)

for i in range(0, 3):
    for j in range(0, 3):
        for k in range(0, 3):
            matriz_resultado[i][j] = matriz_resultado[i][j] + matriz[i][k]*matriz[j][k]

for i in range(0, 3):
    for j in range(0, 3):
        for k in range(0, 3):
            matriz_resultado_2[i][j] = matriz_resultado_2[i][j] + matriz[i][k]*matriz_transposta[k][j]

for i in range(0, 3):
    for j in range(0, 3):
        for k in range(0, 3):
            matriz_resultado_3[i][j] = matriz_resultado_3[i][j] + matriz_transposta[i][k]*matriz[k][j]

for i in range(0, 3):
    for j in range(0, 3):
        for k in range(0, 3):
            matriz_resultado_4[i][j] = matriz_resultado_4[i][j] + matriz[k][i]*matriz[k][j]


    
print(matriz_resultado)
print(matriz_resultado_2)
print(matriz_resultado_3)
print(matriz_resultado_4)