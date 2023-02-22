import pandas as pd
import subprocess

# grupos = {"L3": "\nL3 bandwidth [MBytes/s]", "L2CACHE": "\nL2 miss ratio", "FLOPS_DP": "\nDP MFLOP/s", "\nAVX DP MFLOP/s"}
groups = ["L3", "L2CACHE", "FLOPS_DP"]
keys_groups = ["L3 bandwidth [MBytes/s]", "L2 miss ratio", "DP MFLOP/s", "AVX DP MFLOP/s"]
sizes = [32, 64]
versions = ["v1", "v2"]
# tamanhos = [64, 100, 128, 1024, 2000, 2048]  # sizes of matrix nxn, where size=n

print(f"Pegando valor {keys_groups[0]} do arquivo {sizes[0]}_{groups[0]}_{versions[0]}.csv")
df = pd.read_csv(f"32_L3_v1.csv")
print("Arquivo Lido")
print(df)
valor_1 = df.loc[df['STRUCT'] == keys_groups[0]].values[0] # primeiro parametro
valor_2 = df.loc[df['STRUCT'] == keys_groups[0]].values[1] # segundo parametro
novo_df = pd.DataFrame({ 'parametro': ["conj_grad"], 'valor': [valor_1]}
                       {'parametro': ["resiude"], 'valor': [valor_2]})
novo_df.to_csv('saida.csv', index=False)




# lista = list()

# print("Before the second loop of sizes")
# for tam in sizes:
#     with open(f"{tam}_{grupo[0]}.csv", "r") as arq:
#         f = arq.read()
#         f = f.split(',')
#         lista.append([f[i + 1] for i, elemento in enumerate(f) if elemento == "\nRDTSC Runtime [s]"])
#     lista[len(lista) - 1].insert(0, f"{tam}")
#     with open(f"TIME.csv", 'w') as arq2:
#         arq2.write("Tamanho, conj_grad_pre, calc_residue\n")
#         for m in lista:
#             arq2.write(f"{','.join(m)}\n")
# print("After the second loop of sizes")

# bashCommand = "rm *WITH*.csv" # Bash command
# process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
# output, error = process.communicate()
# print("Output: ", output)
# print("Error: ", error)

# import matplotlib.pyplot as plt
# Generate plots
# arquivos = [("FLOPS_DP", "DP MFLOPS/s"), ("TIME", "TEMPO"), ("L3", "L3 bandwidth [MBytes/s]"), ("L2CACHE", "L2 miss ratio")]
# for arq in arquivos:
#     df = pd.read_csv(f"{arq[0]}.csv")
#     plt.figure(figsize=(10, 3.7))
#     plt.plot(df["Tamanho"], df["t_matmat"], 'r', label="ColunaMatrizxMatriz")
#     plt.plot(df["Tamanho"], df["t_matvet"], 'g', label="ColunaMatrizxVetor")
#     plt.xlabel('Tamanho')
#     plt.ylabel(f'{arq[1]}')
#     plt.legend(loc='upper left', frameon=False)
#     plt.title(f"{arq[0]}", loc='left')
#     plt.savefig(f'{arq[0]}.png')
