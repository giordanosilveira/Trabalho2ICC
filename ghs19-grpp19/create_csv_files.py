import pandas as pd
import numpy as np
import subprocess

groups = ["L3", "L2CACHE", "FLOPS_DP", "FLOPS_AVX"]
params = ["L3 bandwidth [MBytes/s]", "L2 miss ratio", "DP MFLOP/s", "Packed DP MFLOP/s"]
sizes = [32, 64, 128, 256, 512, 1000, 2000, 4000, 6000]
versions = ["v1", "v2"]

print("==> Starting loop.. for L3 and L2CACHE\n")
for group in groups:
    conj_grad_SEM_OTIM = []
    residue_SEM_OTIM = []
    conj_grad_COM_OTIM = []
    residue_COM_OTIM = []
    for size in sizes:
        if group == "L3":
            param = params[0]
        elif group == "L2CACHE":
            param = params[1]
        elif group == "FLOPS_DP":
            param = params[2]
        else:
            param = params[3]
        print(f"==> Group: {group}, Param: {param}.")              
        
        print(f"==> Pegando valor {param} do arquivo {size}_{group}_{versions[0]}.csv")
        df = pd.read_csv(f"{size}_{group}_{versions[0]}.csv")  # Read Pandas File, VERSION 1, not optimized
        print(f"==> File {size}_{group}_{versions[0]}.csv read")
        
        valor_1 = df.loc[ df['STRUCT'] == param ].values[0]  # primeiro parametro, conj_grad
        conj_grad_SEM_OTIM.append(valor_1[1]) # Extract Desired value from numpy array
        print(f"==> List conj_grad_SEM_OTIM: {conj_grad_SEM_OTIM}")
        
        valor_2 = df.loc[ df['STRUCT'] == param ].values[1]  # segundo parametro, residue
        residue_SEM_OTIM.append(valor_2[1])  # Extract Desired value from numpy array
        print(f"==> List residue_SEM_OTIM: {residue_SEM_OTIM}")
        
        print(f"==> Pegando valor {param} do arquivo {size}_{group}_{versions[1]}.csv")
        df2 = pd.read_csv(f"{size}_{group}_{versions[1]}.csv")  # Read Pandas File, VERSION 2, optimized
        print(f"==> File {size}_{group}_{versions[1]}.csv read")
        
        valor_3 = df2.loc[df2['STRUCT'] == param ].values[0]  # primeiro parametro, conj_grad
        conj_grad_COM_OTIM.append(valor_3[1])# Extract Desired value from numpy array
        print(f"==> List conj_grad_COM_OTIM: {conj_grad_COM_OTIM}")
        
        valor_4 = df2.loc[df2['STRUCT'] == param ].values[1]  # segundo parametro, residue
        residue_COM_OTIM.append(valor_4[1])  # Extract Desired value from numpy array
        print(f"==> List residue_COM_OTIM: {residue_COM_OTIM}")

    print(f"\n==> Writing File {group}.csv")
    novo_df = pd.DataFrame({'SIZES': sizes, 'conj_grad_SEM_OTIM': conj_grad_SEM_OTIM, 'residue_SEM_OTIM': residue_SEM_OTIM, 'conj_grad_COM_OTIM': conj_grad_COM_OTIM, 'residue_COM_OTIM': residue_COM_OTIM})
    novo_df.to_csv(f"{group}.csv", index=False)
    print(f"==> File {group}.csv Created Successfully\n")


# Time
conj_grad_SEM_OTIM = []
residue_SEM_OTIM = []
conj_grad_COM_OTIM = []
residue_COM_OTIM = []  
for size in sizes:
    print(f"==> Pegando valor 'Runtime (RDTSC) [s]' do arquivo {size}_L3_{versions[0]}.csv")
    df = pd.read_csv(f"{size}_L3_{versions[0]}.csv") # Example file
    print(f"==> File {size}_L3_{versions[0]}.csv read")
        
    valor_1 = df.loc[ df['STRUCT'] == "Runtime (RDTSC) [s]" ].values[0]  # primeiro parametro, conj_grad
    conj_grad_SEM_OTIM.append(valor_1[1]) # Extract Desired value from numpy array
    print(f"==> List conj_grad_SEM_OTIM: {conj_grad_SEM_OTIM}")
    
    valor_2 = df.loc[ df['STRUCT'] == "Runtime (RDTSC) [s]" ].values[1]  # segundo parametro, residue
    residue_SEM_OTIM.append(valor_2[1])  # Extract Desired value from numpy array
    print(f"==> List residue_SEM_OTIM: {residue_SEM_OTIM}")
    
    print(f"==> Pegando valor 'Runtime (RDTSC) [s]' do arquivo {size}_L3_{versions[1]}.csv")
    df2 = pd.read_csv(f"{size}_L3_{versions[1]}.csv")  # Read Pandas File, VERSION 2, optimized
    print(f"==> File {size}_L3_{versions[1]}.csv read")
    
    valor_3 = df2.loc[df2['STRUCT'] == "Runtime (RDTSC) [s]" ].values[0]  # primeiro parametro, conj_grad
    conj_grad_COM_OTIM.append(valor_3[1])# Extract Desired value from numpy array
    print(f"==> List conj_grad_COM_OTIM: {conj_grad_COM_OTIM}")
    
    valor_4 = df2.loc[df2['STRUCT'] == "Runtime (RDTSC) [s]" ].values[1]  # segundo parametro, residue
    residue_COM_OTIM.append(valor_4[1])  # Extract Desired value from numpy array
    print(f"==> List residue_COM_OTIM: {residue_COM_OTIM}")

# Writing File
print(f"\n==> Writing File TIME.csv")
novo_df = pd.DataFrame({'SIZES': sizes, 'conj_grad_SEM_OTIM': conj_grad_SEM_OTIM, 'residue_SEM_OTIM': residue_SEM_OTIM, 'conj_grad_COM_OTIM': conj_grad_COM_OTIM, 'residue_COM_OTIM': residue_COM_OTIM})
novo_df.to_csv(f"TIME.csv", index=False)
print(f"==> File TIME.csv Created Successfully\n")


import matplotlib.pyplot as plt
# Generate plots
arquivos = [("FLOPS_DP", "DP MFLOPS/s"), ("TIME", "TEMPO"), ("L3", "L3 bandwidth [MBytes/s]"), ("L2CACHE", "L2 miss ratio"), ("FLOPS_AVX", "Packed DP MFLOP/s")]
for arq in arquivos:
    df = pd.read_csv(f"{arq[0]}.csv")
    plt.figure(figsize=(10, 3.7))
    plt.plot(df['SIZES'], df['conj_grad_SEM_OTIM'], 'r', label="conj_grad_SEM_OTIM")
    plt.plot(df['SIZES'], df['residue_SEM_OTIM'], 'g', label="residue_SEM_OTIM")
    plt.plot(df['SIZES'], df['conj_grad_COM_OTIM'], 'b', label="conj_grad_COM_OTIM")
    plt.plot(df['SIZES'], df['residue_COM_OTIM'], 'y', label="residue_COM_OTIM")
    plt.xlabel('SIZE')
    plt.ylabel(f'{arq[1]}')
    plt.yscale('log')
    plt.legend(loc='upper left', frameon=False)
    plt.title(f"{arq[0]}", loc='left')
    plt.savefig(f'{arq[0]}.png')
