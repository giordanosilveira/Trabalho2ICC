#!/bin/bash
METRICA="L3 L2CACHE FLOPS_DP" # Metrics for Likwid Markers

SIZE="32 64" # Size of matrix nxn where SIZE=n

# if [ -f 32_L3_v1_.csv ]; then
#     rm *.csv && echo -e ".csv files cleaned" # Cleaning previous log files
# fi

echo -e "\n> Make purge:" && make purge

echo -e "\n> Compiling..." && make all

echo -e "\n> Running FIRST version, cgsolver_v1:\n"

echo "Compilation Result:"
for k in $METRICA; do
     for s in $SIZE; do
         likwid-perfctr -O -C 3 -g ${k} -m ./cgSolver_v1 -n $s -i 150 -k 7 -p 1 -o saida_v1.txt -d diagonal > ${s}_${k}_v1.csv
         if [ -f ${s}_${k}_v1.csv ]; then
            tail -n +7 ${s}_${k}_v1.csv | grep -v "^TABLE" > tmp && mv tmp ${s}_${k}_v1.csv
           echo -e "=> File ${s}_${k}_v1.csv Created Successfully."
        else
            echo -e "X> Failed to create .csv file"
        fi
     done
done

# sed -i '7 s/CFLAGS=-Wall -I${LIKWID_INCLUDE} -DLIKWID_PERFMON/CFLAGS=-Wall -I${LIKWID_INCLUDE} -DLIKWID_PERFMON -O3 -mavx2 -march=native/' makefile # Modify CFLAGS in makefile to add otimization flags

echo -e "\n> Running SECOND version, cgsolver_v2:\n"

echo "Compilation Result:"
for k in $METRICA; do
    for s in $SIZE; do
        likwid-perfctr -O -C 3 -g ${k} -m ./cgSolver_v2 -n $s -i 150 -k 7 -p 1 -o saida_v2.txt -d diagonal > ${s}_${k}_v2.csv
        if [ -f ${s}_${k}_v2.csv ]; then
            tail -n +7 ${s}_${k}_v2.csv | grep -v "^TABLE" > tmp && mv tmp ${s}_${k}_v2.csv
           echo -e "=> File ${s}_${k}_v2.csv Created Successfully."
        else
            echo -e "X> Failed to create .csv file"
        fi
    done
done
echo -e "\nDone" # Script finished

echo -e "\nCalling create_csv_files.py...\n"
python3 create_csv_files.py
