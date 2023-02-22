CC=gcc -g

PROGV1=cgSolver_v1

PROGV2=cgSolver_v2

CFLAGS=-Wall -I${LIKWID_INCLUDE} -DLIKWID_PERFMON -O3 -mavx2 -march=native
LFLAGS=-lm
LIKWIDFLAGS = -L${LIKWID_LIB} -llikwid

SISLIN=-DSISLIN

OBJSV1=lib_geral_v1.o lib_gradiente_v1.o lib_sislin_v1.o $(PROGV1).o
OBJSV2=lib_geral_v2.o lib_gradiente_v2.o lib_sislin_v2.o $(PROGV2).o

all: $(PROGV1) $(PROGV2)

sislin: CFLAGS += $(SISLIN)
sislin: all

lib_geral_v1.o: lib_geral_v1.c lib_geral_v1.h lib_sislin_v1.h
	$(CC) $(CFLAGS) -c lib_geral_v1.c -lm

lib_gradiente_v1.o: lib_gradiente_v1.c lib_gradiente_v1.h lib_geral_v1.h lib_sislin_v1.h
	$(CC) $(CFLAGS) -c lib_gradiente_v1.c -lm

lib_sislin_v1.o: lib_sislin_v1.c lib_sislin_v1.h
	$(CC) $(CFLAGS) -c lib_sislin_v1.c -lm

$(PROGV1): $(OBJSV1)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS) $(LIKWIDFLAGS)

lib_geral_v2.o: lib_geral_v2.c lib_geral_v2.h lib_sislin_v2.h
	$(CC) $(CFLAGS) -c lib_geral_v2.c -lm

lib_gradiente_v2.o: lib_gradiente_v2.c lib_gradiente_v2.h lib_geral_v2.h lib_sislin_v2.h
	$(CC) $(CFLAGS) -c lib_gradiente_v2.c -lm

lib_sislin_v2.o: lib_sislin_v2.c lib_sislin_v2.h
	$(CC) $(CFLAGS) -c lib_sislin_v2.c -lm

$(PROGV2): $(OBJSV2)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS) $(LIKWIDFLAGS)

clean:
	rm -f *.o *.bak *.csv *.tmp

purge: clean
	rm -f $(PROGV1) $(PROGV2)
 
