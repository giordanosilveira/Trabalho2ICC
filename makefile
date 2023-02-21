CC=gcc
CFLAGS=-Wall -g
SISLIN=-DSISLIN
PROG=cgSolver_v2
OBJS=lib_geral_v2.o lib_gradiente_v2.o lib_sislin_v2.o $(PROG).o

all: $(PROG)

sislin: CFLAGS += $(SISLIN)
sislin: all

lib_geral_v2.o: lib_geral_v2.c lib_geral_v2.h lib_sislin_v2.h
	$(CC) $(CFLAGS) -c lib_geral_v2.c -lm

lib_gradiente_v2.o: lib_gradiente_v2.c lib_gradiente_v2.h lib_geral_v2.h lib_sislin_v2.h
	$(CC) $(CFLAGS) -c lib_gradiente_v2.c -lm

lib_sislin_v2.o: lib_sislin_v2.c lib_sislin_v2.h
	$(CC) $(CFLAGS) -c lib_sislin_v2.c -lm

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -f *.o *.bak

purge: clean
	rm -f $(PROG)
 
