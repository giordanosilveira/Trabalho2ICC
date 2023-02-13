CC=gcc
CFLAGS=-Wall -g
SISLIN=-DSISLIN
PROG=cgSolver
OBJS=lib_geral.o lib_gradiente.o lib_sislin.o $(PROG).o

all: $(PROG)

sislin: CFLAGS += $(SISLIN)
sislin: all

lib_geral.o: lib_geral.c lib_geral.h lib_sislin.h
	$(CC) $(CFLAGS) -c lib_geral.c -lm

lib_gradiente.o: lib_gradiente.c lib_gradiente.h lib_geral.h lib_sislin.h
	$(CC) $(CFLAGS) -c lib_gradiente.c -lm

lib_sislin.o: lib_sislin.c lib_sislin.h
	$(CC) $(CFLAGS) -c lib_sislin.c -lm

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -f *.o *.bak

purge: clean
	rm -f $(PROG)
 
