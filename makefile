SRC=filter_test.c filter.c polynom.c dumper.c
OBJ=$(SRC:.c=.o)
EXE=filter_test

CC=gcc
CFLAGS=-Wall -O3
#LDFLAGS=-static-libgcc -static-libstdc++ -static
LDFLAGS=
RM=rm

all: $(EXE)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

$(EXE): $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $@

clean:
	$(RM) $(OBJ) $(EXE)

depend: $(SRCS)
	rm -f ./.depend
	$(CC) $(CFLAGS) -MM $(SRC) >> ./.depend

-include .depend

.PHONY: all clean
