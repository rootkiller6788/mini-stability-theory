CC=gcc
CFLAGS=-std=c11 -Wall -Wextra -O2 -Iinclude
LIB=libabs.a
SRCS=$(wildcard src/*.c)
OBJS=$(SRCS:.c=.o)
all: $(LIB)
$(LIB): $(OBJS)
	ar rcs $@ $^
src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@
test: $(LIB)
	for t in tests/test_*.c; do $(CC) $(CFLAGS) $$t -L. -labs -lm -o $$t.out && ./$$t.out; rm -f $$t.out; done
clean:
	rm -f src/*.o $(LIB) tests/*.out
