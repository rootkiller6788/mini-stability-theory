CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libfinite_time_stability.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
EXAMPLES = $(wildcard examples/example*.c)
EXE = $(patsubst examples/%.c,build/%.exe,$(EXAMPLES))

.PHONY: all test clean examples

all: builddir $(LIB)

builddir:
	@mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c | builddir
	$(CC) $(CFLAGS) -c $< -o $@

test: $(LIB)
	$(CC) $(CFLAGS) tests/test_fts.c -Lbuild -lfinite_time_stability $(LDFLAGS) -o build/test.exe
	./build/test.exe

examples: $(LIB) $(EXE)
	@for exe in $(EXE); do echo "Running $$exe..."; ./$$exe; done

build/%.exe: examples/%.c | builddir
	$(CC) $(CFLAGS) $< -Lbuild -lfinite_time_stability $(LDFLAGS) -o $@

clean:
	rm -rf build
