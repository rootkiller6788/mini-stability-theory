CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude -D_GNU_SOURCE
LDFLAGS = -lm
LIB = build/liblasalle.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_lasalle.c
TEST_BIN = build/test_lasalle
.PHONY: all builddir clean test
all: builddir $(LIB)
builddir:; @mkdir -p build
$(LIB): $(OBJS); $(AR) rcs $@ $^
build/%.o: src/%.c; $(CC) $(CFLAGS) -c $< -o $@
test: all; $(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -llasalle $(LDFLAGS) && ./$(TEST_BIN)
clean:; rm -rf build
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))
examples: all $(EXAMPLE_BINS); @echo "All examples built."
build/%: examples/%.c $(LIB); $(CC) $(CFLAGS) -o $@ $< -Lbuild -llasalle $(LDFLAGS)
demo: examples; @for e in $(EXAMPLE_BINS); do echo "=== $$(basename $$e) ===" && ./$$e; done
