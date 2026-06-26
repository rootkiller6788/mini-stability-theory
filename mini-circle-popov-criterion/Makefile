CC      = gcc
AR      = ar
CFLAGS  = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB     = build/libcircle.a

SRCS     = $(wildcard src/*.c)
OBJS     = $(patsubst src/%.c,build/%.o,$(SRCS))

TEST_SRCS = $(wildcard tests/test_*.c)
TEST_BINS = $(patsubst tests/%.c,build/%,$(TEST_SRCS))

EXAMPLE_SRCS = $(wildcard examples/example*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))

.PHONY: all builddir test examples clean

all: builddir $(LIB)

builddir:
	@mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: builddir $(LIB) $(TEST_BINS)
	@passed=0; failed=0; \
	for t in $(TEST_BINS); do \
		echo "Running $$t..."; \
		if ./$$t; then passed=$$((passed+1)); \
		else failed=$$((failed+1)); fi; \
	done; \
	echo "=== Test Results: $$passed passed, $$failed failed ==="; \
	if [ $$failed -gt 0 ]; then exit 1; fi

build/test_%: tests/test_%.c builddir $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lcircle $(LDFLAGS)

examples: builddir $(LIB) $(EXAMPLE_BINS)
	@for e in $(EXAMPLE_BINS); do \
		echo "=== Running $$e ==="; ./$$e; \
	done

build/example%: examples/example%.c builddir $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lcircle $(LDFLAGS)

clean:
	rm -rf build
