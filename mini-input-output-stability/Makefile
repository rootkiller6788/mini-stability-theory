CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
AR = ar
LIB = build/libiostability.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_ios.c
TEST_BIN = build/test_ios
.PHONY: all builddir clean test
all: builddir $(LIB)
builddir:
	mkdir -p build
$(LIB): $(OBJS)
	$(AR) rcs $@ $^
build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@
test: all
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -liostability $(LDFLAGS)
	./$(TEST_BIN)
clean:
	rm -rf build
