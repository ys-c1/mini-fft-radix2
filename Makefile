CC ?= cc
CPPFLAGS ?= -Isrc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Wpedantic
LDFLAGS ?=
LDLIBS ?= -lm

BUILD_DIR := build
TARGET := $(BUILD_DIR)/test_fft
SOURCES := src/fft.c tests/test_fft.c
HEADERS := src/fft.h

.PHONY: all test clean

all: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	mkdir -p $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SOURCES) -o $@ $(LDFLAGS) $(LDLIBS)

test: $(TARGET)
	./$(TARGET)

clean:
	rm -rf $(BUILD_DIR)
