# Compiler and flags
CC = gcc
CFLAGS = -O3 -g -fno-omit-frame-pointer 

LDFLAGS = -ljson-c -lhts -lz -llz4 -lbrotlidec -lbrotlienc

# Source and object files
SRCS = main.c reader.c writer.c
OBJS = $(SRCS:.c=.o)

# Output executable
TARGET = main

# Default rule
all: $(TARGET)

# Linking
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Compile each .c to .o
%.o: %.c defs.h meta.h utils.h
	$(CC) $(CFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean