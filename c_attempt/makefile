# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -O0 -g -Wfatal-errors \
-g3 -ggdb -gdwarf-4 -gsplit-dwarf -O0 -fno-inline -fno-omit-frame-pointer \
-fsanitize=address,undefined -fno-eliminate-unused-debug-types

LDFLAGS = -ljson-c -lhts

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