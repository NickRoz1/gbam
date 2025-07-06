# Compiler and flags
CC = gcc
# CFLAGS = -O3 -g -fno-omit-frame-pointer -fsanitize=address -fopenmp

# Default
# LDFLAGS = -ljson-c -lhts -lz -llz4 -lbrotlidec -lbrotlienc -lzstd -fopenmp
# Tell Clang to use libomp's headers and link the runtime correctly

# Needed for mac.
CFLAGS = -O3 -g -fno-omit-frame-pointer -fsanitize=address \
         -Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include

LDFLAGS = -L/usr/local/opt/libomp/lib -lomp -ljson-c -lhts -lz -llz4 -lbrotlidec -lbrotlienc -lzstd

# Needed for linux
CFLAGS = -O3 -g -fno-omit-frame-pointer -I$(HOME)/.local/include
LDFLAGS = -L$(HOME)/.local/lib -ljson-c -lhts -lz -llz4 -lbrotlidec -lbrotlienc -lzstd

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