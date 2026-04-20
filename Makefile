# Makefile

# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -Wextra -std=c99 -lm

# Source files
SRC = $(wildcard *.c)

# Object files
OBJ = $(SRC:.c=.o)

# Executable name
EXEC = my_program

# Default target
all: $(EXEC)

# Linking
$(EXEC): $(OBJ)
	$(CC) $(OBJ) -o $@

# Compilation rule
%.o: %.c
	$(CC) $(CFLAGS) -c $<

# Clean target
clean:
	rm -f $(OBJ) $(EXEC)