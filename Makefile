# Makefile for SIP Solver

# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2
LDFLAGS = -lm

# Source files
SRC = main.c SIP.c heat_equation.c analytical.c comparison.c vtk_output.c
OBJ = $(SRC:.c=.o)

# Executable name
EXEC = sip_solver

# Default target
.PHONY: all clean

all: $(EXEC)

# Linking
$(EXEC): $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $@
	@echo "Build successful: $(EXEC)"

# Compilation rule
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJ) $(EXEC)
	@echo "Clean completed"

# Help target
help:
	@echo "Available targets:"
	@echo "  make all    - Build the SIP solver"
	@echo "  make clean  - Remove build artifacts"
	@echo "  make help   - Show this help message"
