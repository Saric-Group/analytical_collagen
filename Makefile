## Variables
# Compiler
CC = g++

# Options for compiler
CFLAGS = -march=native -O2 -std=c++11 -Wall -c

# cpp files & objects
SRC = $(wildcard src/*.cpp)
OBJ = $(addprefix build/,$(notdir $(SRC:.cpp=.o)))

# Directories
OBJDIR = build
SRCDIR = src


### Make rules

all: main

# Creates the final programm from the object files in build
main: $(OBJ)
	$(CC) -o $@ $^
#-lm -lgsl -lgslcblas -lgd

# Creates all object files in build from the .cpp files in src
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(OBJDIR)/.dirstamp
	$(CC)	$(CFLAGS) -o $@ $<

# Creates the build directory
$(OBJDIR)/.dirstamp:
	mkdir -p $(OBJDIR)
	touch $@

# Cleans the build directory and removes the programm
clean:
	rm -rf $(OBJDIR) main
