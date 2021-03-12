## Variables
# Compiler
CC = g++

# Options for compiler
CFLAGS = -march=native -O2 -I$(SRCDIR)/config4cpp/include/ -std=c++11 -Wall -Wno-deprecated -c

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
	$(CC) -o $@ $^ -L$(SRCDIR)/config4cpp/lib -lconfig4cpp
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
