HOME = ../..

include $(HOME)/Makefile.inc

LIB = libSolvers.a

EXEDIR = ./bin
OUT = interpol_example.exe

OBJDIR  = ./build
LIBDIR  = $(HOME)/lib
INCLUDE = $(HOME)/include

LDFLAGS = -L$(LIBDIR) -lSolvers
IDFLAGS = -I$(INCLUDE)

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=$(OBJDIR)/%.o)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(EXEDIR):
	mkdir -p $(EXEDIR)

$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(IDFLAGS)

main: $(OBJS)
	$(CXX) $(CXXFLAGS) $< -o $(EXEDIR)/$(OUT) $(LDFLAGS)

all: $(OBJDIR) $(EXEDIR) main

exe:
	mkdir -p result
	cd $(EXEDIR); ./$(OUT)

clean:
	rm -rf $(OBJDIR) $(EXEDIR) result

.PHONY: all exe clean