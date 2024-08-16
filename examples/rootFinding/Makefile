HOME = ../..

include $(HOME)/Makefile.inc

LIB = libSolvers.a

EXEDIR = ./bin
OUT = root_example.exe

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

main: $(OBJS) $(LIBDIR)/$(LIB)
	$(CXX) $(CXXFLAGS) $< -o $(EXEDIR)/$(OUT) $(LDFLAGS)

all: $(OBJDIR) $(EXEDIR) main

exe:
	cd $(EXEDIR); ./$(OUT)

clean:
	rm -rf $(OBJDIR) $(EXEDIR)

.PHONY: all exe clean