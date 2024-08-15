PROJECT_ROOT = ../..
include $(PROJECT_ROOT)/makefile.inc

LIB = libSolvers.a

EXEDIR = ./bin
out = interpol.out

OBJDIR = ./obj
LIBDIR = $(PROJECT_ROOT)/lib
INCLUDEDIR = $(PROJECT_ROOT)/include

LDFLAGS = -L$(LIBDIR) -lSolvers
IDFLAGS = -I$(INCLUDEDIR)

SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp, $(OBJDIR)/%.o, $(SRCS))

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(EXEDIR):
	mkdir -p $(EXEDIR)

$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(IDFLAGS)

main: $(OBJS) $(LIBDIR)/$(LIB)
	$(CXX) $(CXXFLAGS) $< -o $(EXEDIR)/$(out) $(LDFLAGS)

all: $(OBJDIR) $(EXEDIR) main

exe:
	mkdir -p result
	cd $(EXEDIR); ./$(out)

clean:
	rm -rf $(OBJDIR) $(EXEDIR) result

.PHONY: all exe clean