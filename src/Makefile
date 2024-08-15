HOME = ..

include $(HOME)/Makefile.inc

LIB = libSolvers.a

OBJDIR  = $(HOME)/build
LIBDIR  = $(HOME)/lib
INCLUDE = $(HOME)/include

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=$(OBJDIR)/%.o)

# Display variable values
$(info Compile with: $(CXX) $(CXXFLAGS))
$(info SRCS = $(SRCS))
$(info OBJS = $(OBJS))
#$(info OBJDIR = $(OBJDIR))

$(OBJDIR):
#   Linux
	mkdir -p $(OBJDIR) 
#	Windows
#	if not exist $(OBJDIR) mkdir $(OBJDIR)

$(LIBDIR):
#   Linux
	mkdir -p $(LIBDIR)  #Linux
#	Windows
#	if not exist $(LIBDIR) mkdir $(LIBDIR) 

$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
all: $(OBJDIR) $(LIBDIR) $(LIB)

$(LIB): $(OBJS)
	ar rcs $(LIBDIR)/$(LIB) $(OBJS)

clean:
#   Linux
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
#   Window
#	if exist $(OBJDIR) rmdir /S /Q $(OBJDIR)
#	if exist $(LIBDIR) rmdir /S /Q $(LIBDIR)

.PHONY: all clean 