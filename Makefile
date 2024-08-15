include ./makefile.inc

OBJDIR = ./build

$(OBJDIR) : 
	if not exist $(OBJDIR) mkdir $(OBJDIR)


lib: $(OBJDIR)
	$(MAKE) -C src all || true

examples:
	$(MAKE) -C examples all || true

root:
	$(MAKE) -C examples root || true

interpolation:
	$(MAKE) -C examples interpolation || true

all: lib examples

exe:
	$(MAKE) -C examples exe

clean:
	$(MAKE) -C src clean || true
	$(MAKE) -C examples clean || true
	rm -rf core.*
	rm -rf run/core.*

.PHONY: lib examples all root interplation clean