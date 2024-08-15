HOME = ..

include $(HOME)/makefile.inc

root:
	$(MAKE) -C rootFinding all || true

interpolation:
	$(MAKE) -C numericalInterpolation all || true

all: root interpolation

exe:
	if [ -f ./rootFinding/bin/root.out ]; then $(MAKE) -C rootFinding exe; fi
	if [ -f ./numericalInterpolation/bin/interpol.out ]; then $(MAKE) -C numericalInterpolation exe; fi

clean:
	$(MAKE) -C rootFinding clean || true
	$(MAKE) -C numericalInterpolation clean || true

.PHONY: root interplation all exe clean