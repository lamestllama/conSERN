#=============================================================================
# Makefile: for Waxman generation code (the stuff I wrote)
#=============================================================================

#-----------------------------------------------------------------------------

# ifeq ($(OSTYPE),linux)
# 	OS = -DLINUX
# 	TARGET = MexFastSERN.mexa64
# endif 
TARGET = conSERN.mexa64
INCLUDES =  -I/usr/local/matlab/extern/include/  -I/opt/local/MATLAB/R2014a/extern/include/

CC=gcc
CFLAGS= -O0  -Ofast   -DMATLAB -fPIC -lm -lpthread  $(INCLUDES) 
MX=mex -cxx -O -largeArrayDims CFLAGS='$$CFLAGS -DMATLAB' $(INCLUDES) 

# -L/usr/local/matlab/bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++
VPATH = ../

HEADERS =  $(shell ls *.h)
CFILES =  $(shell ls *.c)
OBJ =  $(CFILES:.c=.o)

all: $(TARGET) 

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(TARGET): conSERN.o $(OBJ)
	$(MX) $^

FORCE:

# standard cleanup etc
.PHONY : clean
clean: 
	- $(RM) *.mex* *.dll *.oct *.o

.PHONY : distclean
distclean: 
	- $(RM) *.mex* *.dll *.oct *.o

