#=============================================================================
# Makefile: for Waxman generation code (the stuff I wrote)
#=============================================================================

#-----------------------------------------------------------------------------

# ifeq ($(OSTYPE),linux)
# 	OS = -DLINUX
# 	TARGET = MexFastSERN.mexa64
# endif

# -L/usr/local/matlab/bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++
VPATH = ../

# MATLAB = true
ifdef MATLAB
    TARGET = conSERN.mexa64
    # your local matlab install location
    INCLUDES = -I/usr/local/matlab/extern/include/  -I/usr/local/MATLAB/R2016b/extern/include/
    CFLAGS= -O0  -Ofast -DMATLAB -fPIC -lm -lpthread  $(INCLUDES) 
    MX=mex -cxx -O -largeArrayDims CFLAGS='$$CFLAGS -DMATLAB' $(INCLUDES)
    CC=gcc-4.7
    HEADERS =  $(shell ls *.h)
    CFILES =  $(shell ls *.c)
    OBJ =  $(CFILES:.c=.o) 
else
    # just make the non-matlab bits
    TARGET = FastSERN.so
    INCLUDES = 
    CFLAGS= -O0  -Ofast -fPIC -lm -lpthread  $(INCLUDES)
    CC=gcc
    HEADERS = $(shell ls *.h | grep -v conSERN)
    CFILES = $(shell ls *.c | grep -v conSERN)
    OBJ = $(CFILES:.c=.o)
endif 


all: $(TARGET) 

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

conSERN.mexa64: conSERN.o $(OBJ)
	$(MX) $^

FastSERN.so: FastSERN.o $(OBJ)
	$(CC) -shared -o $@ $^ $(CFLAGS)

FORCE:

# standard cleanup etc
.PHONY : clean
clean: 
	- $(RM) *.mex* *.dll *.oct *.o

.PHONY : distclean
distclean: 
	- $(RM) *.mex* *.dll *.oct *.o

