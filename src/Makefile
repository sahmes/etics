METHOD ?= SCF
LMAX ?= 6
NMAX ?= 10

GPUARCH ?= 35
NVCC ?= nvcc


CUDALIB  ?= /usr/local/cuda-6.5/lib64
MPILIB   ?= /usr/local
H5LIB    ?= /usr/lib/x86_64-linux-gnu
MPIINC   ?= /usr/local/include
BOOSTINC ?= /usr/local/sage/local/include
H5INC    ?= /usr/include

H5_CXX_OPTS = -I$(H5INC) -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_BSD_SOURCE -D_FORTIFY_SOURCE=2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security
CUDA_LD_OPTS = -L$(CUDALIB) -lcuda -lcudart
MPI_LD_OPTS = -L$(MPILIB)/lib -lmpicxx -Wl,-rpath -Wl,$(MPILIB)/lib -Wl,--enable-new-dtags -lmpi -lrt -lpthread
H5_LD_OPTS = -L$(H5LIB) $(H5LIB)/libhdf5_hl_cpp.a $(H5LIB)/libhdf5_cpp.a $(H5LIB)/libhdf5_hl.a $(H5LIB)/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lz -ldl -lm -Wl,-rpath -Wl,$(H5LIB)

ifeq ($(BOOSTINC), NOBOOST)
	BOOST = -DNOBOOST
else
	BOOST = -I$(BOOSTINC)
endif

CFLAGS    := -D$(METHOD) -DLMAX=$(LMAX) -DNMAX=$(NMAX) -I$(MPIINC) -O3
LIBS       = -lm $(CUDA_LD_OPTS) $(MPI_LD_OPTS)
CXXFLAGS   = $(CFLAGS) $(BOOST)
CUDAFLAGS  = -arch=sm_$(GPUARCH) $(CFLAGS)

CODELIB = libetics.a

# object files needed for the library and the standalone
CODEOBJS1 = integrate.o mathaux.o scf.o ic.o

# object files only needed for the standalone
CODEOBJS2 = io.o main.o

AR = ar ruv
RANLIB = ranlib

all: $(CODELIB) standalone

clean:
	$(RM) *.o *.a etics

$(CODELIB): $(CODEOBJS1)
	$(NVCC) $(CUDAFLAGS) -dlink $(CODEOBJS1) -o dlink.o
	$(AR) $@ $(CODEOBJS1) dlink.o
	$(RANLIB) $@

library: $(CODELIB)

standalone: $(CODEOBJS1) $(CODEOBJS2)
	$(NVCC) $(CUDAFLAGS) -dlink $(CODEOBJS1) $(CODEOBJS2) -o dlink.o
	$(CXX) -o etics $(CODEOBJS1) $(CODEOBJS2) dlink.o $(LIBS) $(H5_LD_OPTS)

.SUFFIXES: .cu .cpp .o

# only io.cpp should be compiled with all the weird HDF5 flags
io.o:
	$(CXX) $(CXXFLAGS) $(H5_CXX_OPTS) -c -o io.o io.cpp

.cu.o: $<
	$(NVCC) $(CUDAFLAGS) -dc -o $@ $<

.cpp.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<
